import math
from inspect import isfunction
from functools import partial
import pytorch_lightning as pl
from pytorch_lightning import Trainer
import torch
# %matplotlib inline
from tqdm.auto import tqdm
from einops import rearrange, reduce
from einops.layers.torch import Rearrange
import multiprocessing

from torch import nn, einsum
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from torchvision.transforms import Compose, ToTensor, Lambda, ToPILImage, CenterCrop, Resize
import h5py
import math
import os
import numpy as np
import torch._dynamo

def linear_beta_schedule(timesteps):
    beta_start = 0.0001
    beta_end = 0.02
    return torch.linspace(beta_start, beta_end, timesteps)

timesteps = 300

# define beta schedule
betas = linear_beta_schedule(timesteps=timesteps)

# define alphas
alphas = 1. - betas
alphas_cumprod = torch.cumprod(alphas, axis=0)
alphas_cumprod_prev = F.pad(alphas_cumprod[:-1], (1, 0), value=1.0)
sqrt_recip_alphas = torch.sqrt(1.0 / alphas)

# calculations for diffusion q(x_t | x_{t-1}) and others
sqrt_alphas_cumprod = torch.sqrt(alphas_cumprod)
sqrt_one_minus_alphas_cumprod = torch.sqrt(1. - alphas_cumprod)

# calculations for posterior q(x_{t-1} | x_t, x_0)
posterior_variance = betas * (1. - alphas_cumprod_prev) / (1. - alphas_cumprod)

def extract(a, t, x_shape):
    batch_size = t.shape[0]
    out = a.gather(-1, t.cpu())
    return out.reshape(batch_size, *((1,) * (len(x_shape) - 1))).to(t.device)

# forward diffusion (using the nice property)
def q_sample(x_start, t, noise=None):
    if noise is None:
        noise = torch.randn_like(x_start)

    sqrt_alphas_cumprod_t = extract(sqrt_alphas_cumprod, t, x_start.shape)
    sqrt_one_minus_alphas_cumprod_t = extract(
        sqrt_one_minus_alphas_cumprod, t, x_start.shape
    )

    return sqrt_alphas_cumprod_t * x_start + sqrt_one_minus_alphas_cumprod_t * noise

reverse_transform = Compose([
     Lambda(lambda t: (t + 1) / 2),
     Lambda(lambda t: t.permute(1, 2, 0)), # CHW to HWC
     Lambda(lambda t: t * 255.),
     Lambda(lambda t: t.numpy().astype(np.uint8)),
     ToPILImage(),
])

def get_noisy_image(x_start, t):
  # add noise
  x_noisy = q_sample(x_start, t=t)

  # turn back into PIL image
  noisy_image = reverse_transform(x_noisy.squeeze())

  return noisy_image

def p_losses(denoise_model, x_start, t, labels, noise=None, loss_type="l1"):
    if noise is None:
        noise = torch.randn_like(x_start)

    x_noisy = q_sample(x_start=x_start, t=t, noise=noise)
    predicted_noise = denoise_model(x_noisy, t, labels)

    if loss_type == 'l1':
        loss = F.l1_loss(noise, predicted_noise)
    elif loss_type == 'l2':
        loss = F.mse_loss(noise, predicted_noise)
    elif loss_type == "huber":
        loss = F.smooth_l1_loss(noise, predicted_noise)
    else:
        raise NotImplementedError()

    return loss

@torch.no_grad()
def p_sample(model, x, class_vector, t, t_index):
    betas_t = extract(betas, t, x.shape)
    sqrt_one_minus_alphas_cumprod_t = extract(
        sqrt_one_minus_alphas_cumprod, t, x.shape
    )
    sqrt_recip_alphas_t = extract(sqrt_recip_alphas, t, x.shape)

    # Equation 11 in the paper
    # Use our model (noise predictor) to predict the mean
    pred = model(x, t, class_vector)
    model_mean = sqrt_recip_alphas_t * (
        x - betas_t *  pred/ sqrt_one_minus_alphas_cumprod_t
    )

    if t_index == 0:
        return model_mean
    else:
        posterior_variance_t = extract(posterior_variance, t, x.shape)
        noise = torch.randn_like(x)
        # Algorithm 2 line 4:
        return model_mean + torch.sqrt(posterior_variance_t) * noise

# Algorithm 2 (including returning all images)
@torch.no_grad()
def p_sample_loop(model, shape, class_vector):
    device = next(model.parameters()).device

    b = shape[0]
    # start from pure noise (for each example in the batch)
    img = torch.randn(shape, device=device)
    imgs = []

    for i in tqdm(reversed(range(0, timesteps)), desc='sampling loop time step', total=timesteps):
        img = p_sample(model, img, class_vector, torch.full((b,), i, device=device, dtype=torch.long), i)
        imgs.append(img.cpu().numpy())
    return imgs

@torch.no_grad()
def sample(model, image_size, class_vector, batch_size=16, channels=3):
    return p_sample_loop(model, shape=(batch_size, channels, image_size[0], image_size[1]), class_vector=class_vector)


def exists(x):
    return x is not None

def default(val, d):
    if exists(val):
        return val
    return d() if isfunction(d) else d


def num_to_groups(num, divisor):
    groups = num // divisor
    remainder = num % divisor
    arr = [divisor] * groups
    if remainder > 0:
        arr.append(remainder)
    return arr


class Residual(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x, *args, **kwargs):
        return self.fn(x, *args, **kwargs) + x



def Upsample(dim, dim_out=None, normal = True):
    return nn.Sequential(
        nn.Upsample(scale_factor=2 if normal else (1, 2), mode="nearest"),
        nn.Conv2d(dim, default(dim_out, dim), 3, padding=1),
    )


def Downsample(dim, dim_out=None, height=4):
    # No More Strided Convolutions or Pooling
    P1 = 2 if height > 1 else 1
    P2 = 2

    return nn.Sequential(
        Rearrange("b c (h p1) (w p2) -> b (c p1 p2) h w", p1=P1, p2=P2),
        nn.Conv2d(dim * (P1 * P2), default(dim_out, dim), 1),
    )

class SinusoidalPositionEmbeddings(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.dim = dim

    def forward(self, time):
        device = time.device
        half_dim = self.dim // 2
        embeddings = math.log(10000) / (half_dim - 1)
        embeddings = torch.exp(torch.arange(half_dim, device=device) * -embeddings)
        embeddings = time[:, None] * embeddings[None, :]
        embeddings = torch.cat((embeddings.sin(), embeddings.cos()), dim=-1)
        return embeddings

class WeightStandardizedConv2d(nn.Conv2d):
    """
    https://arxiv.org/abs/1903.10520
    weight standardization purportedly works synergistically with group normalization
    """

    def forward(self, x):
        eps = 1e-5 if x.dtype == torch.float32 else 1e-3

        # Calculate mean and variance over the desired dimensions without functools.partial
        weight = self.weight
        mean = weight.mean(dim=(1, 2, 3), keepdim=True)  # Compute mean across spatial dimensions
        var = weight.var(dim=(1, 2, 3), unbiased=False, keepdim=True)  # Compute variance directly

        # Normalize weight
        normalized_weight = (weight - mean) / torch.sqrt(var + eps)

        # Perform the convolution with normalized weights
        return F.conv2d(
            x,
            normalized_weight,
            self.bias,
            self.stride,
            self.padding,
            self.dilation,
            self.groups,
        )


class Block(nn.Module):
    def __init__(self, dim, dim_out, groups=8):
        super().__init__()
        self.proj = WeightStandardizedConv2d(dim, dim_out, 3, padding=1)
        self.norm = nn.GroupNorm(groups, dim_out)
        self.act = nn.SiLU()

    def forward(self, x, scale_shift=None):
        x = self.proj(x)
        x = self.norm(x)

        if exists(scale_shift):
            scale, shift = scale_shift
            x = x * (scale + 1) + shift

        x = self.act(x)
        return x


class ResnetBlock(nn.Module):
    """https://arxiv.org/abs/1512.03385"""

    def __init__(self, dim, dim_out, *, time_emb_dim=None, groups=8):
        super().__init__()
        self.mlp = (
            nn.Sequential(nn.SiLU(), nn.Linear(time_emb_dim, dim_out * 2))
            if exists(time_emb_dim)
            else None
        )

        self.block1 = Block(dim, dim_out, groups=groups)
        self.block2 = Block(dim_out, dim_out, groups=groups)
        self.res_conv = nn.Conv2d(dim, dim_out, 1) if dim != dim_out else nn.Identity()

    def forward(self, x, time_emb=None):
        scale_shift = None
        if exists(self.mlp) and exists(time_emb):
            time_emb = self.mlp(time_emb)
            time_emb = rearrange(time_emb, "b c -> b c 1 1")
            scale_shift = time_emb.chunk(2, dim=1)

        h = self.block1(x, scale_shift=scale_shift)
        h = self.block2(h)
        return h + self.res_conv(x)

class Attention(nn.Module):
    def __init__(self, dim, heads=4, dim_head=32):
        super().__init__()
        self.scale = dim_head**-0.5
        self.heads = heads
        hidden_dim = dim_head * heads
        self.to_qkv = nn.Conv2d(dim, hidden_dim * 3, 1, bias=False)
        self.to_out = nn.Conv2d(hidden_dim, dim, 1)

    def forward(self, x):
        b, c, h, w = x.shape
        qkv = self.to_qkv(x).chunk(3, dim=1)
        q, k, v = map(
            lambda t: rearrange(t, "b (h c) x y -> b h c (x y)", h=self.heads), qkv
        )
        q = q * self.scale

        sim = einsum("b h d i, b h d j -> b h i j", q, k)
        sim = sim - sim.amax(dim=-1, keepdim=True).detach()
        attn = sim.softmax(dim=-1)

        out = einsum("b h i j, b h d j -> b h i d", attn, v)
        out = rearrange(out, "b h (x y) d -> b (h d) x y", x=h, y=w)
        return self.to_out(out)

class LinearAttention(nn.Module):
    def __init__(self, dim, heads=4, dim_head=32):
        super().__init__()
        self.scale = dim_head**-0.5
        self.heads = heads
        hidden_dim = dim_head * heads
        self.to_qkv = nn.Conv2d(dim, hidden_dim * 3, 1, bias=False)

        self.to_out = nn.Sequential(nn.Conv2d(hidden_dim, dim, 1),
                                    nn.GroupNorm(1, dim))

    def forward(self, x):
        b, c, h, w = x.shape
        qkv = self.to_qkv(x).chunk(3, dim=1)
        q, k, v = map(
            lambda t: rearrange(t, "b (h c) x y -> b h c (x y)", h=self.heads), qkv
        )

        q = q.softmax(dim=-2)
        k = k.softmax(dim=-1)

        q = q * self.scale
        context = torch.einsum("b h d n, b h e n -> b h d e", k, v)

        out = torch.einsum("b h d e, b h d n -> b h e n", context, q)
        out = rearrange(out, "b h c (x y) -> b (h c) x y", h=self.heads, x=h, y=w)
        return self.to_out(out)

class PreNorm(nn.Module):
    def __init__(self, dim, fn):
        super().__init__()
        self.fn = fn
        self.norm = nn.GroupNorm(1, dim)

    def forward(self, x):
        x = self.norm(x)
        return self.fn(x)

#change res to tuple
class ClassConditioning(nn.Module):
    def __init__(self, in_dims, out_dims, num_channels=1):
        super().__init__()
        self.block = nn.Sequential(
            nn.Linear(in_dims, out_dims[0] * out_dims[1] * num_channels),
            nn.SiLU(),
            nn.Unflatten(1, (num_channels, out_dims[0], out_dims[1]))
        )
    def forward(self, x):
        return self.block(x)

#Change to tuple
class Unet(pl.LightningModule):
    def __init__(
        self,
        dim,
        num_classes,
        class_embed_dim = 64,
        filters = None,
        init_filters=None,
        out_filters=None,
        filter_mults=(1, 2, 4, 8),
        channels=3,
        self_condition=False,
        resnet_block_groups=4,
    ):
        super().__init__()

        # determine dimensions
        self.dim = dim
        self.channels = channels #1
        self.self_condition = self_condition
        input_channels = channels * (2 if self_condition else 1) #1

        self.class_embeddings = nn.Embedding(num_classes, class_embed_dim) # 10 -> 64

        self.filters = default(filters, int(math.sqrt(dim[0] * dim[1])))
        init_filters = default(init_filters, self.filters) # (4, 200)
        self.init_conv = nn.Conv2d(input_channels, init_filters, 1, padding=0) # int(math.sqrt(4, 200)) instead of init_dim

        dims = [init_filters, *map(lambda m: self.filters * m, filter_mults)] # int(math.sqrt(4, 200)) instead of init_dim and dim
        in_out = list(zip(dims[:-1], dims[1:]))

        block_klass = partial(ResnetBlock, groups=resnet_block_groups)

        # time embeddings
        time_dim = self.filters * 4 #int(math.sqrt(4, 200)) * 4

        self.time_mlp = nn.Sequential(
            SinusoidalPositionEmbeddings(self.filters),
            nn.Linear(self.filters, time_dim),
            nn.GELU(),
            nn.Linear(time_dim, time_dim),
        )

        # layers
        self.downs = nn.ModuleList([])
        self.ups = nn.ModuleList([])
        num_resolutions = len(in_out)
        count = 0
        now_res = dim

        for ind, (filters_in, filters_out) in enumerate(in_out):
            print(now_res)
            print(filters_out)
            is_last = ind >= (num_resolutions - 1)
            if(ind == 3): print(filters_in, filters_out, now_res[0])
            self.downs.append(
                nn.ModuleList(
                    [
                        ClassConditioning(class_embed_dim, now_res),
                        block_klass(filters_in + 1, filters_in, time_emb_dim=time_dim),
                        block_klass(filters_in, filters_in, time_emb_dim=time_dim),
                        Residual(PreNorm(filters_in, LinearAttention(filters_in))),
                        Downsample(filters_in, filters_out, now_res[0])
                        if not is_last
                        else nn.Conv2d(filters_in, filters_out, 3, padding=1),
                    ]
                )
            )
            if not is_last and now_res[0] == 1: count += 1
            if not is_last: now_res = tuple(max(res_dim // 2, 1) for res_dim in now_res)

        mid_dim = dims[-1]
        self.mid_class_conditioning = ClassConditioning(class_embed_dim, now_res)
        self.mid_block1 = block_klass(mid_dim + 1, mid_dim, time_emb_dim=time_dim)
        self.mid_attn = Residual(PreNorm(mid_dim, Attention(mid_dim)))
        self.mid_block2 = block_klass(mid_dim, mid_dim, time_emb_dim=time_dim)
        max_count = count

        for ind, (filters_in, filters_out) in enumerate(reversed(in_out)):
            is_last = ind == (len(in_out) - 1)
            print(now_res)
            print((now_res[0], now_res[1]))
            print(filters_out)
            self.ups.append(
                nn.ModuleList(
                    [
                        ClassConditioning(class_embed_dim, (max(1, now_res[0] // (2**max_count)), now_res[1])),
                        block_klass(filters_out + filters_in + 1, filters_out, time_emb_dim=time_dim),
                        block_klass(filters_out + filters_in, filters_out, time_emb_dim=time_dim),
                        Residual(PreNorm(filters_out, LinearAttention(filters_out))),
                        Upsample(filters_out, filters_in, count <= 0)
                        if not is_last
                        else nn.Conv2d(filters_out, filters_in, 3, padding=1),
                    ]
                )
            )
            count -= 1
            now_res = tuple(res_dim * 2 for res_dim in now_res)

        self.out_dim = default(out_filters, channels)
        self.final_res_block = block_klass(self.filters * 2, self.filters, time_emb_dim=time_dim)
        self.final_conv = nn.Conv2d(self.filters, self.out_dim, 1)

    def forward(self, x, time, class_vector, x_self_cond=None, prints = False):

        if self.self_condition:
            x_self_cond = default(x_self_cond, lambda: torch.zeros_like(x))
            x = torch.cat((x_self_cond, x), dim=1)

        x = self.init_conv(x)
        r = x.clone()

        t = self.time_mlp(time)
        class_vector = self.class_embeddings(class_vector)
        h = []

        for class_conditioning, block1, block2, attn, downsample in self.downs:
            cv = class_conditioning(class_vector)
            x = torch.concat([x, cv], dim = 1)
            x = block1(x, t)
            h.append(x)
            x = block2(x, t)
            x = attn(x)
            h.append(x)
            x = downsample(x)

        cv = self.mid_class_conditioning(class_vector)
        x = torch.concat([x, cv], dim = 1)
        x = self.mid_block1(x, t)
        x = self.mid_attn(x)
        x = self.mid_block2(x, t)

        for class_conditioning, block1, block2, attn, upsample in self.ups:
            cv = class_conditioning(class_vector)
            x = torch.concat([x, cv], dim = 1)
            x = torch.cat((x, h.pop()), dim=1)
            x = block1(x, t)
            x = torch.cat((x, h.pop()), dim=1)
            x = block2(x, t)
            x = attn(x)
            x = upsample(x)


        x = torch.cat((x, r), dim=1)

        x = self.final_res_block(x, t)
        prints = False
        return self.final_conv(x)

    def training_step(self, batch, batch_idx):
      batch_size = batch["pixel_values"].shape[0]
      labels = batch["label"]
      batch = batch["pixel_values"]
      t = torch.randint(0, timesteps, (batch_size,),).long()
      loss = p_losses(self.model, batch, t, labels, loss_type="huber")
      self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True)
      return loss

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=4e-5)
        return optimizer

    def q_sample(self, x_start, t, noise):
        return x_start + noise * torch.sqrt(self.betas[t])

    def extract(self, a, t, shape):
        return a[t].expand(shape)


def generate_promoters(
    model,
    amount,
    class_number,
    sequence_shape,
    device,
):
    return_sequences = []

    one_hot_encode_val = ["A","C","G","T"]
    def find_closest_key(col):
        a = [
            sum(abs(col - np.array([1, 0, 0, 0]))),
            sum(abs(col - np.array([0, 1, 0, 0]))),
            sum(abs(col - np.array([0, 0, 1, 0]))),
            sum(abs(col - np.array([0, 0, 0, 1])))
        ]
        return a.index(min(a))

    def normalize(col):
        total = sum(col)
        new = col/total
        return [round(x) for x in new]

    samples = sample(model=model,\
        image_size=torch.tensor(sequence_shape, device=device), \
        class_vector = torch.tensor(np.array([class_number] * amount), device=device), \
        batch_size=torch.tensor(amount, device=device),\
        channels=torch.tensor(1, device=device)\
        )

    plt.imshow(np.apply_along_axis(normalize, 0, samples[-1][0].reshape(sequence_shape[0], sequence_shape[1]))[:,900:])

    for i in range(amount):
        a = samples[-1][i].reshape(sequence_shape[0], sequence_shape[1])
        a = np.apply_along_axis(normalize, 0, a).tolist()
        result_array = np.apply_along_axis(find_closest_key, 0, a).tolist()
        string = ''.join([one_hot_encode_val[s] for s in result_array])
        return_sequences.append(string)

    return return_sequences

class SequenceDataset(Dataset):
    def __init__(self, data_path):

        with h5py.File(data_path, 'r') as h5f:
            self.pixel_values = h5f['train_data'][:]
            self.labels = h5f['labels'][:]

        self.pixel_values = np.expand_dims(self.pixel_values, axis=1)
        print("Data loaded successfully")

    def __len__(self):
        return len(self.pixel_values)

    def __getitem__(self, idx):
        return {"pixel_values": self.pixel_values[idx], "label": self.labels[idx]}
