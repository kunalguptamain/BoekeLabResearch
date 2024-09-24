from standardized_dna_diffusion_pipeline.model import Unet, SequenceDataset, train
from torch.utils.data import DataLoader
import torch
from torch.optim import Adam

data_path = "C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/processedMouseDNAData/train.h5"

batch_size = 1024
device = "cuda" if torch.cuda.is_available() else "cpu"
image_size = (4, 1024)
channels = 1
num_classes = 4

dataset = SequenceDataset(data_path)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

model = Unet(
    dim=image_size,
    channels=channels,
    num_classes=num_classes,
    filter_mults=(1, 2, 4, 8, 16),
    filters = 64,
)
model.to(device)

optimizer = Adam(model.parameters(), lr=4e-5)

train(
    model=model,
    optimizer=optimizer,
    epochs=300,
    batch_size=1024,
    data_loader=dataloader,
    weight_save_path="C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/weights.npy",
    device=device,
)