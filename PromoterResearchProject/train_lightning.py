from standardized_dna_diffusion_pipeline.model_lightning import Unet, SequenceDataset
from torch.utils.data import DataLoader
import torch
from torch.optim import Adam
from pytorch_lightning import Trainer

data_path = "C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/processedMouseDNAData/train.h5"
weight_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/testweight.pth"

batch_size = 16
image_size = (4, 1024)
channels = 1
num_classes = 4

dataset = SequenceDataset(data_path)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

model = Unet(
    dim=image_size,
    channels=channels,
    num_classes=num_classes,
    filter_mults=(1,),
    filters = 8,
)

trainer = Trainer(limit_train_batches=0.2, max_epochs=1, precision="16-mixed", accelerator="auto", strategy="ddp")
trainer.fit(model, dataloader)

torch.save(model.state_dict(), weight_save_path)