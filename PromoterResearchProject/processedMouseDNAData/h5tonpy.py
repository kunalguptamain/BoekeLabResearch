import h5py
import numpy as np

data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/processedMouseDNAData/train.h5'  # Replace with your actual data path
pixel_values_path = 'pixel_values.npy'
labels_path = 'labels.npy'

with h5py.File(data_path, 'r') as h5f:
    pixel_values = h5f['train_data'][:]
    labels = h5f['labels'][:]

# Save to .npy files
np.save(pixel_values_path, pixel_values)
np.save(labels_path, labels)
