import h5py
import numpy as np
from collections import Counter

data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/processedMouseDNAData/train.h5'  # Replace with your actual data path
pixel_values_path = 'pixel_values.npy'
labels_path = 'labels.npy'

with h5py.File(data_path, 'r') as h5f:
    pixel_values = h5f['train_data'][:]
    labels = h5f['labels'][:]

label_counts = Counter(labels)

# Print the label counts
for label, count in label_counts.items():
    print(f'Label {label}: {count} occurrences')

# Save to .npy files
np.save(pixel_values_path, pixel_values)
np.save(labels_path, labels)
