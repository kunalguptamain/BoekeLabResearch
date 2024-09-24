from pydantic import BaseModel
import numpy as np
from enum import Enum
from typing import List, Optional
import os
import json
import standardized_dna_diffusion_pipeline.helper_functions as bio_methods
import h5py
import gzip

class SequenceSet:

    class Sequence(BaseModel):
        id: int
        dna_type: str
        sequence: str
        sequence_length: int
        one_hot_encoding: List
        is_reversed: bool
        is_generated: bool
        gene_name: Optional[str]

        class Config:
            arbitrary_types_allowed = True

    class SaveData(BaseModel):
        dna_class_types: List[str]
        sequence_padded_length: int
        universal_count: int
        sequences: List['SequenceSet.Sequence']

    def __init__(self, dna_class_types: List[str], sequence_padded_length: int) -> None:
        self.length = 0
        self.universal_count = 0
        self.sequences: List[SequenceSet.Sequence] = []
        self.dna_class_types: List[str] = dna_class_types
        self.sequence_padded_length = sequence_padded_length

    @classmethod
    def from_file(cls, file_path: str) -> 'SequenceSet':
        """Creates a SequenceSet instance from a JSON file."""
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            data = SequenceSet.SaveData(**json.load(f))
        
        instance = cls(
            dna_class_types=data.dna_class_types,
            sequence_padded_length=data.sequence_padded_length,
        )

        instance.universal_count = data.universal_count
        instance.sequences = data.sequences
        instance.length = len(instance.sequences)
        return instance


    def add_sequence(
        self,
        input_dna_type: str,
        sequence: str,
        gene_name: str = None,
        is_reversed: bool = False,
        is_generated: bool = False,
    )-> None:
        
        dna_type = input_dna_type
        if (isinstance(dna_type, Enum)):
            dna_type = dna_type.value
        if dna_type not in self.dna_class_types:
            raise ValueError(f"Invalid DNA type. Allowed types: {self.dna_class_types}")

        sequence = self.Sequence(
            id=self.universal_count,
            dna_type=dna_type,
            sequence=sequence,
            sequence_length=len(sequence),
            one_hot_encoding=bio_methods.one_hot_encode(sequence, self.sequence_padded_length),
            is_generated=is_generated,
            is_reversed=is_reversed,
            gene_name=gene_name,
        )
        self.sequences.append(sequence)
        self.length += 1
        self.universal_count += 1

    def remove_by_id(self, id: int) -> bool:
        """Removes the sequence with the specified ID from the sequences list."""
        for i, seq in enumerate(self.sequences):
            if seq.id == id:
                del self.sequences[i]
                self.length -= 1
                return True
        
        raise Exception(f"id {id} not found")

    def save(self, file_path: str) -> None:
        """Saves the sequences to a JSON file at the specified file path."""
        data = SequenceSet.SaveData(
            dna_class_types=self.dna_class_types,
            sequence_padded_length=self.sequence_padded_length,
            sequences=self.sequences,
            universal_count=self.universal_count
        )
        with gzip.open(file_path, 'wt', encoding='utf-8') as f:
            json.dump(data.model_dump(exclude_none=True), f, indent=4)

    def get_enum(self) -> List[str]:
        return self.dna_class_types
    
    def get_padded_sequence_length(self) -> int:
        return self.sequence_padded_length
    
    def save_train_data(self, file_path):
        # Create a mapping from DNA class types to integer labels
        dna_class_mapping = {dna_class: idx for idx, dna_class in enumerate(self.dna_class_types)}

        # Initialize arrays for storing the one-hot encoded sequences and labels
        train_data = np.array([seq.one_hot_encoding for seq in self.sequences], dtype=np.float32)
        labels = np.array([dna_class_mapping[seq.dna_type] for seq in self.sequences], dtype=np.int32)
        
        assert train_data.shape == (self.length, 4, self.sequence_padded_length)

        with h5py.File(file_path, 'w') as h5f:
            h5f.create_dataset('train_data', data=train_data)
            h5f.create_dataset('labels', data=labels)

        return train_data, labels

