import unittest
import numpy as np
import os
import json
from sequence_set import SequenceSet
import helper_functions as bio_methods
import h5py

class TestSequenceSet(unittest.TestCase):

    def setUp(self):
        self.dna_class_types = ["promoter", "enhancer"]
        self.sequence_padded_length = 10
        self.seq_set = SequenceSet(self.dna_class_types, self.sequence_padded_length)

    def test_init(self):
        self.assertEqual(self.seq_set.length, 0)
        self.assertEqual(self.seq_set.dna_class_types, self.dna_class_types)
        self.assertEqual(self.seq_set.sequence_padded_length, self.sequence_padded_length)
        self.assertEqual(len(self.seq_set.sequences), 0)

    def test_add_sequence_valid(self):
        seq = "ACGT"
        self.seq_set.add_sequence(
            id=1,
            input_dna_type="promoter",
            sequence=seq
        )
        self.assertEqual(self.seq_set.length, 1)
        self.assertEqual(self.seq_set.sequences[0].sequence, seq)
        self.assertEqual(self.seq_set.sequences[0].dna_type, "promoter")
        self.assertEqual(self.seq_set.sequences[0].sequence_length, len(seq))
        # Check one-hot encoding
        expected_one_hot = bio_methods.one_hot_encode(seq, self.sequence_padded_length)
        np.testing.assert_array_equal(self.seq_set.sequences[0].one_hot_encoding, expected_one_hot)

    def test_add_sequence_invalid_dna_type(self):
        with self.assertRaises(ValueError):
            self.seq_set.add_sequence(
                id=2,
                input_dna_type="invalid_type",
                sequence="ACGT"
            )

    def test_remove_sequence_by_id_valid(self):
        self.seq_set.add_sequence(
            id=1,
            input_dna_type="promoter",
            sequence="ACGT"
        )
        result = self.seq_set.remove_by_id(1)
        self.assertTrue(result)
        self.assertEqual(self.seq_set.length, 0)
        self.assertEqual(len(self.seq_set.sequences), 0)

    def test_remove_sequence_by_id_invalid(self):
        with self.assertRaises(Exception):
            self.seq_set.remove_by_id(999)

    def test_save_and_load(self):
        # Add sequences and save them
        seq1 = "ACGT"
        seq2 = "TGCA"
        self.seq_set.add_sequence(
            id=1,
            input_dna_type="promoter",
            sequence=seq1
        )
        self.seq_set.add_sequence(
            id=2,
            input_dna_type="enhancer",
            sequence=seq2
        )

        file_path = "test_sequence_set.json"
        self.seq_set.save(file_path)

        # Load the SequenceSet from the file
        loaded_seq_set = SequenceSet.from_file(file_path)

        self.assertEqual(loaded_seq_set.length, 2)
        self.assertEqual(loaded_seq_set.dna_class_types, self.dna_class_types)
        self.assertEqual(loaded_seq_set.sequence_padded_length, self.sequence_padded_length)

        # Check sequences in the loaded SequenceSet
        self.assertEqual(loaded_seq_set.sequences[0].sequence, seq1)
        self.assertEqual(loaded_seq_set.sequences[1].sequence, seq2)

        # Clean up
        os.remove(file_path)

    def test_get_enum(self):
        result = self.seq_set.get_enum()
        self.assertEqual(result, self.dna_class_types)

    def test_get_padded_sequence_length(self):
        result = self.seq_set.get_padded_sequence_length()
        self.assertEqual(result, self.sequence_padded_length)

    def test_one_hot_encoding(self):
        seq = "ACGT"
        padded_length = 6
        expected_encoding = [
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0]
        ]
        actual_encoding = bio_methods.one_hot_encode(seq, padded_length)
        np.testing.assert_array_equal(actual_encoding, expected_encoding)

    def test_save_train_data(self):
        # Setup the SequenceSet
        dna_class_types = ['type1', 'type2', 'type3']
        sequence_padded_length = 10
        sequence_set = SequenceSet(dna_class_types=dna_class_types, sequence_padded_length=sequence_padded_length)
        
        # Add sequences
        sequence_set.add_sequence(id=1, input_dna_type='type1', sequence="ACGTACGTAC")
        sequence_set.add_sequence(id=2, input_dna_type='type2', sequence="TGCATGCATG")
        
        # Define the file path for saving
        file_path = "test_train_data.h5"

        # Call save_train_data to save the sequences and labels
        sequence_set.save_train_data(file_path)

        # Load the saved data and verify
        with h5py.File(file_path, 'r') as h5f:
            saved_train_data = h5f['train_data'][:]
            saved_labels = h5f['labels'][:]
            
        # Get expected data
        expected_train_data, expected_labels = sequence_set.save_train_data(file_path)
        
        # Check if the saved data matches the expected data
        assert np.array_equal(saved_train_data, expected_train_data), "Train data mismatch"
        assert np.array_equal(saved_labels, expected_labels), "Labels mismatch"
        
        # Cleanup: remove the file after the test
        if os.path.exists(file_path):
            os.remove(file_path)

        print("Test passed successfully!")

if __name__ == "__main__":
    unittest.main()
