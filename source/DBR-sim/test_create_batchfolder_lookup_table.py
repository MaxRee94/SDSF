"""Unit test suite for create_batchfolder_lookup_table module.

This module provides tests for the batch folder lookup table functions.
"""

import unittest
import os
import sys
import tempfile
import shutil
from unittest.mock import patch

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# We'll need to mock the config module to avoid file system dependencies
with patch('create_batchfolder_lookup_table.cfg') as mock_cfg:
    mock_cfg.DATA_OUT_DIR = "/tmp/test_data_out"
    import create_batchfolder_lookup_table as cblt


class TestGetXNewBatchFolders(unittest.TestCase):
    """Test cases for get_x_new_batch_folders function."""

    def test_get_x_new_batch_folders_single(self):
        """Test get_x_new_batch_folders with x=1."""
        with patch.object(cblt.os, 'path') as mock_path:
            # Mock os.path.exists to return False (folders don't exist)
            mock_path.exists.return_value = False
            
            result = cblt.get_x_new_batch_folders(1)
            
            self.assertEqual(len(result), 1)
            self.assertIsInstance(result[0], str)
            self.assertIn("batch_", result[0])

    def test_get_x_new_batch_folders_multiple(self):
        """Test get_x_new_batch_folders with x=5."""
        with patch.object(cblt.os, 'path') as mock_path:
            # Mock os.path.exists to return False (folders don't exist)
            mock_path.exists.return_value = False
            
            result = cblt.get_x_new_batch_folders(5)
            
            self.assertEqual(len(result), 5)
            # All should be strings with batch_ prefix
            for folder in result:
                self.assertIsInstance(folder, str)
                self.assertIn("batch_", folder)

    def test_get_x_new_batch_folders_with_existing(self):
        """Test get_x_new_batch_folders when some folders already exist."""
        with patch.object(cblt.os, 'path') as mock_path:
            # Mock os.path.exists to return True for batch_000001, False for others
            def exists_mock(path):
                return "batch_000001" in path
            
            mock_path.exists.side_effect = exists_mock
            
            result = cblt.get_x_new_batch_folders(2)
            
            self.assertEqual(len(result), 2)
            # Should start from batch_000002
            for folder in result:
                self.assertIn("batch_", folder)

    def test_get_x_new_batch_folders_zero(self):
        """Test get_x_new_batch_folders with x=0."""
        with patch.object(cblt.os, 'path') as mock_path:
            mock_path.exists.return_value = False
            
            result = cblt.get_x_new_batch_folders(0)
            
            self.assertEqual(result, [])


class TestCreateBatchLookupTable(unittest.TestCase):
    """Test cases for create_batch_lookup_table function."""

    def test_create_batch_lookup_table_default(self):
        """Test create_batch_lookup_table with default parameters."""
        with patch.object(cblt, 'get_x_new_batch_folders') as mock_get_folders:
            mock_get_folders.return_value = ["/path/to/batch_000001", "/path/to/batch_000002"]
            
            result = cblt.create_batch_lookup_table()
            
            self.assertIsInstance(result, dict)

    def test_create_batch_lookup_table_with_number(self):
        """Test create_batch_lookup_table with specific number of folders."""
        with patch.object(cblt, 'get_x_new_batch_folders') as mock_get_folders:
            mock_get_folders.return_value = [
                "/path/to/batch_000001", 
                "/path/to/batch_000002",
                "/path/to/batch_000003"
            ]
            
            result = cblt.create_batch_lookup_table(number_of_batchfolders=3)
            
            self.assertIsInstance(result, dict)

    def test_lookup_table_structure(self):
        """Test that the lookup table has the expected structure."""
        with patch.object(cblt, 'get_x_new_batch_folders') as mock_get_folders:
            mock_get_folders.return_value = ["/path/to/batch_000001", "/path/to/batch_000002"]
            
            result = cblt.create_batch_lookup_table(number_of_batchfolders=2)
            
            # Should have keys 0 and 1
            self.assertIn(0, result)
            self.assertIn(1, result)
            
            # Values should be the folder paths
            self.assertEqual(result[0], "/path/to/batch_000001")
            self.assertEqual(result[1], "/path/to/batch_000002")


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)