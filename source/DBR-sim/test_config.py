"""Unit test suite for config module.

This module provides tests for the configuration functions and constants in config.py.
Tests cover the directory derivation and configuration loading functions.
"""

import unittest
import os
import sys
import tempfile
from types import SimpleNamespace
import shutil

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg


class TestConfigConstants(unittest.TestCase):
    """Test cases for configuration constants."""

    def test_constants_defined(self):
        """Test that basic constants are defined."""
        self.assertTrue(hasattr(cfg, 'constants'))
        self.assertTrue(hasattr(cfg, 'cfg'))

    def test_repository_basedir(self):
        """Test that REPOSITORY_BASEDIR is defined."""
        self.assertTrue('REPOSITORY_BASEDIR' in cfg.constants)
        self.assertIsInstance(cfg.constants['REPOSITORY_BASEDIR'], str)

    def test_data_directories_defined(self):
        """Test that data directories are defined."""
        expected_dirs = ['DATA_IN_DIR', 'DATA_OUT_DIR', 'DATA_INTERNAL_DIR', 'BUILD_DIR', 'SOURCE_DIR']
        for dir_name in expected_dirs:
            self.assertTrue(dir_name in cfg.constants, f"{dir_name} not found in constants")


class TestDeriveDirectories(unittest.TestCase):
    """Test cases for directory derivation functions."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a test configuration
        self.test_cfg = SimpleNamespace(
            DATA_OUT_DIR="/tmp/test_out",
            DATA_IN_DIR="/tmp/test_in"
        )

    def test_derive_output_dirs_adds_expected_paths(self):
        """Test that derive_output_dirs adds expected output directory paths."""
        cfg_copy = SimpleNamespace(
            DATA_OUT_DIR="/tmp/test_out"
        )
        
        result = cfg.derive_output_dirs(cfg_copy)
        
        expected_paths = [
            'CPG_OUTPUT_DIR',
            'LEGEND_PATH', 
            'TREE_DBH_FILE',
            'EXPORT_DIR',
            'END2END_TEST_OUTPUT_DIR',
            'END2END_TEST_BENCHMARK_DIR'
        ]
        
        for path_name in expected_paths:
            self.assertTrue(hasattr(result, path_name), f"{path_name} not found in result")

    def test_derive_output_dirs_creates_correct_paths(self):
        """Test that derive_output_dirs creates correct path strings."""
        cfg_copy = SimpleNamespace(
            DATA_OUT_DIR="/tmp/test_out"
        )
        
        result = cfg.derive_output_dirs(cfg_copy)
        
        # Use path normalization to handle different OS path separators
        import os
        self.assertEqual(os.path.normpath(result.CPG_OUTPUT_DIR), os.path.normpath("/tmp/test_out/controlled_pattern_generator"))
        self.assertEqual(os.path.normpath(result.LEGEND_PATH), os.path.normpath("/tmp/test_out/legends"))
        self.assertEqual(os.path.normpath(result.TREE_DBH_FILE), os.path.normpath("/tmp/test_out/state_reports/tree_dbh_values.json"))
        self.assertEqual(os.path.normpath(result.EXPORT_DIR), os.path.normpath("/tmp/test_out/state_data"))

    def test_derive_input_dirs_adds_expected_paths(self):
        """Test that derive_input_dirs adds expected input directory paths."""
        cfg_copy = SimpleNamespace(
            DATA_IN_DIR="/tmp/test_in"
        )
        
        result = cfg.derive_input_dirs(cfg_copy)
        
        expected_paths = [
            'PERLIN_NOISE_DIR',
            'SIMPLE_PATTERNS_DIR',
            'CONTROLLED_PATTERN_DIR',
            'END2END_TESTCASE_DIR'
        ]
        
        for path_name in expected_paths:
            self.assertTrue(hasattr(result, path_name), f"{path_name} not found in result")

    def test_derive_input_dirs_creates_correct_paths(self):
        """Test that derive_input_dirs creates correct path strings."""
        cfg_copy = SimpleNamespace(
            DATA_IN_DIR="/tmp/test_in"
        )
        
        result = cfg.derive_input_dirs(cfg_copy)
        
        self.assertTrue(result.PERLIN_NOISE_DIR.endswith("state_patterns/perlin_noise"))
        self.assertTrue(result.SIMPLE_PATTERNS_DIR.endswith("state_patterns/simple_patterns"))
        self.assertTrue(result.CONTROLLED_PATTERN_DIR.endswith("state_patterns/controlled_patterns"))


class TestApplyLocalOverrides(unittest.TestCase):
    """Test cases for apply_local_overrides function."""

    def test_apply_local_overrides_no_overrides_file(self):
        """Test apply_local_overrides when no overrides file exists."""
        # Create a temporary config
        test_cfg = SimpleNamespace(DATA_IN_DIR="/nonexistent/path")
        
        result = cfg.apply_local_overrides(test_cfg)
        
        # Should return the original config since no overrides file exists
        self.assertEqual(result, test_cfg)

    def test_apply_local_overrides_with_missing_directory(self):
        """Test apply_local_overrides when the local_overrides directory doesn't exist."""
        # Create a temporary directory structure
        with tempfile.TemporaryDirectory() as temp_dir:
            test_cfg = SimpleNamespace(DATA_IN_DIR=temp_dir)
            
            result = cfg.apply_local_overrides(test_cfg)
            
            # Should return the original config since no overrides file exists
            self.assertEqual(result, test_cfg)


class TestConfigIntegration(unittest.TestCase):
    """Test cases for the overall configuration integration."""

    def test_cfg_has_required_attributes(self):
        """Test that the global cfg has required attributes."""
        required_attrs = [
            'REPOSITORY_BASEDIR', 'DATA_IN_DIR', 'DATA_OUT_DIR',
            'DATA_INTERNAL_DIR', 'BUILD_DIR', 'SOURCE_DIR'
        ]
        
        for attr in required_attrs:
            self.assertTrue(hasattr(cfg.cfg, attr), f"cfg missing attribute: {attr}")

    def test_cfg_has_derived_attributes(self):
        """Test that the global cfg has derived output and input directory attributes."""
        # These should have been added by the derivation functions
        expected_output_attrs = [
            'CPG_OUTPUT_DIR', 'LEGEND_PATH', 'TREE_DBH_FILE',
            'EXPORT_DIR', 'END2END_TEST_OUTPUT_DIR', 'END2END_TEST_BENCHMARK_DIR'
        ]
        
        expected_input_attrs = [
            'PERLIN_NOISE_DIR', 'SIMPLE_PATTERNS_DIR',
            'CONTROLLED_PATTERN_DIR', 'END2END_TESTCASE_DIR'
        ]
        
        for attr in expected_output_attrs:
            self.assertTrue(hasattr(cfg.cfg, attr), f"cfg missing output attribute: {attr}")
            
        for attr in expected_input_attrs:
            self.assertTrue(hasattr(cfg.cfg, attr), f"cfg missing input attribute: {attr}")


class TestConfigDefaults(unittest.TestCase):
    """Test cases for configuration defaults (if any exist)."""

    def test_defaults_exist(self):
        """Test that defaults dictionary exists."""
        self.assertTrue(hasattr(cfg, 'defaults'))
        self.assertIsInstance(cfg.defaults, dict)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)