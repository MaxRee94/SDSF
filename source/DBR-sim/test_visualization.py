"""Unit test suite for visualization module.

This module provides comprehensive tests for the visualization functionality.
Tests cover the Visualiser class, image processing, Perlin noise generation, and visualization utilities.
"""

import unittest
import os
import sys
import tempfile
import shutil
import numpy as np
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call
import inspect

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock complex dependencies
sys.modules['cv2'] = MagicMock()
sys.modules['matplotlib'] = MagicMock()
sys.modules['matplotlib.pyplot'] = MagicMock()
sys.modules['colorsys'] = MagicMock()
sys.modules['platform'] = MagicMock()
sys.modules['subprocess'] = MagicMock()
sys.modules['ctypes'] = MagicMock()
sys.modules['re'] = MagicMock()

import visualization as vis_module


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.DATA_IN_DIR = tempfile.mkdtemp()
        self.DATA_OUT_DIR = tempfile.mkdtemp()
        self.LEGEND_PATH = tempfile.mkdtemp()
        self.rng = np.random.default_rng(42)
        self.grid_width = 100


class TestVisualiserClass(unittest.TestCase):
    """Test cases for Visualiser class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = MockConfig()

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        shutil.rmtree(self.cfg.DATA_IN_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.DATA_OUT_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.LEGEND_PATH, ignore_errors=True)

    def test_visualiser_class_exists(self):
        """Test that Visualiser class exists."""
        self.assertTrue(hasattr(vis_module, 'Visualiser'))
        self.assertTrue(callable(vis_module.Visualiser))

    def test_visualiser_initialization(self):
        """Test Visualiser initialization."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            visualiser = vis_module.Visualiser(self.cfg)
            
            self.assertEqual(visualiser.cfg, self.cfg)
            self.assertTrue(hasattr(visualiser, 'dpg'))

    def test_visualiser_init_perlin_noise_global_variables(self):
        """Test that Visualiser initializes Perlin noise variables."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            visualiser = vis_module.Visualiser(self.cfg)
            
            # Should have Perlin noise variables
            self.assertTrue(hasattr(visualiser, 'perm'))
            self.assertTrue(hasattr(visualiser, 'dirs'))


class TestVisualiserImageMethods(unittest.TestCase):
    """Test cases for Visualiser image-related methods."""

    def test_able_to_read_and_resize_image_function_exists(self):
        """Test that able_to_read_and_resize_image function exists."""
        self.assertTrue(hasattr(vis_module.Visualiser, 'able_to_read_and_resize_image'))
        self.assertTrue(callable(vis_module.Visualiser.able_to_read_and_resize_image))

    def test_generate_disk_pattern_function_exists(self):
        """Test that generate_disk_pattern function exists."""
        self.assertTrue(hasattr(vis_module.Visualiser, 'generate_disk_pattern'))
        self.assertTrue(callable(vis_module.Visualiser.generate_disk_pattern))

    def test_generate_perlin_noise_image_function_exists(self):
        """Test that generate_perlin_noise_image function exists."""
        self.assertTrue(hasattr(vis_module.Visualiser, 'generate_perlin_noise_image'))
        self.assertTrue(callable(vis_module.Visualiser.generate_perlin_noise_image))

    def test_get_thresholded_image_function_exists(self):
        """Test that get_thresholded_image function exists."""
        self.assertTrue(hasattr(vis_module.Visualiser, 'get_thresholded_image'))
        self.assertTrue(callable(vis_module.Visualiser.get_thresholded_image))


class TestVisualiserUtilityMethods(unittest.TestCase):
    """Test cases for Visualiser utility methods."""

    def test_visualise_initial_state_function_exists(self):
        """Test that visualise_initial_state function exists."""
        self.assertTrue(hasattr(vis_module, 'visualise_initial_state'))
        self.assertTrue(callable(vis_module.visualise_initial_state))

    def test_create_color_dict_function_exists(self):
        """Test that create_color_dict function exists."""
        self.assertTrue(hasattr(vis_module, 'create_color_dict'))
        self.assertTrue(callable(vis_module.create_color_dict))

    def test_do_visualizations_function_exists(self):
        """Test that do_visualizations function exists."""
        self.assertTrue(hasattr(vis_module, 'do_visualizations'))
        self.assertTrue(callable(vis_module.do_visualizations))


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = ['Visualiser']
        for class_name in expected_classes:
            self.assertTrue(hasattr(vis_module, class_name))
            self.assertTrue(callable(getattr(vis_module, class_name)))

    def test_expected_functions_present(self):
        """Test that all expected functions are present in the module."""
        expected_functions = [
            'visualise_initial_state', 'create_color_dict', 'do_visualizations'
        ]
        for func_name in expected_functions:
            self.assertTrue(hasattr(vis_module, func_name))
            self.assertTrue(callable(getattr(vis_module, func_name)))


class TestVisualiserMethods(unittest.TestCase):
    """Test cases for Visualiser methods by examining source code."""

    def test_visualiser_has_required_methods(self):
        """Test that Visualiser has required methods."""
        required_methods = [
            '__init__', 'init_perlin_noise_global_variables',
            'able_to_read_and_resize_image', 'generate_disk_pattern',
            'generate_perlin_noise_image', 'get_thresholded_image'
        ]
        for method in required_methods:
            self.assertTrue(hasattr(vis_module.Visualiser, method))
            self.assertTrue(callable(getattr(vis_module.Visualiser, method)))

    def test_visualiser_init_perlin_noise_variables(self):
        """Test that init_perlin_noise_global_variables initializes required variables."""
        source = inspect.getsource(vis_module.Visualiser.init_perlin_noise_global_variables)
        
        # Should initialize perm and dirs
        self.assertIn('self.perm', source)
        self.assertIn('self.dirs', source)


class TestModuleImports(unittest.TestCase):
    """Test cases for module imports."""

    def test_disk_pattern_generator_import(self):
        """Test that disk_pattern_generator is imported."""
        source = inspect.getsource(vis_module)
        self.assertIn('import disk_pattern_generator as dpg', source)

    def test_helpers_import(self):
        """Test that helpers is imported."""
        source = inspect.getsource(vis_module)
        self.assertIn('from helpers import *', source)

    def test_file_handling_import(self):
        """Test that file_handling is imported."""
        source = inspect.getsource(vis_module)
        self.assertIn('import file_handling as io', source)

    def test_config_import(self):
        """Test that config is imported."""
        source = inspect.getsource(vis_module)
        self.assertIn('from config import *', source)


class TestPerlinNoiseFunctionality(unittest.TestCase):
    """Test cases for Perlin noise functionality."""

    def test_perlin_noise_initialization(self):
        """Test Perlin noise initialization in Visualiser."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            cfg = MockConfig()
            visualiser = vis_module.Visualiser(cfg)
            
            # Check that Perlin noise variables are initialized
            self.assertTrue(hasattr(visualiser, 'perm'))
            self.assertTrue(hasattr(visualiser, 'dirs'))
            
            # Check that they are lists with expected length
            self.assertIsInstance(visualiser.perm, list)
            self.assertEqual(len(visualiser.perm), 8192)  # Doubled from 4096


class TestModuleConstants(unittest.TestCase):
    """Test cases for module constants and configurations."""

    def test_visualiser_has_cfg_attribute(self):
        """Test that Visualiser has cfg attribute."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            cfg = MockConfig()
            visualiser = vis_module.Visualiser(cfg)
            
            self.assertEqual(visualiser.cfg, cfg)

    def test_visualiser_has_dpg_attribute(self):
        """Test that Visualiser has dpg attribute."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            cfg = MockConfig()
            visualiser = vis_module.Visualiser(cfg)
            
            self.assertEqual(visualiser.dpg, mock_dpg_instance)


class TestVisualiserImageProcessing(unittest.TestCase):
    """Test cases for Visualiser image processing methods."""

    def test_able_to_read_and_resize_image_uses_cv2(self):
        """Test that able_to_read_and_resize_image uses cv2."""
        source = inspect.getsource(vis_module.Visualiser.able_to_read_and_resize_image)
        
        self.assertIn('cv2.imread', source)
        self.assertIn('cv2.resize', source)

    def test_get_thresholded_image_uses_cv2(self):
        """Test that get_thresholded_image uses cv2."""
        source = inspect.getsource(vis_module.Visualiser.get_thresholded_image)
        
        self.assertIn('cv2.threshold', source)


class TestVisualisationFunctions(unittest.TestCase):
    """Test cases for visualization functions."""

    def test_visualise_initial_state_function(self):
        """Test visualise_initial_state function."""
        # This function likely creates visualizations of initial state
        # We can't test the actual behavior without GUI dependencies,
        # but we can verify it exists and has the right signature
        sig = inspect.signature(vis_module.visualise_initial_state)
        params = list(sig.parameters.keys())
        
        self.assertIn('dynamics', params)
        self.assertIn('cfg', params)

    def test_create_color_dict_function(self):
        """Test create_color_dict function."""
        sig = inspect.signature(vis_module.create_color_dict)
        params = list(sig.parameters.keys())
        
        self.assertIn('cfg', params)

    def test_do_visualizations_function(self):
        """Test do_visualizations function."""
        sig = inspect.signature(vis_module.do_visualizations)
        params = list(sig.parameters.keys())
        
        # Should accept many parameters
        self.assertGreater(len(params), 5)


class TestEdgeCases(unittest.TestCase):
    """Test cases for edge cases and error handling."""

    def test_visualiser_with_none_cfg(self):
        """Test Visualiser with None cfg."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            visualiser = vis_module.Visualiser(None)
            self.assertIsNone(visualiser.cfg)

    def test_perlin_noise_with_different_seeds(self):
        """Test that different rng seeds produce different Perlin noise."""
        with patch('visualization.dpg.DiskPatternGenerator') as mock_dpg:
            mock_dpg_instance = MagicMock()
            mock_dpg.return_value = mock_dpg_instance
            
            # Create two visualisers with different RNGs
            cfg1 = MockConfig()
            cfg1.rng = np.random.default_rng(42)
            
            cfg2 = MockConfig()
            cfg2.rng = np.random.default_rng(123)
            
            visualiser1 = vis_module.Visualiser(cfg1)
            visualiser2 = vis_module.Visualiser(cfg2)
            
            # Their perm arrays should be different due to different shuffling
            # (This tests the randomness of the initialization)
            self.assertNotEqual(visualiser1.perm, visualiser2.perm)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)