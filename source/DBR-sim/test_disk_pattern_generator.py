"""Unit test suite for disk_pattern_generator module.

This module provides tests for the DiskPatternGenerator class.
Tests cover the main methods like compute_area and fraction_white_pixels.
"""

import unittest
import numpy as np
from types import SimpleNamespace
import math
import cv2
import sys
import os

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disk_pattern_generator import DiskPatternGenerator


class TestDiskPatternGeneratorInitialization(unittest.TestCase):
    """Test cases for DiskPatternGenerator initialization."""

    def test_initialization(self):
        """Test that DiskPatternGenerator initializes correctly."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        dpg = DiskPatternGenerator(cfg)
        
        self.assertEqual(dpg.cfg, cfg)
        self.assertTrue(hasattr(dpg, 'local_rng_seed_multiplier'))
        self.assertIsInstance(dpg.local_rng_seed_multiplier, (int, np.integer))


class TestComputeArea(unittest.TestCase):
    """Test cases for compute_area method."""

    def setUp(self):
        """Set up test fixtures."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        self.dpg = DiskPatternGenerator(cfg)

    def test_compute_area_too_few_points(self):
        """Test compute_area with too few points."""
        contour = np.array([[0, 0], [1, 1]])  # Only 2 points
        center = (0.5, 0.5)
        
        with self.assertRaises(ValueError):
            self.dpg.compute_area(contour, center)

    def test_compute_area_triangle(self):
        """Test compute_area with a triangle contour."""
        # Create a simple triangle contour
        contour = np.array([[0, 0], [2, 0], [1, 2]], dtype=np.float32)
        center = (1, 1)  # Center of the triangle
        
        area = self.dpg.compute_area(contour, center)
        
        # Area of a triangle with base 2 and height 2 is 2
        expected_area = 2.0
        # Allow some tolerance due to calculation method
        self.assertAlmostEqual(area, expected_area, delta=0.1)

    def test_compute_area_square(self):
        """Test compute_area with a square contour."""
        # Create a square contour (4 points)
        contour = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=np.float32)
        center = (0.5, 0.5)  # Center of the square
        
        area = self.dpg.compute_area(contour, center)
        
        # Area of a 1x1 square is 1
        expected_area = 1.0
        # Allow some tolerance due to calculation method
        self.assertAlmostEqual(area, expected_area, delta=0.1)


class TestFractionWhitePixels(unittest.TestCase):
    """Test cases for fraction_white_pixels method."""

    def setUp(self):
        """Set up test fixtures."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        self.dpg = DiskPatternGenerator(cfg)

    def test_fraction_white_pixels_all_white(self):
        """Test fraction_white_pixels with all white image."""
        img = np.full((10, 10), 255, dtype=np.uint8)
        result = self.dpg.fraction_white_pixels(img)
        self.assertEqual(result, 1.0)

    def test_fraction_white_pixels_all_black(self):
        """Test fraction_white_pixels with all black image."""
        img = np.full((10, 10), 0, dtype=np.uint8)
        result = self.dpg.fraction_white_pixels(img)
        self.assertEqual(result, 0.0)

    def test_fraction_white_pixels_half_white(self):
        """Test fraction_white_pixels with half white image."""
        img = np.zeros((10, 10), dtype=np.uint8)
        img[:5, :] = 255  # Make top half white
        result = self.dpg.fraction_white_pixels(img)
        self.assertEqual(result, 0.5)

    def test_fraction_white_pixels_mixed(self):
        """Test fraction_white_pixels with mixed image."""
        img = np.array([[255, 0], [0, 255]], dtype=np.uint8)
        result = self.dpg.fraction_white_pixels(img)
        self.assertEqual(result, 0.5)  # 2 out of 4 pixels are white


class TestComputeAreaNormalizationFactor(unittest.TestCase):
    """Test cases for compute_area_normalization_factor method."""

    def setUp(self):
        """Set up test fixtures."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        self.dpg = DiskPatternGenerator(cfg)

    def test_compute_area_normalization_factor_with_global_factor(self):
        """Test compute_area_normalization_factor with global factor."""
        # Create a simple contour
        contour = np.array([[0, 0], [2, 0], [1, 2]], dtype=np.float32)
        center = (1, 1)
        base_radius = 1.0
        global_factor = 1.0
        
        result = self.dpg.compute_area_normalization_factor(
            contour, center, base_radius, global_factor
        )
        
        # Should return some normalization factor
        self.assertIsInstance(result, float)
        self.assertGreater(result, 0)

    def test_compute_area_normalization_factor_without_global_factor(self):
        """Test compute_area_normalization_factor without global factor."""
        # Create a simple contour
        contour = np.array([[0, 0], [2, 0], [1, 2]], dtype=np.float32)
        center = (1, 1)
        base_radius = 1.0
        global_factor = None
        
        result = self.dpg.compute_area_normalization_factor(
            contour, center, base_radius, global_factor
        )
        
        # Should return some normalization factor
        self.assertIsInstance(result, float)
        self.assertGreater(result, 0)


class TestDiskPatternGeneratorEdgeCases(unittest.TestCase):
    """Edge case tests for DiskPatternGenerator."""

    def test_fraction_white_pixels_different_dtypes(self):
        """Test fraction_white_pixels with different data types."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        dpg = DiskPatternGenerator(cfg)
        
        # Test with float image
        img_float = np.full((5, 5), 255.0, dtype=np.float32)
        result = dpg.fraction_white_pixels(img_float)
        self.assertEqual(result, 1.0)

    def test_fraction_white_pixels_different_white_values(self):
        """Test fraction_white_pixels with different white value thresholds."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        dpg = DiskPatternGenerator(cfg)
        
        # Test with different white values (though function expects exactly 255)
        img = np.full((5, 5), 254, dtype=np.uint8)  # Close to white but not exactly
        result = dpg.fraction_white_pixels(img)
        self.assertEqual(result, 0.0)  # Should be 0 since it's not exactly 255

    def test_compute_area_different_contour_types(self):
        """Test compute_area with different contour types."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        dpg = DiskPatternGenerator(cfg)
        
        # Test with different contour point counts
        for n in [3, 4, 5, 10]:
            theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
            contour = np.column_stack([np.cos(theta), np.sin(theta)]).astype(np.float32)
            center = (0, 0)
            
            area = dpg.compute_area(contour, center)
            self.assertIsInstance(area, float)
            self.assertGreater(area, 0)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)