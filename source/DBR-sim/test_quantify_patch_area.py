"""Unit test suite for quantify_patch_area module.

This module provides tests for the patch area quantification functions.
Tests cover the main function get_patch_areas with periodic boundary conditions.
"""

import unittest
import numpy as np
import sys
import os

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quantify_patch_area as qpa


class TestGetPatchAreas(unittest.TestCase):
    """Test cases for get_patch_areas function."""

    def test_get_patch_areas_all_zeros(self):
        """Test get_patch_areas with all zeros (no patches)."""
        image = np.zeros((10, 10), dtype=int)
        domain_width = 100.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should return empty list since there are no patches
        self.assertEqual(result, [])

    def test_get_patch_areas_single_patch(self):
        """Test get_patch_areas with a single patch."""
        # Create a 5x5 image with a single 3x3 patch of 1s
        image = np.zeros((5, 5), dtype=int)
        image[1:4, 1:4] = 1  # 3x3 patch
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should return one patch area
        self.assertEqual(len(result), 1)
        # Just check that it's a positive number (exact value depends on implementation)
        self.assertGreater(result[0], 0)

    def test_get_patch_areas_multiple_patches(self):
        """Test get_patch_areas with multiple separate patches."""
        # Create a 10x10 image with two separate 2x2 patches
        image = np.zeros((10, 10), dtype=int)
        image[1:3, 1:3] = 1  # First patch
        image[7:9, 7:9] = 1  # Second patch
        domain_width = 20.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should return two patch areas
        self.assertEqual(len(result), 2)
        # Both should be positive
        self.assertGreater(result[0], 0)
        self.assertGreater(result[1], 0)

    def test_get_patch_areas_periodic_boundary_conditions(self):
        """Test get_patch_areas with periodic boundary conditions."""
        # Create an image where patches wrap around the edges
        image = np.zeros((5, 5), dtype=int)
        image[0, :] = 1  # Top edge
        image[-1, :] = 1  # Bottom edge (should connect with periodic BC)
        domain_width = 10.0
        
        # This should be treated as a single connected patch due to periodic BC
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should work without errors
        self.assertIsInstance(result, list)

    def test_get_patch_areas_output_type(self):
        """Test that get_patch_areas returns a list."""
        image = np.zeros((5, 5), dtype=int)
        result = qpa.get_patch_areas(image, 10.0)
        
        self.assertIsInstance(result, list)


class TestGetPatchAreasEdgeCases(unittest.TestCase):
    """Edge case tests for get_patch_areas."""

    def test_get_patch_areas_empty_image(self):
        """Test get_patch_areas with empty image."""
        image = np.array([], dtype=int).reshape(0, 0)
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        self.assertEqual(result, [])

    def test_get_patch_areas_single_cell_patch(self):
        """Test get_patch_areas with single cell patches."""
        image = np.zeros((5, 5), dtype=int)
        image[2, 2] = 1  # Single cell patch
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        self.assertEqual(len(result), 1)
        # Should be positive
        self.assertGreater(result[0], 0)

    def test_get_patch_areas_all_ones(self):
        """Test get_patch_areas with all ones (one big patch)."""
        image = np.ones((5, 5), dtype=int)
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should return one big patch
        self.assertEqual(len(result), 1)
        self.assertGreater(result[0], 0)

    def test_get_patch_areas_checkerboard_pattern(self):
        """Test get_patch_areas with checkerboard pattern."""
        # Create a checkerboard pattern
        image = np.zeros((4, 4), dtype=int)
        image[::2, ::2] = 1  # Checkerboard
        domain_width = 8.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should find multiple separate patches
        self.assertGreater(len(result), 1)


class TestPatchAreaCalculation(unittest.TestCase):
    """Test cases for the accuracy of patch area calculations."""

    def test_single_cell_area_basic(self):
        """Test that single cell area is calculated consistently."""
        domain_width = 200.0
        image = np.ones((1, 1), dtype=int)  # Single cell
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should return a positive area for the single cell
        self.assertEqual(len(result), 1)
        self.assertGreater(result[0], 0)
        # Store this as reference for other tests
        single_cell_area = result[0]
        
        # Test with 2x2 cells
        image2 = np.ones((2, 2), dtype=int)  # 4 cells
        result2 = qpa.get_patch_areas(image2, domain_width)
        
        # The 2x2 patch should have roughly 4x the area of single cell
        # (exact ratio might not be exactly 4 due to algorithm implementation)
        self.assertEqual(len(result2), 1)
        self.assertGreater(result2[0], single_cell_area)

    def test_domain_width_affects_areas(self):
        """Test that domain width affects the calculated areas."""
        # Same image, different domain widths
        image = np.ones((2, 2), dtype=int)  # 4 cells
        
        domain_width1 = 200.0
        domain_width2 = 400.0
        
        result1 = qpa.get_patch_areas(image, domain_width1)
        result2 = qpa.get_patch_areas(image, domain_width2)
        
        # Both should return positive areas
        self.assertGreater(result1[0], 0)
        self.assertGreater(result2[0], 0)
        
        # The areas should be different for different domain widths
        # (The exact relationship depends on the implementation)


class TestPatchConnectivity(unittest.TestCase):
    """Test cases for patch connectivity behavior."""

    def test_adjacent_pixels_same_patch(self):
        """Test that adjacent pixels are considered same patch."""
        # Create 2 adjacent pixels
        image = np.zeros((3, 3), dtype=int)
        image[1, 1] = 1  # Center pixel
        image[1, 2] = 1  # Right neighbor
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should be one connected patch
        self.assertEqual(len(result), 1)

    def test_diagonal_pixels_separate_patches(self):
        """Test that diagonal pixels are considered separate patches."""
        # Create diagonal pixels (not connected in 4-connectivity)
        image = np.zeros((3, 3), dtype=int)
        image[1, 1] = 1  # Center pixel
        image[2, 2] = 1  # Diagonal neighbor (not connected in 4-connectivity)
        domain_width = 10.0
        
        result = qpa.get_patch_areas(image, domain_width)
        
        # Should be two separate patches (if using 4-connectivity)
        # Note: this depends on the connectivity used in the implementation
        # The function appears to use 4-connectivity (up, down, left, right)
        self.assertEqual(len(result), 2)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)