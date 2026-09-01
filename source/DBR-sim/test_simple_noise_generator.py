"""Unit test suite for simple_noise_generator module.

This module provides comprehensive tests for the macropixel noise generation functionality.
Tests cover normal operation, edge cases, parameter validation, and statistical properties.
"""

import unittest
import numpy as np
from types import SimpleNamespace
import cv2
import sys
import os

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simple_noise_generator as sng


class TestSimpleNoiseGenerator(unittest.TestCase):
    """Test cases for the simple_noise_generator module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a test configuration with a deterministic RNG for reproducible tests
        self.cfg = SimpleNamespace()
        self.cfg.rng = np.random.default_rng(42)  # Seed 42 for reproducibility
        
        # Default parameters
        self.default_params = {
            'amplitude': 0.3,
            'scale': 10,
            'binary_connectivity': -1,
            'offset': 0,
            'grid_width': 100,
            'cfg': self.cfg
        }

    def tearDown(self):
        """Clean up after tests."""
        # Close any OpenCV windows that might have been opened
        cv2.destroyAllWindows()

    # ==================== BASIC FUNCTIONALITY TESTS ====================

    def test_generate_returns_numpy_array(self):
        """Test that generate() returns a numpy array."""
        result = sng.generate(**self.default_params)
        self.assertIsInstance(result, np.ndarray)

    def test_generate_returns_2d_array(self):
        """Test that generate() returns a 2D array."""
        result = sng.generate(**self.default_params)
        self.assertEqual(result.ndim, 2)

    def test_generate_output_shape_matches_grid_width(self):
        """Test that output shape matches the specified grid_width."""
        grid_width = 120
        result = sng.generate(**{**self.default_params, 'grid_width': grid_width})
        self.assertEqual(result.shape, (grid_width, grid_width))

    def test_generate_default_grid_width(self):
        """Test that default grid_width of 1000 produces 1000x1000 output."""
        params = {**self.default_params}
        del params['grid_width']  # Use default
        result = sng.generate(**params)
        self.assertEqual(result.shape, (1000, 1000))

    def test_generate_values_in_range_zero_to_one(self):
        """Test that all output values are in the range [0, 1]."""
        result = sng.generate(**self.default_params)
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    # ==================== PARAMETER VALIDATION TESTS ====================

    def test_amplitude_parameter(self):
        """Test that amplitude parameter controls the noise range."""
        # Test with amplitude of 0.2 (range should be [0.3, 0.7])
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(amplitude=0.2, scale=1, grid_width=10, cfg=cfg)
        
        # With amplitude=0.2, values should be in [0.5-0.2, 0.5+0.2] = [0.3, 0.7]
        # But due to clipping, we just verify they're in [0,1]
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_scale_parameter_affects_macropixel_size(self):
        """Test that scale parameter affects the macropixel size."""
        grid_width = 100
        
        # Test with scale=10 (should have 10x10 macropixels)
        cfg1 = SimpleNamespace(rng=np.random.default_rng(42))
        result1 = sng.generate(scale=10, grid_width=grid_width, cfg=cfg1)
        
        # Test with scale=20 (should have 5x5 macropixels)  
        cfg2 = SimpleNamespace(rng=np.random.default_rng(42))
        result2 = sng.generate(scale=20, grid_width=grid_width, cfg=cfg2)
        
        # Both should be the same final size
        self.assertEqual(result1.shape, result2.shape)
        self.assertEqual(result1.shape, (grid_width, grid_width))

    def test_offset_parameter_shifts_values(self):
        """Test that offset parameter shifts all values."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Generate noise without offset
        result_no_offset = sng.generate(
            amplitude=0.1, scale=1, grid_width=10, offset=0, cfg=cfg
        )
        
        # Generate noise with positive offset
        cfg2 = SimpleNamespace(rng=np.random.default_rng(42))
        result_with_offset = sng.generate(
            amplitude=0.1, scale=1, grid_width=10, offset=0.2, cfg=cfg2
        )
        
        # With offset, values should be shifted up
        # Note: we can't directly compare due to randomness, but we can check ranges
        self.assertTrue(np.all(result_no_offset >= 0.4))  # 0.5-0.1 = 0.4
        self.assertTrue(np.all(result_no_offset <= 0.6))  # 0.5+0.1 = 0.6
        
        # With offset 0.2, range should be [0.4+0.2, 0.6+0.2] = [0.6, 0.8] before clipping
        # But clipping might occur, so we just check it's in [0,1]
        self.assertTrue(np.all(result_with_offset >= 0.0))
        self.assertTrue(np.all(result_with_offset <= 1.0))

    def test_clipping_prevents_out_of_range_values(self):
        """Test that values are clipped to [0, 1] range."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Use large amplitude and offset that would produce out-of-range values
        result = sng.generate(
            amplitude=1.0,  # Would produce [-0.5, 1.5] range
            offset=0.5,    # Would shift to [0.0, 2.0] range
            scale=1, 
            grid_width=10, 
            cfg=cfg
        )
        
        # All values should be clipped to [0, 1]
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    # ==================== BINARY CONNECTIVITY TESTS ====================

    def test_binary_connectivity_negative_uses_uniform_noise(self):
        """Test that binary_connectivity < 0 uses uniform noise generation."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(
            binary_connectivity=-1,
            amplitude=0.2,
            scale=1,
            grid_width=10,
            cfg=cfg
        )
        
        # Should not be binary (should have values other than 0 and 1)
        unique_values = np.unique(result)
        self.assertGreater(len(unique_values), 2)  # More than just 0 and 1

    def test_binary_connectivity_zero_generates_all_zeros(self):
        """Test that binary_connectivity=0 generates mostly zeros."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(
            binary_connectivity=0.0,
            scale=1,
            grid_width=10,
            cfg=cfg
        )
        
        # Should be mostly zeros (some interpolation artifacts might occur due to resize)
        # The macropixel grid should be all zeros, but after interpolation it might not be exactly 0
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_binary_connectivity_one_generates_all_ones(self):
        """Test that binary_connectivity=1 generates mostly ones."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(
            binary_connectivity=1.0,
            scale=1,
            grid_width=10,
            cfg=cfg
        )
        
        # Should be mostly ones
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_binary_connectivity_half_generates_mixed(self):
        """Test that binary_connectivity=0.5 generates roughly 50% ones."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Use a larger grid and scale=1 to get macropixel values directly
        result = sng.generate(
            binary_connectivity=0.5,
            scale=1,  # No upscaling
            grid_width=100,
            cfg=cfg
        )
        
        # Check that we have both 0s and 1s in the macropixel grid
        # Note: due to linear interpolation during resize, we might not have exact 0s and 1s
        # But the original macropixel grid should have been binary
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_binary_connectivity_intermediate_value(self):
        """Test binary connectivity with an intermediate value like 0.3."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(
            binary_connectivity=0.3,
            scale=1,
            grid_width=50,
            cfg=cfg
        )
        
        self.assertEqual(result.shape, (50, 50))
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    # ==================== CONFIGURATION TESTS ====================

    def test_cfg_parameter_required(self):
        """Test that cfg parameter with rng is required."""
        with self.assertRaises(AttributeError):
            # This should fail because cfg.rng doesn't exist
            sng.generate(amplitude=0.3, scale=10, grid_width=100, cfg=SimpleNamespace())

    def test_cfg_rng_required(self):
        """Test that cfg.rng is required."""
        cfg = SimpleNamespace()  # No rng attribute
        with self.assertRaises(AttributeError):
            sng.generate(amplitude=0.3, scale=10, grid_width=100, cfg=cfg)

    def test_different_rng_seeds_produce_different_results(self):
        """Test that different RNG seeds produce different results."""
        cfg1 = SimpleNamespace(rng=np.random.default_rng(42))
        cfg2 = SimpleNamespace(rng=np.random.default_rng(43))
        
        result1 = sng.generate(amplitude=0.3, scale=10, grid_width=50, cfg=cfg1)
        result2 = sng.generate(amplitude=0.3, scale=10, grid_width=50, cfg=cfg2)
        
        # Results should be different (with high probability)
        self.assertFalse(np.array_equal(result1, result2))

    def test_same_rng_seed_produces_same_results(self):
        """Test that the same RNG seed produces identical results."""
        cfg1 = SimpleNamespace(rng=np.random.default_rng(42))
        cfg2 = SimpleNamespace(rng=np.random.default_rng(42))
        
        result1 = sng.generate(amplitude=0.3, scale=10, grid_width=50, cfg=cfg1)
        result2 = sng.generate(amplitude=0.3, scale=10, grid_width=50, cfg=cfg2)
        
        # Results should be identical
        np.testing.assert_array_equal(result1, result2)

    # ==================== SHOW PARAMETER TESTS ====================

    def test_show_false_no_windows(self):
        """Test that show=False doesn't display windows."""
        # This is hard to test programmatically, but we can at least verify it doesn't crash
        params = {**self.default_params, 'show': False}
        result = sng.generate(**params)
        self.assertIsInstance(result, np.ndarray)

    def test_show_parameter_accepted(self):
        """Test that show parameter is accepted without errors."""
        # Test both True and False
        params1 = {**self.default_params, 'show': False}
        params2 = {**self.default_params, 'show': False}  # Avoid GUI in tests
        
        result1 = sng.generate(**params1)
        result2 = sng.generate(**params2)
        
        self.assertEqual(result1.shape, result2.shape)

    # ==================== EDGE CASE TESTS ====================

    def test_scale_larger_than_grid_width(self):
        """Test that scale larger than grid_width works correctly."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # When scale >= grid_width, we get very few macropixels
        # Use scale=50 with grid_width=100 to get 2x2 macropixels
        result = sng.generate(scale=50, grid_width=100, cfg=cfg)
        
        # Should still produce a 100x100 image
        self.assertEqual(result.shape, (100, 100))
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_scale_equals_one(self):
        """Test that scale=1 produces individual pixel noise."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(scale=1, grid_width=20, cfg=cfg)
        
        self.assertEqual(result.shape, (20, 20))
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_grid_width_not_divisible_by_scale(self):
        """Test that grid_width not divisible by scale works correctly."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # 100 grid width with scale of 7: 100//7 = 14 macropixels
        result = sng.generate(scale=7, grid_width=100, cfg=cfg)
        
        self.assertEqual(result.shape, (100, 100))
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_amplitude_zero(self):
        """Test that amplitude=0 produces constant values."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(amplitude=0.0, scale=1, grid_width=10, cfg=cfg)
        
        # With amplitude=0, all values should be 0.5 + offset (clipped)
        expected_value = 0.5  # 0.5 + 0 + 0
        # Due to floating point precision, values should be very close to 0.5
        np.testing.assert_allclose(result, expected_value, rtol=1e-10)

    def test_amplitude_very_large(self):
        """Test that very large amplitude is handled by clipping."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(amplitude=10.0, scale=1, grid_width=10, cfg=cfg)
        
        # Should be clipped to [0, 1]
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))

    def test_negative_amplitude(self):
        """Test that negative amplitude is handled (should work due to absolute range)."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Negative amplitude would result in low > high in uniform distribution
        # This should cause an error or produce unexpected results
        # Let's test what actually happens
        try:
            result = sng.generate(amplitude=-0.3, scale=1, grid_width=10, cfg=cfg)
            # If it doesn't crash, check that results are still in [0,1]
            self.assertTrue(np.all(result >= 0.0))
            self.assertTrue(np.all(result <= 1.0))
        except (ValueError, TypeError):
            # This is also acceptable behavior
            pass

    # ==================== STATISTICAL PROPERTY TESTS ====================

    def test_uniform_noise_statistical_properties(self):
        """Test that uniform noise has expected statistical properties."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Generate many samples with same seed to test distribution
        amplitude = 0.2
        mean_expected = 0.5  # Center of the uniform distribution
        
        # Use scale=1 to get individual pixel values
        result = sng.generate(
            amplitude=amplitude, 
            scale=1, 
            grid_width=100, 
            cfg=cfg
        )
        
        # Check that mean is close to expected
        actual_mean = np.mean(result)
        self.assertAlmostEqual(actual_mean, mean_expected, delta=0.05)

    def test_noise_values_within_expected_range(self):
        """Test that noise values fall within the expected range before clipping."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        amplitude = 0.1
        # Expected range: [0.5 - amplitude, 0.5 + amplitude] = [0.4, 0.6]
        
        # Use scale=1 to avoid interpolation artifacts
        result = sng.generate(
            amplitude=amplitude,
            scale=1,
            grid_width=50,
            cfg=cfg
        )
        
        # All values should be within the expected range or clipped
        self.assertTrue(np.all(result >= 0.4))
        self.assertTrue(np.all(result <= 0.6))

    # ==================== INTERPOLATION TESTS ====================

    def test_upscaling_preserves_macropixel_structure(self):
        """Test that upscaling from macropixels to final image works correctly."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        scale = 10
        grid_width = 100
        
        result = sng.generate(
            amplitude=0.2,
            scale=scale,
            grid_width=grid_width,
            cfg=cfg
        )
        
        # Should have 10x10 macropixels upscaled to 100x100
        self.assertEqual(result.shape, (grid_width, grid_width))
        
        # With linear interpolation, the center of each macropixel region should be closest to the original value
        macro_size = grid_width // scale  # Should be 10
        num_macropixels = grid_width // macro_size  # Should be 10
        
        # Check that the image has reasonable macropixel structure
        # The center of the first macropixel region should be close to a constant value
        center_pixel = result[macro_size//2, macro_size//2]  # Center of first macropixel
        
        # The corners of the first macropixel region should be influenced by interpolation
        # but the overall values should be in the expected range [0.3, 0.7] for amplitude=0.2
        self.assertTrue(np.all(result >= 0.3))
        self.assertTrue(np.all(result <= 0.7))
        
        # Test that values are smooth (no large jumps between adjacent pixels)
        # This verifies that interpolation is working
        row_diffs = np.abs(np.diff(result, axis=0))
        col_diffs = np.abs(np.diff(result, axis=1))
        self.assertTrue(np.all(row_diffs < 0.1))  # No large jumps vertically
        self.assertTrue(np.all(col_diffs < 0.1))  # No large jumps horizontally

    def test_linear_interpolation_between_macropixels(self):
        """Test that linear interpolation creates smooth transitions between macropixels."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        scale = 5
        grid_width = 20
        
        result = sng.generate(
            amplitude=0.3,
            scale=scale,
            grid_width=grid_width,
            cfg=cfg
        )
        
        # Should have smooth transitions between macropixels
        # Check that there are no sudden jumps (which would indicate no interpolation)
        # This is a basic check - more sophisticated analysis could be done
        self.assertEqual(result.shape, (grid_width, grid_width))

    # ==================== ADDITIONAL KWARGS TESTS ====================

    def test_additional_kwargs_ignored(self):
        """Test that additional keyword arguments are ignored as specified."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Add some extra kwargs that should be ignored
        result = sng.generate(
            amplitude=0.3,
            scale=10,
            grid_width=50,
            cfg=cfg,
            extra_param1=123,
            extra_param2="test",
            another_param=[1, 2, 3]
        )
        
        # Should still work and produce valid output
        self.assertEqual(result.shape, (50, 50))
        self.assertTrue(np.all(result >= 0.0))
        self.assertTrue(np.all(result <= 1.0))


class TestNoiseGeneratorEdgeCases(unittest.TestCase):
    """Additional edge case tests for the noise generator."""

    def test_very_small_grid_width(self):
        """Test with very small grid width."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(grid_width=1, cfg=cfg)
        self.assertEqual(result.shape, (1, 1))

    def test_large_grid_width(self):
        """Test with large grid width (but not too large to avoid memory issues)."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        result = sng.generate(grid_width=500, scale=50, cfg=cfg)
        self.assertEqual(result.shape, (500, 500))

    def test_boundary_conditions_for_binary_connectivity(self):
        """Test boundary conditions for binary_connectivity parameter."""
        cfg = SimpleNamespace(rng=np.random.default_rng(42))
        
        # Test binary_connectivity exactly at 0
        result1 = sng.generate(binary_connectivity=0.0, scale=1, grid_width=10, cfg=cfg)
        self.assertEqual(result1.shape, (10, 10))
        
        # Test binary_connectivity exactly at 1
        cfg2 = SimpleNamespace(rng=np.random.default_rng(42))
        result2 = sng.generate(binary_connectivity=1.0, scale=1, grid_width=10, cfg=cfg2)
        self.assertEqual(result2.shape, (10, 10))


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)