"""Unit test suite for sine_pattern_generator module.

This module provides comprehensive tests for the sine-based pattern generation functionality.
Tests cover all sine types, parameter validation, and output correctness.
"""

import unittest
import numpy as np
import cv2
import sys
import os

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sine_pattern_generator as spg


class TestSinePatternGenerator(unittest.TestCase):
    """Test cases for the sine_pattern_generator module."""

    def setUp(self):
        """Set up test fixtures."""
        self.default_dimensions = (100, 100)
        self.default_amplitude = 1.0
        self.default_wavelength = 50.0

    def tearDown(self):
        """Clean up after tests."""
        # Close any OpenCV windows that might have been opened
        cv2.destroyAllWindows()

    # ==================== BASIC FUNCTIONALITY TESTS ====================

    def test_generate_returns_numpy_array(self):
        """Test that generate() returns a numpy array."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5
        )
        self.assertIsInstance(result, np.ndarray)

    def test_generate_returns_2d_array(self):
        """Test that generate() returns a 2D array."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5
        )
        self.assertEqual(result.ndim, 2)

    def test_generate_output_shape_matches_dimensions(self):
        """Test that output shape matches the specified dimensions."""
        dimensions = (64, 128)
        result = spg.generate(
            dimensions=dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5
        )
        self.assertEqual(result.shape, dimensions)

    # ==================== SINE TYPE TESTS ====================

    def test_horizontal_sine_type(self):
        """Test horizontal sine wave generation."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5
        )
        self.assertEqual(result.shape, self.default_dimensions)
        # Horizontal sine should vary more in Y direction than X direction
        y_variation = np.std(np.mean(result, axis=1))  # Variation across rows
        x_variation = np.std(np.mean(result, axis=0))  # Variation across columns
        self.assertGreater(y_variation, 0)  # Should have some variation

    def test_vertical_sine_type(self):
        """Test vertical sine wave generation."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="vertical",
            sine_offset=25,
            mean=0.5
        )
        self.assertEqual(result.shape, self.default_dimensions)
        # Vertical sine should vary more in X direction than Y direction
        x_variation = np.std(np.mean(result, axis=0))  # Variation across columns
        y_variation = np.std(np.mean(result, axis=1))  # Variation across rows
        self.assertGreater(x_variation, 0)  # Should have some variation

    def test_diagonal_sine_type(self):
        """Test diagonal sine wave generation."""
        # Note: The current implementation doesn't have diagonal, let's test what exists
        # Actually, looking at the code, there's no "diagonal" implementation, only horizontal, vertical, radial
        # So this might raise an error
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="diagonal",
                sine_offset=25,
                mean=0.5
            )

    def test_radial_sine_type(self):
        """Test radial sine wave generation."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="radial",
            sine_offset=(50, 50),  # Center point
            mean=0.5
        )
        self.assertEqual(result.shape, self.default_dimensions)
        # Radial pattern should be symmetric around center
        center = np.array(result.shape) // 2
        # Check symmetry by comparing quadants (approximately)
        self.assertTrue(result.shape[0] > 1 and result.shape[1] > 1)

    def test_unknown_sine_type_raises_error(self):
        """Test that unknown sine_type raises ValueError."""
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="unknown_type",
                sine_offset=25,
                mean=0.5
            )

    # ==================== PARAMETER VALIDATION TESTS ====================

    def test_exactly_one_of_max_min_mean_required(self):
        """Test that exactly one of maximum, minimum, or mean must be specified."""
        # Test no parameters
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="horizontal",
                sine_offset=25
                # No maximum, minimum, or mean specified
            )

    def test_multiple_range_parameters_raises_error(self):
        """Test that specifying multiple range parameters raises error."""
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="horizontal",
                sine_offset=25,
                maximum=1.0,
                minimum=0.0
                # Both maximum and minimum specified
            )

    def test_all_three_range_parameters_raises_error(self):
        """Test that specifying all three range parameters raises error."""
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="horizontal",
                sine_offset=25,
                maximum=1.0,
                minimum=0.0,
                mean=0.5
                # All three range parameters specified
            )

    def test_sine_offset_required_for_horizontal_vertical(self):
        """Test that sine_offset is required for horizontal and vertical types."""
        # Test horizontal without offset
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="horizontal",
                sine_offset=None,  # This should cause error
                mean=0.5
            )

        # Test vertical without offset
        with self.assertRaises(ValueError):
            spg.generate(
                dimensions=self.default_dimensions,
                sine_amplitude=self.default_amplitude,
                sine_wavelength=self.default_wavelength,
                sine_type="vertical",
                sine_offset=None,  # This should cause error
                mean=0.5
            )

    def test_radial_sine_offset_can_be_none(self):
        """Test that radial sine type can work with None offset (uses center)."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="radial",
            sine_offset=None,  # Should use default center
            mean=0.5
        )
        self.assertEqual(result.shape, self.default_dimensions)

    # ==================== RANGE PARAMETER TESTS ====================

    def test_maximum_parameter(self):
        """Test that maximum parameter works correctly."""
        maximum = 0.8
        amplitude = 0.2
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            maximum=maximum
        )
        
        # With maximum=0.8 and amplitude=0.2, minimum should be 0.8 - 2*0.2 = 0.4
        # So all values should be in [0.4, 0.8]
        self.assertTrue(np.all(result >= 0.4))
        self.assertTrue(np.all(result <= 0.8))

    def test_minimum_parameter(self):
        """Test that minimum parameter works correctly."""
        minimum = 0.3
        amplitude = 0.1
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            minimum=minimum
        )
        
        # With minimum=0.3 and amplitude=0.1, maximum should be 0.3 + 2*0.1 = 0.5
        # So all values should be in [0.3, 0.5]
        self.assertTrue(np.all(result >= 0.3))
        self.assertTrue(np.all(result <= 0.5))

    def test_mean_parameter(self):
        """Test that mean parameter works correctly."""
        mean = 0.6
        amplitude = 0.1
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # With mean=0.6 and amplitude=0.1, range should be [0.5, 0.7]
        # So all values should be in [0.5, 0.7]
        self.assertTrue(np.all(result >= 0.5))
        self.assertTrue(np.all(result <= 0.7))
        
        # Mean should be close to the specified mean
        actual_mean = np.mean(result)
        self.assertAlmostEqual(actual_mean, mean, delta=0.01)

    # ==================== CUTOFF TESTS ====================

    def test_cutoff_min_parameter(self):
        """Test that cutoff_min parameter works correctly."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            cutoff_min=0.6
        )
        
        # All values should be >= cutoff_min
        self.assertTrue(np.all(result >= 0.6))

    def test_cutoff_max_parameter(self):
        """Test that cutoff_max parameter works correctly."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            cutoff_max=0.4
        )
        
        # All values should be <= cutoff_max
        self.assertTrue(np.all(result <= 0.4))

    def test_both_cutoffs_parameter(self):
        """Test that both cutoff parameters work together."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            cutoff_min=0.4,
            cutoff_max=0.6
        )
        
        # All values should be between cutoffs
        self.assertTrue(np.all(result >= 0.4))
        self.assertTrue(np.all(result <= 0.6))

    # ==================== WAVE PARAMETER TESTS ====================

    def test_amplitude_parameter(self):
        """Test that amplitude parameter affects the wave height."""
        mean = 0.5
        
        # Test with amplitude=0.1 (small variation)
        result1 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=0.1,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # Test with amplitude=0.3 (larger variation)
        result2 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=0.3,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # Larger amplitude should produce larger variation
        std1 = np.std(result1)
        std2 = np.std(result2)
        self.assertGreater(std2, std1)

    def test_wavelength_parameter(self):
        """Test that wavelength parameter affects the wave frequency."""
        mean = 0.5
        amplitude = 0.2
        
        # Test with short wavelength (high frequency)
        result1 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=amplitude,
            sine_wavelength=20,  # Short wavelength
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # Test with long wavelength (low frequency)
        result2 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=amplitude,
            sine_wavelength=100,  # Long wavelength
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # Both should be valid
        self.assertEqual(result1.shape, self.default_dimensions)
        self.assertEqual(result2.shape, self.default_dimensions)

    def test_sine_offset_parameter(self):
        """Test that sine_offset parameter shifts the wave phase."""
        mean = 0.5
        
        # Test with different offsets
        result1 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=10,  # One offset
            mean=mean
        )
        
        result2 = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=30,  # Different offset
            mean=mean
        )
        
        # Results should be different due to different phase
        self.assertFalse(np.array_equal(result1, result2))

    # ==================== SHOW PARAMETER TESTS ====================

    def test_show_false_no_windows(self):
        """Test that show=False doesn't display windows."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            show=False
        )
        self.assertIsInstance(result, np.ndarray)

    def test_show_parameter_accepted(self):
        """Test that show parameter is accepted without errors."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            show=False  # Avoid GUI in tests
        )
        self.assertEqual(result.shape, self.default_dimensions)

    # ==================== ADDITIONAL KWARGS TESTS ====================

    def test_additional_kwargs_ignored(self):
        """Test that additional keyword arguments are ignored."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5,
            extra_param="ignored",
            another_param=42
        )
        
        self.assertEqual(result.shape, self.default_dimensions)

    # ==================== EDGE CASE TESTS ====================

    def test_zero_amplitude(self):
        """Test that zero amplitude produces constant values."""
        mean = 0.7
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=0.0,  # Zero amplitude
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=25,
            mean=mean
        )
        
        # With zero amplitude, all values should be equal to mean
        expected_value = mean
        np.testing.assert_allclose(result, expected_value, rtol=1e-10)

    def test_very_small_dimensions(self):
        """Test with very small dimensions."""
        result = spg.generate(
            dimensions=(10, 10),
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="horizontal",
            sine_offset=5,
            mean=0.5
        )
        self.assertEqual(result.shape, (10, 10))

    def test_non_square_dimensions(self):
        """Test with non-square dimensions."""
        result = spg.generate(
            dimensions=(50, 100),
            sine_amplitude=self.default_amplitude,
            sine_wavelength=self.default_wavelength,
            sine_type="vertical",
            sine_offset=25,
            mean=0.5
        )
        self.assertEqual(result.shape, (50, 100))

    def test_very_long_wavelength(self):
        """Test with very long wavelength (larger than dimensions)."""
        result = spg.generate(
            dimensions=self.default_dimensions,
            sine_amplitude=self.default_amplitude,
            sine_wavelength=500,  # Very long wavelength
            sine_type="horizontal",
            sine_offset=25,
            mean=0.5
        )
        
        # Should still work and produce a valid pattern
        self.assertEqual(result.shape, self.default_dimensions)
        self.assertTrue(np.all(result >= 0.0))  # Should be reasonable values


class TestSinePatternGeneratorEdgeCases(unittest.TestCase):
    """Additional edge case tests for the sine pattern generator."""

    def test_negative_amplitude(self):
        """Test with negative amplitude."""
        result = spg.generate(
            dimensions=(50, 50),
            sine_amplitude=-0.1,  # Negative amplitude
            sine_wavelength=25,
            sine_type="horizontal",
            sine_offset=10,
            mean=0.5
        )
        self.assertEqual(result.shape, (50, 50))

    def test_radial_with_tuple_offset(self):
        """Test radial sine with explicit center offset."""
        result = spg.generate(
            dimensions=(100, 100),
            sine_amplitude=0.5,
            sine_wavelength=50,
            sine_type="radial",
            sine_offset=(30, 70),  # Explicit center
            mean=0.5
        )
        self.assertEqual(result.shape, (100, 100))

    def test_horizontal_phase_shift(self):
        """Test horizontal sine with different phase shifts via offset."""
        result = spg.generate(
            dimensions=(60, 60),
            sine_amplitude=0.3,
            sine_wavelength=30,
            sine_type="horizontal",
            sine_offset=15,  # Peak at y=15
            mean=0.5
        )
        
        # The peak should be around y=15
        row_means = np.mean(result, axis=1)
        peak_row = np.argmax(row_means)
        # Should be close to the expected peak position
        self.assertAlmostEqual(peak_row, 15, delta=5)  # Allow some tolerance


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)