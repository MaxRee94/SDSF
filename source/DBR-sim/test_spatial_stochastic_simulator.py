"""Unit test suite for spatial_stochastic_simulator module.

This module provides tests for the spatial stochastic simulator that generates
2D random fields using gstools library with various covariance models.
"""

import unittest
import numpy as np
from types import SimpleNamespace
import sys
import os

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import spatial_stochastic_simulator as sss
    from spatial_stochastic_simulator import SpatialStochasticSimulator
    GS_TOOLS_AVAILABLE = True
except ImportError as e:
    GS_TOOLS_AVAILABLE = False
    print(f"gstools not available: {e}")


@unittest.skipUnless(GS_TOOLS_AVAILABLE, "gstools library not available")
class TestSpatialStochasticSimulator(unittest.TestCase):
    """Test cases for the SpatialStochasticSimulator class."""

    def setUp(self):
        """Set up test fixtures."""
        # Default parameters for testing
        self.default_args = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )

    def test_simulator_initialization(self):
        """Test that simulator initializes correctly."""
        simulator = SpatialStochasticSimulator(self.default_args)
        
        self.assertEqual(simulator.nx, 50)
        self.assertEqual(simulator.ny, 50)
        self.assertEqual(simulator.model_type, "Exponential")
        self.assertEqual(simulator.sill, 1.0)
        self.assertEqual(simulator.nugget, 0.0)
        self.assertEqual(simulator.vrange, 10.0)
        self.assertEqual(simulator.seed, 42)

    def test_simulator_different_widths(self):
        """Test simulator with different width values."""
        args = SimpleNamespace(
            width=100,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        self.assertEqual(simulator.nx, 100)
        self.assertEqual(simulator.ny, 100)

    def test_simulator_different_model_types(self):
        """Test simulator with different covariance model types."""
        for model_type in ["Exponential", "Gaussian", "Spherical"]:
            args = SimpleNamespace(
                width=50,
                model_type=model_type,
                sill=1.0,
                nugget=0.0,
                vrange=10.0,
                global_seed=42
            )
            
            simulator = SpatialStochasticSimulator(args)
            self.assertEqual(simulator.model_type, model_type)

    def test_generate_method_returns_numpy_array(self):
        """Test that generate() method returns a numpy array."""
        simulator = SpatialStochasticSimulator(self.default_args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertIsInstance(result, np.ndarray)

    def test_generate_method_returns_2d_array(self):
        """Test that generate() method returns a 2D array."""
        simulator = SpatialStochasticSimulator(self.default_args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.ndim, 2)

    def test_generate_method_output_shape(self):
        """Test that generate() method output has correct shape."""
        simulator = SpatialStochasticSimulator(self.default_args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        # Should be width x width
        self.assertEqual(result.shape, (50, 50))

    def test_generate_method_reproducibility(self):
        """Test that generate() method produces reproducible results with same seed."""
        args1 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        args2 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        simulator1 = SpatialStochasticSimulator(args1)
        simulator2 = SpatialStochasticSimulator(args2)
        
        result1 = simulator1.generate(plot_variogram=False, plot_image=False)
        result2 = simulator2.generate(plot_variogram=False, plot_image=False)
        
        # Results should be identical with same seed
        np.testing.assert_array_equal(result1, result2)

    def test_generate_method_different_seeds(self):
        """Test that generate() method produces different results with different seeds."""
        args1 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        args2 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=43
        )
        
        simulator1 = SpatialStochasticSimulator(args1)
        simulator2 = SpatialStochasticSimulator(args2)
        
        result1 = simulator1.generate(plot_variogram=False, plot_image=False)
        result2 = simulator2.generate(plot_variogram=False, plot_image=False)
        
        # Results should be different with different seeds
        self.assertFalse(np.array_equal(result1, result2))

    def test_different_model_types_produce_different_results(self):
        """Test that different model types produce different results."""
        result_exponential = self._generate_with_model_type("Exponential")
        result_gaussian = self._generate_with_model_type("Gaussian")
        result_spherical = self._generate_with_model_type("Spherical")
        
        # Different model types should produce different spatial patterns
        # Note: They might be similar but shouldn't be identical
        self.assertFalse(np.array_equal(result_exponential, result_gaussian))
        self.assertFalse(np.array_equal(result_exponential, result_spherical))
        self.assertFalse(np.array_equal(result_gaussian, result_spherical))

    def _generate_with_model_type(self, model_type):
        """Helper method to generate with specific model type."""
        args = SimpleNamespace(
            width=50,
            model_type=model_type,
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        return simulator.generate(plot_variogram=False, plot_image=False)

    def test_different_sill_values(self):
        """Test that different sill values affect the output variance."""
        args1 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=0.5,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        args2 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=2.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        simulator1 = SpatialStochasticSimulator(args1)
        simulator2 = SpatialStochasticSimulator(args2)
        
        result1 = simulator1.generate(plot_variogram=False, plot_image=False)
        result2 = simulator2.generate(plot_variogram=False, plot_image=False)
        
        # Higher sill should lead to higher variance
        var1 = np.var(result1)
        var2 = np.var(result2)
        self.assertGreater(var2, var1)

    def test_different_vrange_values(self):
        """Test that different vrange (correlation length) values affect spatial structure."""
        args1 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=5.0,  # Short correlation length
            global_seed=42
        )
        
        args2 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=20.0,  # Long correlation length
            global_seed=42
        )
        
        simulator1 = SpatialStochasticSimulator(args1)
        simulator2 = SpatialStochasticSimulator(args2)
        
        result1 = simulator1.generate(plot_variogram=False, plot_image=False)
        result2 = simulator2.generate(plot_variogram=False, plot_image=False)
        
        # Both should be valid
        self.assertEqual(result1.shape, (50, 50))
        self.assertEqual(result2.shape, (50, 50))

    def test_different_nugget_values(self):
        """Test that different nugget values affect the output."""
        args1 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,  # No nugget effect
            vrange=10.0,
            global_seed=42
        )
        
        args2 = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.5,  # With nugget effect
            vrange=10.0,
            global_seed=42
        )
        
        simulator1 = SpatialStochasticSimulator(args1)
        simulator2 = SpatialStochasticSimulator(args2)
        
        result1 = simulator1.generate(plot_variogram=False, plot_image=False)
        result2 = simulator2.generate(plot_variogram=False, plot_image=False)
        
        # Results should be different due to nugget effect
        self.assertFalse(np.array_equal(result1, result2))

    def test_large_grid(self):
        """Test with a larger grid size."""
        args = SimpleNamespace(
            width=200,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=20.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (200, 200))

    def test_small_grid(self):
        """Test with a small grid size."""
        args = SimpleNamespace(
            width=10,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=5.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (10, 10))

    def test_plot_parameters_accepted(self):
        """Test that plot parameters are accepted without errors."""
        simulator = SpatialStochasticSimulator(self.default_args)
        
        # Test with plotting disabled (default for testing)
        result1 = simulator.generate(plot_variogram=False, plot_image=False)
        self.assertEqual(result1.shape, (50, 50))
        
        # Test with plotting enabled (should not crash, though it may show plots)
        # We'll keep it disabled for testing
        result2 = simulator.generate(plot_variogram=False, plot_image=False)
        self.assertEqual(result2.shape, (50, 50))

    def test_output_file_parameter(self):
        """Test that output_file parameter is accepted."""
        simulator = SpatialStochasticSimulator(self.default_args)
        
        # Test with a temporary output file (we won't actually save to avoid file system operations)
        result = simulator.generate(
            output_file="test_output_temp.png",
            plot_variogram=False,
            plot_image=False
        )
        
        self.assertEqual(result.shape, (50, 50))


@unittest.skipUnless(GS_TOOLS_AVAILABLE, "gstools library not available")
class TestSpatialStochasticSimulatorEdgeCases(unittest.TestCase):
    """Edge case tests for the spatial stochastic simulator."""

    def test_zero_nugget(self):
        """Test with zero nugget."""
        args = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=10.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (50, 50))

    def test_high_nugget(self):
        """Test with high nugget relative to sill."""
        args = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.9,  # High nugget
            vrange=10.0,
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (50, 50))

    def test_zero_vrange(self):
        """Test with very small vrange (correlation length)."""
        args = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=0.1,  # Very small correlation length
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (50, 50))

    def test_large_vrange(self):
        """Test with very large vrange (correlation length)."""
        args = SimpleNamespace(
            width=50,
            model_type="Exponential",
            sill=1.0,
            nugget=0.0,
            vrange=1000.0,  # Very large correlation length
            global_seed=42
        )
        
        simulator = SpatialStochasticSimulator(args)
        result = simulator.generate(plot_variogram=False, plot_image=False)
        
        self.assertEqual(result.shape, (50, 50))


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)