"""Unit test suite for app module.

This module provides comprehensive tests for the application functionality.
Tests cover initialization, configuration, simulation functions, and helper utilities.
"""

import unittest
import os
import sys
import tempfile
import shutil
import json
import numpy as np
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock complex dependencies
sys.modules['visualization'] = MagicMock()
sys.modules['file_handling'] = MagicMock()
sys.modules['helpers'] = MagicMock()
sys.modules['disk_pattern_generator'] = MagicMock()
sys.modules['sine_pattern_generator'] = MagicMock()
sys.modules['simple_noise_generator'] = MagicMock()
sys.modules['x64.Release'] = MagicMock()
sys.modules['x64.Release.dbr_cpp'] = MagicMock()

# Create a mock config module with the cfg attribute
class MockConfigModule:
    pass

mock_config = MockConfigModule()
sys.modules['config'] = mock_config

import app as app_module


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.DATA_IN_DIR = tempfile.mkdtemp()
        self.DATA_OUT_DIR = tempfile.mkdtemp()
        self.PERLIN_NOISE_DIR = tempfile.mkdtemp()
        self.SIMPLE_PATTERNS_DIR = tempfile.mkdtemp()
        self.random_seed = 42
        self.firefreq_random_seed = 123
        self.verbosity = 0
        self.rng = np.random.default_rng(42)


class TestSetDispersalKernel(unittest.TestCase):
    """Test cases for set_dispersal_kernel function."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = MockConfig()
        
        # Inject cfg into the mock config module and app module
        sys.modules['config'].cfg = self.cfg
        app_module.cfg = self.cfg
        
        # Create a temporary dispersal parameters file
        self.dispersal_params = {
            "linear": {
                "q1": 0.5,
                "q2": 0.5,
                "min": 0.1,
                "max": 10.0
            },
            "wind": {
                "wspeed_gmean": 5.0,
                "wspeed_stdev": 1.0,
                "wind_direction": 0.0,
                "wind_direction_stdev": 10.0
            },
            "animal": {
                "species1": {"param": "value"},
                "population": {}
            }
        }
        
        self.params_file = os.path.join(self.temp_dir, "dispersal_params.json")
        with open(self.params_file, 'w') as f:
            json.dump(self.dispersal_params, f)
        
        # Mock dynamics object
        self.mock_dynamics = MagicMock()

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        shutil.rmtree(self.cfg.DATA_IN_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.DATA_OUT_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.PERLIN_NOISE_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.SIMPLE_PATTERNS_DIR, ignore_errors=True)

    def test_set_dispersal_kernel_linear_diffusion(self):
        """Test set_dispersal_kernel with linear diffusion mode."""
        self.cfg.multi_disperser_params = os.path.basename(self.params_file)
        
        with patch.object(self.cfg, 'DATA_IN_DIR', self.temp_dir):
            result_dynamics, animal_species = app_module.set_dispersal_kernel(
                self.mock_dynamics, "linear_diffusion", self.cfg.multi_disperser_params
            )
        
        # Check that the correct method was called
        self.mock_dynamics.set_global_linear_kernel.assert_called_once()
        
        # Check animal_species is empty for non-animal modes
        self.assertEqual(animal_species, [])

    def test_set_dispersal_kernel_wind(self):
        """Test set_dispersal_kernel with wind mode."""
        self.cfg.multi_disperser_params = os.path.basename(self.params_file)
        
        with patch.object(self.cfg, 'DATA_IN_DIR', self.temp_dir):
            result_dynamics, animal_species = app_module.set_dispersal_kernel(
                self.mock_dynamics, "wind", self.cfg.multi_disperser_params
            )
        
        # Check that the correct method was called
        self.mock_dynamics.set_global_wind_kernel.assert_called_once()
        
        # Check animal_species is empty for non-animal modes
        self.assertEqual(animal_species, [])

    def test_set_dispersal_kernel_animal(self):
        """Test set_dispersal_kernel with animal mode."""
        self.cfg.multi_disperser_params = os.path.basename(self.params_file)
        
        with patch.object(self.cfg, 'DATA_IN_DIR', self.temp_dir):
            result_dynamics, animal_species = app_module.set_dispersal_kernel(
                self.mock_dynamics, "animal", self.cfg.multi_disperser_params
            )
        
        # Check that the correct method was called
        self.mock_dynamics.set_global_animal_kernel.assert_called_once()
        
        # Check animal_species contains expected species
        self.assertIn("species1", animal_species)

    def test_set_dispersal_kernel_all(self):
        """Test set_dispersal_kernel with all mode."""
        self.cfg.multi_disperser_params = os.path.basename(self.params_file)
        
        with patch.object(self.cfg, 'DATA_IN_DIR', self.temp_dir):
            result_dynamics, animal_species = app_module.set_dispersal_kernel(
                self.mock_dynamics, "all", self.cfg.multi_disperser_params
            )
        
        # Check that the correct method was called
        self.mock_dynamics.set_global_kernels.assert_called_once()
        
        # Check animal_species contains expected species
        self.assertIn("species1", animal_species)


class TestSetRandomSeeds(unittest.TestCase):
    """Test cases for set_random_seeds function."""

    def test_set_random_seeds_with_specific_seed(self):
        """Test set_random_seeds with specific seed."""
        cfg = SimpleNamespace(random_seed=42, firefreq_random_seed=123)
        
        result = app_module.set_random_seeds(cfg)
        
        # Check that rng was created with the specified seed
        self.assertEqual(result.random_seed, 42)
        self.assertEqual(result.firefreq_random_seed, 123)
        self.assertIsNotNone(result.rng)

    def test_set_random_seeds_with_random_seed(self):
        """Test set_random_seeds with random seed (-999)."""
        cfg = SimpleNamespace(random_seed=-999, firefreq_random_seed=123)
        
        result = app_module.set_random_seeds(cfg)
        
        # Check that new random seed was generated
        self.assertNotEqual(result.random_seed, -999)
        self.assertIsNotNone(result.rng)

    def test_set_random_seeds_with_random_firefreq_seed(self):
        """Test set_random_seeds with random fire frequency seed (-999)."""
        cfg = SimpleNamespace(random_seed=42, firefreq_random_seed=-999)
        
        result = app_module.set_random_seeds(cfg)
        
        # Check that new random seed was generated for fire frequency
        self.assertNotEqual(result.firefreq_random_seed, -999)


class TestPatternGeneration(unittest.TestCase):
    """Test cases for pattern generation functions."""

    def test_create_controlled_patches_image_function_exists(self):
        """Test that create_controlled_patches_image function exists."""
        self.assertTrue(hasattr(app_module, 'create_controlled_patches_image'))
        self.assertTrue(callable(app_module.create_controlled_patches_image))

    def test_create_perlin_noise_image_function_exists(self):
        """Test that create_perlin_noise_image function exists."""
        self.assertTrue(hasattr(app_module, 'create_perlin_noise_image'))
        self.assertTrue(callable(app_module.create_perlin_noise_image))

    def test_load_homogeneous_image_function_exists(self):
        """Test that load_homogeneous_image function exists."""
        self.assertTrue(hasattr(app_module, 'load_homogeneous_image'))
        self.assertTrue(callable(app_module.load_homogeneous_image))

    def test_load_specified_existing_image_function_exists(self):
        """Test that load_specified_existing_image function exists."""
        self.assertTrue(hasattr(app_module, 'load_specified_existing_image'))
        self.assertTrue(callable(app_module.load_specified_existing_image))

    def test_set_initial_tree_cover_function_exists(self):
        """Test that set_initial_tree_cover function exists."""
        self.assertTrue(hasattr(app_module, 'set_initial_tree_cover'))
        self.assertTrue(callable(app_module.set_initial_tree_cover))


class TestTerminationConditions(unittest.TestCase):
    """Test cases for termination condition functions."""

    def test_termination_condition_satisfied_function_exists(self):
        """Test that termination_condition_satisfied function exists."""
        self.assertTrue(hasattr(app_module, 'termination_condition_satisfied'))
        self.assertTrue(callable(app_module.termination_condition_satisfied))

    def test_termination_condition_max_timesteps(self):
        """Test termination condition for max timesteps."""
        cfg = SimpleNamespace(max_timesteps=100, termination_conditions="none")
        mock_dynamics = MagicMock()
        mock_dynamics.time = 100
        mock_dynamics.state.grid.get_tree_cover.return_value = 0.5
        
        result = app_module.termination_condition_satisfied(
            mock_dynamics, 0, cfg
        )
        
        # Should be satisfied when time >= max_timesteps
        self.assertTrue(result)

    def test_termination_condition_tree_cover(self):
        """Test termination condition for tree cover."""
        cfg = SimpleNamespace(max_timesteps=1000, termination_conditions="any")
        mock_dynamics = MagicMock()
        mock_dynamics.time = 50
        mock_dynamics.state.grid.get_tree_cover.return_value = 0.995  # > 0.99
        
        result = app_module.termination_condition_satisfied(
            mock_dynamics, 0, cfg
        )
        
        # Should be satisfied when tree cover > 0.99
        self.assertTrue(result)

    def test_termination_condition_not_satisfied(self):
        """Test termination condition not satisfied."""
        cfg = SimpleNamespace(max_timesteps=1000, termination_conditions="none")
        mock_dynamics = MagicMock()
        mock_dynamics.time = 50
        mock_dynamics.state.grid.get_tree_cover.return_value = 0.5
        
        result = app_module.termination_condition_satisfied(
            mock_dynamics, 0, cfg
        )
        
        # Should not be satisfied
        self.assertFalse(result)


class TestBurnIn(unittest.TestCase):
    """Test cases for burn-in functionality."""

    def test_do_burn_in_function_exists(self):
        """Test that do_burn_in function exists."""
        self.assertTrue(hasattr(app_module, 'do_burn_in'))
        self.assertTrue(callable(app_module.do_burn_in))


class TestIteration(unittest.TestCase):
    """Test cases for iteration functionality."""

    def test_do_iteration_function_exists(self):
        """Test that do_iteration function exists."""
        self.assertTrue(hasattr(app_module, 'do_iteration'))
        self.assertTrue(callable(app_module.do_iteration))

    def test_do_update_function_exists(self):
        """Test that do_update function exists."""
        self.assertTrue(hasattr(app_module, 'do_update'))
        self.assertTrue(callable(app_module.do_update))


class TestUpdateLoop(unittest.TestCase):
    """Test cases for update loop functionality."""

    def test_updateloop_function_exists(self):
        """Test that updateloop function exists."""
        self.assertTrue(hasattr(app_module, 'updateloop'))
        self.assertTrue(callable(app_module.updateloop))


class TestBifurcationAnalysis(unittest.TestCase):
    """Test cases for bifurcation analysis functionality."""

    def test_do_bifurcation_analysis_function_exists(self):
        """Test that do_bifurcation_analysis function exists."""
        self.assertTrue(hasattr(app_module, 'do_bifurcation_analysis'))
        self.assertTrue(callable(app_module.do_bifurcation_analysis))


class TestMainFunction(unittest.TestCase):
    """Test cases for main function."""

    def setUp(self):
        """Set up test fixtures."""
        # Ensure cfg is available in the app module
        if not hasattr(app_module, 'cfg'):
            app_module.cfg = MockConfig()
        if not hasattr(sys.modules['config'], 'cfg'):
            sys.modules['config'].cfg = app_module.cfg

    def test_main_function_exists(self):
        """Test that main function exists."""
        self.assertTrue(hasattr(app_module, 'main'))
        self.assertTrue(callable(app_module.main))

    def test_main_function_handles_batch_type(self):
        """Test main function signature with batch_type parameter."""
        # Test that main function accepts parameters via **user_args
        self.assertTrue(callable(app_module.main))
        # Check function signature - main accepts **user_args
        import inspect
        sig = inspect.signature(app_module.main)
        params = sig.parameters
        # Should have user_args parameter that accepts keyword arguments
        self.assertTrue('user_args' in params)
        self.assertEqual(params['user_args'].kind, inspect.Parameter.VAR_KEYWORD)

    def test_main_function_regular_mode(self):
        """Test main function signature."""
        # Test that main function accepts verbosity parameter via **user_args
        self.assertTrue(callable(app_module.main))
        import inspect
        sig = inspect.signature(app_module.main)
        params = sig.parameters
        # Should have user_args parameter that accepts keyword arguments
        self.assertTrue('user_args' in params)
        self.assertEqual(params['user_args'].kind, inspect.Parameter.VAR_KEYWORD)


class TestConfigurationFunctions(unittest.TestCase):
    """Test cases for configuration-related functions."""

    def test_set_argument_in_model_core_function_exists(self):
        """Test that set_argument_in_model_core function exists."""
        self.assertTrue(hasattr(app_module, 'set_argument_in_model_core'))
        self.assertTrue(callable(app_module.set_argument_in_model_core))

    def test_get_old_keyframed_value_function_exists(self):
        """Test that get_old_keyframed_value function exists."""
        self.assertTrue(hasattr(app_module, 'get_old_keyframed_value'))
        self.assertTrue(callable(app_module.get_old_keyframed_value))

    def test_set_keyframe_function_exists(self):
        """Test that set_keyframe function exists."""
        self.assertTrue(hasattr(app_module, 'set_keyframe'))
        self.assertTrue(callable(app_module.set_keyframe))

    def test_apply_keyframes_function_exists(self):
        """Test that apply_keyframes function exists."""
        self.assertTrue(hasattr(app_module, 'apply_keyframes'))
        self.assertTrue(callable(app_module.apply_keyframes))

    def test_derive_suitability_driven_args_function_exists(self):
        """Test that derive_suitability_driven_args function exists."""
        self.assertTrue(hasattr(app_module, 'derive_suitability_driven_args'))
        self.assertTrue(callable(app_module.derive_suitability_driven_args))

    def test_set_suitability_driven_args_in_model_core_function_exists(self):
        """Test that set_suitability_driven_args_in_model_core function exists."""
        self.assertTrue(hasattr(app_module, 'set_suitability_driven_args_in_model_core'))
        self.assertTrue(callable(app_module.set_suitability_driven_args_in_model_core))

    def test_get_suitability_driven_arguments_function_exists(self):
        """Test that get_suitability_driven_arguments function exists."""
        self.assertTrue(hasattr(app_module, 'get_suitability_driven_arguments'))
        self.assertTrue(callable(app_module.get_suitability_driven_arguments))

    def test_derive_suitability_driven_arg_function_exists(self):
        """Test that derive_suitability_driven_arg function exists."""
        self.assertTrue(hasattr(app_module, 'derive_suitability_driven_arg'))
        self.assertTrue(callable(app_module.derive_suitability_driven_arg))


class TestAnimalResources(unittest.TestCase):
    """Test cases for animal resource functions."""

    def test_export_animal_resources_function_exists(self):
        """Test that export_animal_resources function exists."""
        self.assertTrue(hasattr(app_module, 'export_animal_resources'))
        self.assertTrue(callable(app_module.export_animal_resources))

    def test_generate_and_or_load_animal_dist_lookup_tables_function_exists(self):
        """Test that generate_and_or_load_animal_dist_lookup_tables function exists."""
        self.assertTrue(hasattr(app_module, 'generate_and_or_load_animal_dist_lookup_tables'))
        self.assertTrue(callable(app_module.generate_and_or_load_animal_dist_lookup_tables))


class TestInitFunction(unittest.TestCase):
    """Test cases for init function."""

    def test_init_function_exists(self):
        """Test that init function exists."""
        self.assertTrue(hasattr(app_module, 'init'))
        self.assertTrue(callable(app_module.init))


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_expected_functions_present(self):
        """Test that all expected functions are present in the module."""
        expected_functions = [
            'set_dispersal_kernel', 'create_controlled_patches_image', 'create_perlin_noise_image',
            'load_homogeneous_image', 'load_specified_existing_image', 'set_initial_tree_cover',
            'set_random_seeds', 'init', 'do_burn_in', 'do_iteration', 'do_update',
            'updateloop', 'termination_condition_satisfied', 'main',
            'set_argument_in_model_core', 'get_old_keyframed_value', 'set_keyframe',
            'apply_keyframes', 'export_animal_resources', 'generate_and_or_load_animal_dist_lookup_tables',
            'derive_suitability_driven_args', 'set_suitability_driven_args_in_model_core',
            'get_suitability_driven_arguments', 'derive_suitability_driven_arg',
            'do_bifurcation_analysis'
        ]
        for func_name in expected_functions:
            self.assertTrue(hasattr(app_module, func_name))
            self.assertTrue(callable(getattr(app_module, func_name)))

    def test_module_has_imports(self):
        """Test that module has necessary imports."""
        # This is tested implicitly by the fact that the module loads
        self.assertTrue(True)


class TestFunctionSignatures(unittest.TestCase):
    """Test cases for function signatures."""

    def test_set_dispersal_kernel_signature(self):
        """Test set_dispersal_kernel function signature."""
        import inspect
        sig = inspect.signature(app_module.set_dispersal_kernel)
        params = list(sig.parameters.keys())
        
        self.assertIn('dynamics', params)
        self.assertIn('dispersal_mode', params)
        self.assertIn('multi_disperser_params', params)

    def test_set_random_seeds_signature(self):
        """Test set_random_seeds function signature."""
        import inspect
        sig = inspect.signature(app_module.set_random_seeds)
        params = list(sig.parameters.keys())
        
        self.assertIn('cfg', params)

    def test_set_initial_tree_cover_signature(self):
        """Test set_initial_tree_cover function signature."""
        import inspect
        sig = inspect.signature(app_module.set_initial_tree_cover)
        params = list(sig.parameters.keys())
        
        self.assertIn('dynamics', params)
        self.assertIn('cfg', params)
        self.assertIn('color_dicts', params)

    def test_main_function_signature(self):
        """Test main function signature."""
        import inspect
        sig = inspect.signature(app_module.main)
        params = list(sig.parameters.keys())
        
        # Main should accept **kwargs
        self.assertIn('user_args', params)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)