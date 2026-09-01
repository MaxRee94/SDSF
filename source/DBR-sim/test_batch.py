"""Unit test suite for batch module.

This module provides comprehensive tests for the batch processing functionality.
Tests cover the Jobs class, argument parsing, job generation, and batch execution.
"""

import unittest
import os
import sys
import tempfile
import shutil
import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call
import inspect

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock complex dependencies
sys.modules['config'] = MagicMock()
sys.modules['app'] = MagicMock()
sys.modules['file_handling'] = MagicMock()
sys.modules['visualization'] = MagicMock()
sys.modules['helpers'] = MagicMock()

import batch as batch_module


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.DATA_IN_DIR = tempfile.mkdtemp()
        self.csv_parent_dir = tempfile.mkdtemp()


class TestJobsClass(unittest.TestCase):
    """Test cases for the Jobs class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = MockConfig()

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        shutil.rmtree(self.cfg.DATA_IN_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.csv_parent_dir, ignore_errors=True)

    def test_jobs_class_exists(self):
        """Test that Jobs class exists."""
        self.assertTrue(hasattr(batch_module, 'Jobs'))
        self.assertTrue(callable(batch_module.Jobs))

    def test_jobs_initialization(self):
        """Test Jobs class initialization."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "runs": 1}
            arguments = {}
            
            jobs = batch_module.Jobs(runs=1, arguments=arguments, **batch_args)
            
            self.assertEqual(jobs.n_runs, 1)
            self.assertEqual(jobs.arg_changes, arguments)
            self.assertEqual(jobs.batch_type, "regular")

    def test_jobs_initialization_with_defaults(self):
        """Test Jobs initialization with defaults."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {"param1": "default1", "param2": "default2"}
            
            batch_args = {"type": "regular", "runs": 2}
            arguments = {"param1": {"values": [1, 2, 3]}}
            
            jobs = batch_module.Jobs(runs=2, arguments=arguments, **batch_args)
            
            # Should have created defaults
            self.assertIsNotNone(jobs.defaults)
            self.assertIn("param1", jobs.defaults)

    def test_jobs_sampling_mode_default(self):
        """Test Jobs sampling_mode default."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular"}
            jobs = batch_module.Jobs(**batch_args)
            
            self.assertEqual(jobs.sampling_mode, "regular")

    def test_jobs_first_keyframe_with_intersim_variation_default(self):
        """Test Jobs first_keyframe_with_intersim_variation default."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular"}
            jobs = batch_module.Jobs(**batch_args)
            
            self.assertEqual(jobs.first_keyframe_with_intersim_variation, {})

    def test_jobs_parse_method_exists(self):
        """Test that Jobs has parse method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'parse'))
        self.assertTrue(callable(batch_module.Jobs.parse))

    def test_jobs_get_defaults_method_exists(self):
        """Test that Jobs has get_defaults method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'get_defaults'))
        self.assertTrue(callable(batch_module.Jobs.get_defaults))

    def test_jobs_attach_control_variables_method_exists(self):
        """Test that Jobs has attach_control_variables method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'attach_control_variables'))
        self.assertTrue(callable(batch_module.Jobs.attach_control_variables))

    def test_jobs_attach_sim_name_method_exists(self):
        """Test that Jobs has attach_sim_name method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'attach_sim_name'))
        self.assertTrue(callable(batch_module.Jobs.attach_sim_name))

    def test_jobs_get_key_idx_method_exists(self):
        """Test that Jobs has get_key_idx method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'get_key_idx'))
        self.assertTrue(callable(batch_module.Jobs.get_key_idx))

    def test_jobs_get_vec_method_exists(self):
        """Test that Jobs has get_vec method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'get_vec'))
        self.assertTrue(callable(batch_module.Jobs.get_vec))

    def test_jobs_get_keyframe_vec_method_exists(self):
        """Test that Jobs has get_keyframe_vec method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'get_keyframe_vec'))
        self.assertTrue(callable(batch_module.Jobs.get_keyframe_vec))

    def test_jobs_update_default_values_method_exists(self):
        """Test that Jobs has update_default_values method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'update_default_values'))
        self.assertTrue(callable(batch_module.Jobs.update_default_values))

    def test_jobs_generate_and_apply_random_seeds_method_exists(self):
        """Test that Jobs has generate_and_apply_random_seeds method."""
        self.assertTrue(hasattr(batch_module.Jobs, 'generate_and_apply_random_seeds'))
        self.assertTrue(callable(batch_module.Jobs.generate_and_apply_random_seeds))


class TestJobsDefaults(unittest.TestCase):
    """Test cases for Jobs defaults functionality."""

    def test_get_defaults_method(self):
        """Test get_defaults method with mocked config."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {"param1": "default1"}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            arguments = {}
            
            jobs = batch_module.Jobs(runs=1, arguments=arguments, **batch_args)
            defaults = jobs.get_defaults(batch_args, arguments)
            
            # Should have headless=True for batch runs
            self.assertEqual(defaults["headless"], True)
            # Should have verbosity=-1 for batch runs
            self.assertEqual(defaults["verbosity"], -1)
            # Should have EXPORT_DIR set
            self.assertEqual(defaults["EXPORT_DIR"], "/tmp")


class TestJobsHelperMethods(unittest.TestCase):
    """Test cases for Jobs helper methods."""

    def test_get_key_idx_method(self):
        """Test get_key_idx method."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {"param1": "val1", "param2": "val2"}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            arguments = {}
            
            jobs = batch_module.Jobs(runs=1, arguments=arguments, **batch_args)
            
            # Test getting index of existing key
            idx = jobs.get_key_idx("param1")
            self.assertEqual(idx, 0)


class TestJobsRandomSeeds(unittest.TestCase):
    """Test cases for Jobs random seed generation."""

    def test_generate_and_apply_random_seeds_random_seed(self):
        """Test generate_and_apply_random_seeds with random_seed=-999."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(
                runs=1, 
                rng=np.random.default_rng(42),
                firefreq_rng=np.random.default_rng(43),
                **batch_args
            )
            
            # Create a job with random_seed=-999
            job = SimpleNamespace(random_seed=-999, firefreq_random_seed=-999)
            
            result_job = jobs.generate_and_apply_random_seeds(job)
            
            # Should have generated new seeds
            self.assertNotEqual(result_job.random_seed, -999)
            self.assertNotEqual(result_job.firefreq_random_seed, -999)
            
            # Original job should be unchanged
            self.assertEqual(job.random_seed, -999)
            self.assertEqual(job.firefreq_random_seed, -999)

    def test_generate_and_apply_random_seeds_specific_seeds(self):
        """Test generate_and_apply_random_seeds with specific seeds."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(
                runs=1, 
                rng=np.random.default_rng(42),
                firefreq_rng=np.random.default_rng(43),
                **batch_args
            )
            
            # Create a job with specific seeds
            job = SimpleNamespace(random_seed=123, firefreq_random_seed=456)
            
            result_job = jobs.generate_and_apply_random_seeds(job)
            
            # Should keep the specified seeds
            self.assertEqual(result_job.random_seed, 123)
            self.assertEqual(result_job.firefreq_random_seed, 456)


class TestJobsVectorGeneration(unittest.TestCase):
    """Test cases for Jobs vector generation functionality."""

    def test_get_vec_simple_values(self):
        """Test get_vec with simple values."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            arg_cfg = {"values": [1, 2, 3, 4]}
            result = jobs.get_vec(arg_cfg)
            
            self.assertEqual(result, [1, 2, 3, 4])

    def test_get_vec_range_values(self):
        """Test get_vec with range values."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            arg_cfg = {"min": 1, "max": 10, "steps": 5}
            result = jobs.get_vec(arg_cfg)
            
            # Should generate 5 values from 1 to 10
            self.assertEqual(len(result), 5)
            self.assertEqual(result[0], 1)
            self.assertEqual(result[-1], 10)


class TestMainFunction(unittest.TestCase):
    """Test cases for main function."""

    def test_main_function_exists(self):
        """Test that main function exists."""
        self.assertTrue(hasattr(batch_module, 'main'))
        self.assertTrue(callable(batch_module.main))


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = ['Jobs']
        for class_name in expected_classes:
            self.assertTrue(hasattr(batch_module, class_name))
            self.assertTrue(callable(getattr(batch_module, class_name)))

    def test_expected_functions_present(self):
        """Test that all expected functions are present in the module."""
        expected_functions = ['main']
        for func_name in expected_functions:
            self.assertTrue(hasattr(batch_module, func_name))
            self.assertTrue(callable(getattr(batch_module, func_name)))


class TestJobsClassAttributes(unittest.TestCase):
    """Test cases for Jobs class attributes."""

    def test_jobs_has_jobs_attribute(self):
        """Test that Jobs has jobs attribute."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            self.assertTrue(hasattr(jobs, 'jobs'))

    def test_jobs_has_job_indices_attribute(self):
        """Test that Jobs has job_indices attribute."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            self.assertTrue(hasattr(jobs, 'job_indices'))

    def test_jobs_has_default_values_attribute(self):
        """Test that Jobs has default_values attribute."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            self.assertTrue(hasattr(jobs, 'default_values'))


class TestJobsEdgeCases(unittest.TestCase):
    """Test cases for Jobs edge cases."""

    def test_jobs_with_empty_arguments(self):
        """Test Jobs with empty arguments."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=1, arguments={}, **batch_args)
            
            # Should not raise an error
            self.assertIsNotNone(jobs)

    def test_jobs_with_none_runs(self):
        """Test Jobs with None runs."""
        with patch('batch.config.get_all_defaults') as mock_get_defaults:
            mock_get_defaults.return_value = {}
            
            batch_args = {"type": "regular", "csv_parent_dir": "/tmp"}
            jobs = batch_module.Jobs(runs=None, arguments={}, **batch_args)
            
            # Should handle None runs
            self.assertIsNotNone(jobs)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)