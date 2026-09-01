"""Unit test suite for end2end_test module.

This module provides comprehensive tests for the end-to-end testing functionality.
Tests cover the Test class, configuration loading, file management, and integration with the evaluation system.
"""

import unittest
import os
import tempfile
import shutil
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, MagicMock

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import end2end_test as e2e_test


class TestTestClass(unittest.TestCase):
    """Test cases for the Test class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a temporary config file
        self.cfg_file = os.path.join(self.temp_dir, "test_config.json")
        test_config = {
            "random_seed": 42,
            "verbosity": 0,
            "max_timesteps": 10,
            "grid_width": 100,
            "DATA_IN_DIR": self.temp_dir,
            "DATA_OUT_DIR": os.path.join(self.temp_dir, "output")
        }
        with open(self.cfg_file, 'w') as f:
            json.dump(test_config, f)
        
        # Create mock args
        self.args = SimpleNamespace(
            cfg_file=self.cfg_file,
            output_dir=os.path.join(self.temp_dir, "output"),
            mode="evaluate"
        )
        
        # Create a minimal mock config
        self.mock_cfg = SimpleNamespace(
            random_seed=42,
            verbosity=0,
            DATA_IN_DIR=self.temp_dir,
            DATA_OUT_DIR=os.path.join(self.temp_dir, "output")
        )

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    # ==================== INITIALIZATION TESTS ====================

    def test_test_initialization(self):
        """Test that Test class initializes correctly."""
        # Mock the configuration loading functions
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive:
            
            mock_get_defaults.return_value = {}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            # Create a minimal default config
            mock_get_defaults.return_value = {"random_seed": -999}
            mock_load_json.side_effect = lambda x: x
            
            # Mock the check_cfg to avoid the seed assertion
            with patch.object(e2e_test.Test, 'check_cfg'):
                test = e2e_test.Test(self.args)
                
                self.assertEqual(test.output_dir, self.args.output_dir)
                self.assertEqual(test.mode, self.args.mode)
                self.assertEqual(test.name, Path(self.args.cfg_file).stem)

    def test_test_name_extraction(self):
        """Test that test name is correctly extracted from config file path."""
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            expected_name = Path(self.args.cfg_file).stem
            self.assertEqual(test.name, expected_name)

    # ==================== CONFIGURATION LOADING TESTS ====================

    def test_load_casefile_success(self):
        """Test load_casefile successfully loads a JSON config file."""
        result = e2e_test.Test.load_casefile(self.cfg_file)
        
        self.assertIsInstance(result, SimpleNamespace)
        self.assertEqual(result.random_seed, 42)
        self.assertEqual(result.max_timesteps, 10)

    def test_load_casefile_file_not_found(self):
        """Test load_casefile raises when file doesn't exist."""
        with self.assertRaises(FileNotFoundError):
            e2e_test.Test.load_casefile("/nonexistent/file.json")

    def test_load_casefile_invalid_json(self):
        """Test load_casefile raises when JSON is invalid."""
        invalid_json_file = os.path.join(self.temp_dir, "invalid.json")
        with open(invalid_json_file, 'w') as f:
            f.write("invalid json content")
        
        with self.assertRaises(json.JSONDecodeError):
            e2e_test.Test.load_casefile(invalid_json_file)

    def test_apply_case_cfg_basic(self):
        """Test apply_case_cfg applies case config to default config."""
        default_cfg = SimpleNamespace(a=1, b=2, c=3)
        case_cfg = SimpleNamespace(b=10, c=20, d=40)
        
        result = e2e_test.Test.apply_case_cfg(case_cfg, default_cfg)
        
        self.assertEqual(result.a, 1)  # From default
        self.assertEqual(result.b, 10)  # Overridden by case
        self.assertEqual(result.c, 20)  # Overridden by case
        self.assertEqual(result.d, 40)  # Added from case

    def test_apply_case_cfg_preserves_defaults(self):
        """Test apply_case_cfg preserves default values for unspecified keys."""
        default_cfg = SimpleNamespace(a=1, b=2, c=3)
        case_cfg = SimpleNamespace(b=10)
        
        result = e2e_test.Test.apply_case_cfg(case_cfg, default_cfg)
        
        self.assertEqual(result.a, 1)  # Preserved from default
        self.assertEqual(result.b, 10)  # Overridden
        self.assertEqual(result.c, 3)  # Preserved from default

    def test_apply_case_cfg_no_overlap(self):
        """Test apply_case_cfg with no overlapping keys."""
        default_cfg = SimpleNamespace(a=1, b=2)
        case_cfg = SimpleNamespace(c=3, d=4)
        
        result = e2e_test.Test.apply_case_cfg(case_cfg, default_cfg)
        
        self.assertEqual(result.a, 1)
        self.assertEqual(result.b, 2)
        self.assertEqual(result.c, 3)
        self.assertEqual(result.d, 4)

    def test_apply_case_cfg_empty_case_cfg(self):
        """Test apply_case_cfg with empty case config."""
        default_cfg = SimpleNamespace(a=1, b=2)
        case_cfg = SimpleNamespace()
        
        result = e2e_test.Test.apply_case_cfg(case_cfg, default_cfg)
        
        self.assertEqual(result.a, 1)
        self.assertEqual(result.b, 2)

    # ==================== CONFIGURATION VALIDATION TESTS ====================

    def test_check_cfg_valid_seed(self):
        """Test check_cfg passes with valid seed."""
        cfg = SimpleNamespace(random_seed=42)
        
        # Should not raise
        try:
            e2e_test.Test.check_cfg(cfg)
            success = True
        except AssertionError:
            success = False
        
        self.assertTrue(success)

    def test_check_cfg_invalid_seed(self):
        """Test check_cfg fails with invalid seed (-999)."""
        cfg = SimpleNamespace(random_seed=-999)
        
        with self.assertRaises(AssertionError) as context:
            e2e_test.Test.check_cfg(cfg)
        
        self.assertIn("must specify a random seed", str(context.exception))

    # ==================== OUTPUT DIRECTORY OVERRITE TESTS ====================

    def test_overwrite_output_dir(self):
        """Test overwrite_output_dir sets the DATA_OUT_DIR correctly."""
        cfg = SimpleNamespace(
            DATA_OUT_DIR="/original/path",
            OTHER_DIR="/original/other"
        )
        new_output_dir = "/new/output/path"
        
        with patch('end2end_test.derive_output_dirs') as mock_derive:
            mock_derive.side_effect = lambda x: x
            result = e2e_test.Test.overwrite_output_dir(cfg, new_output_dir)
            
            self.assertEqual(result.DATA_OUT_DIR, new_output_dir)

    def test_overwrite_output_dir_preserves_other_attributes(self):
        """Test overwrite_output_dir preserves other configuration attributes."""
        cfg = SimpleNamespace(
            DATA_OUT_DIR="/original/path",
            OTHER_DIR="/original/other",
            some_setting=123
        )
        new_output_dir = "/new/output/path"
        
        with patch('end2end_test.derive_output_dirs') as mock_derive:
            mock_derive.side_effect = lambda x: x
            result = e2e_test.Test.overwrite_output_dir(cfg, new_output_dir)
            
            self.assertEqual(result.DATA_OUT_DIR, new_output_dir)
            self.assertEqual(result.OTHER_DIR, "/original/other")
            self.assertEqual(result.some_setting, 123)

    # ==================== FILE MANAGEMENT TESTS ====================

    def test_copy_cover_image(self):
        """Test copy_cover_image copies the cover image to output directory."""
        # Create a mock Test instance
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            
            # Create a source cover image
            src_cover = os.path.join(self.temp_dir, "cover.png")
            with open(src_cover, 'wb') as f:
                f.write(b"fake image data")
            test.cfg.cover_img_path = src_cover
            
            # Copy the cover image
            test.copy_cover_image()
            
            # Check that file was copied
            dst_cover = os.path.join(test.output_dir, "cover.png")
            self.assertTrue(os.path.exists(dst_cover))

    def test_rename_state_data_file(self):
        """Test rename_state_data_file renames the state data file correctly."""
        # Create a mock Test instance
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            
            # Create a state data file
            state_dir = os.path.join(test.output_dir, "state_data")
            os.makedirs(state_dir, exist_ok=True)
            original_file = os.path.join(state_dir, "Simulation_run_001.csv")
            with open(original_file, 'w') as f:
                f.write("fake csv data")
            
            # Rename the file
            test.rename_state_data_file()
            
            # Check that file was renamed
            expected_file = os.path.join(state_dir, "Simulation.csv")
            self.assertTrue(os.path.exists(expected_file))
            self.assertFalse(os.path.exists(original_file))

    def test_rename_state_data_file_no_matching_files(self):
        """Test rename_state_data_file when no matching files are found."""
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            
            # Create state_data directory but no matching files
            state_dir = os.path.join(test.output_dir, "state_data")
            os.makedirs(state_dir, exist_ok=True)
            
            # Should raise IndexError when glob returns empty list
            with self.assertRaises(IndexError):
                test.rename_state_data_file()

    def test_do_post_test_filemanagement(self):
        """Test do_post_test_filemanagement executes both file operations."""
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'), \
             patch.object(e2e_test.Test, 'copy_cover_image'), \
             patch.object(e2e_test.Test, 'rename_state_data_file'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            test.do_post_test_filemanagement()
            
            # Check that both methods were called
            test.copy_cover_image.assert_called_once()
            test.rename_state_data_file.assert_called_once()

    # ==================== RUN METHOD TESTS ====================

    def test_run_method_integration(self):
        """Test the run method calls app.main and post-test file management."""
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'), \
             patch('end2end_test.app.main') as mock_app_main, \
             patch.object(e2e_test.Test, 'do_post_test_filemanagement'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            mock_app_main.return_value = (None, self.mock_cfg)
            
            test = e2e_test.Test(self.args)
            test.run()
            
            # Check that app.main was called
            mock_app_main.assert_called_once()
            
            # Check that post-test file management was called
            test.do_post_test_filemanagement.assert_called_once()

    # ==================== REPLACE BENCHMARK FUNCTION TESTS ====================

    def test_replace_benchmark_function(self):
        """Test the replace_benchmark function."""
        # Create temporary benchmark directory
        benchmark_dir = os.path.join(self.temp_dir, "benchmark")
        os.makedirs(benchmark_dir, exist_ok=True)
        
        # Create output directory with some content
        output_dir = os.path.join(self.temp_dir, "output")
        output_subdir = os.path.join(output_dir, "state_data")
        os.makedirs(output_subdir, exist_ok=True)
        
        test_file = os.path.join(output_subdir, "test.csv")
        with open(test_file, 'w') as f:
            f.write("test data")
        
        # Create mock objects
        mock_test = MagicMock()
        mock_test.output_dir = output_dir
        
        mock_evaluator = MagicMock()
        mock_evaluator.benchmark_dir = benchmark_dir
        
        mock_args = MagicMock()
        
        # Mock the file handling functions
        with patch('end2end_test.io.remove_dir_contents'), \
             patch('end2end_test.io.copy_tree_if_needed') as mock_copy:
            
            e2e_test.replace_benchmark(mock_test, mock_evaluator, mock_args)
            
            # Check that functions were called
            mock_copy.assert_called_once_with(mock_test.output_dir, mock_evaluator.benchmark_dir)

    # ==================== MAIN FUNCTION TESTS ====================

    def test_main_function_generate_benchmark_mode(self):
        """Test main function in generate_benchmark mode."""
        with patch('end2end_test.ArgumentParser') as mock_parser, \
             patch('end2end_test.Test') as mock_test_class, \
             patch('end2end_test.Evaluator') as mock_evaluator_class, \
             patch('end2end_test.replace_benchmark') as mock_replace:
            
            # Setup mocks
            mock_args = MagicMock()
            mock_args.cfg_file = self.cfg_file
            mock_args.output_dir = os.path.join(self.temp_dir, "output")
            mock_args.mode = "generate_benchmark"
            
            mock_parser_instance = MagicMock()
            mock_parser_instance.parse_args.return_value = mock_args
            mock_parser.return_value = mock_parser_instance
            
            mock_test_instance = MagicMock()
            mock_test_class.return_value = mock_test_instance
            
            mock_evaluator_instance = MagicMock()
            mock_evaluator_class.return_value = mock_evaluator_instance
            
            # Call main
            e2e_test.main()
            
            # Check that replace_benchmark was called
            mock_replace.assert_called_once()

    def test_main_function_evaluate_mode(self):
        """Test main function in evaluate mode."""
        with patch('end2end_test.ArgumentParser') as mock_parser, \
             patch('end2end_test.Test') as mock_test_class, \
             patch('end2end_test.Evaluator') as mock_evaluator_class, \
             patch('end2end_test.replace_benchmark'), \
             patch.object(e2e_test.Evaluator, 'evaluate') as mock_evaluate:
            
            # Setup mocks
            mock_args = MagicMock()
            mock_args.cfg_file = self.cfg_file
            mock_args.output_dir = os.path.join(self.temp_dir, "output")
            mock_args.mode = "evaluate"
            
            mock_parser_instance = MagicMock()
            mock_parser_instance.parse_args.return_value = mock_args
            mock_parser.return_value = mock_parser_instance
            
            mock_test_instance = MagicMock()
            mock_test_class.return_value = mock_test_instance
            
            mock_evaluator_instance = MagicMock()
            mock_evaluator_class.return_value = mock_evaluator_instance
            
            # Call main
            e2e_test.main()
            
            # Check that evaluate was called
            mock_evaluate.assert_called_once_with(mock_test_instance)

    def test_main_function_unknown_mode(self):
        """Test main function with unknown mode raises ValueError."""
        with patch('end2end_test.ArgumentParser') as mock_parser:
            # Setup mocks
            mock_args = MagicMock()
            mock_args.cfg_file = self.cfg_file
            mock_args.output_dir = os.path.join(self.temp_dir, "output")
            mock_args.mode = "unknown_mode"
            
            mock_parser_instance = MagicMock()
            mock_parser_instance.parse_args.return_value = mock_args
            mock_parser.return_value = mock_parser_instance
            
            with self.assertRaises(ValueError) as context:
                e2e_test.main()
            
            self.assertIn("Unknown mode", str(context.exception))

    # ==================== EDGE CASE TESTS ====================

    def test_load_casefile_empty_file(self):
        """Test load_casefile with empty JSON file."""
        empty_json_file = os.path.join(self.temp_dir, "empty.json")
        with open(empty_json_file, 'w') as f:
            f.write("")
        
        with self.assertRaises(json.JSONDecodeError):
            e2e_test.Test.load_casefile(empty_json_file)

    def test_apply_case_cfg_none_values(self):
        """Test apply_case_cfg with None values."""
        default_cfg = SimpleNamespace(a=1, b=2)
        case_cfg = SimpleNamespace(a=None)
        
        result = e2e_test.Test.apply_case_cfg(case_cfg, default_cfg)
        
        self.assertIsNone(result.a)
        self.assertEqual(result.b, 2)

    def test_overwrite_output_dir_empty_string(self):
        """Test overwrite_output_dir with empty string."""
        cfg = SimpleNamespace(DATA_OUT_DIR="/original/path")
        
        with patch('end2end_test.derive_output_dirs') as mock_derive:
            mock_derive.side_effect = lambda x: x
            result = e2e_test.Test.overwrite_output_dir(cfg, "")
            
            self.assertEqual(result.DATA_OUT_DIR, "")

    def test_copy_cover_image_file_not_exists(self):
        """Test copy_cover_image when source file doesn't exist."""
        with patch('end2end_test.get_all_defaults') as mock_get_defaults, \
             patch('end2end_test.load_json_strings_if_any') as mock_load_json, \
             patch('end2end_test.derive_output_dirs') as mock_derive, \
             patch.object(e2e_test.Test, 'check_cfg'):
            
            mock_get_defaults.return_value = {"random_seed": 42}
            mock_load_json.side_effect = lambda x: x
            mock_derive.side_effect = lambda x: x
            
            test = e2e_test.Test(self.args)
            test.cfg.cover_img_path = "/nonexistent/cover.png"
            
            with self.assertRaises(FileNotFoundError):
                test.copy_cover_image()


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)