"""Unit test suite for end2end_test_administrator module.

This module provides comprehensive tests for the end-to-end test administration functionality.
Tests cover the OutputTests class, test file management, subprocess execution, and reporting.
"""

import unittest
import os
import sys
import tempfile
import shutil
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call
import inspect

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock complex dependencies
sys.modules['visualization'] = MagicMock()
sys.modules['file_handling'] = MagicMock()
sys.modules['helpers'] = MagicMock()
sys.modules['config'] = MagicMock()
sys.modules['cv2'] = MagicMock()
sys.modules['numpy'] = MagicMock()
sys.modules['matplotlib'] = MagicMock()
sys.modules['matplotlib.pyplot'] = MagicMock()
sys.modules['scipy'] = MagicMock()
sys.modules['scipy.stats'] = MagicMock()
sys.modules['app'] = MagicMock()

import end2end_test_administrator as e2e_admin


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.END2END_TESTCASE_DIR = tempfile.mkdtemp()
        self.END2END_TEST_OUTPUT_DIR = tempfile.mkdtemp()


class TestOutputTestsClass(unittest.TestCase):
    """Test cases for OutputTests class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = MockConfig()
        
        # Create a mock test case file
        self.test_case_dir = self.cfg.END2END_TESTCASE_DIR
        os.makedirs(self.test_case_dir, exist_ok=True)
        
        self.test_file = os.path.join(self.test_case_dir, "test_case.json")
        with open(self.test_file, 'w') as f:
            json.dump({"param": "value"}, f)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        shutil.rmtree(self.cfg.END2END_TESTCASE_DIR, ignore_errors=True)
        shutil.rmtree(self.cfg.END2END_TEST_OUTPUT_DIR, ignore_errors=True)

    def test_output_tests_class_exists(self):
        """Test that OutputTests class exists."""
        self.assertTrue(hasattr(e2e_admin, 'OutputTests'))
        self.assertTrue(callable(e2e_admin.OutputTests))

    def test_output_tests_initialization_evaluate_mode(self):
        """Test OutputTests initialization in evaluate mode."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists'), \
             patch('end2end_test_administrator.glob.glob') as mock_glob:
            
            mock_glob.return_value = [self.test_file]
            
            tests = e2e_admin.OutputTests("evaluate")
            
            self.assertEqual(tests.mode, "evaluate")
            self.assertEqual(tests.test_files, [self.test_file])

    def test_output_tests_initialization_reset_mode(self):
        """Test OutputTests initialization in reset mode."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists'), \
             patch('end2end_test_administrator.glob.glob') as mock_glob:
            
            mock_glob.return_value = [self.test_file]
            
            tests = e2e_admin.OutputTests("end2end_reset")
            
            # Should be converted to generate_benchmark mode
            self.assertEqual(tests.mode, "generate_benchmark")

    def test_output_tests_initialization_benchmark_mode(self):
        """Test OutputTests initialization in benchmark mode."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists'), \
             patch('end2end_test_administrator.glob.glob') as mock_glob:
            
            mock_glob.return_value = [self.test_file]
            
            tests = e2e_admin.OutputTests("generate_benchmark")
            
            self.assertEqual(tests.mode, "generate_benchmark")


class TestGetTestFiles(unittest.TestCase):
    """Test cases for get_test_files method."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_case_dir = os.path.join(self.temp_dir, "test_cases")
        os.makedirs(self.test_case_dir, exist_ok=True)
        
        # Create some test files
        for i in range(3):
            with open(os.path.join(self.test_case_dir, f"test_{i}.json"), 'w') as f:
                json.dump({"test": i}, f)
        
        self.cfg = MockConfig()
        self.cfg.END2END_TESTCASE_DIR = self.test_case_dir

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_get_test_files_returns_json_files(self):
        """Test that get_test_files returns JSON files from test case directory."""
        with patch('end2end_test_administrator.cfg', self.cfg):
            result = e2e_admin.OutputTests.get_test_files()
            
            self.assertEqual(len(result), 3)
            for file_path in result:
                self.assertTrue(file_path.endswith('.json'))
                self.assertIn('test_cases', file_path)


class TestCreateMainOutputDir(unittest.TestCase):
    """Test cases for create_main_output_dir method."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = MockConfig()
        self.cfg.END2END_TEST_OUTPUT_DIR = self.temp_dir

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_create_main_output_dir_creates_directory(self):
        """Test that create_main_output_dir creates a timestamped directory."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists') as mock_create, \
             patch('end2end_test_administrator.time.strftime') as mock_time:
            
            mock_time.return_value = "2023_01_01-12,00,00"
            
            tests = e2e_admin.OutputTests("evaluate")
            output_dir = tests.create_main_output_dir()
            
            expected_dir = os.path.join(self.temp_dir, "2023_01_01-12,00,00")
            self.assertEqual(output_dir, expected_dir)
            mock_create.assert_called_once_with(expected_dir)


class TestRunTestProcess(unittest.TestCase):
    """Test cases for run_test_process method."""

    def test_run_test_process_function_exists(self):
        """Test that run_test_process method exists."""
        self.assertTrue(hasattr(e2e_admin.OutputTests, 'run_test_process'))
        self.assertTrue(callable(e2e_admin.OutputTests.run_test_process))

    def test_run_test_process_uses_subprocess(self):
        """Test that run_test_process uses subprocess to run tests."""
        source = inspect.getsource(e2e_admin.OutputTests.run_test_process)
        
        self.assertIn('subprocess.run', source)
        self.assertIn('sys.executable', source)
        self.assertIn('end2end_test.py', source)

    def test_run_test_process_returns_output(self):
        """Test that run_test_process returns stdout, stderr, and returncode."""
        source = inspect.getsource(e2e_admin.OutputTests.run_test_process)
        
        self.assertIn('p.stdout', source)
        self.assertIn('p.stderr', source)
        self.assertIn('p.returncode', source)


class TestWriteCmdlineOutput(unittest.TestCase):
    """Test cases for write_cmdline_output_to_file method."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_write_cmdline_output_to_file_function_exists(self):
        """Test that write_cmdline_output_to_file method exists."""
        self.assertTrue(hasattr(e2e_admin.OutputTests, 'write_cmdline_output_to_file'))
        self.assertTrue(callable(e2e_admin.OutputTests.write_cmdline_output_to_file))

    def test_write_cmdline_output_to_file_creates_file(self):
        """Test that write_cmdline_output_to_file creates output file."""
        output_content = "Test output content"
        output_dir = self.temp_dir
        filename = "test_output.txt"
        
        e2e_admin.OutputTests.write_cmdline_output_to_file(
            output_content, output_dir, filename
        )
        
        expected_path = os.path.join(output_dir, filename)
        self.assertTrue(os.path.exists(expected_path))
        
        with open(expected_path, 'r') as f:
            content = f.read()
        self.assertEqual(content, output_content)


class TestPrepareTest(unittest.TestCase):
    """Test cases for prepare_test method."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.temp_dir, "test.json")
        with open(self.test_file, 'w') as f:
            json.dump({"test": "data"}, f)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_prepare_test_creates_output_dirs(self):
        """Test that prepare_test creates output directories."""
        with patch('end2end_test_administrator.io.create_output_dirs') as mock_create_dirs:
            tests = e2e_admin.OutputTests("evaluate")
            tests.main_output_dir = self.temp_dir
            
            test_name, output_dir, message = tests.prepare_test(0, self.test_file)
            
            expected_name = "test"
            expected_output_dir = os.path.join(self.temp_dir, expected_name)
            
            self.assertEqual(test_name, expected_name)
            self.assertEqual(output_dir, expected_output_dir)
            self.assertIn("Running", message)
            mock_create_dirs.assert_called_once_with(expected_output_dir)

    def test_prepare_test_message_format(self):
        """Test that prepare_test creates proper message."""
        with patch('end2end_test_administrator.io.create_output_dirs'):
            tests = e2e_admin.OutputTests("evaluate")
            tests.main_output_dir = self.temp_dir
            tests.test_files = [self.test_file]
            
            test_name, output_dir, message = tests.prepare_test(0, self.test_file)
            
            self.assertIn("Running", message)
            self.assertIn("1/1", message)
            self.assertIn("test", message)


class TestProduceTestReport(unittest.TestCase):
    """Test cases for produce_test_report method."""

    def test_produce_test_report_success_evaluate_mode(self):
        """Test produce_test_report for successful test in evaluate mode."""
        tests = e2e_admin.OutputTests("evaluate")
        
        # This method just prints, so we can't easily test the output
        # But we can test that it doesn't crash
        try:
            tests.produce_test_report("Running test", True)
            success = True
        except Exception:
            success = False
        
        self.assertTrue(success)

    def test_produce_test_report_failure_evaluate_mode(self):
        """Test produce_test_report for failed test in evaluate mode."""
        tests = e2e_admin.OutputTests("evaluate")
        
        try:
            tests.produce_test_report("Running test", False)
            success = True
        except Exception:
            success = False
        
        self.assertTrue(success)

    def test_produce_test_report_benchmark_mode(self):
        """Test produce_test_report in benchmark mode."""
        tests = e2e_admin.OutputTests("generate_benchmark")
        
        try:
            tests.produce_test_report("Generating benchmark", True)
            success = True
        except Exception:
            success = False
        
        self.assertTrue(success)


class TestRunAll(unittest.TestCase):
    """Test cases for run_all method."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_case_dir = os.path.join(self.temp_dir, "test_cases")
        os.makedirs(self.test_case_dir, exist_ok=True)
        
        # Create a test file
        self.test_file = os.path.join(self.test_case_dir, "test.json")
        with open(self.test_file, 'w') as f:
            json.dump({"test": "data"}, f)
        
        self.cfg = MockConfig()
        self.cfg.END2END_TESTCASE_DIR = self.test_case_dir
        self.cfg.END2END_TEST_OUTPUT_DIR = os.path.join(self.temp_dir, "output")

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_run_all_returns_counts(self):
        """Test that run_all returns failure and total counts."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists'), \
             patch('end2end_test_administrator.io.create_output_dirs'), \
             patch('end2end_test_administrator.glob.glob') as mock_glob, \
             patch('end2end_test_administrator.OutputTests.run_test') as mock_run_test:
            
            mock_glob.return_value = [self.test_file]
            mock_run_test.return_value = True
            
            tests = e2e_admin.OutputTests("evaluate")
            n_failed, n_total = tests.run_all()
            
            self.assertEqual(n_total, 1)
            self.assertEqual(n_failed, 0)  # No failures

    def test_run_all_counts_failures(self):
        """Test that run_all correctly counts failures."""
        with patch('end2end_test_administrator.cfg', self.cfg), \
             patch('end2end_test_administrator.io.create_directory_if_not_exists'), \
             patch('end2end_test_administrator.io.create_output_dirs'), \
             patch('end2end_test_administrator.glob.glob') as mock_glob, \
             patch('end2end_test_administrator.OutputTests.run_test') as mock_run_test:
            
            mock_glob.return_value = [self.test_file]
            mock_run_test.return_value = False  # Simulate failure
            
            tests = e2e_admin.OutputTests("evaluate")
            n_failed, n_total = tests.run_all()
            
            self.assertEqual(n_total, 1)
            self.assertEqual(n_failed, 1)  # One failure


class TestMoveSubfoldersToOld(unittest.TestCase):
    """Test cases for move_subfolders_to_old function."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create directory structure
        self.root_dir = os.path.join(self.temp_dir, "root")
        os.makedirs(self.root_dir, exist_ok=True)
        
        # Create old directory
        self.old_dir = os.path.join(self.root_dir, "old")
        os.makedirs(self.old_dir, exist_ok=True)
        
        # Create some subdirectories to move
        for i in range(3):
            subdir = os.path.join(self.root_dir, f"run_{i:03d}")
            os.makedirs(subdir, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_move_subfolders_to_old_function_exists(self):
        """Test that move_subfolders_to_old function exists."""
        self.assertTrue(hasattr(e2e_admin, 'move_subfolders_to_old'))
        self.assertTrue(callable(e2e_admin.move_subfolders_to_old))

    def test_move_subfolders_to_old_moves_directories(self):
        """Test that move_subfolders_to_old moves subdirectories to old directory."""
        e2e_admin.move_subfolders_to_old(self.root_dir, "old")
        
        # Check that subdirectories were moved
        for i in range(3):
            old_subdir = os.path.join(self.old_dir, f"run_{i:03d}")
            self.assertTrue(os.path.exists(old_subdir))
            
            # Check that original directories are gone
            original_subdir = os.path.join(self.root_dir, f"run_{i:03d}")
            self.assertFalse(os.path.exists(original_subdir))

    def test_move_subfolders_to_old_preserves_old_dir(self):
        """Test that move_subfolders_to_old preserves the old directory itself."""
        e2e_admin.move_subfolders_to_old(self.root_dir, "old")
        
        # Old directory should still exist
        self.assertTrue(os.path.exists(self.old_dir))
        self.assertTrue(os.path.isdir(self.old_dir))

    def test_move_subfolders_to_old_nonexistent_old_dir(self):
        """Test that move_subfolders_to_old raises when old directory doesn't exist."""
        with self.assertRaises(FileNotFoundError):
            e2e_admin.move_subfolders_to_old(self.root_dir, "nonexistent")

    def test_move_subfolders_to_old_skips_non_directories(self):
        """Test that move_subfolders_to_old skips non-directory files."""
        # Create a file in root directory
        test_file = os.path.join(self.root_dir, "test_file.txt")
        with open(test_file, 'w') as f:
            f.write("test content")
        
        e2e_admin.move_subfolders_to_old(self.root_dir, "old")
        
        # File should still exist in root directory
        self.assertTrue(os.path.exists(test_file))


class TestPrepareFunction(unittest.TestCase):
    """Test cases for prepare function."""

    def test_prepare_function_exists(self):
        """Test that prepare function exists."""
        self.assertTrue(hasattr(e2e_admin, 'prepare'))
        self.assertTrue(callable(e2e_admin.prepare))

    def test_prepare_calls_move_subfolders_to_old(self):
        """Test that prepare calls move_subfolders_to_old."""
        with patch('end2end_test_administrator.move_subfolders_to_old') as mock_move, \
             patch('end2end_test_administrator.cfg') as mock_cfg:
            
            mock_cfg.END2END_TEST_OUTPUT_DIR = "/tmp/output"
            
            e2e_admin.prepare(mock_cfg)
            
            mock_move.assert_called_once_with("/tmp/output")


class TestMainFunction(unittest.TestCase):
    """Test cases for main function."""

    def test_main_function_exists(self):
        """Test that main function exists."""
        self.assertTrue(hasattr(e2e_admin, 'main'))
        self.assertTrue(callable(e2e_admin.main))

    def test_main_function_evaluate_mode(self):
        """Test main function in evaluate mode."""
        with patch('end2end_test_administrator.prepare'), \
             patch('end2end_test_administrator.OutputTests') as mock_output_tests, \
             patch('end2end_test_administrator.cfg') as mock_cfg:
            
            mock_tests = MagicMock()
            mock_tests.run_all.return_value = (0, 1)  # 0 failed, 1 total
            mock_output_tests.return_value = mock_tests
            
            e2e_admin.main("evaluate")
            
            mock_output_tests.assert_called_once_with("evaluate")
            mock_tests.run_all.assert_called_once()

    def test_main_function_generate_mode(self):
        """Test main function in generate mode."""
        with patch('end2end_test_administrator.prepare'), \
             patch('end2end_test_administrator.OutputTests') as mock_output_tests, \
             patch('end2end_test_administrator.cfg') as mock_cfg:
            
            mock_tests = MagicMock()
            mock_tests.run_all.return_value = (0, 1)
            mock_output_tests.return_value = mock_tests
            
            e2e_admin.main("generate_benchmark")
            
            mock_output_tests.assert_called_once_with("generate_benchmark")


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = ['OutputTests']
        for class_name in expected_classes:
            self.assertTrue(hasattr(e2e_admin, class_name))
            self.assertTrue(callable(getattr(e2e_admin, class_name)))

    def test_expected_functions_present(self):
        """Test that all expected functions are present in the module."""
        expected_functions = [
            'move_subfolders_to_old', 'prepare', 'main'
        ]
        for func_name in expected_functions:
            self.assertTrue(hasattr(e2e_admin, func_name))
            self.assertTrue(callable(getattr(e2e_admin, func_name)))


class TestModuleDocumentation(unittest.TestCase):
    """Test cases for module and function documentation."""

    def test_move_subfolders_to_old_has_docstring(self):
        """Test that move_subfolders_to_old has docstring."""
        self.assertIsNotNone(e2e_admin.move_subfolders_to_old.__doc__)
        self.assertGreater(len(e2e_admin.move_subfolders_to_old.__doc__), 0)

    def test_move_subfolders_to_old_docstring_contains_example(self):
        """Test that move_subfolders_to_old docstring contains example."""
        docstring = e2e_admin.move_subfolders_to_old.__doc__
        self.assertIn("Example layout", docstring)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)