"""Unit test suite for end2end_evaluation module.

This module provides comprehensive tests for the end-to-end evaluation functionality.
Tests cover the Evaluator class, directory validation, file comparison, and error handling.
"""

import unittest
import os
import tempfile
import shutil
import sys

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import end2end_evaluation as e2e_eval
import file_handling as io


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self, benchmark_dir):
        self.END2END_TEST_BENCHMARK_DIR = benchmark_dir


class MockTest:
    """Mock test object for testing."""
    def __init__(self, output_dir):
        self.output_dir = output_dir


class TestEvaluator(unittest.TestCase):
    """Test cases for the Evaluator class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directories for testing
        self.temp_dir = tempfile.mkdtemp()
        self.benchmark_dir = os.path.join(self.temp_dir, "benchmarks")
        self.output_dir = os.path.join(self.temp_dir, "output")
        
        os.makedirs(self.benchmark_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create mock config
        self.cfg = MockConfig(self.benchmark_dir)
        
        # Create required subdirectories
        for subdir in e2e_eval.Evaluator.REQUIRED_SUBDIRS:
            bench_subdir = os.path.join(self.benchmark_dir, subdir)
            output_subdir = os.path.join(self.output_dir, subdir)
            os.makedirs(bench_subdir, exist_ok=True)
            os.makedirs(output_subdir, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        # Remove temporary directory
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    # ==================== INITIALIZATION TESTS ====================

    def test_evaluator_initialization(self):
        """Test that Evaluator initializes correctly."""
        testcase = "test_case_1"
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        
        self.assertEqual(evaluator.testcase, testcase)
        self.assertEqual(evaluator.benchmark_dir, os.path.join(self.benchmark_dir, testcase))

    def test_get_benchmark_dir_creates_directory(self):
        """Test that get_benchmark_dir creates directory if it doesn't exist."""
        testcase = "new_test_case"
        new_benchmark_path = os.path.join(self.benchmark_dir, testcase)
        
        # Ensure directory doesn't exist initially
        if os.path.exists(new_benchmark_path):
            shutil.rmtree(new_benchmark_path)
        
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        
        # Check that directory was created
        self.assertTrue(os.path.isdir(evaluator.benchmark_dir))

    def test_get_benchmark_dir_existing_directory(self):
        """Test that get_benchmark_dir works with existing directory."""
        testcase = "existing_case"
        existing_path = os.path.join(self.benchmark_dir, testcase)
        os.makedirs(existing_path, exist_ok=True)
        
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        self.assertEqual(evaluator.benchmark_dir, existing_path)

    # ==================== PATH HELPER TESTS ====================

    def test_out_path_helper(self):
        """Test the _out_path helper method."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        result = evaluator._out_path("/base/path", "subdir/file.txt")
        expected = "/base/path/subdir/file.txt"
        
        # Use os.path.normpath to handle platform-specific path separators
        self.assertEqual(os.path.normpath(result), os.path.normpath(expected))

    def test_bench_path_helper(self):
        """Test the _bench_path helper method."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        result = evaluator._bench_path("/base/path", "subdir/file.txt")
        expected = "/base/path/subdir/file.txt"
        
        # Use os.path.normpath to handle platform-specific path separators
        self.assertEqual(os.path.normpath(result), os.path.normpath(expected))

    # ==================== DIRECTORY VALIDATION TESTS ====================

    def test_assert_is_dir_valid_directory(self):
        """Test that _assert_is_dir works with valid directories."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Should not raise for valid directory
        try:
            evaluator._assert_is_dir(self.benchmark_dir)
            success = True
        except AssertionError:
            success = False
        
        self.assertTrue(success)

    def test_assert_is_dir_nonexistent_directory(self):
        """Test that _assert_is_dir raises for non-existent directories."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_is_dir("/nonexistent/path/12345")
        
        self.assertIn("Directory missing or not a directory", str(context.exception))

    def test_assert_is_dir_file_instead_of_directory(self):
        """Test that _assert_is_dir raises when given a file instead of directory."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create a file
        test_file = os.path.join(self.temp_dir, "test_file.txt")
        with open(test_file, 'w') as f:
            f.write("test content")
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_is_dir(test_file)
        
        self.assertIn("Directory missing or not a directory", str(context.exception))

    def test_assert_required_dirs_exist_all_present(self):
        """Test _assert_required_dirs_exist when all required directories exist."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Should not raise when all directories exist
        try:
            evaluator._assert_required_dirs_exist(self.output_dir, self.benchmark_dir)
            success = True
        except AssertionError:
            success = False
        
        self.assertTrue(success)

    def test_assert_required_dirs_exist_missing_output_root(self):
        """Test _assert_required_dirs_exist when output root is missing."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_required_dirs_exist("/nonexistent/output", self.benchmark_dir)
        
        self.assertIn("Directory missing or not a directory", str(context.exception))

    def test_assert_required_dirs_exist_missing_benchmark_root(self):
        """Test _assert_required_dirs_exist when benchmark root is missing."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_required_dirs_exist(self.output_dir, "/nonexistent/benchmark")
        
        self.assertIn("Directory missing or not a directory", str(context.exception))

    def test_assert_required_dirs_exist_missing_required_subdir(self):
        """Test _assert_required_dirs_exist when a required subdirectory is missing."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Remove one required subdirectory from output
        missing_subdir = os.path.join(self.output_dir, e2e_eval.Evaluator.REQUIRED_SUBDIRS[0])
        shutil.rmtree(missing_subdir)
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_required_dirs_exist(self.output_dir, self.benchmark_dir)
        
        self.assertIn("Directory missing or not a directory", str(context.exception))

    # ==================== FILE ITERATION TESTS ====================

    def test_iter_files_relative_to_empty_directory(self):
        """Test _iter_files_relative_to with empty directory."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create an empty subdirectory
        empty_subdir = os.path.join(self.benchmark_dir, "empty_subdir")
        os.makedirs(empty_subdir, exist_ok=True)
        
        files = list(evaluator._iter_files_relative_to(empty_subdir))
        self.assertEqual(len(files), 0)

    def test_iter_files_relative_to_with_files(self):
        """Test _iter_files_relative_to finds files in directory."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create some test files
        test_subdir = os.path.join(self.benchmark_dir, "test_subdir")
        os.makedirs(test_subdir, exist_ok=True)
        
        with open(os.path.join(test_subdir, "file1.txt"), 'w') as f:
            f.write("content1")
        with open(os.path.join(test_subdir, "file2.txt"), 'w') as f:
            f.write("content2")
        
        files = list(evaluator._iter_files_relative_to(self.benchmark_dir))
        
        # Should find the files (relative paths)
        self.assertGreater(len(files), 0)
        self.assertTrue(any("file1.txt" in f for f in files))
        self.assertTrue(any("file2.txt" in f for f in files))

    def test_iter_files_relative_to_nested_directories(self):
        """Test _iter_files_relative_to with nested directory structure."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create nested directory structure
        nested_dir = os.path.join(self.benchmark_dir, "level1", "level2")
        os.makedirs(nested_dir, exist_ok=True)
        
        with open(os.path.join(nested_dir, "nested_file.txt"), 'w') as f:
            f.write("nested content")
        
        files = list(evaluator._iter_files_relative_to(self.benchmark_dir))
        
        # Should find the nested file
        nested_file_found = any("level1" in f and "level2" in f and "nested_file.txt" in f for f in files)
        self.assertTrue(nested_file_found)

    # ==================== FILE COMPARISON TESTS ====================

    def test_assert_output_file_matches_identical_files(self):
        """Test _assert_output_file_matches with identical files."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create identical files in both benchmark and output
        rel_path = "test_identical.txt"
        bench_file = os.path.join(self.benchmark_dir, rel_path)
        output_file = os.path.join(self.output_dir, rel_path)
        
        content = "This is identical content for testing"
        with open(bench_file, 'w') as f:
            f.write(content)
        with open(output_file, 'w') as f:
            f.write(content)
        
        # Should not raise
        try:
            evaluator._assert_output_file_matches(self.benchmark_dir, self.output_dir, rel_path)
            success = True
        except AssertionError:
            success = False
        
        self.assertTrue(success)

    def test_assert_output_file_matches_different_content(self):
        """Test _assert_output_file_matches with different file content."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create files with different content
        rel_path = "test_different.txt"
        bench_file = os.path.join(self.benchmark_dir, rel_path)
        output_file = os.path.join(self.output_dir, rel_path)
        
        with open(bench_file, 'w') as f:
            f.write("benchmark content")
        with open(output_file, 'w') as f:
            f.write("different output content")
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_output_file_matches(self.benchmark_dir, self.output_dir, rel_path)
        
        self.assertIn("Output differs from benchmark", str(context.exception))

    def test_assert_output_file_matches_missing_output_file(self):
        """Test _assert_output_file_matches when output file is missing."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create only benchmark file
        rel_path = "test_missing.txt"
        bench_file = os.path.join(self.benchmark_dir, rel_path)
        
        with open(bench_file, 'w') as f:
            f.write("benchmark content")
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_output_file_matches(self.benchmark_dir, self.output_dir, rel_path)
        
        self.assertIn("Missing output file", str(context.exception))

    def test_assert_output_file_matches_output_is_directory(self):
        """Test _assert_output_file_matches when output path is a directory."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create benchmark file and output directory with same name
        rel_path = "test_dir.txt"
        bench_file = os.path.join(self.benchmark_dir, rel_path)
        output_dir_path = os.path.join(self.output_dir, rel_path)
        
        with open(bench_file, 'w') as f:
            f.write("benchmark content")
        os.makedirs(output_dir_path, exist_ok=True)
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_output_file_matches(self.benchmark_dir, self.output_dir, rel_path)
        
        self.assertIn("Output path is not a file", str(context.exception))

    # ==================== EXTRA FILE DETECTION TESTS ====================

    def test_assert_no_extra_output_files_none_extra(self):
        """Test _assert_no_extra_output_files when there are no extra files."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create identical files in both benchmark and output
        for subdir in e2e_eval.Evaluator.REQUIRED_SUBDIRS:
            bench_subdir = os.path.join(self.benchmark_dir, subdir)
            output_subdir = os.path.join(self.output_dir, subdir)
            
            test_file = os.path.join(bench_subdir, "common_file.txt")
            with open(test_file, 'w') as f:
                f.write("common content")
            
            test_file = os.path.join(output_subdir, "common_file.txt")
            with open(test_file, 'w') as f:
                f.write("common content")
        
        # Should not raise
        try:
            evaluator._assert_no_extra_output_files(self.benchmark_dir, self.output_dir)
            success = True
        except AssertionError:
            success = False
        
        self.assertTrue(success)

    def test_assert_no_extra_output_files_with_extra_file(self):
        """Test _assert_no_extra_output_files when there are extra files."""
        evaluator = e2e_eval.Evaluator("test", self.cfg)
        
        # Create an extra file in output that doesn't exist in benchmark
        extra_file = os.path.join(self.output_dir, "extra_file.txt")
        with open(extra_file, 'w') as f:
            f.write("extra content")
        
        with self.assertRaises(AssertionError) as context:
            evaluator._assert_no_extra_output_files(self.benchmark_dir, self.output_dir)
        
        self.assertIn("Unexpected extra output file", str(context.exception))

    # ==================== INTEGRATION TESTS ====================

    def test_evaluate_successful_case(self):
        """Test the full evaluate() method with a successful case."""
        testcase = "test_case"
        # Create benchmark directory with required subdirectories
        benchmark_path = os.path.join(self.benchmark_dir, testcase)
        os.makedirs(benchmark_path, exist_ok=True)
        
        # Create required subdirectories in both benchmark and output
        for subdir in e2e_eval.Evaluator.REQUIRED_SUBDIRS:
            bench_subdir = os.path.join(benchmark_path, subdir)
            output_subdir = os.path.join(self.output_dir, subdir)
            os.makedirs(bench_subdir, exist_ok=True)
            os.makedirs(output_subdir, exist_ok=True)
            
            # Create a test file in each
            test_file_rel = "test_file.txt"
            with open(os.path.join(bench_subdir, test_file_rel), 'w') as f:
                f.write("benchmark content")
            with open(os.path.join(output_subdir, test_file_rel), 'w') as f:
                f.write("benchmark content")  # Same content
        
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        mock_test = MockTest(self.output_dir)
        
        # Should return True
        result = evaluator.evaluate(mock_test)
        self.assertTrue(result)

    def test_evaluate_fails_missing_directory(self):
        """Test the evaluate() method fails when directories are missing."""
        testcase = "test_case"
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        mock_test = MockTest("/nonexistent/output")
        
        with self.assertRaises(AssertionError):
            evaluator.evaluate(mock_test)

    # ==================== EDGE CASE TESTS ====================

    def test_empty_testcase_name(self):
        """Test Evaluator with empty testcase name."""
        evaluator = e2e_eval.Evaluator("", self.cfg)
        self.assertEqual(evaluator.testcase, "")
        self.assertEqual(evaluator.benchmark_dir, os.path.join(self.benchmark_dir, ""))

    def test_testcase_with_special_characters(self):
        """Test Evaluator with special characters in testcase name."""
        testcase = "test-case_with.special+chars"
        evaluator = e2e_eval.Evaluator(testcase, self.cfg)
        
        self.assertEqual(evaluator.testcase, testcase)
        self.assertEqual(evaluator.benchmark_dir, os.path.join(self.benchmark_dir, testcase))

    def test_required_subdirs_constant(self):
        """Test that REQUIRED_SUBDIRS is properly defined."""
        self.assertIsInstance(e2e_eval.Evaluator.REQUIRED_SUBDIRS, tuple)
        self.assertGreater(len(e2e_eval.Evaluator.REQUIRED_SUBDIRS), 0)
        
        # Check that all subdirs are strings
        for subdir in e2e_eval.Evaluator.REQUIRED_SUBDIRS:
            self.assertIsInstance(subdir, str)
            self.assertGreater(len(subdir), 0)


class TestFileHandlingIntegration(unittest.TestCase):
    """Test integration between end2end_evaluation and file_handling modules."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.benchmark_dir = os.path.join(self.temp_dir, "benchmarks")
        os.makedirs(self.benchmark_dir, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_create_directory_if_not_exists_integration(self):
        """Test that file_handling.create_directory_if_not_exists is used correctly."""
        testcase = "integration_test"
        benchmark_path = os.path.join(self.benchmark_dir, testcase)
        
        # Directory should not exist initially
        self.assertFalse(os.path.exists(benchmark_path))
        
        # Create mock config
        cfg = MockConfig(self.benchmark_dir)
        
        # This should create the directory
        evaluator = e2e_eval.Evaluator(testcase, cfg)
        
        # Directory should now exist
        self.assertTrue(os.path.exists(benchmark_path))
        self.assertTrue(os.path.isdir(benchmark_path))


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)