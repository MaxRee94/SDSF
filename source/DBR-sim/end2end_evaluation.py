import filecmp
import os

import file_handling as io


class Evaluator:
    REQUIRED_SUBDIRS = ("image_timeseries", "state_data")

    def __init__(self, testcase, cfg):
        self.testcase = testcase
        self.benchmark_dir = self.get_benchmark_dir(cfg)

    def get_benchmark_dir(self, cfg):
        benchmark_dir = os.path.join(cfg.END2END_TEST_BENCHMARK_DIR, self.testcase)
        io.create_directory_if_not_exists(benchmark_dir)
        return benchmark_dir

    # -------------------------
    # Public API
    # -------------------------
    def evaluate(self, test):
        self._assert_required_dirs_exist(test.output_dir, self.benchmark_dir)
        self._compare_all_benchmark_files(self.benchmark_dir, test.output_dir)
        self._assert_no_extra_output_files(self.benchmark_dir, test.output_dir)

        return True

    # -------------------------
    # Path helpers
    # -------------------------
    def _out_path(self, out_root: str, rel_path: str) -> str:
        return os.path.join(out_root, rel_path)

    def _bench_path(self, bench_root: str, rel_path: str) -> str:
        return os.path.join(bench_root, rel_path)

    # -------------------------
    # Assertions / comparisons
    # -------------------------
    def _assert_required_dirs_exist(self, out_root: str, bench_root: str) -> None:
        # Root must exist
        self._assert_is_dir(out_root)
        self._assert_is_dir(bench_root)  # print output-root on failure

        # Required subdirs must exist in both
        for sub in self.REQUIRED_SUBDIRS:
            out_sub = self._out_path(out_root, sub)
            bench_sub = self._bench_path(bench_root, sub)
            self._assert_is_dir(out_sub)
            self._assert_is_dir(bench_sub)  # print output-sub on failure

    def _assert_is_dir(self, path: str) -> None:
        if not os.path.isdir(path):
            raise AssertionError(f"Directory missing or not a directory: {path}")

    def _compare_all_benchmark_files(self, bench_root: str, out_root: str) -> None:
        """
        Walk benchmark tree; for every file found, ensure output tree contains the
        same relative file with identical bytes.
        """
        for rel_path in self._iter_files_relative_to(bench_root):
            self._assert_output_file_matches(bench_root, out_root, rel_path)

    def _iter_files_relative_to(self, root: str):
        for dirpath, _, filenames in os.walk(root):
            for fname in filenames:
                full = os.path.join(dirpath, fname)
                yield os.path.relpath(full, root)

    def _assert_output_file_matches(self, bench_root: str, out_root: str, rel_path: str) -> None:
        bench_path = self._bench_path(bench_root, rel_path)
        out_path = self._out_path(out_root, rel_path)

        if not os.path.exists(out_path):
            print(out_path)
            raise AssertionError(f"Missing output file: {out_path}")

        if not os.path.isfile(out_path):
            print(out_path)
            raise AssertionError(f"Output path is not a file: {out_path}")

        if not filecmp.cmp(bench_path, out_path, shallow=False):
            print(out_path)
            raise AssertionError(f"Output differs from benchmark: {out_path}")

    def _assert_no_extra_output_files(self, bench_root: str, out_root: str) -> None:
        """
        Fail if output contains files that benchmark doesn't.
        """
        for rel_path in self._iter_files_relative_to(out_root):
            print("rel_path:", rel_path)
            bench_path = self._bench_path(bench_root, rel_path)
            if not os.path.exists(bench_path):
                out_path = self._out_path(out_root, rel_path)
                print(out_path)
                raise AssertionError(f"Unexpected extra output file: {out_path}") 

