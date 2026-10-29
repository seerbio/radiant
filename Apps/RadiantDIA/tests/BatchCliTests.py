#!/usr/bin/env python3
"""Exercise batch orchestration without running an expensive spectrum search."""

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class BatchCliTests(unittest.TestCase):
    binary = None
    fixtures = None

    def setUp(self):
        temporary_directory = tempfile.TemporaryDirectory(prefix="radiant-batch-cli-")
        self.addCleanup(temporary_directory.cleanup)
        self.work_directory = Path(temporary_directory.name)

    def run_batch(self, inputs, output_folder=None, library_path=None):
        if library_path is None:
            library_path = self.fixtures / "FragLibReaderTests.fragLibFF"
        command = [
            str(self.binary),
            str(library_path),
            str(self.fixtures / "human_plasma_entrapment_super_trunc.fasta"),
            str(self.fixtures / "test_params_wide_window.pythiaConfig"),
            *[str(path) for path in inputs],
        ]
        if output_folder is not None:
            command.extend(["--output-folder", str(output_folder)])
        return subprocess.run(
            command,
            cwd=self.work_directory,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8",
            errors="replace",
            timeout=60,
            check=False,
        )

    def assert_shared_argument_rejected_before_library_load(self, result):
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertNotIn("Reading library", result.stdout)
        self.assertNotIn("Processing file", result.stdout)
        self.assertNotIn("Batch summary:", result.stdout)
        self.assertIn("Usage:", result.stdout)
        self.assertEqual(list(self.work_directory.rglob("*.radiantDIA")), [])

    def assert_all_inputs_failed_in_order(self, result, inputs):
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertEqual(result.stdout.count("Reading library"), 1, result.stdout)
        previous_position = -1
        for index, path in enumerate(inputs, start=1):
            marker = f"Processing file {index}/{len(inputs)}: {path}"
            self.assertIn(marker, result.stdout)
            position = result.stdout.index(marker)
            self.assertGreater(position, previous_position, result.stdout)
            previous_position = position
            self.assertIn(f"FAILED: {path}", result.stdout)
        self.assertIn(
            f"Batch summary: 0 completed, {len(inputs)} failed, {len(inputs)} total",
            result.stdout,
        )
        self.assertNotIn("COMPLETED:", result.stdout)
        self.assertEqual(list(self.work_directory.rglob("*.radiantDIA")), [])

    @staticmethod
    def write_corrupt_mzml(path):
        # The malformed spectrum attribute causes an XML reader error before scoring.
        path.write_text(
            '<mzML>\n<spectrum id="scan=1" index="0" invalid>\n</mzML>\n',
            encoding="utf-8",
        )

    def test_unsupported_and_missing_inputs_both_report_failure(self):
        unsupported = self.work_directory / "unsupported.txt"
        unsupported.write_text("unsupported input", encoding="utf-8")
        missing = self.work_directory / "missing.mzML"
        inputs = [unsupported, missing]

        result = self.run_batch(inputs)

        self.assert_all_inputs_failed_in_order(result, inputs)
        self.assertIn("expected mzML, Thermo raw, Bruker d, or prqFF input", result.stdout)

    def test_missing_data_argument_returns_failure_before_library_load(self):
        result = self.run_batch([])

        self.assert_shared_argument_rejected_before_library_load(result)
        self.assertIn(
            "Expected library, FASTA, config, and at least one data path",
            result.stdout,
        )

    def test_wrong_library_suffix_returns_failure_before_library_load(self):
        wrong_library = self.work_directory / "library.txt"
        wrong_library.write_text("wrong library format", encoding="utf-8")

        result = self.run_batch(
            [self.work_directory / "missing.mzML"],
            library_path=wrong_library,
        )

        self.assert_shared_argument_rejected_before_library_load(result)
        self.assertIn("First command line argument must be a library file", result.stdout)

    def test_corrupt_xml_does_not_prevent_attempting_later_input(self):
        corrupt = self.work_directory / "corrupt.mzML"
        self.write_corrupt_mzml(corrupt)
        inputs = [corrupt, self.work_directory / "later-missing.mzML"]

        result = self.run_batch(inputs)

        self.assert_all_inputs_failed_in_order(result, inputs)
        self.assertIn("XML error:", result.stdout)

    def test_duplicate_output_names_rejected_before_loading_library(self):
        inputs = []
        for directory_name in ["first", "second"]:
            directory = self.work_directory / directory_name
            directory.mkdir()
            path = directory / "sample.mzML"
            self.write_corrupt_mzml(path)
            inputs.append(path)
        output_folder = self.work_directory / "results"

        result = self.run_batch(inputs, output_folder=output_folder)

        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("Batch inputs would overwrite the same output:", result.stdout)
        self.assertNotIn("Reading library", result.stdout)
        self.assertNotIn("Processing file", result.stdout)
        self.assertEqual(list(output_folder.iterdir()), [])

    def test_file_symlink_alias_with_distinct_name_does_not_collide(self):
        original = self.work_directory / "original.mzML"
        self.write_corrupt_mzml(original)
        alias = self.work_directory / "alias.mzML"
        try:
            alias.symlink_to(original)
        except (OSError, NotImplementedError) as error:
            self.skipTest(f"File symlinks are unavailable: {error}")
        inputs = [original, alias]

        result = self.run_batch(inputs)

        self.assert_all_inputs_failed_in_order(result, inputs)
        self.assertNotIn("Batch inputs would overwrite the same output:", result.stdout)

    def test_all_seventeen_missing_inputs_are_attempted(self):
        inputs = [
            self.work_directory / f"missing-{index}.mzML"
            for index in range(17)
        ]

        result = self.run_batch(inputs)

        self.assert_all_inputs_failed_in_order(result, inputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--fixtures", required=True, type=Path)
    arguments, unittest_arguments = parser.parse_known_args()
    BatchCliTests.binary = arguments.binary.resolve()
    BatchCliTests.fixtures = arguments.fixtures.resolve()
    if not BatchCliTests.binary.is_file():
        parser.error(f"RadiantDIA executable does not exist: {BatchCliTests.binary}")
    for filename in [
        "FragLibReaderTests.fragLibFF",
        "human_plasma_entrapment_super_trunc.fasta",
        "test_params_wide_window.pythiaConfig",
    ]:
        fixture = BatchCliTests.fixtures / filename
        if not fixture.is_file():
            parser.error(f"Required fixture does not exist: {fixture}")
    unittest.main(argv=[sys.argv[0], *unittest_arguments])


if __name__ == "__main__":
    main()
