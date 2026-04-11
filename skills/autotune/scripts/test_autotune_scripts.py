from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ADD_LOG_ENTRY = SCRIPT_DIR / "add_log_entry.py"
INIT_EXPERIMENT = SCRIPT_DIR / "init_experiment.py"
RECORD_RESULT = SCRIPT_DIR / "record_result.py"


class AutotuneScriptTests(unittest.TestCase):
    def test_add_log_entry_inserts_a_separator_for_existing_logs_without_newline(self) -> None:
        today = date.today().isoformat()

        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "log.md"
            log_file.write_text(f"## {today}", encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(ADD_LOG_ENTRY), str(log_file), "New finding"],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(log_file.read_text(encoding="utf-8"), f"## {today}\n- New finding\n")

    def test_init_experiment_rejects_empty_slugs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "docs" / "autotune"

            result = subprocess.run(
                [
                    sys.executable,
                    str(INIT_EXPERIMENT),
                    "!!!",
                    "--root",
                    str(root),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(root.exists())
            self.assertIn("error", result.stderr.lower())

    def test_record_result_writes_metric_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                    "--metric",
                    "score=1.5",
                    "--metric",
                    "time=2",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                metric_file.read_text(encoding="utf-8"),
                (
                    "[[metric]]\n"
                    '"commit" = "abc123"\n'
                    '"status" = "keep"\n'
                    '"description" = "baseline run"\n'
                    '"score" = 1.5\n'
                    '"time" = 2.0\n\n'
                ),
            )

    def test_record_result_separates_existing_content_without_trailing_newline(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"
            metric_file.write_text(
                (
                    "[[metric]]\n"
                    '"commit" = "old"\n'
                    '"status" = "keep"\n'
                    '"description" = "previous run"'
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                    "--metric",
                    "score=1.5",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                metric_file.read_text(encoding="utf-8"),
                (
                    "[[metric]]\n"
                    '"commit" = "old"\n'
                    '"status" = "keep"\n'
                    '"description" = "previous run"\n'
                    "[[metric]]\n"
                    '"commit" = "abc123"\n'
                    '"status" = "keep"\n'
                    '"description" = "baseline run"\n'
                    '"score" = 1.5\n\n'
                ),
            )

    def test_record_result_rejects_missing_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(metric_file.exists())
            self.assertIn("error", result.stderr.lower())

    def test_record_result_rejects_malformed_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                    "--metric",
                    "invalid",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(metric_file.exists())
            self.assertIn("error", result.stderr.lower())

    def test_record_result_rejects_reserved_metric_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                    "--metric",
                    "commit=1",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(metric_file.exists())
            self.assertIn("reserved metric key", result.stderr.lower())

    def test_record_result_rejects_duplicate_metric_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            metric_file = Path(tmpdir) / "metric.toml"

            result = subprocess.run(
                [
                    sys.executable,
                    str(RECORD_RESULT),
                    str(metric_file),
                    "--commit",
                    "abc123",
                    "--status",
                    "keep",
                    "--description",
                    "baseline run",
                    "--metric",
                    "score=1.5",
                    "--metric",
                    "score=2.5",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(metric_file.exists())
            self.assertIn("duplicate metric key", result.stderr.lower())

    def test_record_result_rejects_non_finite_metrics(self) -> None:
        for metric_value in ("nan", "inf", "-inf"):
            with self.subTest(metric_value=metric_value):
                with tempfile.TemporaryDirectory() as tmpdir:
                    metric_file = Path(tmpdir) / "metric.toml"

                    result = subprocess.run(
                        [
                            sys.executable,
                            str(RECORD_RESULT),
                            str(metric_file),
                            "--commit",
                            "abc123",
                            "--status",
                            "keep",
                            "--description",
                            "baseline run",
                            "--metric",
                            f"score={metric_value}",
                        ],
                        capture_output=True,
                        text=True,
                        check=False,
                    )

                    self.assertNotEqual(result.returncode, 0)
                    self.assertFalse(metric_file.exists())
                    self.assertIn("finite numeric value", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
