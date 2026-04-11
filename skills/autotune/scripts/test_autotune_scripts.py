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


if __name__ == "__main__":
    unittest.main()
