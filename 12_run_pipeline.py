"""Root-level launcher for the corrected candidate-gene pipeline.

This file is intentionally kept in the repository root so the pipeline can be
started from Windows with:

    python 12_run_pipeline.py --help

The implementation delegates to modules/12_run_pipeline.py and ensures that
modules/ is importable regardless of the current working directory.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODULES_DIR = ROOT / "modules"
RUNNER = MODULES_DIR / "12_run_pipeline.py"

if not RUNNER.exists():
    raise FileNotFoundError(
        f"Corrected pipeline runner was not found: {RUNNER}"
    )

if str(MODULES_DIR) not in sys.path:
    sys.path.insert(0, str(MODULES_DIR))

spec = importlib.util.spec_from_file_location(
    "candidate_gene_pipeline_runner",
    RUNNER,
)

if spec is None or spec.loader is None:
    raise ImportError(f"Could not load pipeline runner: {RUNNER}")

runner = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = runner
spec.loader.exec_module(runner)


if __name__ == "__main__":
    runner.main()
