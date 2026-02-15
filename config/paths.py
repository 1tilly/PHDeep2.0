import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    data_root: Path
    output_root: Path
    reference_root: Path
    logs_root: Path


def load_paths() -> ProjectPaths:
    base = Path(os.getenv("PHDEEP_ROOT", Path.cwd()))
    data_root = Path(os.getenv("PHDEEP_DATA_ROOT", base / "data"))
    output_root = Path(os.getenv("PHDEEP_OUTPUT_ROOT", base / "results"))
    reference_root = Path(os.getenv("PHDEEP_REFERENCE_ROOT", data_root / "reference_genome"))
    logs_root = Path(os.getenv("PHDEEP_LOGS_ROOT", output_root / "logs"))

    return ProjectPaths(
        data_root=data_root,
        output_root=output_root,
        reference_root=reference_root,
        logs_root=logs_root,
    )
