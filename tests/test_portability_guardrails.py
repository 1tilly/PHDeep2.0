from pathlib import Path

FORBIDDEN_PATH_SNIPPETS = ("/home/", "/rds/")


def _python_sources(repo_root: Path) -> list[Path]:
    candidates: list[Path] = []
    for root in (repo_root / "src", repo_root / "config"):
        candidates.extend(sorted(root.rglob("*.py")))
    candidates.append(repo_root / "main.py")
    return candidates


def test_no_hardcoded_site_paths_in_core_code():
    repo_root = Path(__file__).resolve().parents[1]
    offenders: list[str] = []

    for path in _python_sources(repo_root):
        text = path.read_text(encoding="utf-8")
        for snippet in FORBIDDEN_PATH_SNIPPETS:
            if snippet in text:
                offenders.append(f"{path.relative_to(repo_root)} contains '{snippet}'")

    assert not offenders, "Hardcoded site-specific paths found:\n" + "\n".join(offenders)
