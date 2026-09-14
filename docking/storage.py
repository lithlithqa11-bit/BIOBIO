# Storage and Manifest Management
import os
import json
import uuid
import hashlib
import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VINA_PATH = PROJECT_ROOT / "tools" / "vina" / "vina.exe"
BASE_RUNS_DIR = PROJECT_ROOT / "docking_runs"

DISCLAIMER_TEXT = (
    "تنبيه علمي: هذه النتائج تمثل فرضية ارتباط حاسوبية (In-Silico Docking Hypothesis) "
    "تهدف لترتيب الوضعيات والمركبات ضمن هذا البروتوكول المحدد. "
    "لا تُعد هذه الدرجات قياساً معملياً لطاقة الارتباط الحرة أو دليلاً على الفاعلية السريرية أو الجرعة العلاجية."
)

def get_vina_path() -> Path:
    return VINA_PATH

def create_run_directory(prefix="run") -> tuple[str, Path]:
    BASE_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    uid = uuid.uuid4().hex[:6]
    run_id = f"{prefix}_{ts}_{uid}"
    run_dir = BASE_RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_id, run_dir

def write_initial_manifest(run_dir: Path, data: dict) -> Path:
    initial_path = run_dir / "manifest.initial.json"
    with open(initial_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return initial_path

def finalize_manifest(run_dir: Path, data: dict) -> Path:
    manifest_path = run_dir / "manifest.json"
    temp_path = run_dir / "manifest.json.tmp"
    with open(temp_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    temp_path.replace(manifest_path)
    return manifest_path

def load_manifest(manifest_path: Path) -> dict:
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def compute_file_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum for manifest provenance tracking."""
    if not filepath.exists() or not filepath.is_file():
        return ""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def build_and_save_manifest(
    run_dir: Path,
    run_id: str,
    status: str,
    target_meta: dict,
    binding_site_meta: dict,
    ligand_meta: dict,
    prep_meta: dict,
    engine_meta: dict,
    grid_meta: dict,
    results_meta: dict,
    warnings: list[str] = None
) -> Path:
    """Build and write the final manifest for a docking run."""
    run_dir = Path(run_dir)

    valid_statuses = [
        "validation passed",
        "validation passed (Top-10 recovery)",
        "validation failed",
        "validation unavailable",
        "exploratory—no reference validation",
        "exploratory—insufficient independent seeds",
        "run failed"
    ]
    if status not in valid_statuses:
        status = "run failed" if "fail" in status.lower() else "exploratory—no reference validation"

    created_at_utc = target_meta.get("created_at_utc") or ""

    manifest_data = {
        "schema_version": 1,
        "run_id": run_id,
        "status": status,
        "created_at_utc": created_at_utc,
        "target": target_meta,
        "binding_site": binding_site_meta,
        "ligand": ligand_meta,
        "receptor_preparation": prep_meta,
        "engine": engine_meta,
        "grid": grid_meta,
        "inputs": {
            "sha256": {
                p.name: compute_file_sha256(p)
                for p in sorted(run_dir.iterdir())
                if p.is_file() and not p.name.startswith("manifest")
            }
        },
        "results": results_meta,
        "warnings": warnings or [],
        "scientific_disclaimer": DISCLAIMER_TEXT
    }
    return finalize_manifest(run_dir, manifest_data)
