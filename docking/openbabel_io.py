from __future__ import annotations

from pathlib import Path

from openbabel import openbabel as ob


def _conversion(input_format: str, output_format: str) -> ob.OBConversion:
    conversion = ob.OBConversion()
    if not conversion.SetInFormat(input_format):
        raise RuntimeError(f"Open Babel does not provide input format: {input_format}")
    if not conversion.SetOutFormat(output_format):
        raise RuntimeError(f"Open Babel does not provide output format: {output_format}")
    return conversion


def require_openbabel_formats() -> list[str]:
    """Fail early with a clear deployment error if required formats are absent."""
    required = (("pdb", "pdbqt"), ("mol2", "pdbqt"), ("sdf", "pdbqt"), ("pdbqt", "pdb"), ("pdbqt", "mol"))
    verified = set()
    for input_format, output_format in required:
        _conversion(input_format, output_format)
        verified.add(input_format)
        verified.add(output_format)
    return sorted(verified)


def read_file(path: str | Path, input_format: str) -> ob.OBMol:
    path = Path(path)
    conversion = _conversion(input_format, "mol")
    molecule = ob.OBMol()
    if not conversion.ReadFile(molecule, str(path)) or molecule.NumAtoms() == 0:
        raise RuntimeError(
            f"Open Babel could not read {path.name} as {input_format}; "
            f"atoms read: {molecule.NumAtoms()}."
        )
    return molecule


def read_text(text: str, input_format: str) -> ob.OBMol:
    conversion = _conversion(input_format, "mol")
    molecule = ob.OBMol()
    if not conversion.ReadString(molecule, text) or molecule.NumAtoms() == 0:
        raise RuntimeError(
            f"Open Babel could not parse in-memory {input_format}; "
            f"atoms read: {molecule.NumAtoms()}."
        )
    return molecule


def write_file(
    molecule: ob.OBMol,
    path: str | Path,
    output_format: str,
    *,
    rigid_receptor: bool = False,
) -> Path:
    path = Path(path)
    conversion = _conversion("mol", output_format)
    if rigid_receptor:
        conversion.AddOption("r", ob.OBConversion.OUTOPTIONS)
    if not conversion.WriteFile(molecule, str(path)):
        raise RuntimeError(f"Open Babel could not write {path.name} as {output_format}.")
    conversion.CloseOutFile()
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Open Babel wrote no usable output file: {path.name}.")
    return path


def write_text(molecule: ob.OBMol, output_format: str) -> str:
    conversion = _conversion("mol", output_format)
    text = conversion.WriteString(molecule)
    if not text.strip():
        raise RuntimeError(f"Open Babel produced empty {output_format} text.")
    return text
