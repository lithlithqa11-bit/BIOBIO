# Binding site definition and grid calculations
import re
import numpy as np
from io import StringIO
from Bio.PDB import PDBParser
from docking.models import GridBox

def extract_reference_ligand(
    pdb_string: str,
    ligand_resname: str,
    chain_id: str | None = None,
    resseq: str | None = None,
) -> tuple[np.ndarray, str]:
    """Extract coordinates and PDB lines for a specific co-crystallised reference ligand.
    If multiple instances exist across chains or residue numbers and none is explicitly specified,
    the first complete instance (chain, resseq) is selected to represent the single binding site."""
    ligand_resname = (ligand_resname or "").strip().upper()
    if not ligand_resname:
        raise ValueError("No reference ligand code provided.")

    coords = []
    ligand_lines = []
    instances = []
    lines_by_instance = {}
    coords_by_instance = {}

    for line in pdb_string.splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            rname = line[17:20].strip().upper()
            if rname == ligand_resname:
                ligand_lines.append(line)
                ch = line[21:22]
                seq = line[22:26].strip()
                inst_key = (ch, seq)
                if inst_key not in lines_by_instance:
                    instances.append(inst_key)
                    lines_by_instance[inst_key] = []
                    coords_by_instance[inst_key] = []

                lines_by_instance[inst_key].append(line)
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    coords_by_instance[inst_key].append([x, y, z])
                except ValueError:
                    pass

    if not instances:
        raise ValueError(f"Reference ligand '{ligand_resname}' was not found in PDB coordinates.")

    chosen_key = None
    if chain_id is not None or resseq is not None:
        matching = []
        for k in instances:
            ch_match = (chain_id is None) or (k[0].upper() == str(chain_id).strip().upper())
            seq_match = (resseq is None) or (k[1] == str(resseq).strip())
            if ch_match and seq_match:
                matching.append(k)
        if not matching:
            raise ValueError(
                f"Reference ligand '{ligand_resname}' matching chain={chain_id}, resseq={resseq} not found."
            )
        if len(matching) > 1:
            choices = ", ".join(f"{chain or '_'}:{seq}" for chain, seq in matching)
            raise ValueError(
                f"Reference ligand {ligand_resname!r} has {len(matching)} matching instances ({choices}). "
                "Choose its chain and residue number explicitly."
            )
        chosen_key = matching[0]
    else:
        if len(instances) != 1:
            choices = ", ".join(f"{chain or '_'}:{seq}" for chain, seq in instances)
            raise ValueError(
                f"Reference ligand {ligand_resname!r} has {len(instances)} instances ({choices}). "
                "Choose its chain and residue number explicitly."
            )
        chosen_key = instances[0]

    coords = coords_by_instance[chosen_key]
    ligand_lines = lines_by_instance[chosen_key]

    if not coords:
        raise ValueError(f"No coordinates found for reference ligand '{ligand_resname}'.")

    coords_arr = np.array(coords, dtype=np.float32)
    return coords_arr, "\n".join(ligand_lines)

def define_grid_from_ligand(
    pdb_string: str,
    ligand_resname: str,
    padding: float = 8.0,
    chain_id: str | None = None,
    resseq: str | None = None,
) -> tuple[GridBox, np.ndarray, str]:
    """Define a finite 3D docking search box centered on the crystal reference ligand."""
    coords_arr, ligand_pdb_text = extract_reference_ligand(
        pdb_string, ligand_resname, chain_id=chain_id, resseq=resseq
    )
    center = coords_arr.mean(axis=0)

    ranges = coords_arr.max(axis=0) - coords_arr.min(axis=0)
    sizes = np.maximum(ranges + 2 * padding, 15.0)

    grid = GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(sizes[0]), 1),
        size_y=round(float(sizes[1]), 1),
        size_z=round(float(sizes[2]), 1)
    )
    return grid, coords_arr, ligand_pdb_text

def _parse_residue_identifier(raw: str) -> tuple[int, str]:
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", raw.strip())
    if not match:
        raise ValueError(
            f"Invalid residue identifier {raw!r}. Use a number optionally followed "
            "by one insertion code, for example A:100 or A:100A."
        )
    return int(match.group(1)), match.group(2).upper()


def _residue_identifier(residue) -> str:
    insertion_code = str(residue.id[2]).strip().upper()
    return f"{residue.id[1]}{insertion_code}"


def define_grid_from_residues(pdb_string: str, residue_specs: list[str], padding: float = 10.0) -> GridBox:
    if not residue_specs:
        raise ValueError("No binding-site residues provided. Docking requires a validated pocket.")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", StringIO(pdb_string))
    model = structure[0]

    parsed_specs: list[tuple[str | None, str]] = []
    for raw_spec in residue_specs:
        raw_spec = str(raw_spec).strip()
        if not raw_spec:
            continue
        if ":" in raw_spec:
            chain_id, residue_number = (part.strip() for part in raw_spec.split(":", 1))
            if not chain_id or not residue_number:
                raise ValueError(f"Invalid residue specification: {raw_spec!r}. Use CHAIN:RESIDUE, e.g. A:790.")
            parsed_specs.append((chain_id, residue_number))
        else:
            parsed_specs.append((None, raw_spec))

    coordinates = []
    for chain_id, residue_number in parsed_specs:
        wanted_number, wanted_insertion_code = _parse_residue_identifier(residue_number)
        matches = []
        for chain in model:
            if chain_id is not None and chain.id != chain_id:
                continue
            for residue in chain:
                if (
                    residue.id[1] == wanted_number
                    and str(residue.id[2]).strip().upper() == wanted_insertion_code
                ):
                    matches.append(residue)

        if not matches:
            label = f"{chain_id}:{residue_number}" if chain_id else residue_number
            raise ValueError(f"Binding-site residue {label} was not found in the receptor.")
        if chain_id is None and len(matches) > 1:
            matching_labels = ", ".join(f"{res.get_parent().id}:{_residue_identifier(res)}" for res in matches)
            raise ValueError(
                f"Residue {residue_number} occurs in more than one chain ({matching_labels}). "
                "Use explicit CHAIN:RESIDUE values, for example A:790, A:858."
            )
        for residue in matches:
            coordinates.extend(atom.get_coord() for atom in residue.get_atoms())

    coords_arr = np.asarray(coordinates, dtype=np.float32)
    if coords_arr.size == 0:
        raise ValueError("No atom coordinates were found for the selected binding-site residues.")

    center = coords_arr.mean(axis=0)
    ranges = coords_arr.max(axis=0) - coords_arr.min(axis=0)
    sizes = np.maximum(ranges + 2 * padding, 18.0)

    return GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(sizes[0]), 1),
        size_y=round(float(sizes[1]), 1),
        size_z=round(float(sizes[2]), 1)
    )
