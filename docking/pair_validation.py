"""Pair validation and homologous residue mapping for Matched Docking comparisons.

Ensures that matched comparisons are only conducted between scientifically validated
homologous chains (e.g. wild-type and its genuine variant), preventing misleading delta scores
on unrelated proteins.
"""

from io import StringIO
from Bio.PDB import PDBParser
from analysis.constants import AA_3TO1
from analysis.alignment import get_alignment
from docking.binding_site import _residue_identifier

MIN_MATCHED_PAIR_IDENTITY: float = 90.0


def extract_chain_sequence_with_ids(pdb_text: str, chain_id: str) -> list[dict]:
    """Extract standard amino acid residues for a given chain with exact residue numbering and insertion codes."""
    if not pdb_text or not chain_id:
        return []

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("protein", StringIO(pdb_text))
    model = struct[0]

    residues = []
    if chain_id in [c.id for c in model]:
        for res in model[chain_id]:
            # Standard amino acids have hetflag ' '
            if res.id[0] == ' ':
                resname = res.get_resname().strip().upper()
                one_letter = AA_3TO1.get(resname, 'X')
                resseq = res.id[1]
                icode = str(res.id[2]).strip().upper()
                res_id_str = f"{resseq}{icode}"
                residues.append({
                    "chain": chain_id,
                    "res_num": resseq,
                    "icode": icode,
                    "res_id": res_id_str,
                    "full_id": f"{chain_id}:{res_id_str}",
                    "res_name": resname,
                    "one_letter": one_letter,
                })
    return residues


def validate_matched_pair(
    healthy_pdb_text: str,
    mutant_pdb_text: str,
    healthy_chain: str,
    mutant_chain: str,
    min_identity: float = MIN_MATCHED_PAIR_IDENTITY,
) -> dict:
    """Validate that the selected healthy and mutant chains form a homologous WT/variant pair.

    Performs global sequence alignment, computes sequence identity strictly over aligned non-gap positions,
    and constructs an explicit healthy-to-mutant residue mapping.
    """
    if not healthy_pdb_text or not mutant_pdb_text:
        return {
            "is_valid": False,
            "identity_pct": 0.0,
            "min_identity": min_identity,
            "error": "Missing coordinates for healthy or mutant receptor.",
            "mapping": {},
            "aligned_residue_count": 0,
            "mutation_count": 0,
            "mutations": [],
        }

    if not healthy_chain or not mutant_chain:
        return {
            "is_valid": False,
            "identity_pct": 0.0,
            "min_identity": min_identity,
            "error": "Healthy or mutant chain is not selected.",
            "mapping": {},
            "aligned_residue_count": 0,
            "mutation_count": 0,
            "mutations": [],
        }

    h_residues = extract_chain_sequence_with_ids(healthy_pdb_text, healthy_chain)
    m_residues = extract_chain_sequence_with_ids(mutant_pdb_text, mutant_chain)

    if not h_residues or not m_residues:
        return {
            "is_valid": False,
            "identity_pct": 0.0,
            "min_identity": min_identity,
            "error": f"Could not extract amino acid sequence for chain {healthy_chain} or {mutant_chain}.",
            "mapping": {},
            "aligned_residue_count": 0,
            "mutation_count": 0,
            "mutations": [],
        }

    seq_h = "".join(r["one_letter"] for r in h_residues)
    seq_m = "".join(r["one_letter"] for r in m_residues)

    _, score, aligned_h, aligned_m = get_alignment(seq_h, seq_m, mode="global")
    if not aligned_h or not aligned_m:
        return {
            "is_valid": False,
            "identity_pct": 0.0,
            "min_identity": min_identity,
            "error": "Global sequence alignment failed.",
            "mapping": {},
            "aligned_residue_count": 0,
            "mutation_count": 0,
            "mutations": [],
        }

    # Walk alignment and build healthy -> mutant residue mapping
    h_ptr = 0
    m_ptr = 0
    mapping = {}
    mutations = []
    aligned_non_gap = 0
    identical_count = 0

    for ch_h, ch_m in zip(aligned_h, aligned_m):
        curr_h = h_residues[h_ptr] if ch_h != '-' else None
        curr_m = m_residues[m_ptr] if ch_m != '-' else None

        if curr_h is not None and curr_m is not None:
            aligned_non_gap += 1
            if ch_h == ch_m:
                identical_count += 1
            else:
                mutations.append({
                    "healthy_res": curr_h["full_id"],
                    "healthy_aa": curr_h["res_name"],
                    "mutant_res": curr_m["full_id"],
                    "mutant_aa": curr_m["res_name"],
                })

            # Map both full spec (e.g. "B:6") and bare spec (e.g. "6")
            mapping[curr_h["full_id"]] = curr_m["full_id"]
            mapping[curr_h["res_id"]] = curr_m["full_id"]

        if ch_h != '-':
            h_ptr += 1
        if ch_m != '-':
            m_ptr += 1

    identity_pct = (identical_count / aligned_non_gap * 100.0) if aligned_non_gap > 0 else 0.0
    identity_pct = round(identity_pct, 2)
    is_valid = identity_pct >= min_identity

    error_msg = None
    if not is_valid:
        error_msg = (
            f"Matched comparison was not run: selected chains are not a validated WT/variant pair "
            f"(sequence identity {identity_pct:.1f}%, required >= {min_identity:.1f}%). "
            f"Use Exploratory Cross-Protein Docking instead."
        )

    return {
        "is_valid": is_valid,
        "identity_pct": identity_pct,
        "min_identity": min_identity,
        "aligned_residue_count": aligned_non_gap,
        "mutation_count": len(mutations),
        "mutations": mutations,
        "healthy_chain": healthy_chain,
        "mutant_chain": mutant_chain,
        "mapping": mapping,
        "alignment_score": score,
        "error": error_msg,
    }


def extract_protein_chains(pdb_text: str) -> list[str]:
    """Extract list of chain IDs that contain standard amino acid residues."""
    if not pdb_text:
        return []
    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("protein", StringIO(pdb_text))
        chains = []
        for c in struct[0]:
            if any(r.id[0] == ' ' for r in c):
                chains.append(c.id)
        return chains or ["A"]
    except Exception:
        return ["A"]


def map_healthy_pocket_residues(
    healthy_pocket_residues: list[str],
    mapping: dict[str, str],
    mutant_pdb_text: str,
    mutant_chain: str,
    healthy_chain: str | None = None,
) -> tuple[list[str], list[str], list[str]]:
    """Map healthy pocket residue specifications to equivalent mutant residue specifications.

    Returns:
        (mapped_mutant_residues, unmapped_healthy_residues, validation_errors)
    """
    mapped_mutant = []
    unmapped = []
    errors = []

    # Parse mutant structure to verify mapped residues actually exist
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("mutant", StringIO(mutant_pdb_text))
    model = struct[0]

    mutant_existing = set()
    if mutant_chain in [c.id for c in model]:
        for res in model[mutant_chain]:
            if res.id[0] == ' ':
                icode = str(res.id[2]).strip().upper()
                mutant_existing.add(f"{mutant_chain}:{res.id[1]}{icode}")
                mutant_existing.add(f"{res.id[1]}{icode}")

    for raw_spec in healthy_pocket_residues:
        spec = raw_spec.strip()
        if not spec:
            continue

        target_m = mapping.get(spec)
        if not target_m and ":" not in spec and healthy_chain:
            target_m = mapping.get(f"{healthy_chain}:{spec}")

        if target_m:
            # Verify mapped residue exists in mutant
            if target_m in mutant_existing:
                mapped_mutant.append(target_m)
            else:
                errors.append(f"Mapped residue {target_m} was not found in mutant chain {mutant_chain}.")
        else:
            unmapped.append(spec)
            errors.append(f"Healthy pocket residue {spec} does not map to any homologous position in mutant chain {mutant_chain}.")

    return mapped_mutant, unmapped, errors
