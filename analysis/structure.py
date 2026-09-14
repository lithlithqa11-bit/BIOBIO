# Biopython structure parsing, chain handling, KDTree contact analysis, and SASA maps
from io import StringIO
import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley
import streamlit as st

from analysis.constants import AA_3TO1, AA_PROPS

def structure_to_pdb_str(structure):
    """تحويل كائن Structure الخاص بـ Biopython إلى نص PDB صالح للعرض والتحليل."""
    try:
        io = PDBIO()
        io.set_structure(structure)
        out = StringIO()
        io.save(out)
        return out.getvalue()
    except Exception:
        return ""

def process_protein_structure(pdb_string: str, pdb_id: str):
    """Parse one already-selected coordinate set and calculate residue SASA.

    `pdb_string` must be either the deposited coordinates or an official RCSB
    biological-assembly file. This function must never delete protein chains.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, StringIO(pdb_string))

        # Remove only water. Keep protein chains and functional ligands such as
        # heme (HEM), because deleting HEM makes hemoglobin visualization and
        # contact analysis less meaningful.
        for model in structure:
            for chain in model:
                water_ids = [residue.id for residue in chain if residue.id[0] == 'W']
                for residue_id in water_ids:
                    chain.detach_child(residue_id)

        sr = ShrakeRupley()
        sr.compute(structure, level='R')
        return structure
    except Exception as error:
        st.error(f"Could not read structure {pdb_id}: {error}")
        return None

def get_all_chains(structure):
    """استخراج جميع معرفات السلاسل (Chains) الموجودة في بنية البروتين."""
    try:
        return [chain.id for chain in structure[0]]
    except Exception as ex:
        st.error(f"فشل قراءة السلاسل: {ex}")
        return []

def get_protein_sequence(structure, chain_id):
    """استخراج تسلسل الأحماض الأمينية لسلسلة محددة مع أرقامها التسلسلية."""
    sequence = []
    try:
        model = structure[0]
        if chain_id in [c.id for c in model]:
            for residue in model[chain_id]:
                if residue.id[0] == ' ':  # التأكد من أنه حمض أميني حقيقي
                    sequence.append({
                        'res_num': residue.id[1],
                        'res_name': residue.get_resname()
                    })
    except Exception as ex:
        st.error(f"فشل قراءة التسلسل للسلسلة {chain_id}: {ex}")
    return sequence

def sequence_to_fasta(structure, chain_id, protein_name="protein"):
    """تحويل تسلسل الأحماض الأمينية إلى تنسيق FASTA القياسي."""
    seq_data = get_protein_sequence(structure, chain_id)
    if not seq_data:
        return None
    one_letter = ''.join([AA_3TO1.get(r['res_name'], 'X') for r in seq_data])
    lines = [one_letter[i:i+80] for i in range(0, len(one_letter), 80)]
    header = f">{protein_name}|Chain_{chain_id}|length={len(one_letter)}\n"
    return header + '\n'.join(lines) + '\n'

def calculate_all_distances(structure, chain_id: str, radius: float = 5.0):
    """Calculate residue contacts inside the exact structure shown to the user."""
    if not structure:
        return []

    try:
        model = structure[0]
        if chain_id not in [chain.id for chain in model]:
            return []

        all_atoms = list(model.get_atoms())
        all_coords = np.array([atom.get_coord() for atom in all_atoms], dtype=np.float32)
        tree = KDTree(all_coords)

        atom_info = []
        for atom in all_atoms:
            residue = atom.get_parent()
            atom_info.append((residue.get_parent().id, residue.id[1], residue.id[2]))

        results = []
        for residue in model[chain_id]:
            if residue.id[0] != ' ':
                continue

            target_id = (chain_id, residue.id[1], residue.id[2])
            target_coords = np.array(
                [atom.get_coord() for atom in residue.get_atoms()], dtype=np.float32
            )
            nearby_indices = set()
            for index_list in tree.query_ball_point(target_coords, radius):
                for index in index_list:
                    other_chain, other_resnum, _ = atom_info[index]
                    # Exclude self and immediate covalent sequence neighbors (|i - j| <= 1)
                    if other_chain == chain_id and abs(other_resnum - residue.id[1]) <= 1:
                        continue
                    nearby_indices.add(index)

            min_dist = "-"
            if nearby_indices:
                neighbour_coords = all_coords[list(nearby_indices)]
                min_dist = round(float(cdist(target_coords, neighbour_coords).min()), 2)

            sasa = getattr(residue, 'sasa', '-')
            if isinstance(sasa, (float, int)):
                sasa = round(sasa, 2)

            residue_name = residue.get_resname()
            results.append({
                'chain': chain_id,
                'res_num': residue.id[1],
                'res_name': residue_name,
                'one_letter': AA_3TO1.get(residue_name, 'X'),
                'class': AA_PROPS.get(residue_name, '-'),
                'min_dist': min_dist,
                'sasa': sasa,
            })

        return results
    except Exception as error:
        st.error(f"Distance analysis failed: {error}")
        return []

def calculate_sasa_map(structure, chain_id):
    """حساب مساحة السطح المعرضة للمذيب (SASA) وإنشاء خريطة تربط رقم الحمض بقيمته."""
    try:
        sasa_map = {}
        for res in structure[0][chain_id]:
            if res.id[0] == ' ':
                sasa_map[res.id[1]] = round(getattr(res, 'sasa', 0), 2)
        return sasa_map
    except Exception:
        return {}
