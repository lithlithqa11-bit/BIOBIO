# Constants and biochemical properties for amino acids and mutations
"""
Constants module:
- AA_3TO1: Mapping from 3-letter to 1-letter amino acid codes.
- AA_PROPS: Physicochemical classification of amino acids.
- REMOTE_JSON_URL: GitHub repository URL for documented mutations.
"""

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

AA_PROPS = {
    'ALA': 'Non-polar',  'GLY': 'Non-polar', 'ILE': 'Non-polar', 
    'LEU': 'Non-polar',  'MET': 'Non-polar', 'PRO': 'Non-polar',
    'VAL': 'Non-polar',  'ASN': 'Polar',     'CYS': 'Disulfide-forming',   
    'GLN': 'Polar',      'SER': 'Polar',     'THR': 'Polar',
    'ASP': 'Acidic (-)', 'GLU': 'Acidic (-)',
    'ARG': 'Basic (+)',  'HIS': 'Basic (+)',  'LYS': 'Basic (+)',
    'PHE': 'Aromatic',   'TRP': 'Aromatic',   'TYR': 'Aromatic'
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}

# Van der Waals amino acid volumes in Å³ (Zamyatnin / Chothia standard scale)
AA_VOLUMES = {
    'GLY': 60.1,  'ALA': 88.6,  'SER': 89.0,  'CYS': 108.5, 'ASP': 111.1,
    'PRO': 112.7, 'ASN': 114.1, 'THR': 116.1, 'GLU': 138.4, 'VAL': 140.0,
    'GLN': 143.8, 'HIS': 153.2, 'MET': 162.9, 'ILE': 166.7, 'LEU': 166.7,
    'LYS': 168.6, 'ARG': 173.4, 'PHE': 189.9, 'TYR': 193.6, 'TRP': 227.8
}

# Remote repository URL for documented known mutations
REMOTE_JSON_URL = "https://raw.githubusercontent.com/ha9160034-create/Bio/main/mutation.json"
