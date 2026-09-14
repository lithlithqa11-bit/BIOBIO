# PubChem Integration Module for Drug Lookup
from __future__ import annotations
import json
import urllib.parse
import urllib.request
from typing import Optional

# Built-in offline fallback cache for key pharmacological agents
_OFFLINE_DRUG_CACHE: dict[str, dict] = {
    "gefitinib": {
        "CID": 123631,
        "Title": "Gefitinib",
        "CanonicalSMILES": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
        "MolecularWeight": "446.9",
        "IUPACName": "N-(3-chloro-4-fluorophenyl)-7-methoxy-6-(3-morpholin-4-ylpropoxy)quinazolin-4-amine"
    },
    "erlotinib": {
        "CID": 2826731,
        "Title": "Erlotinib",
        "CanonicalSMILES": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
        "MolecularWeight": "393.4",
        "IUPACName": "N-(3-ethynylphenyl)-6,7-bis(2-methoxyethoxy)quinazolin-4-amine"
    },
    "osimertinib": {
        "CID": 71496458,
        "Title": "Osimertinib",
        "CanonicalSMILES": "CN1CCN(CC1)C2=CC(=C(C=C2)NC3=NC=CC(=N3)C4=CN(C5=CC=CC=C54)C)NC(=O)C=C",
        "MolecularWeight": "499.6",
        "IUPACName": "N-(2-{[2-(dimethylamino)ethyl](methyl)amino}-4-methoxy-5-{[4-(1-methyl-1H-indol-3-yl)pyrimidin-2-yl]amino}phenyl)acrylamide"
    },
    "imatinib": {
        "CID": 5291,
        "Title": "Imatinib",
        "CanonicalSMILES": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        "MolecularWeight": "493.6",
        "IUPACName": "4-[(4-methylpiperazin-1-yl)methyl]-N-[4-methyl-3-[(4-pyridin-3-ylpyrimidin-2-yl)amino]phenyl]benzamide"
    },
    "aspirin": {
        "CID": 2244,
        "Title": "Aspirin",
        "CanonicalSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "MolecularWeight": "180.16",
        "IUPACName": "2-acetyloxybenzoic acid"
    },
    "paracetamol": {
        "CID": 1983,
        "Title": "Acetaminophen",
        "CanonicalSMILES": "CC(=O)NC1=CC=C(C=C1)O",
        "MolecularWeight": "151.16",
        "IUPACName": "N-(4-hydroxyphenyl)acetamide"
    },
    "remdesivir": {
        "CID": 121304016,
        "Title": "Remdesivir",
        "CanonicalSMILES": "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4",
        "MolecularWeight": "602.6",
        "IUPACName": "2-ethylbutyl (2S)-2-[[[(2R,3S,4R,5R)-5-(4-aminopyrrolo[2,1-f][1,2,4]triazin-7-yl)-5-cyano-3,4-dihydroxyoxolan-2-yl]methoxy-phenoxyphosphoryl]amino]propanoate"
    },
    "vemurafenib": {
        "CID": 42611257,
        "Title": "Vemurafenib",
        "CanonicalSMILES": "CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F",
        "MolecularWeight": "489.9",
        "IUPACName": "N-[3-[5-(4-chlorophenyl)-1H-pyrrolo[2,3-b]pyridine-3-carbonyl]-2,4-difluorophenyl]propane-1-sulfonamide"
    },
    "voxelotor": {
        "CID": 89851851,
        "Title": "Voxelotor",
        "CanonicalSMILES": "CC1=NN(C=C1)CC2=CC=CC(=C2O)C(=O)OCC3=CC=CC(=C3)C#N",
        "MolecularWeight": "337.3",
        "IUPACName": "2-hydroxy-6-[(1-methylpyrazol-4-yl)methoxy]benzaldehyde"
    }
}

def fetch_drug_from_pubchem(drug_name: str, timeout: int = 5) -> dict:
    """
    Look up a drug compound by generic or trade name via NCBI PubChem PUG-REST API.
    Returns a dictionary with success status, canonical SMILES, MW, IUPAC name, and CID.
    Falls back to offline curated knowledge for core oncology/clinical agents if network fails.
    """
    clean_name = drug_name.strip()
    if not clean_name:
        return {"success": False, "error": "يرجى إدخال اسم الدواء للبحث."}

    norm_key = clean_name.lower()

    # Try live query first
    try:
        encoded = urllib.parse.quote(clean_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES,ConnectivitySMILES,MolecularWeight,IUPACName,Title/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "BioImpactAnalyzer/2.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            if response.status == 200:
                data = json.loads(response.read().decode("utf-8"))
                props = data.get("PropertyTable", {}).get("Properties", [])
                if props:
                    p = props[0]
                    smiles = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES", "")
                    return {
                        "success": True,
                        "title": p.get("Title", clean_name),
                        "smiles": smiles,
                        "molecular_weight": float(p.get("MolecularWeight", 0.0)),
                        "iupac_name": p.get("IUPACName", ""),
                        "cid": p.get("CID", None),
                        "source": "PubChem API (Live)"
                    }
    except Exception as e:
        # Check offline cache before failing
        if norm_key in _OFFLINE_DRUG_CACHE:
            p = _OFFLINE_DRUG_CACHE[norm_key]
            return {
                "success": True,
                "title": p.get("Title", clean_name),
                "smiles": p.get("CanonicalSMILES", ""),
                "molecular_weight": float(p.get("MolecularWeight", 0.0)),
                "iupac_name": p.get("IUPACName", ""),
                "cid": p.get("CID", None),
                "source": "Offline Clinical Pharmacopoeia"
            }
        return {
            "success": False,
            "error": f"تعذر جلب بيانات الدواء '{clean_name}' من PubChem ({e}). تأكد من صحة الاسم باللغة الإنجليزية أو الاتصال بالإنترنت."
        }

    # Fallback to cache if empty live response
    if norm_key in _OFFLINE_DRUG_CACHE:
        p = _OFFLINE_DRUG_CACHE[norm_key]
        return {
            "success": True,
            "title": p.get("Title", clean_name),
            "smiles": p.get("CanonicalSMILES", ""),
            "molecular_weight": float(p.get("MolecularWeight", 0.0)),
            "iupac_name": p.get("IUPACName", ""),
            "cid": p.get("CID", None),
            "source": "Offline Clinical Pharmacopoeia"
        }

    return {"success": False, "error": f"لم يتم العثور على مركب كيميائي مسجل باسم '{clean_name}' في PubChem."}
