# RCSB PDB coordinate downloads and mutation dataset fetching
import gzip
import json
import requests
from pathlib import Path
import streamlit as st
from analysis.constants import REMOTE_JSON_URL

@st.cache_data(ttl=3600)
def fetch_deposited_pdb(pdb_id: str):
    """Download the coordinates deposited in the PDB entry (the asymmetric unit)."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        return response.text
    except requests.RequestException as error:
        st.error(f"Could not download deposited coordinates for {pdb_id}: {error}")
        return None

@st.cache_data(ttl=3600)
def fetch_biological_assembly(pdb_id: str, assembly_id: int = 1):
    """Download the official RCSB biological assembly coordinate file.

    This file already contains the correct symmetry-generated oligomer. It must
    be used instead of guessing assemblies from chain names or sequences.
    """
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_id}.gz"
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        return gzip.decompress(response.content).decode("utf-8")
    except (requests.RequestException, OSError, UnicodeDecodeError) as error:
        st.warning(
            f"Biological Assembly {assembly_id} is unavailable for {pdb_id}. "
            f"The deposited coordinates will be used unchanged. Details: {error}"
        )
        return None

@st.cache_data(ttl=60)
def load_mutation_db():
    """تحميل قاعدة بيانات الطفرات الموثقة لربط البروتينات السليمة بالمصابة تلقائياً (محلياً أو من GitHub)."""
    base_dir = Path(__file__).resolve().parents[1]
    
    # 1. البحث في مجلد data/ أولاً
    local_data_file = base_dir / "data" / "mutation.json"
    if local_data_file.exists():
        try:
            with open(local_data_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass

    # 2. البحث في الجذر كاحتياط
    root_file = base_dir / "mutation.json"
    if root_file.exists():
        try:
            with open(root_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass

    # 3. الجلب السحابي من GitHub
    url = REMOTE_JSON_URL
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            return response.json()
    except Exception:
        pass

    return {}
