# Protein-Ligand Interaction and Pharmacological Profiler
"""
وحدة تحليل التفاعلات الكيميائية بين الدواء والبروتين، وكشف الروابط الهيدروجينية،
وتقييم الخصائص الصيدلانية ومعايير ليبينسكي (Lipinski's Rule of 5)، وتفسير مقاومة الأدوية.
"""

from io import StringIO
import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def calculate_drug_properties(smiles: str) -> dict:
    """
    حساب الخصائص الصيدلانية والكيميائية للدواء وفق قواعد ليبينسكي (Lipinski Rule of Five)
    لتقييم مدى ملائمة الجزيء كدواء فموي (Drug-likeness / ADME).
    """
    if not smiles or not smiles.strip():
        return {}
    
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return {"error": "صيغة SMILES غير صالحة كيميائياً"}

    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # التحقق من قواعد ليبينسكي الخمس (Lipinski Violations)
    violations = 0
    if mw > 500: violations += 1
    if logp > 5.0: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    lipinski_pass = (violations <= 1)

    return {
        "formula": formula,
        "molecular_weight": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "tpsa": tpsa,
        "lipinski_violations": violations,
        "lipinski_pass": lipinski_pass,
        "drug_likeness": "ممتازة (مطابق لقواعد ليبينسكي)" if lipinski_pass else f"تنبيه: يوجد {violations} تجاوزات لقواعد ليبينسكي"
    }

def _parse_atoms_from_block(block_text: str, is_ligand: bool = False) -> list[dict]:
    """استخراج إحداثيات وخصائص الذرات من نصوص PDB أو PDBQT."""
    atoms = []
    if not block_text:
        return atoms

    for line in block_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 54:
                continue
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip() or "A"
                res_num = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                # استخراج العنصر الكيميائي
                element = ""
                if len(line) >= 78:
                    element = line[76:78].strip().upper()
                if not element:
                    # تخمين العنصر من الحرف الأول لاسم الذرة
                    element = "".join([c for c in atom_name if c.isalpha()][:1]).upper()

                atoms.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain": chain_id,
                    "res_num": res_num,
                    "coord": np.array([x, y, z], dtype=np.float32),
                    "element": element,
                    "is_ligand": is_ligand
                })
            except (ValueError, IndexError):
                continue

    return atoms

def analyze_protein_ligand_interactions(
    receptor_pdb_text: str,
    ligand_pdbqt_block: str,
    hbond_cutoff: float = 3.5,
    contact_cutoff: float = 4.0
) -> dict:
    """
    تحليل دقيق للتفاعلات الكيميائية بين وضعية الدواء المستقبلة والمستقبل البروتيني:
    1. الروابط الهيدروجينية (Hydrogen Bonds): ذرات O/N بمسافة <= 3.5 Å.
    2. التماسات الكارهة للماء (Hydrophobic Contacts): ذرات كربون بمسافة <= 4.0 Å.
    3. شبكة الأحماض المحيطة بجيب الارتباط (Interacting Pocket Residues).
    """
    rec_atoms = _parse_atoms_from_block(receptor_pdb_text, is_ligand=False)
    lig_atoms = _parse_atoms_from_block(ligand_pdbqt_block, is_ligand=True)

    if not rec_atoms or not lig_atoms:
        return {
            "hbonds": [],
            "hydrophobic": [],
            "interacting_residues": [],
            "hbond_count": 0,
            "hydrophobic_count": 0
        }

    rec_coords = np.array([a["coord"] for a in rec_atoms])
    lig_coords = np.array([a["coord"] for a in lig_atoms])

    dist_matrix = cdist(rec_coords, lig_coords)

    hbonds = []
    hydrophobic = []
    interacting_res_map = {}

    hbond_elements = {"O", "N", "F"}

    for r_idx, r_atom in enumerate(rec_atoms):
        for l_idx, l_atom in enumerate(lig_atoms):
            dist = float(dist_matrix[r_idx, l_idx])

            # تسجيل الأحماض المتفاعلة بشكل عام ضمن النطاق القريب
            if dist <= contact_cutoff:
                res_key = f"{r_atom['res_name']} {r_atom['res_num']} (Chain {r_atom['chain']})"
                if res_key not in interacting_res_map or dist < interacting_res_map[res_key]["min_distance"]:
                    interacting_res_map[res_key] = {
                        "res_name": r_atom["res_name"],
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "min_distance": round(dist, 2)
                    }

            # 1. كشف الروابط الهيدروجينية (H-Bonds)
            if dist <= hbond_cutoff:
                if r_atom["element"] in hbond_elements and l_atom["element"] in hbond_elements:
                    hbonds.append({
                        "residue": f"{r_atom['res_name']}{r_atom['res_num']}",
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "rec_atom": r_atom["atom_name"],
                        "lig_atom": l_atom["atom_name"],
                        "distance_angstrom": round(dist, 2)
                    })

            # 2. كشف التماسات الكارهة للماء (Carbon-Carbon)
            elif dist <= contact_cutoff:
                if r_atom["element"] == "C" and l_atom["element"] == "C":
                    hydrophobic.append({
                        "residue": f"{r_atom['res_name']}{r_atom['res_num']}",
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "distance_angstrom": round(dist, 2)
                    })

    # إزالة التكرارات للروابط لنفس زوج الذرات
    unique_hbonds = []
    seen_hb = set()
    for hb in hbonds:
        key = (hb["residue"], hb["rec_atom"], hb["lig_atom"])
        if key not in seen_hb:
            seen_hb.add(key)
            unique_hbonds.append(hb)

    unique_interacting_res = sorted(
        list(interacting_res_map.values()),
        key=lambda x: x["min_distance"]
    )

    return {
        "hbonds": unique_hbonds,
        "hydrophobic": hydrophobic[:25],
        "interacting_residues": unique_interacting_res,
        "hbond_count": len(unique_hbonds),
        "hydrophobic_count": len(hydrophobic)
    }

def compare_interactions(
    healthy_interactions: dict,
    mutant_interactions: dict,
    mutation_residue_nums: list = None
) -> dict:
    """
    مقارنة التفاعلات الكيميائية بين السليم والمصاب لتحديد:
    - الروابط الهيدروجينية المفقودة (Lost H-Bonds) بسبب الطفرة.
    - الروابط الهيدروجينية الجديدة (Gained H-Bonds).
    - هل يدخل حمض الطفرة في تفاعل مباشر مع الدواء؟
    """
    h_hb_res = {hb["residue"] for hb in healthy_interactions.get("hbonds", [])}
    m_hb_res = {hb["residue"] for hb in mutant_interactions.get("hbonds", [])}

    lost_hb_residues = list(h_hb_res - m_hb_res)
    gained_hb_residues = list(m_hb_res - h_hb_res)
    conserved_hb_residues = list(h_hb_res.intersection(m_hb_res))

    h_all_res_nums = {r["res_num"] for r in healthy_interactions.get("interacting_residues", [])}
    m_all_res_nums = {r["res_num"] for r in mutant_interactions.get("interacting_residues", [])}

    mut_nums_clean = [str(n).strip() for n in (mutation_residue_nums or [])]

    mutation_engaged_healthy = any(n in h_all_res_nums for n in mut_nums_clean)
    mutation_engaged_mutant = any(n in m_all_res_nums for n in mut_nums_clean)

    return {
        "lost_hb_residues": lost_hb_residues,
        "gained_hb_residues": gained_hb_residues,
        "conserved_hb_residues": conserved_hb_residues,
        "healthy_hbond_count": healthy_interactions.get("hbond_count", 0),
        "mutant_hbond_count": mutant_interactions.get("hbond_count", 0),
        "mutation_engaged_in_binding": mutation_engaged_healthy or mutation_engaged_mutant,
        "mutation_residues_checked": mut_nums_clean
    }

def interpret_resistance_risk(
    delta_score: float,
    lost_hbonds_count: int,
    mutation_engaged: bool
) -> dict:
    """
    تفسير صيدلاني سريري لنتائج فارق طاقة الارتباط ΔScore ومخاطر مقاومة الدواء:
    ΔScore = Score_mutant - Score_healthy (تذكر: الدرجة الأقل بالسالب تعني ارتباط أقوى)
    إذا كانت ΔScore موجبة بشكل كبير (+1.5 kcal/mol فأكثر) -> ضعف شديد في الارتباط (مقاومة دوائية).
    """
    if delta_score >= 1.5:
        risk_level = "🔴 خطر مقاومة دوائية مرتفع (High Resistance Risk)"
        risk_class = "danger"
        summary = (
            f"أظهرت المحاكاة ضعفاً كبيراً في تقارب الدواء للبروتين المصاب (فارق الطاقة ΔScore = +{delta_score:.2f} kcal/mol). "
            f"يرتبط هذا التراجع عادةً بفشل العلاج السريري بسبب الإعاقة الفراغية أو فقدان روابط تثبيت حيوية."
        )
    elif delta_score >= 0.6:
        risk_level = "🟠 خطر مقاومة متوسط / انخفاض في الفعالية (Moderate Resistance Risk)"
        risk_class = "warning"
        summary = (
            f"هناك انخفاض ملحوظ في تقارب الدواء (ΔScore = +{delta_score:.2f} kcal/mol). "
            f"قد يتطلب هذا رفع الجرعة أو استبدال الدواء بمثبط من جيل أحدث."
        )
    elif delta_score <= -1.5:
        risk_level = "🟢 ارتباط انتقائي بالطافر (Mutant-Selective Binding)"
        risk_class = "success"
        summary = (
            f"أظهر الدواء تفضيلاً كبيراً للارتباط بالبروتين المصاب (ΔScore = {delta_score:.2f} kcal/mol). "
            f"هذا السلوك مثالي للأدوية الموجهة للطفرات (مثل مثبطات الجيل الثالث) لتقليل السمية على الخلايا السليمة."
        )
    else:
        risk_level = "⚪ استجابة متقاربة / فعالية محفوظة (Preserved Binding / Minor Impact)"
        risk_class = "info"
        summary = (
            f"فارق طاقة الارتباط طفيف جداً (ΔScore = {delta_score:+.2f} kcal/mol). "
            f"من المتوقع احتفاظ الدواء بفعاليته الأساسية تجاه البروتين الطافر."
        )

    clinical_pointers = []
    if lost_hbonds_count > 0:
        clinical_pointers.append(f"تم فقدان عدد ({lost_hbonds_count}) من الروابط الهيدروجينية المحورية في موقع الارتباط.")
    if mutation_engaged:
        clinical_pointers.append("حمض الطفرة يتواجد في مسافة تماس مباشر مع جزيء الدواء (Direct Binding Site Mutation).")
    else:
        clinical_pointers.append("حمض الطفرة يقع على أطراف الجيب أو يؤثر بشكل خيفي غير مباشر (Allosteric / Peripheral Effect).")

    return {
        "risk_level": risk_level,
        "risk_class": risk_class,
        "summary": summary,
        "clinical_pointers": clinical_pointers
    }

def generate_2d_depiction_svg(smiles: str, width: int = 340, height: int = 220) -> str:
    """
    Generate an SVG string for 2D depiction of a molecule from SMILES string using RDKit.
    Returns empty string if invalid SMILES.
    """
    clean_s = smiles.strip()
    if not clean_s:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(clean_s)
        if not mol:
            return ""
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ""

