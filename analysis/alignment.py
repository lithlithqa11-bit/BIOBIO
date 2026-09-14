# Sequence alignment and mutation biochemical impact classification
from Bio.Align import PairwiseAligner
import streamlit as st
from analysis.constants import AA_PROPS, AA_1TO3, AA_VOLUMES

def analyze_impact(h_res, m_res, h_sasa, m_sasa):
    """تحليل التأثير العلمي للطفرة بناءً على الخواص الكيميائية ومساحة السطح (SASA)."""
    # التعامل مع حالات الحذف والإضافة (Indels) بشكل صحيح بيولوجياً
    if h_res == '-' and m_res != '-':
        return "Residue Insertion (Critical Structural Shift)"
    if h_res != '-' and m_res == '-':
        return "Residue Deletion (Critical Structural Shift)"
    
    # تحويل الأكواد الأحادية إلى ثلاثية إن وجدت
    h_code = AA_1TO3.get(h_res.upper(), h_res.upper()) if len(h_res) == 1 else h_res.upper()
    m_code = AA_1TO3.get(m_res.upper(), m_res.upper()) if len(m_res) == 1 else m_res.upper()

    # إذا كان الحَمضان متطابقين تماماً فلا توجد طفرة
    if h_code == m_code:
        return "Identical (No Mutation)"
    
    h_type = AA_PROPS.get(h_code, 'Unknown')
    m_type = AA_PROPS.get(m_code, 'Unknown')
    
    impacts = []
    
    # 1. تحليل الشحنة الكهربائية والتغير الكيميائي
    is_h_acidic = "Acidic" in h_type
    is_m_acidic = "Acidic" in m_type
    is_h_basic = "Basic" in h_type
    is_m_basic = "Basic" in m_type
    is_h_charged = is_h_acidic or is_h_basic
    is_m_charged = is_m_acidic or is_m_basic

    if (is_h_acidic and is_m_basic) or (is_h_basic and is_m_acidic):
        impacts.append("Charge Flip (Critical)")
    elif is_h_charged and not is_m_charged:
        impacts.append("Charge Loss (Electrostatic Disruption)")
    elif not is_h_charged and is_m_charged:
        impacts.append("Charge Gain (Electrostatic Shift)")
    elif h_type != m_type:
        impacts.append(f"Chem-Class Change ({h_type} ➔ {m_type})")

    # 2. تغير الأروماتية (Aromaticity Change: Pi-Stacking)
    is_h_arom = "Aromatic" in h_type
    is_m_arom = "Aromatic" in m_type
    if is_h_arom != is_m_arom:
        impacts.append("Aromaticity Change (Pi-Stacking Risk)")
    
    # 3. روابط ثنائي الكبريت (Disulfide Bond Risk)
    if h_code == 'CYS' or m_code == 'CYS':
        impacts.append("Disulfide Bond Risk (Critical)")
    
    # 4. تأثير البرولين وتكسير الحلزونات (Helix Breaker)
    if (h_code == 'PRO') != (m_code == 'PRO'):
        impacts.append("Proline Substitution (Helix Breaker)")
    
    # 5. مرونة الجلايسين (Glycine Flexibility)
    if (h_code == 'GLY') != (m_code == 'GLY'):
        impacts.append("Glycine Flexibility Shift")

    # 6. التغير الحجمي الفراغي (Steric Volume Change)
    h_vol = AA_VOLUMES.get(h_code)
    m_vol = AA_VOLUMES.get(m_code)
    if h_vol is not None and m_vol is not None:
        delta_vol = m_vol - h_vol
        if delta_vol >= 55.0:
            impacts.append(f"Steric Bulk Increase (+{delta_vol:.0f} Å³ Clash Risk)")
        elif delta_vol <= -55.0:
            impacts.append(f"Steric Void Creation ({delta_vol:.0f} Å³)")

    # 7. تحليل مساحة السطح والموقع الهيكلي (SASA & Structural Environment)
    try:
        diff_sasa = m_sasa - h_sasa
        if h_sasa < 15.0:
            # الحمض مدفون في لب البروتين (Core residue)
            if is_m_charged or (h_vol and m_vol and abs(m_vol - h_vol) >= 50):
                impacts.append("Buried Core Mutation (High Destabilization Risk)")
            else:
                impacts.append("Core Mutation")

        if diff_sasa >= 15.0:
            impacts.append("Increased Exposure (+ΔSASA)")
        elif diff_sasa <= -15.0:
            impacts.append("Increased Burial (-ΔSASA)")
    except Exception:
        pass
        
    if not impacts:
        if h_type == m_type:
            return "Conservative Substitution (Same Chemical Class)"
        return "Minor Physicochemical Change"

    return " | ".join(impacts)

def get_alignment(seq1, seq2, mode='global'):
    """إجراء محاذاة تسلسلية (Sequence Alignment) بين بروتينين وحساب النتيجة باستخدام خوارزميات Biopython مع مصفوفة BLOSUM62."""
    aligner = PairwiseAligner()
    aligner.mode = mode
    try:
        from Bio.Align import substitution_matrices
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
    except Exception:
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -1.0

    try:
        best_aln = aligner.align(seq1, seq2)[0]
        # استخراج سلاسل المحاذاة مع الفجوات (-) بشكل دقيق من الكائن مباشرة
        aligned_h = best_aln[0, :]
        aligned_m = best_aln[1, :]
        return str(best_aln), best_aln.score, aligned_h, aligned_m
    except Exception as error:
        st.warning(f"فشل إجراء المحاذاة التسلسلية: {error}")
        return "", 0, "", ""
