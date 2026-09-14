# Human Pharmaceutical & Clinical Presets Library
"""
قاعدة بيانات الحالات السريرية البشرية الموثقة في تصميم الأدوية ومقاومة الطفرات.
مخصصة للصيادلة والباحثين في الطب البشري وعلم الأورام وأمراض الدم والفيروسات.
"""

HUMAN_CLINICAL_PRESETS = {
    "egfr_erlotinib_t790m": {
        "title": "🎗️ سرطان الرئة: مستقبل EGFR مع عقار Erlotinib (طفرة المقاومة T790M)",
        "category": "الأورام (Oncology / NSCLC)",
        "target_name": "EGFR Tyrosine Kinase Domain",
        "healthy_id": "1M17",
        "mutant_id": "2JIT",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Erlotinib (Tarceva)",
        "drug_smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
        "ref_ligand": "AQ4",
        "pocket_residues": "790, 793, 858",
        "target_domain": "Kinase Domain",
        "clinical_note": (
            "طفرة T790M (Gatekeeper Mutation) تُغير الثريونين الصغير إلى ميثيونين أكبر حجماً، "
            "مما يُحدث إعاقة فراغية (Steric Hindrance) تمنع أدوية الجيل الأول مثل Erlotinib من الارتباط، "
            "مسببة مقاومة سريرية لدى مرضى سرطان الرئة غير صغير الخلايا (NSCLC)."
        )
    },
    "egfr_osimertinib_t790m": {
        "title": "💊 سرطان الرئة: مستقبل EGFR مع عقار Osimertinib (الجيل الثالث الموجه للطفرة)",
        "category": "الأورام (Oncology / NSCLC)",
        "target_name": "EGFR T790M Active Conformation",
        "healthy_id": "1M17",
        "mutant_id": "2JIT",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Osimertinib (Tagrisso)",
        "drug_smiles": "CC(=O)NC1=C(C=CC(=C1)NC2=NC=C(C(=N2)C3=CN(C4=CC=CC=C34)C)N(C)CCN(C)C)NC(=O)C=C",
        "ref_ligand": "AQ4",
        "pocket_residues": "790, 793, 797, 858",
        "target_domain": "Kinase Domain",
        "clinical_note": (
            "عقار Osimertinib صُمم خصيصاً كعلاج موجه ليتفادى الإعاقة الفراغية لطفرة T790M "
            "ويرتبط تساهمياً بالحامض Cys797، ليعيد التحكم بالورم بعد فشل أدوية الجيل الأول."
        )
    },
    "sickle_cell_voxelotor": {
        "title": "🩸 فقر الدم المنجلي: الهيموجلوبين السليم مقابل المنجلي مع دواء Voxelotor",
        "category": "أمراض الدم (Hematology)",
        "target_name": "Human Hemoglobin Tetramer (HbA vs HbS)",
        "healthy_id": "2HHB",
        "mutant_id": "2HBS",
        "healthy_chain": "B",
        "mutant_chain": "B",
        "drug_name": "Voxelotor (Oxbryta)",
        "drug_smiles": "CC1=C(C(=NO1)C)C2=NC=C(C=C2)OCC3=C(C(=CC=C3)O)C=O",
        "ref_ligand": "",
        "pocket_residues": "1, 2, 3, 6",
        "target_domain": "Hemoglobin Beta Subunit",
        "clinical_note": (
            "طفرة Glu6Val في سلسلة بيتا تجعل الهيموجلوبين يتبلمر في حالات نقص الأكسجين مشكلاً خلايا منجلية. "
            "عقار Voxelotor يرتبط برأس السلسلة ويزيد تقارب الهيموجلوبين للأكسجين مانعاً البلمرة والتمنجل."
        )
    },
    "cml_bcr_abl_imatinib": {
        "title": "🎗️ سرطان الدم CML: إنزيم BCR-ABL1 مع عقار Imatinib (طفرة المقاومة T315I)",
        "category": "أمراض الدم والأورام (Hematology / Oncology)",
        "target_name": "BCR-ABL1 Tyrosine Kinase",
        "healthy_id": "1IEP",
        "mutant_id": "2GQG",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Imatinib (Gleevec)",
        "drug_smiles": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        "ref_ligand": "STI",
        "pocket_residues": "315, 317, 381, 382",
        "target_domain": "Abl Kinase Domain",
        "clinical_note": (
            "طفرة T315I تُعرف بطفرة بوابة الكيناز (Gatekeeper)، حيث يُستبدل الثريونين بالآيزولوسين "
            "مما يؤدي إلى فقدان رابطة هيدروجينية أساسية وحدوث إعاقة فراغية تطرد عقار Imatinib."
        )
    },
    "melanoma_braf_vemurafenib": {
        "title": "☀️ سرطان الجلد الميلانوما: إنزيم BRAF مع عقار Vemurafenib (طفرة V600E)",
        "category": "الأورام (Oncology / Melanoma)",
        "target_name": "BRAF Kinase Domain",
        "healthy_id": "4MNE",
        "mutant_id": "3OG7",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Vemurafenib (Zelboraf)",
        "drug_smiles": "CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F",
        "ref_ligand": "032",
        "pocket_residues": "509, 595, 600",
        "target_domain": "BRAF Kinase Domain",
        "clinical_note": (
            "طفرة V600E تُحاكي حالة الفسفرة النشطة لإنزيم BRAF وتجعله نشطاً بشكل دائم محفزاً نمو الورم. "
            "عقار Vemurafenib يثبط الإنزيم النشط انتقائياً ويحقق استجابة سريرية سريعة لدى مرضى الميلانوما."
        )
    },
    "covid_mpro_paxlovid": {
        "title": "🦠 الفيروسات التنفسية: بروتياز كورونا مع دواء Paxlovid (Nirmatrelvir)",
        "category": "الأمراض المعدية والفيروسات (Infectious Diseases)",
        "target_name": "SARS-CoV-2 Main Protease (Mpro)",
        "healthy_id": "6LU7",
        "mutant_id": "8D4N",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Nirmatrelvir (Paxlovid)",
        "drug_smiles": "CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C",
        "ref_ligand": "N3",
        "pocket_residues": "140, 142, 144, 145, 163, 166",
        "target_domain": "Chymotrypsin-like Protease",
        "clinical_note": (
            "إنزيم Mpro أساسي لتضاعف الفيروس عبر تقطيع متعدد البروتين. "
            "يرتبط عقار Nirmatrelvir بـ Cys145 و His41 مثبطاً عمل الإنزيم بصورة فعالة لتقليل شدة المرض."
        )
    }
}
