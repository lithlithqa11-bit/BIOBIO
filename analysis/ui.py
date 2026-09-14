# Streamlit User Interface for Structural Analysis and Mutation Comparison
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd

from analysis.constants import AA_3TO1
from analysis.fetch import fetch_deposited_pdb, fetch_biological_assembly, load_mutation_db
from analysis.structure import (
    process_protein_structure,
    structure_to_pdb_str,
    get_all_chains,
    get_protein_sequence,
    sequence_to_fasta,
    calculate_all_distances,
    calculate_sasa_map,
)
from analysis.alignment import get_alignment, analyze_impact
from analysis.visualization import render_protein_3d, build_sasa_figure

def initialize_session_state():
    """تهيئة متغيرات الجلسة (Session State) لضمان بقاء البيانات أثناء التفاعل مع التطبيق."""
    defaults = {
        'h_pdb': None, 'm_pdb': None,
        'h_assembly_pdb': None, 'm_assembly_pdb': None,
        'h_id': '', 'm_id': '',
        'h_results': None, 'm_results': None,
        'h_selected_chain': None, 'm_selected_chain': None,
        'h_id_in_val': '', 'm_id_in_val': ''
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val

def reset_protein_state(prefix: str, preserve_keys: list = None):
    """تنظيف متغيرات الحالة القديمة لبروتين محدد لبدء تحليل جديد ونظيف ومنع تضارب البيانات."""
    preserve = set(preserve_keys or [])
    default_preserved = {
        f"{prefix}_pdb", f"{prefix}_assembly_pdb", f"{prefix}_id",
        f"{prefix}_src", f"{prefix}_id_in_val", f"{prefix}_id_widget", f"{prefix}_up"
    }
    keep = preserve.union(default_preserved)
    for k in list(st.session_state.keys()):
        if k.startswith(f"{prefix}_") and k not in keep:
            st.session_state.pop(k, None)

def protein_ui_panel(p):
    """واجهة المستخدم لتحميل ملفات البروتين أو إدخال أكواد PDB وعرض الخيارات الخاصة بكل بروتين."""
    with p["col"]:
        icon = '🟢' if p['prefix'] == 'h' else '🔴'
        st.header(f"{icon} {p['label']}")
        source = st.radio("المصدر:", ["PDB ID", "رفع ملف"], key=f"{p['prefix']}_src", horizontal=True)

        if source == "PDB ID":
            pdb_input = st.text_input(
                "كود PDB",
                value=st.session_state[f"{p['prefix']}_id_in_val"],
                key=f"{p['prefix']}_id_widget"
            ).strip().upper()
            st.session_state[f"{p['prefix']}_id_in_val"] = pdb_input

            if st.button(f"تحميل {p['label']}", key=f"btn_{p['prefix']}"):
                with st.spinner('جاري التحميل...'):
                    # ميزة الذكاء: إذا تم تحميل البروتين المصاب، يتم جلب السليم المقابل له تلقائياً
                    if p['prefix'] == 'm':
                        mdb = load_mutation_db()
                        if pdb_input in mdb:
                            h_id = mdb[pdb_input]
                            st.session_state['h_id_in_val'] = h_id
                            st.session_state['h_id_widget'] = h_id
                            h_data = fetch_deposited_pdb(h_id)
                            if h_data:
                                st.session_state["h_pdb"] = h_data
                                st.session_state["h_assembly_pdb"] = fetch_biological_assembly(h_id, assembly_id=1)
                                st.session_state["h_id"] = h_id
                                reset_protein_state('h')
                                st.info(f"✅ تم تحميل البروتين السليم تلقائياً: {h_id}")
                    deposited_data = fetch_deposited_pdb(pdb_input)
                    assembly_data = fetch_biological_assembly(pdb_input, assembly_id=1)
                    if deposited_data:
                        st.session_state[f"{p['prefix']}_pdb"] = deposited_data
                        st.session_state[f"{p['prefix']}_assembly_pdb"] = assembly_data
                        st.session_state[f"{p['prefix']}_id"] = pdb_input
                        reset_protein_state(p['prefix'])
                        st.rerun()
                    else:
                        st.error(f"لم يتم العثور على البروتين بالكود: {pdb_input}")
        else:
            file = st.file_uploader(f"ارفع ملف {p['label']} (.pdb):", type=["pdb"], key=f"{p['prefix']}_up")
            if file:
                st.session_state[f"{p['prefix']}_pdb"] = file.getvalue().decode("utf-8")
                st.session_state[f"{p['prefix']}_assembly_pdb"] = None
                st.session_state[f"{p['prefix']}_id"] = file.name
                reset_protein_state(p['prefix'])

def render_analysis_tab():
    """عرض تبويب التحليل الهيكلي ومقارنة الطفرات بالكامل."""
    # ── إعدادات الشريط الجانبي (Sidebar) ──
    st.sidebar.header("⚙️ التحليل والإعدادات")
    with st.sidebar.expander("🎨 خيارات العرض", expanded=True):
        search_radius = st.slider("🔍 نصف قطر البحث (Å)", 3.0, 12.0, 5.0)
        view_style = st.selectbox("نمط العرض", ["cartoon", "stick", "sphere"])
        show_ligands = st.checkbox("إظهار الليجاند (Show Ligand)", value=True)
        if show_ligands:
            ligand_style = st.selectbox(
                "نمط تمثيل الليجاند (Ligand Style):",
                ["أعواد ملونة مع كرات (Sticks & Spheres)", "أعواد فقط (Sticks Only)", "كرات فضاء (Spheres)"],
                index=0
            )
        else:
            ligand_style = "إخفاء"
        show_surface = st.checkbox("إظهار السطح (Surface)", value=False)
        surface_type = st.selectbox(
            "نوع السطح الجزيئي:",
            ["MS (Molecular Surface)", "SAS (Solvent Accessible)", "VDW (Van der Waals)"],
            index=0
        )
        surface_opacity = st.slider("شفافية السطح", 0.0, 1.0, 0.3)
    
    with st.sidebar.expander("🧬 خيارات نمط الإحداثيات والتجمع الحيوي", expanded=True):
        assembly_mode_display = st.radio(
            "نمط الإحداثيات:",
            [
                "الوحدة الحيوية الوظيفية (Biological Assembly 1)",
                "الوحدة البلورية المودعة (Deposited Asymmetric Unit)"
            ],
            index=0,
            help=(
                "الوحدة الحيوية الوظيفية: الهيكل البيولوجي الفعلي المعتمد وظيفياً من RCSB (الخيار الأنسب للمقارنة).\n"
                "الوحدة البلورية المودعة: ملف الإحداثيات الأصلي المودع للبلورة (قد يشمل نسخ تبلور متعددة لدراسة البلمرة)."
            )
        )
        use_biological_assembly = "الوحدة الحيوية الوظيفية" in assembly_mode_display
        current_coordinate_mode = (
            "biological_assembly" if use_biological_assembly else "deposited_asymmetric_unit"
        )

        if st.session_state.get("last_coordinate_mode") != current_coordinate_mode:
            st.session_state["h_results"] = None
            st.session_state["m_results"] = None
            st.session_state["last_coordinate_mode"] = current_coordinate_mode

        show_mutations = st.checkbox("تلوين الطفرات", value=True)
        zoom_mutations = st.checkbox("تركيز العرض على الطفرات", value=False)
        alignment_mode = st.selectbox("نوع المحاذاة (Alignment)", ["global", "local"])

    st.sidebar.info("التحليل الهيكلي يشمل البحث في جميع سلاسل البروتين المختارة.")

    # ── الخطوة 1: لوحات إدخال البروتينات ──
    col1, col2 = st.columns(2)
    proteins = [
        {"label": "المصاب", "prefix": "m", "col": col2, "bg": "#1E0D0D"},
        {"label": "السليم", "prefix": "h", "col": col1, "bg": "#0D1B1E"}
    ]

    for p in proteins:
        protein_ui_panel(p)

    st.divider()

    # ── الخطوة 2: المعالجة واختيار السلاسل ──
    v_col1, v_col2 = st.columns(2)
    structures = {}
    
    for p in proteins:
        current_col = v_col1 if p['prefix'] == 'h' else v_col2
        with current_col:
            prefix = p['prefix']
            deposited_pdb = st.session_state.get(f"{prefix}_pdb")
            assembly_pdb = st.session_state.get(f"{prefix}_assembly_pdb")

            if not deposited_pdb and not assembly_pdb:
                continue

            if use_biological_assembly and assembly_pdb:
                active_pdb = assembly_pdb
                active_source = "الوحدة الحيوية الوظيفية (Biological Assembly 1)"
            elif use_biological_assembly and deposited_pdb:
                active_pdb = deposited_pdb
                active_source = "الوحدة البلورية المودعة (غير موثق حيوياً)"
                st.warning(
                    f"{p['label']}: لا يتوفر ملف وحدة حيوية رسمي من RCSB لهذا الإدخال. تم استخدام الإحداثيات المودعة كاملة دون حذف أي سلسلة."
                )
            else:
                active_pdb = deposited_pdb
                active_source = "الوحدة البلورية المودعة (Deposited Asymmetric Unit)"

            st.subheader(f"هيكل {p['label']}")
            st.caption(f"📁 نمط الإحداثيات النشط: **{active_source}**")
            struct = process_protein_structure(active_pdb, prefix)
            structures[prefix] = struct

            if struct:
                chains = get_all_chains(struct)
                if chains:
                    selected_chain = st.selectbox(f"اختر السلسلة - {p['label']}", options=chains, key=f"{prefix}_chain_sel")
                    if st.session_state.get(f"{prefix}_selected_chain") != selected_chain:
                        st.session_state[f"{prefix}_results"] = None
                    st.session_state[f"{prefix}_selected_chain"] = selected_chain
                else:
                    st.warning("لم يتم العثور على سلاسل في هذا الهيكل.")

    # ── الخطوة 3: حساب الطفرات تلقائياً عبر مقارنة التسلسلات ──
    highlight_map = {'h': None, 'm': None}
    h_chain = st.session_state.get('h_selected_chain')
    m_chain = st.session_state.get('m_selected_chain')
    
    alignment_data = None
    if structures.get('h') and structures.get('m') and h_chain and m_chain:
        healthy_sequence = get_protein_sequence(structures['h'], h_chain)
        mutant_sequence = get_protein_sequence(structures['m'], m_chain)
        
        if healthy_sequence and mutant_sequence:
            healthy_str = "".join([AA_3TO1.get(r['res_name'], 'X') for r in healthy_sequence])
            mutant_str = "".join([AA_3TO1.get(r['res_name'], 'X') for r in mutant_sequence])
            
            alignment_text, score, aligned_healthy, aligned_mutant = get_alignment(healthy_str, mutant_str, alignment_mode)
            alignment_data = {
                'text': alignment_text, 'score': score, 
                'aligned_healthy': aligned_healthy, 'aligned_mutant': aligned_mutant,
                'healthy_seq': healthy_sequence, 'mutant_seq': mutant_sequence,
                'healthy_str': healthy_str, 'mutant_str': mutant_str
            }
            
            mutations_healthy, mutations_mutant = [], []
            healthy_ptr, mutant_ptr = 0, 0
            
            for char_h, char_m in zip(aligned_healthy, aligned_mutant):
                h_res = healthy_sequence[healthy_ptr] if char_h != '-' else None
                m_res = mutant_sequence[mutant_ptr] if char_m != '-' else None
                
                if char_h != char_m:
                    if m_res:
                        mutations_mutant.append({'resi': str(m_res['res_num']), 'chain': str(m_chain)})
                    if h_res:
                        mutations_healthy.append({'resi': str(h_res['res_num']), 'chain': str(h_chain)})
                
                if char_h != '-':
                    healthy_ptr += 1
                if char_m != '-':
                    mutant_ptr += 1
                
            highlight_map['m'] = mutations_mutant if mutations_mutant else None
            highlight_map['h'] = mutations_healthy if mutations_healthy else None

    # ── الخطوة 4: العرض التفاعلي ثلاثي الأبعاد والتحليل الرقمي ──
    for p in proteins:
        current_col = v_col1 if p['prefix'] == 'h' else v_col2
        with current_col:
            prefix = p['prefix']
            struct = structures.get(prefix)
            selected_chain = st.session_state.get(f"{prefix}_selected_chain")

            if not struct or not selected_chain:
                continue

            highlight = highlight_map[prefix]
            focus_mut = None
            if highlight:
                mut_opts = ["الكل"] + [f"Residue {m['resi']} (Chain {m['chain']})" for m in highlight]
                sel_mut = st.selectbox(f"🔍 التركيز على طفرة - {p['label']}", mut_opts, key=f"focus_{prefix}")
                if sel_mut != "الكل":
                    sel_idx = mut_opts.index(sel_mut) - 1
                    focus_mut = highlight[sel_idx]

            pdb_for_view = structure_to_pdb_str(struct)
            surf_code = surface_type.split()[0]
            view_html = render_protein_3d(
                pdb_for_view, bg_color=p['bg'], style_type=view_style,
                show_surface=show_surface, surface_opacity=surface_opacity,
                surface_type=surf_code,
                mutations=highlight if show_mutations else None,
                mut_color='#F44336' if prefix == 'm' else '#4CAF50',
                zoom_to_mutations=zoom_mutations, focus_mut=focus_mut,
                show_ligands=show_ligands, ligand_style=ligand_style
            )
            components.html(view_html, height=460)

            try:
                total_res = sum(1 for r in struct[0].get_residues() if r.id[0] == ' ')
            except Exception:
                total_res = 0
            metric_col1, metric_col2, metric_col3 = st.columns(3)
            metric_col1.metric("عدد السلاسل", len(get_all_chains(struct)))
            metric_col2.metric("إجمالي الأحماض", total_res)
            metric_col3.metric("السلسلة الحالية", selected_chain)

            fasta = sequence_to_fasta(struct, selected_chain, st.session_state.get(f"{prefix}_id", p['label']))
            if fasta:
                with st.expander(f"🧬 تنزيل FASTA - {p['label']}"):
                    st.download_button(
                        "⬇️ تنزيل FASTA", fasta,
                        f"{st.session_state.get(f'{prefix}_id', 'protein')}_{selected_chain}.fasta",
                        "text/plain", key=f"dl_f_{prefix}"
                    )

            st.divider()
            if st.button(f"🔬 تحليل كامل لـ {p['label']}", key=f"analyze_btn_{prefix}", type="primary"):
                with st.spinner("جاري التحليل..."):
                    results = calculate_all_distances(struct, selected_chain, radius=search_radius)
                    st.session_state[f"{prefix}_results"] = results
            
            if st.session_state.get(f"{prefix}_results"):
                res = st.session_state[f"{prefix}_results"]
                with st.expander("📊 نتائج التحليل"):
                    df = pd.DataFrame(res)
                    df_display = df.rename(columns={
                        'res_num': 'رقم البقية', 'res_name': 'الحمض', 'class': 'الفئة الكيميائية',
                        'min_dist': 'أقرب مسافة', 'sasa': 'SASA'
                    })
                    st.dataframe(df_display, use_container_width=True, hide_index=True)

    # ── الخطوة 5: المقارنة المباشرة بين السليم والمصاب ──
    if structures.get('h') and structures.get('m') and h_chain and m_chain:
        st.divider()
        st.header("📋 مقارنة السلسلة (Comparison)")
        
        if alignment_data:
            healthy_sasa_map = calculate_sasa_map(structures['h'], h_chain)
            mutant_sasa_map = calculate_sasa_map(structures['m'], m_chain)

            rows = []
            healthy_pos, mutant_pos = 0, 0
            for char_h, char_m in zip(alignment_data['aligned_healthy'], alignment_data['aligned_mutant']):
                h_res = alignment_data['healthy_seq'][healthy_pos] if char_h != '-' else None
                m_res = alignment_data['mutant_seq'][mutant_pos] if char_m != '-' else None
                
                h_name = h_res['res_name'] if h_res else '-'
                m_name = m_res['res_name'] if m_res else '-'
                h_res_num = h_res['res_num'] if h_res else '-'
                m_res_num = m_res['res_num'] if m_res else '-'
                
                res_num = m_res['res_num'] if m_res else (h_res['res_num'] if h_res else 0)
                res_num_label = f"H:{h_res_num} | M:{m_res_num}" if (h_res and m_res) else (f"H:{h_res_num}" if h_res else f"M:{m_res_num}")
                
                sasa_healthy = healthy_sasa_map.get(h_res['res_num'], 0) if h_res else 0
                sasa_mutant = mutant_sasa_map.get(m_res['res_num'], 0) if m_res else 0
                
                rows.append({
                    'res_num': res_num,
                    'res_num_h': h_res_num,
                    'res_num_m': m_res_num,
                    'res_num_label': str(res_num),
                    'السليم': h_name,
                    'المصاب': m_name,
                    'SASA_H': sasa_healthy,
                    'SASA_M': sasa_mutant,
                    'SASA_Delta': round(sasa_mutant - sasa_healthy, 2),
                    'الحالة': '🔴 طفرة' if h_name != m_name else '🟢 محافظ',
                    'Impact': analyze_impact(h_name, m_name, sasa_healthy, sasa_mutant)
                })
                
                if char_h != '-':
                    healthy_pos += 1
                if char_m != '-':
                    mutant_pos += 1
            
            comparison_df = pd.DataFrame(rows)
            
            with st.expander("جدول المقارنة المتقدم (Structural & Chemical Impact)"):
                display_df = comparison_df.rename(columns={
                    'res_num_h': 'رقم السليم',
                    'res_num_m': 'رقم المصاب',
                    'SASA_H': 'SASA H',
                    'SASA_M': 'SASA M',
                    'SASA_Delta': 'ΔSASA',
                    'Impact': 'نوع التأثير العلمي'
                })[['رقم السليم', 'رقم المصاب', 'السليم', 'المصاب', 'SASA H', 'SASA M', 'ΔSASA', 'الحالة', 'نوع التأثير العلمي']]
                st.dataframe(
                    display_df.style.apply(
                        lambda r: ['background-color: #3e2723' if r['السليم'] != r['المصاب'] else ''] * len(r),
                        axis=1
                    ),
                    use_container_width=True,
                    hide_index=True
                )

            # ── رسم بياني لمقارنة SASA المتقدمة ──
            st.subheader("📊 مقارنة SASA المتقدمة")
            comparison_df['plot_index'] = range(1, len(comparison_df) + 1)
            sasa_figure = build_sasa_figure(comparison_df)
            st.plotly_chart(sasa_figure, use_container_width=True)

            # ── نتيجة المحاذاة (Alignment) ──
            st.header(f"🧬 المحاذاة التسلسلية (Alignment: {alignment_mode.capitalize()})")
            aligned_len = len(alignment_data['aligned_healthy'])
            identity_pct = (sum(a == b and a != '-' for a, b in zip(alignment_data['aligned_healthy'], alignment_data['aligned_mutant'])) / aligned_len * 100) if aligned_len > 0 else 0.0
            mutations_count = sum(a != b and a != '-' and b != '-' for a, b in zip(alignment_data['aligned_healthy'], alignment_data['aligned_mutant']))
            indels_count = sum(a != b and (a == '-' or b == '-') for a, b in zip(alignment_data['aligned_healthy'], alignment_data['aligned_mutant']))

            score_c1, score_c2, score_c3, score_c4 = st.columns(4)
            score_c1.metric("درجة المحاذاة (Score)", round(alignment_data['score'], 1))
            score_c2.metric("نسبة التطابق (Identity)", f"{identity_pct:.1f}%")
            score_c3.metric("عدد الطفرات (Mutations)", f"{mutations_count}")
            score_c4.metric("الفجوات والإزاحات (Indels)", f"{indels_count}")
            st.code(alignment_data['text'], language='text')
