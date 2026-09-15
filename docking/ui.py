# Docking Streamlit UI Module
import os
import datetime
from pathlib import Path
import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import py3Dmol

from docking.models import GridBox
from docking.storage import (
    create_run_directory,
    write_initial_manifest,
    finalize_manifest,
    get_vina_path,
    build_and_save_manifest,
    DISCLAIMER_TEXT,
)
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.binding_site import define_grid_from_ligand, define_grid_from_residues
from docking.engine import run_vina_multi_seeds, get_vina_version
from docking.validation import calculate_mapped_redocking_rmsd, cluster_poses_across_seeds, validate_redocking_top_n
from docking.rcsb import fetch_deposited_pdb_text, fetch_biological_assembly_text
from docking.comparison import run_matched_docking_comparison, COMPARISON_LABEL
from docking.interactions import (
    calculate_drug_properties,
    analyze_protein_ligand_interactions,
    compare_interactions,
    interpret_resistance_risk,
    generate_2d_depiction_svg
)
from docking.pubchem import fetch_drug_from_pubchem
import plotly.graph_objects as go
import numpy as np

def extract_chains_from_pdb(pdb_text: str) -> list[str]:
    chains = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
            ch = line[21].strip()
            if ch:
                chains.add(ch)
    return sorted(list(chains))

def render_docked_complex_3d(receptor_pdb_path: Path, pose_pdbqt_block: str, ref_ligand_pdb_text: str = "",
                              show_surface: bool = False, surface_opacity: float = 0.3, surface_type: str = "MS",
                              pocket_residues: list = None, hbonds: list = None,
                              show_ligand: bool = True, ligand_style: str = "أعواد (Sticks)",
                              grid_box: GridBox = None, show_grid_box: bool = False) -> str:
    view = py3Dmol.view(width="100%", height=480)
    if receptor_pdb_path.exists():
        with open(receptor_pdb_path, "r", encoding="utf-8") as f:
            view.addModel(f.read(), "pdb")
        view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
        if show_surface:
            surf_kind = py3Dmol.MS if surface_type == "MS" else (py3Dmol.VDW if surface_type == "VDW" else py3Dmol.SAS)
            view.addSurface(surf_kind, {"opacity": surface_opacity, "color": "#ECEFF1"}, {"model": 0})
        
        # تمييز أحماض جيب الارتباط بالأعواد (Sticks)
        if pocket_residues:
            for r in pocket_residues:
                r_clean = str(r).strip()
                if r_clean:
                    view.addStyle({"model": 0, "resi": r_clean}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.25}})

        # تمييز الأحماض المشاركة في الروابط الهيدروجينية
        if hbonds:
            for hb in hbonds:
                r_num = str(hb.get("res_num", "")).strip()
                if r_num:
                    view.addStyle({"model": 0, "resi": r_num}, {"stick": {"colorscheme": "magentaCarbon", "radius": 0.32}})

    if ref_ligand_pdb_text and show_ligand:
        view.addModel(ref_ligand_pdb_text, "pdb")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "grayCarbon", "radius": 0.25}})
    if pose_pdbqt_block and show_ligand:
        view.addModel(pose_pdbqt_block, "pdb")
        if "كرات" in ligand_style and "أعواد" not in ligand_style:
            l_style = {"sphere": {"colorscheme": "greenCarbon", "radius": 0.8}}
        elif "أعواد وكرات" in ligand_style or "Ball" in ligand_style:
            l_style = {"stick": {"colorscheme": "greenCarbon", "radius": 0.35}, "sphere": {"colorscheme": "greenCarbon", "radius": 0.5}}
        else:
            l_style = {"stick": {"colorscheme": "greenCarbon", "radius": 0.35}}
        view.setStyle({"model": -1}, l_style)
        view.zoomTo({"model": -1})
    else:
        view.zoomTo()

    # 📦 إظهار صندوق حيز الارتباط ثلاثي الأبعاد (3D Grid Box Visualizer)
    if grid_box and show_grid_box:
        try:
            view.addBox({
                "center": {"x": float(grid_box.center_x), "y": float(grid_box.center_y), "z": float(grid_box.center_z)},
                "dimensions": {"w": float(grid_box.size_x), "h": float(grid_box.size_y), "d": float(grid_box.size_z)},
                "color": "#FFD700",
                "opacity": 0.25
            })
        except Exception:
            pass

    return view._make_html()



def sync_target_binding_site_state(current_target_identity: str, session_state_dict: dict) -> bool:
    """
    Syncs and resets target-specific binding site inputs when the target protein changes.
    Returns True if a reset occurred, False otherwise.
    """
    previous_target_identity = session_state_dict.get("dock_previous_target_identity")
    if previous_target_identity != current_target_identity:
        session_state_dict["dock_ref_lig_input"] = (
            "AQ4" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_res_input"] = (
            "790, 858" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_padding"] = 8.0
        session_state_dict["dock_previous_target_identity"] = current_target_identity
        return True
    return False

def find_recovery_pose_index(poses: list, validation_details: dict | None) -> int:
    """
    Finds the index of the validated recovery pose from validation details.
    Matches both seed and rank. Returns 0 if missing, not found, or invalid.
    """
    if not poses or not validation_details:
        return 0
    recovered = validation_details.get("recovery", {})
    recovered_seed = recovered.get("seed")
    recovered_rank = recovered.get("rank")
    if recovered_seed is None or recovered_rank is None:
        return 0
    for index, pose in enumerate(poses):
        if getattr(pose, "seed", None) == recovered_seed and getattr(pose, "rank", None) == recovered_rank:
            return index
    return 0

def render_docking_tab():
    st.header("🔬 محاكاة الارتباط الجزيئي (Docking)")
    vina_ver = get_vina_version()
    st.caption(f"⚙️ المحرك: {vina_ver} | إعدادات الجزيئات: RDKit & OpenBabel")

    vina_exec = get_vina_path()
    if not vina_exec.exists():
        st.error(f"❌ لم يتم العثور على ملف المحرك Vina في المسار: {vina_exec}")
        return



    with st.container(border=True):
        tb_col1, tb_col2 = st.columns(2)
        with tb_col1:
            workflow_mode = st.radio(
                "🎯 نوع التحليل المطلوب:",
                ["single_docking", "matched_healthy_vs_mutant"],
                format_func=lambda m: "إرساء فردي (Single Docking)" if m == "single_docking" else "مقارنة السليم والمصاب (Matched Comparison)",
                horizontal=True,
                key="docking_workflow_mode"
            )
        with tb_col2:
            ui_mode = st.radio(
                "⚡ نمط الواجهة والتجربة:",
                ["⚡ النمط السريري السريع (Quick Clinical Mode)", "🔬 النمط البحثي المتقدم (Advanced Research Mode)"],
                horizontal=True,
                key="docking_interface_mode"
            )
            is_advanced_mode = "المتقدم" in ui_mode

        if is_advanced_mode:
            st.info("🔬 **النمط البحثي المتقدم مفعّل:** يتيح لك التحكم الكامل في معاملات محرك AutoDock Vina (دقة البحث، البذور، عدد الوضعيات، ونطاق الطاقة) مع استعراض تحليلي إحصائي موسّع للعناقيد والتحقق من الوضعيات.")
        else:
            st.success("⚡ **النمط السريري السريع مفعّل:** إعدادات معيارية محددة مسبقاً (دقة = 8، 3 بذور مستقلة). واجهة سلسة وموجزة تركز على التقييم الصيدلاني وقوة الارتباط المباشرة دون إرباك بالبيانات الحسابية المعقدة.")

    if workflow_mode == "single_docking":
        _render_single_docking_ui(vina_ver, is_advanced_mode=is_advanced_mode)
    else:
        _render_matched_comparison_ui(vina_ver, is_advanced_mode=is_advanced_mode)

def _render_single_docking_ui(vina_ver: str, is_advanced_mode: bool = False):
    svg_2d = ""
    d_props = None

    col_t1, col_t2 = st.columns(2)
    with col_t1:
        with st.container(border=True):
            st.markdown("#### 🏛️ 1️⃣ اختيار البروتين وجيب الارتباط")
            target_mode = st.radio(
                "مصدر البروتين:",
                ["custom_pdb", "upload_pdb", "from_tab1"],
                format_func=lambda m: {
                    "custom_pdb": "رمز PDB",
                    "upload_pdb": "رفع ملف PDB",
                    "from_tab1": "من التحليل الهيكلي"
                }[m],
                horizontal=True,
                key="dock_target_source_mode"
            )

            custom_target_id = ""
            uploaded_pdb_file = None
            selected_target_label = ""
            target_key = ""

            if target_mode == "custom_pdb":
                custom_target_id = st.text_input(
                    "كود PDB:",
                    value=st.session_state.get("dock_custom_pdb_val", ""),
                    placeholder="مثال: 1M17",
                    key="dock_custom_pdb_input"
                ).strip().upper()
                st.session_state["dock_custom_pdb_val"] = custom_target_id
                selected_target_label = custom_target_id if custom_target_id else "CUSTOM"
                target_key = "custom_pdb"
            elif target_mode == "upload_pdb":
                uploaded_pdb_file = st.file_uploader("ملف PDB:", type=["pdb"], key="dock_upload_pdb")
                custom_target_id = st.text_input("اسم البروتين:", value="CUSTOM_RECEPTOR", key="dock_upload_name").strip()
                selected_target_label = custom_target_id
                target_key = "upload_pdb"
            elif target_mode == "from_tab1":
                target_options = []
                if st.session_state.get("h_id"):
                    target_options.append(("السليم: " + st.session_state["h_id"], "h"))
                if st.session_state.get("m_id"):
                    target_options.append(("المصاب: " + st.session_state["m_id"], "m"))
                if target_options:
                    selected_target_label, target_key = st.selectbox(
                        "اختر البروتين:",
                        target_options,
                        format_func=lambda x: x[0],
                        key="dock_tab1_target_select"
                    )
                else:
                    st.info("لم يتم تحميل بروتينات في التبويب الأول بعد.")
                    target_key = "none"
                    selected_target_label = "NONE"

            # Stable identity for the selected target
            if target_mode == "custom_pdb":
                current_target_identity = f"custom:{custom_target_id}"
            elif target_mode == "upload_pdb":
                current_target_identity = (
                    f"upload:{uploaded_pdb_file.name}"
                    if uploaded_pdb_file is not None
                    else "upload:none"
                )
            elif target_mode == "from_tab1":
                current_target_identity = (
                    f"tab1:{target_key}:{st.session_state.get(f'{target_key}_id', '')}"
                )
            else:
                current_target_identity = "custom:none"

            previous_target_identity = st.session_state.get("dock_previous_target_identity")
            was_reset = sync_target_binding_site_state(current_target_identity, st.session_state)
            if was_reset and previous_target_identity is not None:
                st.info("تمت إعادة ضبط إعدادات جيب الارتباط لتغيير البروتين الهدف.")

            st.markdown("---")
            site_mode = st.radio(
                "تحديد جيب الارتباط:",
                ["co_crystal_ligand", "residue_defined"],
                format_func=lambda m: "ربيطة مرجعية (Ligand)" if m == "co_crystal_ligand" else "أحماض أمينية (Residues)",
                horizontal=True,
                key="dock_site_mode_radio"
            )

            ref_ligand_code = ""
            pocket_residues_input = ""
            is_1m17 = (target_mode == "custom_pdb" and custom_target_id == "1M17")

            if site_mode == "co_crystal_ligand":
                default_ref_lig = st.session_state.get("dock_ref_lig_input", "") or ("AQ4" if is_1m17 else st.session_state.get("dock_ref_lig_val", ""))
                ref_ligand_code = st.text_input(
                    "كود الربيطة (Ligand):",
                    value=default_ref_lig,
                    placeholder="مثال: AQ4",
                    key="dock_ref_lig_input"
                ).strip().upper()
                if not is_1m17:
                    st.session_state["dock_ref_lig_val"] = ref_ligand_code
            else:
                default_pocket_res = st.session_state.get("dock_pocket_res_input", "") or ("790, 858" if is_1m17 else st.session_state.get("dock_pocket_res_val", ""))
                pocket_residues_input = st.text_input(
                    "أرقام الأحماض الأمينية:",
                    value=default_pocket_res,
                    placeholder="مثال: 790, 858",
                    key="dock_pocket_res_input"
                )
                if not is_1m17:
                    st.session_state["dock_pocket_res_val"] = pocket_residues_input

            padding_kwargs = {
                "min_value": 4.0,
                "max_value": 14.0,
                "step": 0.5,
                "format": "%.1f Å",
                "key": "dock_pocket_padding",
            }
            if "dock_pocket_padding" not in st.session_state:
                padding_kwargs["value"] = 8.0
            pocket_padding = st.slider("توسيع الجيب (Å Padding):", **padding_kwargs)

    with col_t2:
        with st.container(border=True):
            st.markdown("#### 💊 2️⃣ مواصفات الدواء ومعاملات المحاكاة")
            # 🔍 البحث المباشر في PubChem بالاسم العلمي أو التجاري
            col_pc1, col_pc2 = st.columns([2.2, 1.0])
            with col_pc1:
                pubchem_query = st.text_input(
                    "🔍 البحث عن دواء في PubChem:",
                    placeholder="مثال: Gefitinib, Aspirin, Imatinib",
                    key="dock_pubchem_search_input"
                )
            with col_pc2:
                st.write("")
                st.write("")
                if st.button("🌐 جلب الدواء", key="btn_pubchem_fetch", use_container_width=True):
                    if pubchem_query.strip():
                        with st.spinner(f"جارٍ البحث عن '{pubchem_query}' في PubChem..."):
                            pc_res = fetch_drug_from_pubchem(pubchem_query)
                            if pc_res.get("success"):
                                st.session_state["dock_lig_name_val"] = pc_res["title"]
                                st.session_state["dock_lig_smiles_val"] = pc_res["smiles"]
                                st.session_state["dock_lig_name_input"] = pc_res["title"]
                                st.session_state["dock_lig_smiles_input"] = pc_res["smiles"]
                                st.success(f"✅ تم جلب {pc_res['title']} من {pc_res.get('source')}!")
                                st.rerun()
                            else:
                                st.error(pc_res.get("error", "فشل البحث في PubChem."))
                    else:
                        st.warning("يرجى إدخال اسم الدواء أولاً.")

            col_dn1, col_dn2 = st.columns([1, 1.6])
            with col_dn1:
                default_lig_name = "Erlotinib" if is_1m17 else st.session_state.get("dock_lig_name_val", "")
                ligand_name = st.text_input(
                    "اسم الدواء:",
                    value=default_lig_name,
                    placeholder="مثال: Erlotinib",
                    key="dock_lig_name_input"
                )
                if not is_1m17 and ligand_name:
                    st.session_state["dock_lig_name_val"] = ligand_name

            with col_dn2:
                default_smiles = (
                    "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
                    if is_1m17
                    else st.session_state.get("dock_lig_smiles_val", "")
                )
                ligand_smiles = st.text_input(
                    "صيغة SMILES:",
                    value=default_smiles,
                    placeholder="أدخل صيغة SMILES الكيميائية للدواء",
                    key="dock_lig_smiles_input"
                )
                if not is_1m17 and ligand_smiles:
                    st.session_state["dock_lig_smiles_val"] = ligand_smiles

            # 🧪 المعاينة الكيميائية ثنائية الأبعاد والخصائص الصيدلانية
            if ligand_smiles.strip():
                tab_struct, tab_lip = st.tabs(["🧪 التركيب الكيميائي 2D", "📊 معايير ليبينسكي (ADME)"])
                with tab_struct:
                    svg_2d = generate_2d_depiction_svg(ligand_smiles, width=320, height=160)
                    if svg_2d:
                        try:
                            st.image(svg_2d, caption="التركيب الكيميائي ثنائي الأبعاد (2D Depiction)", use_container_width=True)
                        except Exception:
                            st.caption("تعذر عرض المخطط ثنائي الأبعاد.")
                    else:
                        st.caption("تعذر رسم الصيغة الكيميائية ثنائية الأبعاد.")

                with tab_lip:
                    d_props = calculate_drug_properties(ligand_smiles)
                    if d_props and "error" not in d_props:
                        dp1, dp2, dp3 = st.columns(3)
                        dp1.metric("الوزن الجزيئي", f"{d_props['molecular_weight']} Da")
                        dp2.metric("LogP", f"{d_props['logp']}")
                        dp3.metric("TPSA", f"{d_props['tpsa']} Å²")
                        dp4, dp5, dp6 = st.columns(3)
                        dp4.metric("مانح H-Bond", d_props['hbd'])
                        dp5.metric("مستقبل H-Bond", d_props['hba'])
                        dp6.metric("الروابط الدوارة", d_props['rotatable_bonds'])
                        if d_props['lipinski_pass']:
                            st.caption("✅ الجزيء يتوافق تماماً مع معايير ليبينسكي للدواء الفموي (Drug-like).")
                        else:
                            st.caption(f"⚠️ {d_props['drug_likeness']}")

            st.markdown("---")
            # دقة البحث والبذور
            if is_advanced_mode:
                st.markdown("##### ⚙️ إعدادات محرك Vina المتقدمة:")
                col_cfg1, col_cfg2 = st.columns(2)
                with col_cfg1:
                    exhaustiveness = st.select_slider(
                        "دقة البحث (Exhaustiveness):",
                        options=[4, 8, 16, 32, 64],
                        value=8,
                        help="كلما زادت القيمة زادت دقة استكشاف الجيب والارتباط على حساب وقت المعالجة.",
                        key="dock_single_exh"
                    )
                    num_modes = st.slider(
                        "أقصى عدد للوضعيات (num_modes):",
                        min_value=1,
                        max_value=20,
                        value=9,
                        help="عدد أوضاع الارتباط الناتجة لكل بذرة.",
                        key="dock_single_num_modes"
                    )
                with col_cfg2:
                    seed_choice = st.selectbox(
                        "استراتيجية البذور (Seeds):",
                        ["three_seeds", "single_seed"],
                        format_func=lambda x: "3 بذور [42, 101, 2024] (تحقق إحصائي)" if x == "three_seeds" else "بذرة واحدة [42] (اختبار سريع)",
                        key="dock_single_seeds"
                    )
                    energy_range = st.slider(
                        "نطاق طاقة الوضعيات (Energy Range kcal/mol):",
                        min_value=1.0,
                        max_value=5.0,
                        value=3.0,
                        step=0.5,
                        help="الحد الأقصى لفارق الطاقة عن الوضعية الأفضل لاعتماد الوضعية.",
                        key="dock_single_energy_range"
                    )
            else:
                exhaustiveness = 8
                seed_choice = "three_seeds"
                num_modes = 9
                energy_range = 3.0
                st.caption("⚡ **إعدادات سريعة قياسية:** دقة البحث = 8 | 3 بذور عشوائية مستقلة [42, 101, 2024] لضمان التكرارية | 9 وضعيات.")

    if st.button("🚀 بدء الإرساء (Run Docking)", type="primary"):
        if site_mode == "co_crystal_ligand" and not ref_ligand_code:
            st.error("يرجى إدخال كود الربيطة المرجعية لتحديد جيب الارتباط.")
            return
        if site_mode == "residue_defined" and not pocket_residues_input.strip():
            st.error("يرجى إدخال أرقام الأحماض الأمينية لتحديد جيب الارتباط.")
            return
        if not ligand_smiles.strip():
            st.error("يرجى إدخال صيغة SMILES للدواء.")
            return

        seeds_to_run = [42, 101, 2024] if seed_choice == "three_seeds" else [42]
        created_at_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

        # 1. Coordinate retrieval and exact provenance
        if target_key == "demo_1M17":
            target_id = "1M17"
            assembly_text = fetch_biological_assembly_text("1M17", 1)
            if assembly_text:
                pdb_data = assembly_text
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = "https://files.rcsb.org/download/1M17.pdb1.gz"
            else:
                pdb_data = fetch_deposited_pdb_text("1M17")
                coord_source = "deposited"
                assembly_id = None
                source_url = "https://files.rcsb.org/download/1M17.pdb"
        elif target_key == "custom_pdb":
            if not custom_target_id:
                st.error("يرجى إدخال رمز PDB صالح.")
                return
            target_id = custom_target_id
            with st.spinner(f"جارٍ جلب بنية البروتين {target_id} من RCSB PDB..."):
                assembly_text = fetch_biological_assembly_text(target_id, 1)
                if assembly_text:
                    pdb_data = assembly_text
                    coord_source = "rcsb_biological_assembly"
                    assembly_id = 1
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
                else:
                    pdb_data = fetch_deposited_pdb_text(target_id)
                    coord_source = "deposited"
                    assembly_id = None
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            if not pdb_data or len(pdb_data.strip()) == 0:
                st.error(f"❌ تعذر جلب ملف PDB للرمز {target_id}. يرجى التحقق من صحة الرمز أو الاتصال بالإنترنت.")
                return
        elif target_key == "upload_pdb":
            if not uploaded_pdb_file:
                st.error("يرجى رفع ملف PDB من جهازك أولاً.")
                return
            target_id = custom_target_id if custom_target_id else "UPLOADED_PROTEIN"
            pdb_data = uploaded_pdb_file.getvalue().decode("utf-8", errors="ignore")
            coord_source = "user_upload"
            assembly_id = None
            source_url = uploaded_pdb_file.name
        else:
            target_id = st.session_state.get(f"{target_key}_id", "PROTEIN")
            if st.session_state.get(f"{target_key}_assembly_pdb"):
                pdb_data = st.session_state[f"{target_key}_assembly_pdb"]
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
            elif st.session_state.get(f"{target_key}_pdb"):
                pdb_data = st.session_state[f"{target_key}_pdb"]
                coord_source = "deposited"
                assembly_id = None
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            else:
                st.error("لم يتم العثور على إحداثيات للبروتين المختار.")
                return

        chains = extract_chains_from_pdb(pdb_data)
        run_id, run_dir = create_run_directory(f"dock_{target_id}")

        # Save immutable source PDB
        source_pdb_path = run_dir / "receptor_source.pdb"
        with open(source_pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_data)

        target_meta = {
            "pdb_id": target_id,
            "coordinate_source": coord_source,
            "assembly_id": assembly_id,
            "source_url": source_url,
            "chains": chains,
            "created_at_utc": created_at_utc,
            "mutation_label": st.session_state.get(f"{target_key}_mutation_label", "")
        }
        binding_site_meta = {
            "method": site_mode,
            "reference_ligand": ref_ligand_code if site_mode == "co_crystal_ligand" else None,
            "residues": pocket_residues_input if site_mode == "residue_defined" else None,
            "padding": pocket_padding
        }
        ligand_meta = {
            "name": ligand_name,
            "input_smiles": ligand_smiles,
            "protonation_assumed_ph": 7.4
        }
        engine_meta = {
            "name": "AutoDock Vina",
            "version": vina_ver,
            "seeds": seeds_to_run,
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range
        }

        # 2. Binding site definition (Prohibits origin fallback)
        ref_ligand_pdb_text = ""
        try:
            if site_mode == "co_crystal_ligand":
                grid, ref_coords, ref_ligand_pdb_text = define_grid_from_ligand(
                    pdb_data, ref_ligand_code, padding=pocket_padding
                )
                with open(run_dir / "reference_ligand.pdb", "w", encoding="utf-8") as f:
                    f.write(ref_ligand_pdb_text + "\n")
            else:
                res_specs = [r.strip() for r in pocket_residues_input.split(",") if r.strip()]
                grid = define_grid_from_residues(pdb_data, res_specs, padding=pocket_padding)
        except Exception as e:
            manifest_fail = {
                "schema_version": 1,
                "run_id": run_id,
                "status": "run failed",
                "created_at_utc": created_at_utc,
                "error": str(e),
                "target": target_meta,
                "binding_site": binding_site_meta,
                "ligand": ligand_meta,
                "engine": engine_meta,
                "warnings": [f"Binding site error: {e}"]
            }
            finalize_manifest(run_dir, manifest_fail)
            st.error(f"Docking was not run: define a validated binding site. خطأ: {e}")
            return

        grid_meta = {
            "center": [grid.center_x, grid.center_y, grid.center_z],
            "size": [grid.size_x, grid.size_y, grid.size_z]
        }

        # 3. Write initial manifest before running executables
        initial_data = {
            "schema_version": 1,
            "run_id": run_id,
            "status": "running",
            "created_at_utc": created_at_utc,
            "target": target_meta,
            "binding_site": binding_site_meta,
            "ligand": ligand_meta,
            "grid": grid_meta,
            "engine": engine_meta
        }
        write_initial_manifest(run_dir, initial_data)

        # 4. Receptor & Ligand Preparation
        with st.spinner("جاري إعداد المستقبل والدواء وتشغيل محاكي الارتباط..."):
            try:
                clean_pdb, rec_pdbqt, prep_rep = prepare_receptor(
                    pdb_data,
                    run_dir,
                    keep_waters=False,
                    reference_ligand_resname=ref_ligand_code if site_mode == "co_crystal_ligand" else None
                )
                sdf_path, lig_pdbqt, lig_rep = prepare_ligand_from_smiles(
                    ligand_smiles, ligand_name, run_dir, ph=7.4
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=ligand_meta,
                    prep_meta={},
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Preparation failure: {e}"]
                )
                st.error(f"فشل إعداد المدخلات: {e}")
                return

            # 5. Execute Multi-Seed Vina runs
            try:
                poses, runs_meta = run_vina_multi_seeds(
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=lig_pdbqt,
                    grid=grid,
                    output_dir=run_dir,
                    seeds=seeds_to_run,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    energy_range=energy_range
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=lig_rep,
                    prep_meta=prep_rep,
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Engine failure: {e}"]
                )
                st.error(f"فشل تشغيل المحرك Vina: {e}")
                return

            # 6. Cluster poses across seeds
            clusters = cluster_poses_across_seeds(poses, rmsd_threshold=2.0)

            # 7. Redocking validation (Top-N validation per REDOCKING_VALIDATION_TOP_N_REPAIR)
            VALIDATION_TOP_N = 10
            MIN_VALIDATION_SEEDS = 3
            MIN_CLUSTER_SEED_SUPPORT = 2
            RMSD_CUTOFF_ANGSTROM = 2.0

            redocking_rmsd = None
            validation_status = "exploratory—no reference validation"
            warnings = []
            validated_pose_info = {}
            mapping_details = {}
            validation_details = None

            if len(seeds_to_run) < MIN_VALIDATION_SEEDS:
                validation_status = "exploratory—insufficient independent seeds"
                warnings.append("Single-seed/insufficient seeds: pose reproducibility across independent runs cannot be verified.")

            if site_mode == "co_crystal_ligand" and ref_ligand_pdb_text:
                try:
                    validation_details = validate_redocking_top_n(
                        ref_ligand_pdb_text,
                        poses,
                        clusters,
                        top_n=VALIDATION_TOP_N,
                        rmsd_cutoff=RMSD_CUTOFF_ANGSTROM,
                        min_seed_support=MIN_CLUSTER_SEED_SUPPORT,
                    )

                    top1_info = validation_details["top1"]
                    rec_info = validation_details["recovery"]

                    redocking_rmsd = top1_info["rmsd_angstrom"]
                    mapping_details = top1_info["mapping_details"]

                    validated_pose_info = {
                        "seed": rec_info["seed"],
                        "rank": rec_info["rank"],
                        "cluster_id": rec_info["cluster_id"],
                        "output_file": f"out_seed_{rec_info['seed']}.pdbqt",
                        "score": rec_info["score"],
                        "rmsd_angstrom": rec_info["rmsd_angstrom"],
                        "cluster_seed_count": rec_info["cluster_seed_count"],
                        "supporting_seeds": rec_info["supporting_seeds"]
                    }

                    if len(seeds_to_run) < MIN_VALIDATION_SEEDS:
                        validation_status = "exploratory—insufficient independent seeds"
                    elif validation_details["recovery_rmsd_pass"] and validation_details["recovery_cluster_pass"]:
                        validation_status = "validation passed (Top-10 recovery)"
                        if not validation_details["top1_rmsd_pass"]:
                            warnings.append("Crystal pose recovered reproducibly in Top-10, but Vina did not rank it first.")
                    else:
                        validation_status = "validation failed"
                        if not validation_details["recovery_rmsd_pass"]:
                            warnings.append(f"No pose within Top-{VALIDATION_TOP_N} achieved RMSD <= {RMSD_CUTOFF_ANGSTROM} Å (best recovery: {rec_info['rmsd_angstrom']} Å).")
                        elif not validation_details["recovery_cluster_pass"]:
                            warnings.append(f"Gate A Warning: Pose RMSD <= {RMSD_CUTOFF_ANGSTROM} Å, but cluster lacks support from at least {MIN_CLUSTER_SEED_SUPPORT} independent seeds.")
                except Exception as ex:
                    validation_status = "validation unavailable"
                    warnings.append(f"Redocking calculation unavailable: {ex}")

            # Calculate global top_score across every pose from every seed
            global_top_score = min((p.score for p in poses), default=None)

            redocking_manifest_data = {
                "reference_ligand_file": "reference_ligand.pdb" if ref_ligand_pdb_text else None,
                "validated_pose": validated_pose_info,
                "rmsd_angstrom": redocking_rmsd,
                "mapping_details": mapping_details
            }
            if validation_details is not None:
                redocking_manifest_data.update({
                    "criterion": {
                        "top_n": VALIDATION_TOP_N,
                        "rmsd_cutoff_angstrom": RMSD_CUTOFF_ANGSTROM,
                        "minimum_independent_seeds": MIN_VALIDATION_SEEDS,
                        "minimum_cluster_seed_support": MIN_CLUSTER_SEED_SUPPORT
                    },
                    "top1_scored_pose": {
                        "seed": validation_details["top1"]["seed"],
                        "rank": validation_details["top1"]["rank"],
                        "score": validation_details["top1"]["score"],
                        "cluster_id": validation_details["top1"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['top1']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["top1"]["rmsd_angstrom"],
                        "mapping_details": validation_details["top1"]["mapping_details"]
                    },
                    "best_recovered_pose_in_top_n": {
                        "seed": validation_details["recovery"]["seed"],
                        "rank": validation_details["recovery"]["rank"],
                        "score": validation_details["recovery"]["score"],
                        "cluster_id": validation_details["recovery"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['recovery']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["recovery"]["rmsd_angstrom"],
                        "cluster_seed_count": validation_details["recovery"]["cluster_seed_count"],
                        "supporting_seeds": validation_details["recovery"]["supporting_seeds"],
                        "mapping_details": validation_details["recovery"]["mapping_details"]
                    },
                    "top1_rmsd_pass": validation_details["top1_rmsd_pass"],
                    "recovery_rmsd_pass": validation_details["recovery_rmsd_pass"],
                    "recovery_cluster_pass": validation_details["recovery_cluster_pass"]
                })

            results_meta = {
                "top_score": global_top_score,
                "validation_status": validation_status,
                "per_seed": runs_meta,
                "clusters": [
                    {
                        "cluster_id": cl.cluster_id,
                        "supporting_seeds": cl.supporting_seeds,
                        "seed_count": cl.seed_count,
                        "top_score": cl.top_score,
                        "median_score": cl.median_score,
                        "representative_pose_source": cl.representative_pose_source
                    }
                    for cl in clusters
                ],
                "redocking": redocking_manifest_data
            }

            # 8. Finalize Manifest atomically
            final_manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status=validation_status,
                target_meta=target_meta,
                binding_site_meta=binding_site_meta,
                ligand_meta=lig_rep,
                prep_meta=prep_rep,
                engine_meta=engine_meta,
                grid_meta=grid_meta,
                results_meta=results_meta,
                warnings=warnings
            )

            st.session_state["current_docking_result"] = {
                "run_id": run_id,
                "target_id": target_id,
                "ligand_name": ligand_name,
                "clean_pdb": clean_pdb,
                "prep_report": prep_rep,
                "poses": poses,
                "clusters": clusters,
                "grid": grid,
                "ref_ligand_pdb": ref_ligand_pdb_text,
                "redocking_rmsd": redocking_rmsd,
                "validation_status": validation_status,
                "validation_details": validation_details,
                "warnings": warnings,
                "run_dir": run_dir
            }
            st.success(f"✅ اكتمل الإرساء بنجاح! ({run_id})")

    # Display results
    res = st.session_state.get("current_docking_result")
    if res:
        st.markdown("---")
        st.subheader(f"📊 نتائج الإرساء: {res['ligand_name']} ⚯ {res['target_id']}")

        val_stat = res["validation_status"]
        val_det = res.get("validation_details")
        best_affinity = f"{res['poses'][0].score:.2f} kcal/mol" if res.get("poses") else "—"
        total_poses_count = len(res.get("poses", []))
        clusters_count = len(res.get("clusters", []))

        # Top KPI Metric Cards
        with st.container(border=True):
            kpi_col1, kpi_col2, kpi_col3, kpi_col4 = st.columns(4)
            if "passed" in val_stat:
                gate_badge = "🟢 تم التحقق بنجاح"
            elif "failed" in val_stat:
                gate_badge = "🔴 لم يجتز العتبة"
            elif "insufficient" in val_stat:
                gate_badge = "🟡 استكشافي (بذور غير كافية)"
            else:
                gate_badge = "⚪ غير متاح (دواء جديد)"

            kpi_col1.metric("🎯 حالة التحقق (Gate A)", gate_badge)
            kpi_col2.metric("⚡ أفضل طاقة ارتباط", best_affinity)
            kpi_col3.metric("🔢 عدد الوضعيات المحسوبة", f"{total_poses_count} وضعية")
            kpi_col4.metric("👥 مجموعات العناقيد", f"{clusters_count} عناقيد")

        # Clinical Summary for Quick Mode
        if not is_advanced_mode:
            best_score_val = res['poses'][0].score if res.get('poses') else 0.0
            if best_score_val <= -8.0:
                strength_label = "🟢 ارتباط قوي جداً (High Binding Affinity)"
                strength_desc = "طاقة الارتباط منخفضة جداً، مما يشير إلى ثبات واستقرار عالٍ للدواء داخل جيب البروتين وتثبيط فعال."
            elif best_score_val <= -6.0:
                strength_label = "🟡 ارتباط متوسط إلى جيد (Moderate Affinity)"
                strength_desc = "طاقة الارتباط تقع ضمن النطاق الدوائي النموذجي لمعظم المثبطات المعتمدة سريرياً."
            else:
                strength_label = "🔴 ارتباط ضعيف (Weak Affinity)"
                strength_desc = "طاقة الارتباط مرتفعة نسبياً، مما قد يشير إلى ضعف الحساسية أو مقاومة دوائية محتملة."

            with st.container(border=True):
                st.markdown(f"#### 🩺 التقييم السريري المباشر: {strength_label}")
                st.markdown(f"- **طاقة الارتباط الأفضل:** `{best_score_val:.2f} kcal/mol` — {strength_desc}")
                if d_props and "lipinski_pass" in d_props:
                    lip_msg = "✅ مطابق لقواعد ليبينسكي للدواء الفموي" if d_props["lipinski_pass"] else f"⚠️ {d_props.get('drug_likeness', 'غير مطابق لمعايير الدواء الفموي')}"
                    st.markdown(f"- **الخواص الصيدلانية الفموية:** {lip_msg}")

        # Status badge & warnings
        if val_stat == "validation passed (Top-10 recovery)":
            rec_r = val_det['recovery']['rmsd_angstrom'] if val_det else res['redocking_rmsd']
            st.success(f"🎯 Gate A Status: {val_stat} (Best Top-10 recovery RMSD = {rec_r} Å <= 2.0 Å)")
        elif val_stat == "validation passed":
            st.success(f"🎯 Gate A Status: validation passed (RMSD = {res['redocking_rmsd']} Å <= 2.0 Å)")
        elif val_stat == "validation failed":
            st.error("❌ Gate A Status: validation failed")
        elif val_stat == "validation unavailable":
            st.warning("⚠️ Gate A Status: validation unavailable (لا يتطابق عدد الذرات أو الهوية الكيميائية مع المرجع).")
        else:
            st.info(f"ℹ️ Status: {val_stat}")

        # Exact user-facing metrics: shown directly in Advanced Mode, folded in Quick Mode
        if val_det:
            t1 = val_det["top1"]
            rec = val_det["recovery"]
            val_summary_md = (
                f"- **Top-1 Vina pose RMSD:** `{t1['rmsd_angstrom']:.3f} Å` (Score: `{t1['score']} kcal/mol`)\n"
                f"- **Best Top-10 recovery RMSD:** `{rec['rmsd_angstrom']:.3f} Å` (Seed `{rec['seed']}`, Rank `{rec['rank']}`, Score: `{rec['score']} kcal/mol`)\n"
                f"- **Recovery cluster support:** `{rec['cluster_seed_count']}` independent seeds `{rec['supporting_seeds']}`"
            )
            if is_advanced_mode:
                st.markdown(val_summary_md)
            else:
                with st.expander("🔬 تفاصيل التحقق الإحصائي والحسابي المتقدم (Gate A Metrics)"):
                    st.markdown(val_summary_md)

        for w in res.get("warnings", []):
            st.warning(f"⚠️ {w}")

        # Organized Result Tabs
        tab_3d, tab_poses = st.tabs([
            "🔮 العرض ثلاثي الأبعاد والتفاعلات (3D & Interactions)",
            "📋 جدول الوضعيات والعناقيد (Poses & Clusters)"
        ])

        with tab_3d:
            # 3D Toolbar controls
            with st.container(border=True):
                c_row1, c_row2, c_row3, c_row4 = st.columns([1.5, 1.2, 1.0, 1.0])
                with c_row1:
                    default_pose_idx = find_recovery_pose_index(res["poses"], res.get("validation_details"))
                    sel_pose_idx = st.selectbox(
                        "اختر الوضعية (Pose) للعرض:",
                        range(len(res["poses"])),
                        index=default_pose_idx,
                        key=f"dock_pose_selector_{res['run_id']}",
                        format_func=lambda i: (
                            f"Pose #{res['poses'][i].rank} [Seed {res['poses'][i].seed}] "
                            f"({res['poses'][i].score} kcal/mol)"
                        ),
                    )
                with c_row2:
                    dock_show_ligand = st.checkbox("إظهار الليجاند (Show Ligand)", value=True, key="dock_single_show_ligand")
                    if dock_show_ligand:
                        dock_lig_style = st.selectbox("نمط الليجاند:", ["أعواد (Sticks)", "أعواد وكرات (Ball & Stick)", "كرات (Spheres)"], key="dock_single_lig_style")
                    else:
                        dock_lig_style = "أعواد (Sticks)"
                with c_row3:
                    dock_show_grid = st.checkbox("📦 صندوق البحث (Grid Box)", value=False, key="dock_single_show_grid")
                    dock_show_surface = st.checkbox("إظهار السطح (Surface)", value=False, key="dock_single_show_surface")
                with c_row4:
                    if dock_show_surface:
                        dock_surf_type = st.selectbox("نوع السطح:", ["MS (Molecular Surface)", "SAS", "VDW"], key="dock_single_surf_type")
                        dock_surf_opacity = st.slider("شفافية السطح:", 0.0, 1.0, 0.3, key="dock_single_surf_opacity")
                    else:
                        dock_surf_type = "MS"
                        dock_surf_opacity = 0.3

            validation_details = res.get("validation_details")
            if validation_details and "recovery" in validation_details:
                recovered = validation_details["recovery"]
                sel_p = res["poses"][sel_pose_idx]
                if getattr(sel_p, "seed", None) == recovered.get("seed") and getattr(sel_p, "rank", None) == recovered.get("rank"):
                    st.caption(
                        "🎯 الوضعية المعروضة: وضعية التحقق المسترجعة ضمن Top-10 "
                        f"(RMSD {recovered['rmsd_angstrom']:.3f} Å)."
                    )
                else:
                    st.caption("ℹ️ الوضعية المعروضة باختيار المستخدم؛ وليست بالضرورة وضعية التحقق المسترجعة.")

            # Analyze interactions for selected pose
            active_pose = res["poses"][sel_pose_idx]
            try:
                with open(res["clean_pdb"], "r", encoding="utf-8") as f:
                    rec_clean_text = f.read()
                single_inter = analyze_protein_ligand_interactions(
                    receptor_pdb_text=rec_clean_text,
                    ligand_pdbqt_block=active_pose.pdbqt_block
                )
            except Exception:
                single_inter = {"hbonds": [], "hydrophobic": [], "interacting_residues": [], "hbond_count": 0, "hydrophobic_count": 0}

            # Layout: 60% 3D Viewer, 40% Interactions Details
            col_v3d, col_vint = st.columns([1.5, 1.0])
            with col_v3d:
                st.caption("🟢 الأخضر: وضعية الدواء المحسوبة | ⚪ الرمادي: موضع الربيطة البلورية المرجعية | 🟣 الوردي: أحماض الروابط الهيدروجينية")
                active_pose_pdbqt = active_pose.pdbqt_block
                view_html = render_docked_complex_3d(
                    res["clean_pdb"],
                    active_pose_pdbqt,
                    res["ref_ligand_pdb"],
                    show_surface=dock_show_surface,
                    surface_opacity=dock_surf_opacity,
                    surface_type=dock_surf_type.split()[0],
                    hbonds=single_inter.get("hbonds"),
                    show_ligand=dock_show_ligand,
                    ligand_style=dock_lig_style,
                    grid_box=res.get("grid"),
                    show_grid_box=dock_show_grid
                )
                components.html(view_html, height=520)

            with col_vint:
                with st.container(border=True):
                    st.markdown("##### 🔗 التفاعلات الدوائية للوضعية")
                    si_col1, si_col2 = st.columns(2)
                    si_col1.metric("الروابط الهيدروجينية", single_inter["hbond_count"])
                    si_col2.metric("التماسات الكارهة للماء", single_inter["hydrophobic_count"])

                    if single_inter["hbonds"]:
                        st.markdown("###### 💧 الروابط الهيدروجينية (H-Bonds <= 3.5 Å):")
                        df_hb = pd.DataFrame(single_inter["hbonds"]).rename(columns={
                            "residue": "الحمض",
                            "chain": "السلسلة",
                            "rec_atom": "ذرة المستقبل",
                            "lig_atom": "ذرة الدواء",
                            "distance_angstrom": "المسافة (Å)"
                        })[["الحمض", "السلسلة", "ذرة المستقبل", "ذرة الدواء", "المسافة (Å)"]]
                        st.dataframe(df_hb, use_container_width=True, hide_index=True, height=160)
                    else:
                        st.info("لا توجد روابط هيدروجينية مباشرة مسافتها <= 3.5 Å في هذه الوضعية.")

                    if single_inter["interacting_residues"]:
                        st.markdown("###### 🧱 أحماض الجيب المحيطة (Contacts <= 4.0 Å):")
                        df_ir = pd.DataFrame(single_inter["interacting_residues"]).rename(columns={
                            "res_name": "الحمض",
                            "res_num": "الرقم",
                            "chain": "السلسلة",
                            "min_distance": "أقرب مسافة (Å)"
                        })[["الحمض", "الرقم", "السلسلة", "أقرب مسافة (Å)"]]
                        st.dataframe(df_ir, use_container_width=True, hide_index=True, height=160)

        with tab_poses:
            poses_data = [
                {
                    "الترتيب (Rank)": p.rank,
                    "درجة الإرساء (Vina، kcal/mol)": p.score,
                    "البذرة (Seed)": p.seed,
                    "المجموعة (Cluster)": p.cluster_id,
                    "RMSD l.b.": p.rmsd_lb,
                    "RMSD u.b.": p.rmsd_ub,
                    "RMSD للمرجع (Å)": f"{p.rmsd_to_reference:.3f}" if p.rmsd_to_reference is not None else "—"
                }
                for p in res["poses"]
            ]
            df_poses = pd.DataFrame(poses_data)

            st.markdown("#### 📋 جدول كافة الوضعيات المحسوبة (Calculated Poses)")
            st.dataframe(df_poses, use_container_width=True, hide_index=True)

            if res.get("clusters"):
                st.markdown("#### 👥 مجموعات الوضعيات والعناقيد (Pose Clusters)")
                cl_data = [
                    {
                        "المجموعة (Cluster)": cl.cluster_id,
                        "عدد البذور الداعمة": cl.seed_count,
                        "البذور": str(cl.supporting_seeds),
                        "أفضل درجة (kcal/mol)": cl.top_score,
                        "الدرجة الوسيطة (kcal/mol)": cl.median_score
                    }
                    for cl in res["clusters"]
                ]
                st.dataframe(pd.DataFrame(cl_data), use_container_width=True, hide_index=True)



def _render_matched_comparison_ui(vina_ver: str, is_advanced_mode: bool = False):
    st.subheader("⚖️ مقارنة الإرساء بين السليم والمصاب (Matched WT vs Mutant Docking)")

    if not st.session_state.get("h_pdb") or not st.session_state.get("m_pdb"):
        st.warning("يرجى تحميل بيانات البروتين السليم والمصاب أولاً من التبويب الأول (التحليل الهيكلي).")
        return

    # 1. بطاقة بروتينات المقارنة وجيب الارتباط المشترك
    with st.container(border=True):
        st.markdown("#### 🏛️ 1️⃣ المستقبلات المستهدفة وجيب الارتباط المشترك")
        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        with col_m1:
            st.metric("البروتين السليم (WT)", st.session_state.get("h_id", "WT"))
        with col_m2:
            st.metric("البروتين المصاب (Mutant)", st.session_state.get("m_id", "MUTANT"))
        with col_m3:
            default_pocket = st.session_state.get("comp_pocket_val") or st.session_state.get("dock_pocket_res_val", "790, 858")
            pocket_res = st.text_input("أحماض جيب الارتباط:", value=default_pocket, placeholder="مثال: 790, 858", key="comp_pocket_input")
        with col_m4:
            default_domain = st.session_state.get("comp_domain_val") or st.session_state.get("dock_target_domain", "Kinase Domain")
            target_domain = st.text_input("النطاق الوظيفي (Domain):", value=default_domain, key="comp_domain_input")

    # 2. تبويبات سير العمل لمقارنة الإرساء
    tab_single_comp, tab_panel_screen = st.tabs([
        "🔬 مقارنة دواء محدد (Single Drug Analysis)",
        "💊 الفرز الافتراضي لباقة أدوية (Multi-Drug Panel Screening)"
    ])

    with tab_single_comp:
        comp_svg_2d = ""
        comp_drug_props = None

        with st.container(border=True):
            st.markdown("#### 💊 2️⃣ مواصفات الدواء ومعاملات المحاكاة")
            # 🔍 البحث المباشر في PubChem
            col_cpc1, col_cpc2 = st.columns([2.2, 1.0])
            with col_cpc1:
                comp_pubchem_query = st.text_input(
                    "🔍 البحث عن دواء في PubChem:",
                    placeholder="مثال: Gefitinib, Osimertinib, Erlotinib",
                    key="comp_pubchem_search_input"
                )
            with col_cpc2:
                st.write("")
                st.write("")
                if st.button("🌐 جلب الدواء", key="btn_comp_pubchem_fetch", use_container_width=True):
                    if comp_pubchem_query.strip():
                        with st.spinner(f"جارٍ البحث عن '{comp_pubchem_query}' في PubChem..."):
                            pc_res = fetch_drug_from_pubchem(comp_pubchem_query)
                            if pc_res.get("success"):
                                st.session_state["comp_lig_name_val"] = pc_res["title"]
                                st.session_state["comp_smiles_val"] = pc_res["smiles"]
                                st.session_state["comp_lig_name_input"] = pc_res["title"]
                                st.session_state["comp_smiles_input"] = pc_res["smiles"]
                                st.success(f"✅ تم جلب {pc_res['title']}!")
                                st.rerun()
                            else:
                                st.error(pc_res.get("error", "فشل البحث في PubChem."))
                    else:
                        st.warning("يرجى إدخال اسم الدواء أولاً.")

            col_cdn1, col_cdn2 = st.columns([1, 1.6])
            with col_cdn1:
                default_l_name = st.session_state.get("comp_lig_name_val") or st.session_state.get("dock_lig_name_val", "Erlotinib")
                lig_name = st.text_input("اسم الدواء:", value=default_l_name, key="comp_lig_name_input")
            with col_cdn2:
                default_smiles = st.session_state.get("comp_smiles_val") or st.session_state.get("dock_lig_smiles_val", "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC")
                lig_smiles = st.text_input("صيغة SMILES:", value=default_smiles, key="comp_smiles_input")

            # 🧪 المعاينة الكيميائية ثنائية الأبعاد والخصائص الصيدلانية
            if lig_smiles.strip():
                tab_c_struct, tab_c_lip = st.tabs(["🧪 التركيب الكيميائي 2D", "📊 خصائص ليبينسكي (ADME)"])
                with tab_c_struct:
                    comp_svg_2d = generate_2d_depiction_svg(lig_smiles, width=320, height=160)
                    if comp_svg_2d:
                        try:
                            st.image(comp_svg_2d, caption="التركيب الكيميائي ثنائي الأبعاد (2D Depiction)", use_container_width=True)
                        except Exception:
                            st.caption("تعذر عرض المخطط ثنائي الأبعاد.")
                    else:
                        st.caption("تعذر إنشاء المخطط ثنائي الأبعاد لصيغة SMILES الحالية.")
                with tab_c_lip:
                    comp_drug_props = calculate_drug_properties(lig_smiles)
                    if comp_drug_props and "error" not in comp_drug_props:
                        cdp1, cdp2, cdp3 = st.columns(3)
                        cdp1.metric("الوزن الجزيئي", f"{comp_drug_props['molecular_weight']} Da")
                        cdp2.metric("LogP", f"{comp_drug_props['logp']}")
                        cdp3.metric("TPSA", f"{comp_drug_props['tpsa']} Å²")
                        cdp4, cdp5, cdp6 = st.columns(3)
                        cdp4.metric("مانح H-Bond", comp_drug_props['hbd'])
                        cdp5.metric("مستقبل H-Bond", comp_drug_props['hba'])
                        cdp6.metric("الروابط الدوارة", comp_drug_props['rotatable_bonds'])
                        if comp_drug_props['lipinski_pass']:
                            st.caption("✅ الجزيء يتوافق مع معايير ليبينسكي للدواء الفموي البشري.")
                        else:
                            st.caption(f"⚠️ {comp_drug_props['drug_likeness']}")

            st.markdown("---")
            if is_advanced_mode:
                st.markdown("##### ⚙️ إعدادات المقارنة المتقدمة:")
                col_mc1, col_mc2 = st.columns(2)
                with col_mc1:
                    exh = st.select_slider("دقة البحث (Exhaustiveness):", options=[4, 8, 16, 32], value=8, key="comp_exh")
                with col_mc2:
                    st.caption("🔬 يتم تطبيق نفس الإعدادات الصارمة على البروتينين (السليم والمصاب) بالتوازي عبر 3 بذور مستقلة لضمان التكافؤ الإحصائي التام.")
            else:
                exh = 8
                st.caption("⚡ **النمط السريع:** دقة البحث = 8 | فحص السليم والمصاب بالتوازي عبر 3 بذور مستقلة.")

        if st.button("🚀 تشغيل المقارنة للدواء المحدد (Run Comparison)", type="primary"):
            res_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
            if not res_list:
                st.error("Docking was not run: define a validated binding site.")
                return

            with st.spinner("جاري تنفيذ الإرساء المتطابق للسليم والمصاب عبر 3 بذور مستقلة..."):
                h_pdb = st.session_state.get("h_assembly_pdb") or st.session_state.get("h_pdb")
                m_pdb = st.session_state.get("m_assembly_pdb") or st.session_state.get("m_pdb")
                h_meta = {"pdb_id": st.session_state.get("h_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
                m_meta = {"pdb_id": st.session_state.get("m_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}

                try:
                    comp_result = run_matched_docking_comparison(
                        healthy_pdb_text=h_pdb,
                        mutant_pdb_text=m_pdb,
                        healthy_target_meta=h_meta,
                        mutant_target_meta=m_meta,
                        ligand_smiles=lig_smiles,
                        ligand_name=lig_name,
                        pocket_residues=res_list,
                        target_domain=target_domain,
                        exhaustiveness=exh
                    )
                    st.session_state["matched_comparison_result"] = comp_result
                    st.success("✅ اكتملت المقارنة بنجاح!")
                except Exception as e:
                    st.error(f"❌ فشلت المقارنة: {e}")

        comp = st.session_state.get("matched_comparison_result")
        if comp:
            st.markdown("---")
            st.subheader(f"📊 نتائج المقارنة الإرسائية: {comp['comparison_label']}")

            # ── حساب التفاعلات الصيدلانية ومقاومة الدواء ──
            h_clean_text = ""
            m_clean_text = ""
            try:
                with open(comp["healthy"]["clean_pdb"], "r", encoding="utf-8") as f:
                    h_clean_text = f.read()
                with open(comp["mutant"]["clean_pdb"], "r", encoding="utf-8") as f:
                    m_clean_text = f.read()
            except Exception:
                pass

            h_top_pose = comp["healthy"]["clusters"][0].representative_pose.pdbqt_block if comp["healthy"].get("clusters") else ""
            m_top_pose = comp["mutant"]["clusters"][0].representative_pose.pdbqt_block if comp["mutant"].get("clusters") else ""

            h_inter = analyze_protein_ligand_interactions(h_clean_text, h_top_pose)
            m_inter = analyze_protein_ligand_interactions(m_clean_text, m_top_pose)

            pocket_nums = [r.strip() for r in pocket_res.split(",") if r.strip()] if pocket_res else []
            comp_inter = compare_interactions(h_inter, m_inter, mutation_residue_nums=pocket_nums)
            resistance_info = interpret_resistance_risk(
                delta_score=comp["delta_score"],
                lost_hbonds_count=len(comp_inter["lost_hb_residues"]),
                mutation_engaged=comp_inter["mutation_engaged_in_binding"]
            )

            # Top KPI Metric Cards
            with st.container(border=True):
                kpi_c1, kpi_c2, kpi_c3, kpi_c4 = st.columns(4)
                kpi_c1.metric("السليم (Healthy WT)", f"{comp['healthy']['dominant_cluster_median']} kcal/mol")
                kpi_c2.metric("المصاب (Mutant MT)", f"{comp['mutant']['dominant_cluster_median']} kcal/mol")
                kpi_c3.metric("فارق الدرجة (ΔScore)", f"{comp['delta_score']:+.2f} kcal/mol")
                kpi_c4.metric("تقييم المقاومة", resistance_info["risk_level"].split()[0] + " " + resistance_info["risk_level"].split()[1] if len(resistance_info["risk_level"].split()) > 1 else resistance_info["risk_level"])

            col_s1, col_s2 = st.columns(2)
            with col_s1:
                st.write(f"**بذور السليم الداعمة:** {comp['healthy']['dominant_cluster_seeds']} ({comp['healthy']['seed_count']})")
            with col_s2:
                st.write(f"**بذور المصاب الداعمة:** {comp['mutant']['dominant_cluster_seeds']} ({comp['mutant']['seed_count']})")

            st.caption("ملاحظة: ΔScore = الدرجة الوسيطة للمصاب - الدرجة الوسيطة للسليم (القيم السالبة الأكبر تعني تقارباً أقوى).")

            # تفاصيل نتائج المقارنة في تبويبات منظمة
            comp_tab_3d, comp_tab_resist = st.tabs([
                "🔮 المقارنة ثلاثية الأبعاد المتزامنة (WT vs MT 3D)",
                "🔬 تقييم مقاومة الدواء والتفاعلات (Resistance & H-Bonds)"
            ])

            with comp_tab_3d:
                with st.container(border=True):
                    c_col1, c_col2, c_col3 = st.columns(3)
                    with c_col1:
                        comp_show_ligand = st.checkbox("إظهار الليجاند (Show Ligand)", value=True, key="dock_comp_show_ligand")
                        if comp_show_ligand:
                            comp_lig_style = st.selectbox("نمط الليجاند:", ["أعواد (Sticks)", "أعواد وكرات (Ball & Stick)", "كرات (Spheres)"], key="dock_comp_lig_style")
                        else:
                            comp_lig_style = "أعواد (Sticks)"
                    with c_col2:
                        comp_show_grid = st.checkbox("📦 صندوق البحث (Grid Box)", value=False, key="dock_comp_show_grid")
                        comp_show_surface = st.checkbox("إظهار السطح (Surface)", value=False, key="dock_comp_show_surface")
                    with c_col3:
                        if comp_show_surface:
                            comp_surf_type = st.selectbox("نوع السطح:", ["MS (Molecular Surface)", "SAS", "VDW"], key="dock_comp_surf_type")
                            comp_surf_opacity = st.slider("شفافية السطح:", 0.0, 1.0, 0.3, key="dock_comp_surf_opacity")
                        else:
                            comp_surf_type = "MS"
                            comp_surf_opacity = 0.3

                v_col1, v_col2 = st.columns(2)
                with v_col1:
                    st.markdown(f"**🟢 السليم (Healthy): {comp['healthy']['target_id']}**")
                    st.caption("🟣 الأحماض الوردية: أحماض الروابط الهيدروجينية المكتشفة بالسليم")
                    h_view_html = render_docked_complex_3d(
                        comp['healthy']['clean_pdb'],
                        h_top_pose,
                        "",
                        show_surface=comp_show_surface,
                        surface_opacity=comp_surf_opacity,
                        surface_type=comp_surf_type.split()[0],
                        pocket_residues=pocket_nums,
                        hbonds=h_inter.get("hbonds"),
                        show_ligand=comp_show_ligand,
                        ligand_style=comp_lig_style,
                        grid_box=comp['healthy'].get('grid'),
                        show_grid_box=comp_show_grid
                    )
                    components.html(h_view_html, height=480)

                with v_col2:
                    st.markdown(f"**🔴 المصاب (Mutant): {comp['mutant']['target_id']}**")
                    st.caption("🟣 الأحماض الوردية: أحماض الروابط الهيدروجينية بالمصاب | 🟠 البرتقالية: موقع أحماض الجيب")
                    m_view_html = render_docked_complex_3d(
                        comp['mutant']['clean_pdb'],
                        m_top_pose,
                        "",
                        show_surface=comp_show_surface,
                        surface_opacity=comp_surf_opacity,
                        surface_type=comp_surf_type.split()[0],
                        pocket_residues=pocket_nums,
                        hbonds=m_inter.get("hbonds"),
                        show_ligand=comp_show_ligand,
                        ligand_style=comp_lig_style,
                        grid_box=comp['mutant'].get('grid'),
                        show_grid_box=comp_show_grid
                    )
                    components.html(m_view_html, height=480)

            with comp_tab_resist:
                if resistance_info["risk_class"] == "danger":
                    st.error(f"### {resistance_info['risk_level']}")
                elif resistance_info["risk_class"] == "warning":
                    st.warning(f"### {resistance_info['risk_level']}")
                elif resistance_info["risk_class"] == "success":
                    st.success(f"### {resistance_info['risk_level']}")
                else:
                    st.info(f"### {resistance_info['risk_level']}")

                st.write(resistance_info["summary"])
                for ptr in resistance_info["clinical_pointers"]:
                    st.markdown(f"- 📌 {ptr}")

                with st.container(border=True):
                    rc_col1, rc_col2, rc_col3 = st.columns(3)
                    rc_col1.metric("روابط السليم الهيدروجينية", comp_inter["healthy_hbond_count"])
                    rc_col2.metric("روابط المصاب الهيدروجينية", comp_inter["mutant_hbond_count"])
                    rc_col3.metric("الروابط المفقودة في المصاب", len(comp_inter["lost_hb_residues"]))

                    if comp_inter["lost_hb_residues"]:
                        st.error(f"⚠️ **الروابط الهيدروجينية المفقودة في البروتين الطافر:** {', '.join(comp_inter['lost_hb_residues'])}")
                    if comp_inter["gained_hb_residues"]:
                        st.success(f"✨ **الروابط الهيدروجينية الجديدة في البروتين الطافر:** {', '.join(comp_inter['gained_hb_residues'])}")
                    if comp_inter["conserved_hb_residues"]:
                        st.info(f"🛡️ **الروابط الهيدروجينية المحفوظة:** {', '.join(comp_inter['conserved_hb_residues'])}")

                    col_tb1, col_tb2 = st.columns(2)
                    with col_tb1:
                        st.markdown("**🟢 روابط السليم الهيدروجينية (WT H-Bonds):**")
                        if h_inter["hbonds"]:
                            df_h_hb = pd.DataFrame(h_inter["hbonds"]).rename(columns={
                                "residue": "الحمض", "chain": "السلسلة", "rec_atom": "ذرة المستقبل", "lig_atom": "ذرة الدواء", "distance_angstrom": "المسافة (Å)"
                            })[["الحمض", "السلسلة", "ذرة المستقبل", "ذرة الدواء", "المسافة (Å)"]]
                            st.dataframe(df_h_hb, use_container_width=True, hide_index=True)
                        else:
                            st.caption("لا توجد روابط هيدروجينية بمسافة <= 3.5 Å.")

                    with col_tb2:
                        st.markdown("**🔴 روابط المصاب الهيدروجينية (Mutant H-Bonds):**")
                        if m_inter["hbonds"]:
                            df_m_hb = pd.DataFrame(m_inter["hbonds"]).rename(columns={
                                "residue": "الحمض", "chain": "السلسلة", "rec_atom": "ذرة المستقبل", "lig_atom": "ذرة الدواء", "distance_angstrom": "المسافة (Å)"
                            })[["الحمض", "السلسلة", "ذرة المستقبل", "ذرة الدواء", "المسافة (Å)"]]
                            st.dataframe(df_m_hb, use_container_width=True, hide_index=True)
                        else:
                            st.caption("لا توجد روابط هيدروجينية بمسافة <= 3.5 Å.")



    # 3. تبويب الفرز الافتراضي متعدد الأدوية
    with tab_panel_screen:
        with st.container(border=True):
            st.markdown("#### 💊 الفرز الافتراضي والمقارنة متعددة الأدوية (Multi-Drug Virtual Screening)")
            st.markdown("مقارنة عدة أدوية سريرية ضد الطفرة في جولة فحص واحدة لتحديد الدواء الأكثر استقراراً والأقل عرضة للمقاومة:")
            screening_panels = {
                "EGFR TKI Panel (أدوية سرطان الرئة)": [
                    ("Erlotinib (جيل 1)", "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"),
                    ("Gefitinib (جيل 1)", "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"),
                    ("Osimertinib (جيل 3)", "CN1CCN(CC1)C2=CC(=C(C=C2)NC3=NC=CC(=N3)C4=CN(C5=CC=CC=C54)C)NC(=O)C=C")
                ],
                "BCR-ABL Inhibitors (أدوية سرطان الدم CML)": [
                    ("Imatinib (جيل 1)", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"),
                    ("Dasatinib (جيل 2)", "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO")
                ],
                "Clinical Antivirals (مضادات الفيروسات)": [
                    ("Nirmatrelvir (Paxlovid)", "CC(C)(C)C(C(=O)NC(C1CC1)C#N)NC(=O)C2CC3(C2)C(C(O3)(F)F)C(=O)O"),
                    ("Remdesivir", "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4")
                ]
            }
            chosen_panel_name = st.selectbox("اختر باقة الفرز الصيدلانية:", list(screening_panels.keys()), key="screen_panel_select")
            chosen_drugs = screening_panels[chosen_panel_name]

            if st.button("🚀 تشغيل الفرز الافتراضي للباقة (Screen Drug Panel)", key="btn_run_panel_screen", type="primary"):
                screen_pocket_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
                if not screen_pocket_list:
                    st.error("يرجى تحديد أحماض جيب الارتباط أولاً في البطاقة العلوية.")
                else:
                    progress_bar = st.progress(0.0)
                    screen_results = []
                    h_p_data = st.session_state.get("h_assembly_pdb") or st.session_state.get("h_pdb")
                    m_p_data = st.session_state.get("m_assembly_pdb") or st.session_state.get("m_pdb")
                    h_meta_s = {"pdb_id": st.session_state.get("h_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
                    m_meta_s = {"pdb_id": st.session_state.get("m_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}

                    for idx, (d_name, d_smiles) in enumerate(chosen_drugs):
                        try:
                            cr = run_matched_docking_comparison(
                                healthy_pdb_text=h_p_data,
                                mutant_pdb_text=m_p_data,
                                healthy_target_meta=h_meta_s,
                                mutant_target_meta=m_meta_s,
                                ligand_smiles=d_smiles,
                                ligand_name=d_name,
                                pocket_residues=screen_pocket_list,
                                target_domain=target_domain,
                                exhaustiveness=4
                            )
                            wt_s = cr['healthy']['dominant_cluster_median']
                            mt_s = cr['mutant']['dominant_cluster_median']
                            delta_s = cr['delta_score']
                            # Classification
                            if delta_s >= 1.5:
                                r_tag = "🔴 مقاومة مرتفعة (Resistant)"
                            elif delta_s >= 0.6:
                                r_tag = "🟠 تراجع الفعالية (Moderate)"
                            elif delta_s <= -1.5:
                                r_tag = "🟢 ارتباط انتقائي بالطافر (Highly Effective)"
                            else:
                                r_tag = "⚪ فعالية محفوظة (Neutral)"

                            screen_results.append({
                                "الدواء": d_name,
                                "السليم (WT) kcal/mol": wt_s,
                                "المصاب (Mutant) kcal/mol": mt_s,
                                "فارق الطاقة (ΔScore)": delta_s,
                                "التقييم السريري": r_tag
                            })
                        except Exception as ex:
                            st.warning(f"تعذر فرز {d_name}: {ex}")
                        progress_bar.progress((idx + 1) / len(chosen_drugs))

                    if screen_results:
                        st.success("✅ اكتمل الفرز الافتراضي بنجاح!")
                        df_screen = pd.DataFrame(screen_results)
                        st.dataframe(df_screen, use_container_width=True, hide_index=True)

                        # Plotly Grouped Bar Chart
                        fig_screen = go.Figure()
                        d_names = [r["الدواء"] for r in screen_results]
                        wt_vals = [r["السليم (WT) kcal/mol"] for r in screen_results]
                        mt_vals = [r["المصاب (Mutant) kcal/mol"] for r in screen_results]

                        fig_screen.add_trace(go.Bar(
                            x=d_names, y=wt_vals, name="السليم (WT)",
                            marker_color="#2ea043"
                        ))
                        fig_screen.add_trace(go.Bar(
                            x=d_names, y=mt_vals, name="المصاب (Mutant)",
                            marker_color="#da3633"
                        ))
                        fig_screen.update_layout(
                            title="مقارنة طاقة الارتباط الحاسوبية عبر باقة الأدوية (kcal/mol - الأقل أفضل)",
                            barmode="group",
                            template="plotly_dark",
                            height=380,
                            margin=dict(l=20, r=20, t=40, b=20)
                        )
                        st.plotly_chart(fig_screen, use_container_width=True)
