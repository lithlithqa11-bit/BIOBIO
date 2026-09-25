"""
Bio-Impact Analyzer - Main Application Entrypoint
=================================================
A comprehensive computational platform for:
1. Structural analysis, sequence alignment, SASA profiling, and 3D visualization.
2. Molecular docking simulation using AutoDock Vina, pocket mapping, and interaction profiling.
"""

import streamlit as st
from analysis.ui import render_analysis_tab, initialize_session_state
from docking.ui import render_docking_tab

def main():
    st.set_page_config(
        page_title="Bio-Impact Analyzer",
        page_icon="🧬",
        layout="wide"
    )
    initialize_session_state()
    st.title("🧬 Bio-Impact Analyzer")

    tab1, tab2 = st.tabs([
        "🧬 التحليل الهيكلي والمقارنة",
        "💊 محاكاة الارتباط (Docking)"
    ])

    with tab1:
        render_analysis_tab()

    with tab2:
        render_docking_tab()

if __name__ == "__main__":
    main()
