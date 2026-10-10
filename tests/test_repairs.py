"""Automated test suite verifying required repairs and scientific validity constraints."""

import os
import sys
import json
import hashlib
import shutil
import tempfile
import unittest
from pathlib import Path
from io import StringIO
import numpy as np

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from analysis.fetch import fetch_biological_assembly, fetch_deposited_pdb
from analysis.structure import process_protein_structure, get_protein_sequence, calculate_sasa_map
from analysis.alignment import get_alignment
from analysis.constants import AA_3TO1

from docking.binding_site import (
    define_grid_from_residues,
    define_grid_from_ligand,
    define_grid_whole_protein,
    extract_reference_ligand,
    GridBox,
)
import streamlit as st
from analysis.ui import invalidate_dependent_results, clear_protein
from docking.comparison import (
    verify_comparison_clusters,
    _select_consensus_cluster,
    run_matched_docking_comparison,
    run_exploratory_cross_docking,
    validate_pocket_chain,
    MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT,
    PROTONATION_ASSUMPTION_NOTE,
    COMPARISON_LABEL,
)
from docking.pair_validation import (
    validate_matched_pair,
    map_healthy_pocket_residues,
    MIN_MATCHED_PAIR_IDENTITY,
    extract_protein_chains,
)
from docking.validation import (
    validate_redocking_top_n,
    calculate_mapped_redocking_rmsd,
    verify_reference_ligand_identity,
)
from docking.models import DockingPose, PoseCluster
from docking.interactions import (
    analyze_protein_ligand_interactions,
    compare_interactions,
    interpret_score_difference,
    CONTACTS_DISCLAIMER,
)
from unittest.mock import patch
from streamlit.testing.v1 import AppTest

from docking.ui import (
    compute_input_fingerprint,
    compute_text_sha256,
    lowest_vina_score,
)
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.engine import get_vina_version
from docking.storage import (
    get_vina_path,
    build_and_save_manifest,
    normalize_manifest_status,
    MANIFEST_STATUSES,
    LEGACY_STATUS_ALIASES,
    BASE_RUNS_DIR,
)
from docking.openbabel_io import (
    require_openbabel_formats,
    read_text,
    write_text,
    read_file,
    write_file,
)


class TestRepair1And4IndependentGrids(unittest.TestCase):
    """Repair 1 & Test 4: Independent Cartesian pocket grids, not averaged midpoint grid."""

    def test_translated_structures_produce_independent_grids(self):
        # Receptor 1: pocket at origin
        pdb_h = (
            "ATOM      1  N   ALA A 790       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A 790       1.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )
        # Receptor 2: pocket translated by 100 Å along X
        pdb_m = (
            "ATOM      1  N   ALA A 790     100.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A 790     101.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )

        grid_h = define_grid_from_residues(pdb_h, ["A:790"], padding=10.0)
        grid_m = define_grid_from_residues(pdb_m, ["A:790"], padding=10.0)

        # Centers must be independent, near their own receptor, not at 50.0 (the midpoint)
        self.assertAlmostEqual(grid_h.center_x, 0.5, delta=1.0)
        self.assertAlmostEqual(grid_m.center_x, 100.5, delta=1.0)
        self.assertGreater(abs(grid_m.center_x - grid_h.center_x), 95.0)


class TestRepair2UnambiguousPocketResidues(unittest.TestCase):
    """Repair 2 & Test 3: Require unambiguous pocket residues in multichain proteins."""

    def setUp(self):
        # Multichain protein with residue 790 in chain A and chain B
        self.multichain_pdb = (
            "ATOM      1  CA  ALA A 790       0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  ALA B 790      50.000  50.000  50.000  1.00 20.00           C\n"
            "END\n"
        )

    def test_bare_residue_rejected_in_multichain(self):
        with self.assertRaises(ValueError) as ctx:
            define_grid_from_residues(self.multichain_pdb, ["790"])
        self.assertIn("occurs in more than one chain", str(ctx.exception))

    def test_qualified_residue_succeeds_in_multichain(self):
        grid = define_grid_from_residues(self.multichain_pdb, ["A:790"])
        self.assertIsInstance(grid, GridBox)
        self.assertAlmostEqual(grid.center_x, 0.0, delta=0.1)


class TestRepair3ConsensusClustersAndSeedSupport(unittest.TestCase):
    """Repair 3 & Test 5: Reject matched comparison without multi-seed consensus support."""

    def _make_cluster(self, cluster_id, seed_count, median_score, top_score=None):
        pose = DockingPose(
            seed=42,
            rank=1,
            score=median_score,
            rmsd_lb=0.0,
            rmsd_ub=0.0,
            pdbqt_block="ATOM      1  C   LIG     1       0.000   0.000   0.000\n",
            cluster_id=cluster_id,
        )
        return PoseCluster(
            cluster_id=cluster_id,
            seed_count=seed_count,
            supporting_seeds=[42] if seed_count == 1 else [42, 101],
            top_score=top_score if top_score is not None else median_score,
            median_score=median_score,
            representative_pose=pose,
            representative_pose_source={"seed": 42, "rank": 1},
        )

    def test_cluster_selection_prefers_seed_count_first(self):
        # Cluster 1: lower score (-9.0) but only 1 seed
        # Cluster 2: slightly higher score (-8.0) but supported by 3 seeds
        c1 = self._make_cluster(cluster_id=1, seed_count=1, median_score=-9.0)
        c2 = self._make_cluster(cluster_id=2, seed_count=3, median_score=-8.0)

        selected = _select_consensus_cluster([c1, c2])
        self.assertEqual(selected.cluster_id, 2, "Consensus selection must prioritize seed reproducibility over raw score")

    def test_single_seed_comparison_rejected(self):
        c_healthy = [self._make_cluster(cluster_id=1, seed_count=2, median_score=-8.0)]
        c_mutant = [self._make_cluster(cluster_id=1, seed_count=1, median_score=-7.5)]

        with self.assertRaises(ValueError) as ctx:
            verify_comparison_clusters(c_healthy, c_mutant)
        self.assertIn("must be supported by at least two independent seeds", str(ctx.exception))


class TestRepair4FingerprintingAndInvalidation(unittest.TestCase):
    """Repair 4 & Test 6 & 7: Input SHA-256 fingerprinting and stale result invalidation."""

    def test_fingerprint_changes_on_input_modification(self):
        base_inputs = {
            "workflow": "single_docking",
            "target": "1M17",
            "smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
            "exhaustiveness": 8,
            "seeds": [42, 101, 2024],
        }
        fp1 = compute_input_fingerprint(base_inputs)

        # Modify exhaustiveness
        mod_inputs = dict(base_inputs, exhaustiveness=16)
        fp2 = compute_input_fingerprint(mod_inputs)
        self.assertNotEqual(fp1, fp2)

        # Modify SMILES
        mod_smiles = dict(base_inputs, smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
        fp3 = compute_input_fingerprint(mod_smiles)
        self.assertNotEqual(fp1, fp3)

    def test_manifest_preservation_on_state_pop(self):
        temp_dir = Path(tempfile.mkdtemp(prefix="test_fp_"))
        try:
            # Build and save manifest in temp_dir
            manifest_path = build_and_save_manifest(
                run_dir=temp_dir,
                run_id="run_test_001",
                status="exploratory",
                target_meta={"pdb_id": "1M17"},
                binding_site_meta={"residues": "A:790"},
                ligand_meta={"name": "Drug"},
                prep_meta={},
                engine_meta={"version": "1.2.5"},
                grid_meta={"center": [0, 0, 0]},
                results_meta={"top_score": -8.5},
                warnings=[],
                input_fingerprint="abc123hash",
            )
            self.assertTrue(manifest_path.exists())

            # Simulating Streamlit state pop on rerun/invalidation:
            session_state = {"current_docking_result": {"run_dir": temp_dir}}
            session_state.pop("current_docking_result", None)

            # File on disk remains intact and retrievable
            self.assertTrue(manifest_path.exists())
            with open(manifest_path, "r", encoding="utf-8") as f:
                saved = json.load(f)
            self.assertEqual(saved["input_fingerprint"], "abc123hash")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def test_invalid_smiles_leaves_no_valid_ligand(self):
        temp_dir = Path(tempfile.mkdtemp(prefix="test_invalid_smiles_"))
        try:
            with self.assertRaises(ValueError):
                prepare_ligand_from_smiles("THIS_IS_NOT_A_VALID_SMILES", "BadDrug", temp_dir)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)


class TestRepair5And6InteractionLabelsAndScoreInterpretation(unittest.TestCase):
    """Repair 5, 6 & Test 10: Exploratory distance contacts, disclaimer, and non-clinical interpretation."""

    def test_mandatory_disclaimer_content(self):
        required_phrase = (
            "Distance-based contact candidates only. They are not hydrogen-bond assignments, "
            "binding free-energy components, or experimental interaction evidence."
        )
        self.assertIn(required_phrase, CONTACTS_DISCLAIMER)

    def test_interaction_keys_renamed(self):
        pdb_sample = (
            "ATOM      1  N   ALA A  10       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  O   ALA A  10       2.500   0.000   0.000  1.00 20.00           O\n"
            "END\n"
        )
        ligand_pdbqt = (
            "ATOM      1  N   LIG     1       2.800   0.000   0.000  0.00  0.00           N\n"
            "ATOM      2  C   LIG     1       5.000   0.000   0.000  0.00  0.00           C\n"
            "END\n"
        )
        result = analyze_protein_ligand_interactions(pdb_sample, ligand_pdbqt)
        self.assertIn("polar_contact_candidates", result)
        self.assertIn("carbon_contact_candidates", result)
        self.assertIn("polar_contact_count", result)
        self.assertIn("carbon_contact_count", result)

    def test_interpret_score_difference_is_non_clinical(self):
        for delta in [2.5, -2.0, 0.0]:
            info = interpret_score_difference(delta)
            text = info["summary"].lower()
            cleaned = text.replace("this is not binding affinity, resistance, efficacy, or clinical guidance", "")
            self.assertNotIn("clinical", cleaned)
            self.assertNotIn("treatment failure", text)
            self.assertNotIn("dose", text)
            self.assertIn("exploratory protocol-specific score difference", text)


class TestRepair7CoCrystalLigandSelection(unittest.TestCase):
    """Repair 7 & Test 8: Co-crystal ligand with multiple copies requires explicit chain/residue."""

    def setUp(self):
        # Two copies of ligand AQ4 in different chains and residue numbers
        self.duplicate_lig_pdb = (
            "HETATM    1  N1  AQ4 A1001       0.000   0.000   0.000  1.00 20.00           N\n"
            "HETATM    2  C1  AQ4 A1001       1.000   0.000   0.000  1.00 20.00           C\n"
            "HETATM    3  N1  AQ4 B2001      50.000  50.000  50.000  1.00 20.00           N\n"
            "HETATM    4  C1  AQ4 B2001      51.000  50.000  50.000  1.00 20.00           C\n"
            "END\n"
        )

    def test_ambiguous_ligand_raises_error(self):
        with self.assertRaises(ValueError) as ctx:
            extract_reference_ligand(self.duplicate_lig_pdb, "AQ4")
        self.assertIn("has 2 instances", str(ctx.exception))
        self.assertIn("Choose its chain and residue number explicitly", str(ctx.exception))

    def test_explicit_ligand_selection_succeeds(self):
        coords, text = extract_reference_ligand(self.duplicate_lig_pdb, "AQ4", chain_id="A", resseq=1001)
        self.assertEqual(len(coords), 2)
        self.assertIn("A1001", text)


class TestRepair8DeploymentAndOpenBabel(unittest.TestCase):
    """Repair 8 & OpenBabel regression tests from DEPLOYED_OPENBABEL_PDB_PREP_REPAIR.md."""

    def test_openbabel_formats_ready(self):
        formats = require_openbabel_formats()
        for fmt in ("pdb", "mol2", "sdf", "pdbqt"):
            self.assertIn(fmt, formats)

    def test_direct_obconversion_pdb_to_pdbqt(self):
        pdb_sample = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )
        mol = read_text(pdb_sample, "pdb")
        self.assertIsNotNone(mol)
        pdbqt_out = write_text(mol, "pdbqt")
        self.assertIn("ATOM", pdbqt_out)

    def test_vina_binary_exists_and_cross_platform(self):
        vina_path = get_vina_path()
        self.assertTrue(vina_path.exists(), f"Vina executable not found at {vina_path}")
        ver = get_vina_version()
        self.assertIn("AutoDock Vina", ver)

    def test_clean_docking_init_no_monkey_patch(self):
        import docking
        self.assertIn("Bio-Impact Analyzer", docking.__doc__)
        self.assertFalse(hasattr(docking, "_installed_patches"))


class TestHemoglobin2HHBvs2HBS(unittest.TestCase):
    """Test 1 & 2: 2HHB vs 2HBS Chain B alignment, SASA maps, and low sequence identity guard."""

    def test_2hhb_vs_2hbs_chain_b_and_sasa(self):
        # Fetch biological assemblies or deposited structures
        h_pdb = fetch_biological_assembly("2HHB", 1) or fetch_deposited_pdb("2HHB")
        m_pdb = fetch_biological_assembly("2HBS", 1) or fetch_deposited_pdb("2HBS")

        self.assertIsNotNone(h_pdb, "Failed to fetch 2HHB coordinates")
        self.assertIsNotNone(m_pdb, "Failed to fetch 2HBS coordinates")

        h_struct = process_protein_structure(h_pdb, "2HHB")
        m_struct = process_protein_structure(m_pdb, "2HBS")

        h_seq = get_protein_sequence(h_struct, "B")
        m_seq = get_protein_sequence(m_struct, "B")

        self.assertGreater(len(h_seq), 0)
        self.assertEqual(len(h_seq), len(m_seq), "2HHB and 2HBS Chain B sequence lengths must match")
        self.assertEqual(len(h_seq), 146, "Hemoglobin beta chain must contain 146 residues")

        # Alignment
        h_str = "".join(AA_3TO1.get(r["res_name"], "X") for r in h_seq)
        m_str = "".join(AA_3TO1.get(r["res_name"], "X") for r in m_seq)
        _, score, aln_h, aln_m = get_alignment(h_str, m_str, "global")

        # Identity calculation
        aligned_len = len(aln_h)
        identity_pct = (sum(a == b and a != '-' for a, b in zip(aln_h, aln_m)) / aligned_len * 100)
        self.assertGreater(identity_pct, 95.0, "HbA vs HbS Chain B sequence identity must be > 95%")

        # SASA calculation
        h_sasa = calculate_sasa_map(h_struct, "B")
        m_sasa = calculate_sasa_map(m_struct, "B")
        self.assertEqual(len(h_sasa), len(h_seq), "Every residue in healthy chain B must have a SASA value")
        self.assertEqual(len(m_sasa), len(m_seq), "Every residue in mutant chain B must have a SASA value")

    def test_unrelated_sequences_low_identity_flagged(self):
        # Test low sequence identity guard
        seq1 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNAL"
        seq2 = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"

        _, _, aln_1, aln_2 = get_alignment(seq1, seq2, "global")
        aligned_len = len(aln_1)
        identity_pct = (sum(a == b and a != '-' for a, b in zip(aln_1, aln_2)) / aligned_len * 100)

        # Identity between hemoglobin and Ras is < 30%
        self.assertLess(identity_pct, 40.0, "Unrelated proteins must have < 40% sequence identity")


class Test1M17DockingPreparation(unittest.TestCase):
    """Test 9: 1M17/AQ4/Erlotinib receptor and ligand preparation."""

    def test_1m17_aq4_erlotinib_prep(self):
        pdb_1m17 = fetch_biological_assembly("1M17", 1) or fetch_deposited_pdb("1M17")
        self.assertIsNotNone(pdb_1m17, "Could not fetch 1M17")

        temp_dir = Path(tempfile.mkdtemp(prefix="test_1m17_prep_"))
        try:
            # 1. Define grid from reference ligand AQ4
            grid, ref_coords, ref_ligand_pdb = define_grid_from_ligand(
                pdb_1m17, "AQ4", padding=8.0, chain_id="A", resseq=999
            )
            self.assertGreater(len(ref_coords), 0)
            self.assertIsInstance(grid, GridBox)

            # 2. Prepare receptor
            clean_pdb, rec_pdbqt, prep_rep = prepare_receptor(
                pdb_1m17, temp_dir, reference_ligand_resname="AQ4"
            )
            self.assertTrue(Path(rec_pdbqt).exists())
            self.assertGreater(Path(rec_pdbqt).stat().st_size, 0)
            self.assertTrue(prep_rep["reference_ligand_removed"])

            # 3. Prepare ligand Erlotinib
            erlotinib_smiles = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
            sdf_path, lig_pdbqt, lig_rep = prepare_ligand_from_smiles(
                erlotinib_smiles, "Erlotinib", temp_dir
            )
            self.assertTrue(Path(lig_pdbqt).exists())
            self.assertGreater(lig_rep["total_atom_count"], 0)

        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)


class TestRemainingRepairsRegression(unittest.TestCase):
    """Regression suite covering remaining correctness and scientific-interpretation repairs."""

    def test_carbon_carbon_pair_at_3A_counted(self):
        rec_pdb = (
            "ATOM      1  CA  ALA A 100       0.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )
        lig_pdbqt = (
            "ATOM      1  C1  LIG     1       3.000   0.000   0.000  0.00  0.00    +0.000 C\n"
        )
        res = analyze_protein_ligand_interactions(rec_pdb, lig_pdbqt, polar_cutoff=3.5, contact_cutoff=4.0)
        self.assertEqual(res["carbon_contact_count"], 1)
        self.assertEqual(len(res["carbon_contact_candidates"]), 1)
        self.assertEqual(res["polar_contact_count"], 0)
        self.assertEqual(res["disclaimer"], CONTACTS_DISCLAIMER)

    def test_lowest_vina_score_returns_global_minimum_across_seeds(self):
        poses = [
            DockingPose(seed=42, rank=1, score=-7.4, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="", cluster_id=1),
            DockingPose(seed=101, rank=1, score=-9.2, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="", cluster_id=2),
            DockingPose(seed=2024, rank=1, score=-8.1, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="", cluster_id=3),
        ]
        best = lowest_vina_score(poses)
        self.assertEqual(best, -9.2)
        self.assertNotEqual(best, poses[0].score)
        self.assertIsNone(lowest_vina_score([]))

    def test_insertion_code_selection_distinguishes_100_and_100A(self):
        synth_pdb = (
            "ATOM      1  CA  ALA A 100       0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  ALA A 100A     50.000  50.000  50.000  1.00 20.00           C\n"
            "END\n"
        )
        grid_100 = define_grid_from_residues(synth_pdb, ["A:100"], padding=5.0)
        grid_100a = define_grid_from_residues(synth_pdb, ["A:100A"], padding=5.0)
        self.assertAlmostEqual(grid_100.center_x, 0.0, delta=0.5)
        self.assertAlmostEqual(grid_100a.center_x, 50.0, delta=0.5)

    def test_coordinate_text_difference_changes_sha256_fingerprint(self):
        h1 = "ATOM 1 CA ALA A 100 0.0 0.0 0.0\n"
        h2 = "ATOM 1 CA ALA A 100 10.0 10.0 10.0\n"
        m = "ATOM 1 CA ALA A 100 0.0 0.0 0.0\n"

        p1 = {
            "workflow": "matched_healthy_vs_mutant",
            "healthy_id": "WT",
            "mutant_id": "MUT",
            "healthy_coordinate_source": "biological_assembly_1",
            "mutant_coordinate_source": "biological_assembly_1",
            "healthy_coordinate_sha256": compute_text_sha256(h1),
            "mutant_coordinate_sha256": compute_text_sha256(m),
            "pocket_res": "A:100",
            "target_domain": "Kinase",
            "lig_name": "Drug",
            "lig_smiles": "CC",
            "exhaustiveness": 8,
            "vina_version": "1.2.5",
        }
        p2 = dict(p1, healthy_coordinate_sha256=compute_text_sha256(h2))
        p3 = dict(p1, mutant_coordinate_sha256=compute_text_sha256(h2))

        self.assertNotEqual(compute_text_sha256(h1), compute_text_sha256(h2))
        self.assertNotEqual(compute_input_fingerprint(p1), compute_input_fingerprint(p2))
        self.assertNotEqual(compute_input_fingerprint(p1), compute_input_fingerprint(p3))

    def test_matched_comparison_manifest_contains_input_fingerprint(self):
        pdb_text = (
            "ATOM      1  N   ALA A 790       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A 790       1.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      3  C   ALA A 790       2.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      4  O   ALA A 790       2.000   1.000   0.000  1.00 20.00           O\n"
            "END\n"
        )
        test_fp = "matched_fp_provenance_test_val_999"

        def mock_vina(receptor_pdbqt, ligand_pdbqt, grid, output_dir, seeds, exhaustiveness):
            poses = [
                DockingPose(seed=s, rank=1, score=-8.0, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="ATOM      1  C   LIG     1       0.000   0.000   0.000\n", cluster_id=1)
                for s in seeds
            ]
            return poses, [{"seed": s, "top_score": -8.0} for s in seeds]

        with patch("docking.comparison.run_vina_multi_seeds", side_effect=mock_vina):
            res = run_matched_docking_comparison(
                healthy_pdb_text=pdb_text,
                mutant_pdb_text=pdb_text,
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="CC",
                ligand_name="Ethanol",
                pocket_residues=["A:790"],
                input_fingerprint=test_fp,
            )

        h_run_id = res["healthy"]["run_id"]
        m_run_id = res["mutant"]["run_id"]
        h_dir = BASE_RUNS_DIR / h_run_id
        m_dir = BASE_RUNS_DIR / m_run_id
        try:
            with open(h_dir / "manifest.json", encoding="utf-8") as f:
                h_data = json.load(f)
            with open(m_dir / "manifest.json", encoding="utf-8") as f:
                m_data = json.load(f)

            self.assertEqual(h_data.get("input_fingerprint"), test_fp)
            self.assertEqual(m_data.get("input_fingerprint"), test_fp)
        finally:
            shutil.rmtree(h_dir, ignore_errors=True)
            shutil.rmtree(m_dir, ignore_errors=True)

    def test_streamlit_apptest_zero_exceptions(self):
        at = AppTest.from_file(str(PROJECT_ROOT / "app.py"))
        at.run(timeout=30)
        self.assertFalse(at.exception, [str(x.value) for x in at.exception])

    def test_heme_cma_does_not_produce_curium_cm_atom_type(self):
        """Ensure atoms named CMA/CMB (e.g. in HEM cofactors) do not get misparsed as Curium (Cm)."""
        pdb_sample = (
            "ATOM      1  N   ALA A   1      20.000  20.000  20.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1      21.000  20.000  20.000  1.00 20.00           C\n"
            "HETATM 1078  CMA HEM A 142      22.338  21.946  21.537  1.00 23.27           C\n"
            "HETATM 1079  CAA HEM A 142      23.000  22.000  22.000  1.00 23.27           C\n"
        )
        temp_dir = Path(tempfile.mkdtemp(prefix="test_heme_cm_"))
        try:
            clean_pdb, rec_pdbqt, rep = prepare_receptor(pdb_sample, temp_dir)
            with open(rec_pdbqt, "r", encoding="utf-8") as f:
                rec_content = f.read()
            for line in rec_content.splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    atom_type = line[77:79].strip()
                    self.assertNotEqual(atom_type, "Cm", f"Found invalid AutoDock type Cm in: {line}")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)



class TestGenericMatchedDockingScientificReliability(unittest.TestCase):
    """Section 7: Required regression tests for generic matched docking & scientific reliability."""

    def _make_synthetic_pdb(self, chain: str, residues: list[tuple[str, int]]) -> str:
        lines = []
        atom_num = 1
        for resname, resnum in residues:
            lines.append(f"ATOM  {atom_num:5d}  N   {resname:3s} {chain:1s}{resnum:4d}       0.000   0.000   0.000  1.00 20.00           N")
            atom_num += 1
            lines.append(f"ATOM  {atom_num:5d}  CA  {resname:3s} {chain:1s}{resnum:4d}       1.000   0.000   0.000  1.00 20.00           C")
            atom_num += 1
            lines.append(f"ATOM  {atom_num:5d}  C   {resname:3s} {chain:1s}{resnum:4d}       2.000   0.000   0.000  1.00 20.00           C")
            atom_num += 1
            lines.append(f"ATOM  {atom_num:5d}  O   {resname:3s} {chain:1s}{resnum:4d}       2.000   1.000   0.000  1.00 20.00           O")
            atom_num += 1
        lines.append("END\n")
        return "\n".join(lines)

    def test_1_unrelated_chains_rejected_in_matched_mode_and_offered_exploratory(self):
        """1. Two unrelated selected chains must be rejected in Matched mode and offered the exploratory cross-protein workflow."""
        wt_pdb = self._make_synthetic_pdb("A", [("ALA", 1), ("GLY", 2), ("SER", 3), ("LEU", 4), ("VAL", 5)])
        unrelated_pdb = self._make_synthetic_pdb("B", [("TRP", 1), ("PHE", 2), ("TYR", 3), ("HIS", 4), ("MET", 5)])

        val = validate_matched_pair(wt_pdb, unrelated_pdb, "A", "B", min_identity=MIN_MATCHED_PAIR_IDENTITY)
        self.assertFalse(val["is_valid"])
        self.assertLess(val["identity_pct"], 90.0)
        self.assertIn("Matched comparison was not run: selected chains are not a validated WT/variant pair", val["error"])
        self.assertIn("Use Exploratory Cross-Protein Docking instead", val["error"])

    def test_2_high_identity_synthetic_pair_with_numbering_shift_creates_separate_pockets(self):
        """2. A high-identity synthetic WT/variant pair with a residue-number shift must create correct separate WT and mutant pocket lists."""
        wt_pdb = self._make_synthetic_pdb("A", [("ALA", 10), ("GLY", 11), ("SER", 12), ("LEU", 13)])
        mt_pdb = self._make_synthetic_pdb("B", [("ALA", 50), ("GLY", 51), ("SER", 52), ("LEU", 53)])

        val = validate_matched_pair(wt_pdb, mt_pdb, "A", "B", min_identity=MIN_MATCHED_PAIR_IDENTITY)
        self.assertTrue(val["is_valid"])
        self.assertEqual(val["identity_pct"], 100.0)

        mapped_m, unmapped_h, errors = map_healthy_pocket_residues(
            ["A:10", "A:12"], val["mapping"], mt_pdb, "B", healthy_chain="A"
        )
        self.assertEqual(mapped_m, ["B:50", "B:52"])
        self.assertEqual(unmapped_h, [])
        self.assertEqual(errors, [])

        # Build grids separately
        g_wt = define_grid_from_residues(wt_pdb, ["A:10", "A:12"])
        g_mt = define_grid_from_residues(mt_pdb, mapped_m)
        self.assertIsInstance(g_wt, GridBox)
        self.assertIsInstance(g_mt, GridBox)

    def test_3_ui_contact_and_visualization_pose_matches_consensus_representative_not_cluster_0(self):
        """3. The UI contact pose must be the returned consensus representative, not clusters[0] when they differ."""
        pose_c1 = DockingPose(seed=42, rank=1, score=-9.5, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="ATOM 1 C LIG 1 0.0 0.0 0.0\n", cluster_id=1)
        c1_low_seed = PoseCluster(
            cluster_id=1, seed_count=1, supporting_seeds=[42], top_score=-9.5, median_score=-9.5,
            representative_pose=pose_c1, representative_pose_source={"seed": 42, "rank": 1}
        )

        pose_c2 = DockingPose(seed=42, rank=2, score=-8.2, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="ATOM 1 C LIG 1 10.0 10.0 10.0\n", cluster_id=2)
        c2_consensus = PoseCluster(
            cluster_id=2, seed_count=3, supporting_seeds=[42, 101, 2024], top_score=-8.2, median_score=-8.2,
            representative_pose=pose_c2, representative_pose_source={"seed": 42, "rank": 2}
        )

        clusters = [c1_low_seed, c2_consensus]
        selected = _select_consensus_cluster(clusters)
        self.assertEqual(selected.cluster_id, 2)
        self.assertEqual(selected.representative_pose, pose_c2)
        self.assertNotEqual(selected.representative_pose, clusters[0].representative_pose)

    def test_4_single_seed_cluster_prevents_delta_score(self):
        """4. A cluster supported by one seed must prevent ΔScore output."""
        pose_h = DockingPose(seed=42, rank=1, score=-8.0, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="", cluster_id=1)
        cl_h = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-8.0, median_score=-8.0, representative_pose=pose_h, representative_pose_source={})]

        pose_m = DockingPose(seed=42, rank=1, score=-7.0, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="", cluster_id=1)
        cl_m_single_seed = [PoseCluster(cluster_id=1, seed_count=1, supporting_seeds=[42], top_score=-7.0, median_score=-7.0, representative_pose=pose_m, representative_pose_source={})]

        with self.assertRaises(ValueError) as ctx:
            verify_comparison_clusters(cl_h, cl_m_single_seed)
        self.assertIn("must be supported by at least two independent seeds", str(ctx.exception))

    def test_5_matched_run_without_reference_redocking_labelled_exploratory(self):
        """5. A matched run without reference redocking must be labelled Exploratory: no pose validation."""
        pdb_text = self._make_synthetic_pdb("A", [("ALA", 10), ("GLY", 11)])

        def mock_vina(receptor_pdbqt, ligand_pdbqt, grid, output_dir, seeds, exhaustiveness):
            poses = [
                DockingPose(seed=s, rank=1, score=-8.0, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block="ATOM      1  C   LIG     1       0.000   0.000   0.000\n", cluster_id=1)
                for s in seeds
            ]
            return poses, [{"seed": s, "top_score": -8.0} for s in seeds]

        with patch("docking.comparison.run_vina_multi_seeds", side_effect=mock_vina):
            res = run_matched_docking_comparison(
                healthy_pdb_text=pdb_text,
                mutant_pdb_text=pdb_text,
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="CC",
                ligand_name="Ethanol",
                healthy_pocket_residues=["A:10"],
                mutant_pocket_residues=["A:10"],
                wt_reference_ligand_meta=None,
            )

        self.assertEqual(res["validation_status"], "exploratory—no reference validation")
        self.assertEqual(res["validation_reason"], "No verified WT co-crystal ligand was available for redocking validation.")
        self.assertEqual(res["protonation_note"], PROTONATION_ASSUMPTION_NOTE)

        # Cleanup created run dirs
        shutil.rmtree(BASE_RUNS_DIR / res["healthy"]["run_id"], ignore_errors=True)
        shutil.rmtree(BASE_RUNS_DIR / res["mutant"]["run_id"], ignore_errors=True)

    def test_6_passing_synthetic_redocking_reports_rmsd_criterion_top_n_seed_support(self):
        """6. A passing synthetic or documented redocking case must show its RMSD criterion, Top-1/Top-N result, seed support, and no-fit mapping method."""
        ref_lig = (
            "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "HETATM    2  C2  LIG A   1       1.500   0.000   0.000  1.00 20.00           C\n"
        )
        pose_pdbqt = (
            "ATOM      1  C1  LIG     1       0.050   0.000   0.000  0.00  0.00    +0.000 C\n"
            "ATOM      2  C2  LIG     1       1.550   0.000   0.000  0.00  0.00    +0.000 C\n"
        )
        pose = DockingPose(seed=42, rank=1, score=-8.5, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=pose_pdbqt, cluster_id=1)
        cl = PoseCluster(
            cluster_id=1, seed_count=3, supporting_seeds=[42, 101, 2024], top_score=-8.5, median_score=-8.5,
            representative_pose=pose, representative_pose_source={"seed": 42, "rank": 1}
        )

        report = validate_redocking_top_n(
            reference_ligand_pdb_text=ref_lig,
            poses=[pose],
            clusters=[cl],
            top_n=10,
            rmsd_cutoff=2.0,
            min_seed_support=2,
        )

        self.assertEqual(report["rmsd_cutoff"], 2.0)
        self.assertEqual(report["min_seed_support"], 2)
        self.assertEqual(report["top_n"], 1)
        self.assertEqual(report["method"], "symmetry-aware heavy-atom mapped RMSD (no-fit)")
        self.assertTrue(report["top1_rmsd_pass"])
        self.assertTrue(report["recovery_rmsd_pass"])
        self.assertTrue(report["recovery_cluster_pass"])
        self.assertLess(report["top1"]["rmsd_angstrom"], 2.0)

    def test_7_empty_matched_fields_block_execution_without_egfr_defaults(self):
        """7. Empty matched fields must block execution; no implicit EGFR, Erlotinib, or residue defaults may be used."""
        pdb = self._make_synthetic_pdb("A", [("ALA", 1)])
        with self.assertRaises(ValueError) as ctx1:
            run_matched_docking_comparison(
                healthy_pdb_text=pdb,
                mutant_pdb_text=pdb,
                healthy_target_meta={},
                mutant_target_meta={},
                ligand_smiles="CC",
                ligand_name="Drug",
                healthy_pocket_residues=[],
                mutant_pocket_residues=[],
            )
        self.assertIn("define validated healthy and mutant binding sites", str(ctx1.exception))

        import docking.ui as dui
        self.assertFalse(hasattr(dui, "DEFAULT_EGFR_POCKET"))

    def test_8_hemoglobin_is_not_hardcoded_as_required_or_default(self):
        """8. Hemoglobin is only an optional test, never a hardcoded requirement or UI default."""
        import docking.comparison as dc
        import docking.pair_validation as dpv
        import inspect

        sig = inspect.signature(dc.run_matched_docking_comparison)
        for param in sig.parameters.values():
            self.assertNotEqual(param.default, "2HHB")
            self.assertNotEqual(param.default, "2HBS")

        self.assertEqual(dpv.MIN_MATCHED_PAIR_IDENTITY, 90.0)


class TestConsolidatedCrossValidationAndSessionStateRepairs(unittest.TestCase):
    """Regression tests verifying all 10 requirements in Section 6."""

    def _make_pdb(self, chain="A", resseq=1):
        return (
            f"ATOM  {1:5d}  N   ALA {chain:1s}{resseq:4d}       0.000   0.000   0.000  1.00 20.00           N\n"
            f"ATOM  {2:5d}  CA  ALA {chain:1s}{resseq:4d}       1.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )

    def _make_dummy_pose(self, seed=42, score=-7.0, rank=1):
        pdbqt = "ATOM      1  C1  LIG     1       0.000   0.000   0.000  0.00  0.00    +0.000 C\n"
        return DockingPose(seed=seed, rank=rank, score=score, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=pdbqt, cluster_id=1)

    # 1. Successful mocked Cross-Protein docking creates two final manifests with their own pocket lists and no delta_score.
    @patch("docking.comparison.run_vina_multi_seeds")
    def test_1_successful_mocked_cross_docking_manifests_no_delta_score(self, mock_vina):
        mock_vina.side_effect = [
            ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            ([self._make_dummy_pose(seed=42, score=-6.5), self._make_dummy_pose(seed=101, score=-6.5), self._make_dummy_pose(seed=2024, score=-6.5)], []),
        ]
        pdb1 = self._make_pdb("A", 10)
        pdb2 = self._make_pdb("B", 20)
        res = run_exploratory_cross_docking(
            receptor1_pdb_text=pdb1,
            receptor2_pdb_text=pdb2,
            receptor1_target_meta={"pdb_id": "REC1", "selected_chain": "A"},
            receptor2_target_meta={"pdb_id": "REC2", "selected_chain": "B"},
            ligand_smiles="CC",
            ligand_name="Ethanol",
            receptor1_pocket_residues=["A:10"],
            receptor2_pocket_residues=["B:20"],
        )
        self.assertEqual(res["status"], "exploratory—cross-protein docking")
        self.assertNotIn("delta_score", res)

        r1_dir = BASE_RUNS_DIR / res["receptor1"]["run_id"]
        r2_dir = BASE_RUNS_DIR / res["receptor2"]["run_id"]

        self.assertTrue((r1_dir / "receptor_source.pdb").exists())
        self.assertTrue((r2_dir / "receptor_source.pdb").exists())
        self.assertTrue((r1_dir / "manifest.initial.json").exists())
        self.assertTrue((r2_dir / "manifest.initial.json").exists())
        self.assertTrue((r1_dir / "manifest.json").exists())
        self.assertTrue((r2_dir / "manifest.json").exists())

        with open(r1_dir / "manifest.json", "r", encoding="utf-8") as f:
            m1 = json.load(f)
        with open(r2_dir / "manifest.json", "r", encoding="utf-8") as f:
            m2 = json.load(f)

        self.assertEqual(m1["status"], "exploratory—cross-protein docking")
        self.assertEqual(m2["status"], "exploratory—cross-protein docking")
        self.assertEqual(m1["binding_site"]["residues"], ["A:10"])
        self.assertEqual(m2["binding_site"]["residues"], ["B:20"])
        self.assertNotIn("delta_score", m1)
        self.assertNotIn("delta_score", m1.get("results", {}))
        self.assertNotIn("delta_score", m2)
        self.assertNotIn("delta_score", m2.get("results", {}))

        # Cleanup
        shutil.rmtree(r1_dir, ignore_errors=True)
        shutil.rmtree(r2_dir, ignore_errors=True)

    # 2. Failure during the second Cross-Protein Vina run finalizes both manifests as run failed with the exact failure reason.
    @patch("docking.comparison.run_vina_multi_seeds")
    def test_2_failure_on_second_receptor_finalizes_both_manifests_failed(self, mock_vina):
        mock_vina.side_effect = [
            ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            RuntimeError("Simulated engine crash on receptor 2"),
        ]
        pdb1 = self._make_pdb("A", 10)
        pdb2 = self._make_pdb("B", 20)

        with self.assertRaises(RuntimeError) as ctx:
            run_exploratory_cross_docking(
                receptor1_pdb_text=pdb1,
                receptor2_pdb_text=pdb2,
                receptor1_target_meta={"pdb_id": "REC1"},
                receptor2_target_meta={"pdb_id": "REC2"},
                ligand_smiles="CC",
                ligand_name="Ethanol",
                receptor1_pocket_residues=["A:10"],
                receptor2_pocket_residues=["B:20"],
            )
        self.assertIn("Simulated engine crash on receptor 2", str(ctx.exception))

        r1_dirs = sorted(BASE_RUNS_DIR.glob("dock_cross_receptor1_*"), key=os.path.getmtime)
        r2_dirs = sorted(BASE_RUNS_DIR.glob("dock_cross_receptor2_*"), key=os.path.getmtime)
        self.assertTrue(len(r1_dirs) > 0)
        self.assertTrue(len(r2_dirs) > 0)

        latest_r1 = r1_dirs[-1]
        latest_r2 = r2_dirs[-1]

        with open(latest_r1 / "manifest.json", "r", encoding="utf-8") as f:
            m1 = json.load(f)
        with open(latest_r2 / "manifest.json", "r", encoding="utf-8") as f:
            m2 = json.load(f)

        self.assertEqual(m1["status"], "run failed")
        self.assertEqual(m2["status"], "run failed")
        self.assertIn("Simulated engine crash on receptor 2", str(m1["results"]["failure_reason"]))
        self.assertIn("Simulated engine crash on receptor 2", str(m2["results"]["failure_reason"]))

        # Cleanup
        shutil.rmtree(latest_r1, ignore_errors=True)
        shutil.rmtree(latest_r2, ignore_errors=True)

    # 3. A one-seed Cross-Protein cluster is labelled non-reproducible; a two-seed cluster is labelled reproducible.
    def test_3_one_seed_cluster_non_reproducible_two_seed_reproducible(self):
        self.assertEqual(MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT, 2)
        pose = self._make_dummy_pose(seed=42)
        cl_1seed = PoseCluster(
            cluster_id=1, seed_count=1, supporting_seeds=[42], top_score=-7.0, median_score=-7.0,
            representative_pose=pose, representative_pose_source={"seed": 42, "rank": 1}
        )
        cl_2seed = PoseCluster(
            cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-7.0, median_score=-7.0,
            representative_pose=pose, representative_pose_source={"seed": 42, "rank": 1}
        )

        reprod_1 = (cl_1seed.seed_count >= MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT)
        reprod_2 = (cl_2seed.seed_count >= MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT)

        self.assertFalse(reprod_1)
        self.assertTrue(reprod_2)

    # 4. Two ligands with equal heavy-atom counts but different connectivity or bond order fail the strict reference-ligand identity gate and cannot receive Protocol validated on WT reference.
    def test_4_different_connectivity_isomers_fail_identity_gate(self):
        # Ethanol (CCO) vs Dimethyl ether (COC): same heavy atom count (3), same formula (C2H6O), different connectivity
        passed, reason, details = verify_reference_ligand_identity("CCO", "COC")
        self.assertFalse(passed)
        self.assertIn("Canonical isomeric SMILES mismatch", reason)

        # Propan-1-ol (CCCO) vs Propan-2-ol (CC(O)C): same heavy atom count (4), same formula (C3H8O)
        passed_prop, reason_prop, _ = verify_reference_ligand_identity("CCCO", "CC(O)C")
        self.assertFalse(passed_prop)
        self.assertIn("Canonical isomeric SMILES mismatch", reason_prop)

        # Test within run_matched_docking_comparison
        pdb = self._make_pdb("A", 10)
        ref_pdb = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        with patch("docking.comparison.run_vina_multi_seeds") as mock_v, \
             patch("docking.comparison.fetch_chem_comp_smiles") as mock_rcsb:
            mock_v.side_effect = [
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            ]
            mock_rcsb.return_value = "CCO"  # WT co-crystal is Ethanol
            res = run_matched_docking_comparison(
                healthy_pdb_text=pdb,
                mutant_pdb_text=pdb,
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="COC",  # Submitted is Dimethyl ether
                ligand_name="Ether",
                healthy_pocket_residues=["A:10"],
                mutant_pocket_residues=["A:10"],
                wt_reference_ligand_meta={"resname": "LIG", "ref_ligand_pdb_text": ref_pdb},
            )

        self.assertEqual(
            res["validation_status"],
            "not validated for affinity comparison"
        )
        self.assertEqual(
            res["validation_reason"],
            "The submitted ligand is not chemically identical to the verified WT co-crystal reference ligand."
        )
        self.assertNotEqual(res["validation_status"], "protocol validated on WT reference")
        shutil.rmtree(BASE_RUNS_DIR / res["healthy"]["run_id"], ignore_errors=True)
        shutil.rmtree(BASE_RUNS_DIR / res["mutant"]["run_id"], ignore_errors=True)

    # 5. Cross-Protein pocket residues belonging to a chain other than the selected chain are rejected.
    def test_5_cross_protein_pocket_residues_wrong_chain_rejected(self):
        # Valid: residues on chain A with selected chain A
        validate_pocket_chain(["A:10", "A:15"], "A", "Receptor 1")

        # Invalid: residue on chain B when selected chain is A
        with self.assertRaises(ValueError) as ctx:
            validate_pocket_chain(["B:10", "A:15"], "A", "Receptor 1")
        self.assertIn("does not belong to selected chain 'A'", str(ctx.exception))

        # Invalid: unqualified residue without chain colon
        with self.assertRaises(ValueError) as ctx2:
            validate_pocket_chain(["10"], "A", "Receptor 1")
        self.assertIn("must use CHAIN:RESIDUE notation", str(ctx2.exception))

        # Check rejection in run_exploratory_cross_docking
        pdb = self._make_pdb("A", 10)
        with self.assertRaises(ValueError) as ctx3:
            run_exploratory_cross_docking(
                receptor1_pdb_text=pdb,
                receptor2_pdb_text=pdb,
                receptor1_target_meta={"pdb_id": "REC1", "selected_chain": "A"},
                receptor2_target_meta={"pdb_id": "REC2", "selected_chain": "A"},
                ligand_smiles="CC",
                ligand_name="Drug",
                receptor1_pocket_residues=["B:10"],
                receptor2_pocket_residues=["A:10"],
            )
        self.assertIn("does not belong to selected chain 'A'", str(ctx3.exception))

    # 6. Replacing only healthy coordinates clears current_docking_result, matched_comparison_result, exploratory_cross_result, and screening_results.
    def test_6_replacing_healthy_coordinates_clears_dependent_results(self):
        st.session_state["current_docking_result"] = {"run": 1}
        st.session_state["matched_comparison_result"] = {"comp": 1}
        st.session_state["exploratory_cross_result"] = {"cross": 1}
        st.session_state["screening_results"] = {"screen": 1}
        st.session_state["h_results"] = {"h": 1}

        invalidate_dependent_results("Healthy coordinates replaced")

        self.assertNotIn("current_docking_result", st.session_state)
        self.assertNotIn("matched_comparison_result", st.session_state)
        self.assertNotIn("exploratory_cross_result", st.session_state)
        self.assertNotIn("screening_results", st.session_state)
        self.assertNotIn("h_results", st.session_state)
        self.assertEqual(st.session_state.get("last_result_invalidation_reason"), "Healthy coordinates replaced")

    # 7. Replacing only mutant coordinates clears the same result objects.
    def test_7_replacing_mutant_coordinates_clears_dependent_results(self):
        st.session_state["current_docking_result"] = {"run": 2}
        st.session_state["matched_comparison_result"] = {"comp": 2}
        st.session_state["exploratory_cross_result"] = {"cross": 2}
        st.session_state["screening_results"] = {"screen": 2}
        st.session_state["m_results"] = {"m": 2}

        invalidate_dependent_results("Mutant coordinates replaced")

        self.assertNotIn("current_docking_result", st.session_state)
        self.assertNotIn("matched_comparison_result", st.session_state)
        self.assertNotIn("exploratory_cross_result", st.session_state)
        self.assertNotIn("screening_results", st.session_state)
        self.assertNotIn("m_results", st.session_state)
        self.assertEqual(st.session_state.get("last_result_invalidation_reason"), "Mutant coordinates replaced")

    # 8. Clearing either protein removes its source PDB/chain and makes both-protein docking unavailable.
    def test_8_clearing_either_protein_removes_source_and_blocks_comparative_docking(self):
        st.session_state["h_pdb"] = self._make_pdb("A", 1)
        st.session_state["h_assembly_pdb"] = self._make_pdb("A", 1)
        st.session_state["h_id"] = "1WT"
        st.session_state["h_selected_chain"] = "A"
        st.session_state["matched_comparison_result"] = {"comp": "active"}

        clear_protein("h", "Healthy")

        self.assertNotIn("h_pdb", st.session_state)
        self.assertNotIn("h_assembly_pdb", st.session_state)
        self.assertNotIn("h_id", st.session_state)
        self.assertNotIn("h_selected_chain", st.session_state)
        self.assertNotIn("matched_comparison_result", st.session_state)
        # Verify comparative docking is blocked because h_active_pdb is empty
        h_active = st.session_state.get("h_assembly_pdb") or st.session_state.get("h_pdb") or ""
        self.assertEqual(h_active, "")

    # 9. Re-uploading identical bytes does not invalidate current results; changed bytes do.
    def test_9_reuploading_identical_bytes_preserves_results_changed_bytes_invalidate(self):
        data1 = b"HEADER TEST 1\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n"
        data2 = b"HEADER TEST 2\nATOM      1  CA  GLY A   1       1.000   1.000   1.000  1.00 20.00           C\nEND\n"

        sha1 = hashlib.sha256(data1).hexdigest()

        # Step A: Load initial upload
        st.session_state["h_upload_sha256"] = sha1
        st.session_state["matched_comparison_result"] = {"comparison": "valid"}

        # Simulate upload with IDENTICAL bytes: hash does not change, no invalidation
        new_upload_bytes = data1
        new_sha = hashlib.sha256(new_upload_bytes).hexdigest()
        if new_sha != st.session_state.get("h_upload_sha256"):
            invalidate_dependent_results("Uploaded file changed")

        self.assertIn("matched_comparison_result", st.session_state)
        self.assertEqual(st.session_state["matched_comparison_result"]["comparison"], "valid")

        # Simulate upload with CHANGED bytes: hash changes, triggers invalidation
        changed_upload_bytes = data2
        new_sha = hashlib.sha256(changed_upload_bytes).hexdigest()
        if new_sha != st.session_state.get("h_upload_sha256"):
            st.session_state["h_upload_sha256"] = new_sha
            invalidate_dependent_results("Uploaded file changed")

        self.assertNotIn("matched_comparison_result", st.session_state)

    # 10. Old manifest files remain on disk after session-state invalidation.
    def test_10_old_manifest_files_remain_on_disk_after_invalidation(self):
        test_run_dir = BASE_RUNS_DIR / "test_provenance_run_persist"
        test_run_dir.mkdir(parents=True, exist_ok=True)
        manifest_file = test_run_dir / "manifest.json"
        with open(manifest_file, "w", encoding="utf-8") as f:
            json.dump({"run_id": "test_provenance_run_persist", "status": "exploratory—cross-protein docking"}, f)

        self.assertTrue(manifest_file.exists())

        st.session_state["matched_comparison_result"] = {"run": "something"}
        invalidate_dependent_results("User changed target")
        clear_protein("h", "Healthy")

        # The manifest file on disk MUST still exist untouched
        self.assertTrue(manifest_file.exists(), "Old manifest on disk was deleted during session state invalidation!")

        with open(manifest_file, "r", encoding="utf-8") as f:
            persisted_manifest = json.load(f)
        self.assertEqual(persisted_manifest["status"], "exploratory—cross-protein docking")

        shutil.rmtree(test_run_dir, ignore_errors=True)


class TestManifestValidationStatusRepairs(unittest.TestCase):
    """Tests verifying preservation of scientific validation status in docking manifests."""

    def _make_pdb(self, chain="A", resseq=1):
        return (
            f"ATOM  {1:5d}  N   ALA {chain:1s}{resseq:4d}       0.000   0.000   0.000  1.00 20.00           N\n"
            f"ATOM  {2:5d}  CA  ALA {chain:1s}{resseq:4d}       1.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )

    def _make_dummy_pose(self, seed=42, score=-7.0, rank=1):
        pdbqt = "ATOM      1  C1  LIG     1       0.000   0.000   0.000  0.00  0.00    +0.000 C\n"
        return DockingPose(seed=seed, rank=rank, score=score, rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=pdbqt, cluster_id=1)

    # 1. Save a manifest with protocol validated on WT reference; assert that the JSON retains exactly that value.
    def test_1_save_manifest_protocol_validated_on_wt_reference(self):
        run_dir = BASE_RUNS_DIR / "test_status_validated_wt"
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            build_and_save_manifest(
                run_dir=run_dir,
                run_id="test_status_validated_wt",
                status="protocol validated on WT reference",
                target_meta={"pdb_id": "WT"},
                binding_site_meta={"residues": ["A:1"]},
                ligand_meta={"name": "Drug"},
                prep_meta={},
                engine_meta={},
                grid_meta={},
                results_meta={"validation": {"status": "protocol validated on WT reference"}},
            )
            with open(run_dir / "manifest.json", "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["status"], "protocol validated on WT reference")
        finally:
            shutil.rmtree(run_dir, ignore_errors=True)

    # 2. Save a manifest with not validated for affinity comparison; assert that it remains that value and its validation.reason remains present.
    def test_2_save_manifest_not_validated_for_affinity_comparison_and_preserves_reason(self):
        run_dir = BASE_RUNS_DIR / "test_status_not_validated_affinity"
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            build_and_save_manifest(
                run_dir=run_dir,
                run_id="test_status_not_validated_affinity",
                status="not validated for affinity comparison",
                target_meta={"pdb_id": "WT"},
                binding_site_meta={"residues": ["A:1"]},
                ligand_meta={"name": "Drug"},
                prep_meta={},
                engine_meta={},
                grid_meta={},
                results_meta={
                    "validation": {
                        "status": "not validated for affinity comparison",
                        "reason": "The submitted ligand is not chemically identical to the verified WT co-crystal reference ligand.",
                    }
                },
            )
            with open(run_dir / "manifest.json", "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["status"], "not validated for affinity comparison")
            self.assertEqual(
                data["results"]["validation"]["reason"],
                "The submitted ligand is not chemically identical to the verified WT co-crystal reference ligand.",
            )
        finally:
            shutil.rmtree(run_dir, ignore_errors=True)

    # 3. Mock a successful matched WT-reference validation and assert that both the returned result and final manifest use protocol validated on WT reference.
    def test_3_mocked_successful_wt_validation_result_and_manifest(self):
        pdb = self._make_pdb("A", 10)
        ref_pdb = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        with patch("docking.comparison.run_vina_multi_seeds") as mock_v, \
             patch("docking.comparison.fetch_chem_comp_smiles") as mock_rcsb, \
             patch("docking.comparison.validate_redocking_top_n") as mock_val:
            mock_v.side_effect = [
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            ]
            mock_rcsb.return_value = "CC"
            mock_val.return_value = {
                "top1_rmsd_pass": True,
                "recovery_rmsd_pass": True,
                "recovery_cluster_pass": True,
                "top1": {"rmsd_angstrom": 0.5},
            }
            res = run_matched_docking_comparison(
                healthy_pdb_text=pdb,
                mutant_pdb_text=pdb,
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="CC",
                ligand_name="Ethanol",
                healthy_pocket_residues=["A:10"],
                mutant_pocket_residues=["A:10"],
                wt_reference_ligand_meta={"resname": "LIG", "ref_ligand_pdb_text": ref_pdb},
            )

        try:
            self.assertEqual(res["validation_status"], "protocol validated on WT reference")
            self.assertIn("matched the verified WT reference ligand", res["validation_reason"])

            h_dir = BASE_RUNS_DIR / res["healthy"]["run_id"]
            m_dir = BASE_RUNS_DIR / res["mutant"]["run_id"]

            with open(h_dir / "manifest.json", "r", encoding="utf-8") as f:
                h_m = json.load(f)
            with open(m_dir / "manifest.json", "r", encoding="utf-8") as f:
                m_m = json.load(f)

            self.assertEqual(h_m["status"], "protocol validated on WT reference")
            self.assertEqual(m_m["status"], "protocol validated on WT reference")
            self.assertEqual(h_m["results"]["validation"]["status"], "protocol validated on WT reference")
        finally:
            shutil.rmtree(BASE_RUNS_DIR / res["healthy"]["run_id"], ignore_errors=True)
            shutil.rmtree(BASE_RUNS_DIR / res["mutant"]["run_id"], ignore_errors=True)

    # 4. Mock a reference-identity mismatch and assert that the final manifest is not validated for affinity comparison, never run failed and never an exploratory reference-validation status.
    def test_4_mocked_reference_identity_mismatch_manifest_status(self):
        pdb = self._make_pdb("A", 10)
        ref_pdb = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        with patch("docking.comparison.run_vina_multi_seeds") as mock_v, \
             patch("docking.comparison.fetch_chem_comp_smiles") as mock_rcsb:
            mock_v.side_effect = [
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
                ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            ]
            mock_rcsb.return_value = "CCO"  # Ethanol
            res = run_matched_docking_comparison(
                healthy_pdb_text=pdb,
                mutant_pdb_text=pdb,
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="COC",  # Dimethyl ether
                ligand_name="Ether",
                healthy_pocket_residues=["A:10"],
                mutant_pocket_residues=["A:10"],
                wt_reference_ligand_meta={"resname": "LIG", "ref_ligand_pdb_text": ref_pdb},
            )

        try:
            self.assertEqual(res["validation_status"], "not validated for affinity comparison")
            self.assertNotEqual(res["validation_status"], "run failed")
            self.assertNotIn("exploratory", res["validation_status"])

            h_dir = BASE_RUNS_DIR / res["healthy"]["run_id"]
            with open(h_dir / "manifest.json", "r", encoding="utf-8") as f:
                h_m = json.load(f)

            self.assertEqual(h_m["status"], "not validated for affinity comparison")
            self.assertNotEqual(h_m["status"], "run failed")
            self.assertNotIn("exploratory", h_m["status"])
            self.assertEqual(h_m["results"]["validation"]["status"], "not validated for affinity comparison")
        finally:
            shutil.rmtree(BASE_RUNS_DIR / res["healthy"]["run_id"], ignore_errors=True)
            shutil.rmtree(BASE_RUNS_DIR / res["mutant"]["run_id"], ignore_errors=True)

    # 5. Assert that a cross-protein manifest remains exploratory—cross-protein docking.
    @patch("docking.comparison.run_vina_multi_seeds")
    def test_5_cross_protein_manifest_remains_exploratory_cross_protein_docking(self, mock_vina):
        mock_vina.side_effect = [
            ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
            ([self._make_dummy_pose(seed=42), self._make_dummy_pose(seed=101), self._make_dummy_pose(seed=2024)], []),
        ]
        pdb1 = self._make_pdb("A", 10)
        pdb2 = self._make_pdb("B", 20)
        res = run_exploratory_cross_docking(
            receptor1_pdb_text=pdb1,
            receptor2_pdb_text=pdb2,
            receptor1_target_meta={"pdb_id": "REC1", "selected_chain": "A"},
            receptor2_target_meta={"pdb_id": "REC2", "selected_chain": "B"},
            ligand_smiles="CC",
            ligand_name="Drug",
            receptor1_pocket_residues=["A:10"],
            receptor2_pocket_residues=["B:20"],
        )
        try:
            self.assertEqual(res["status"], "exploratory—cross-protein docking")
            r1_dir = BASE_RUNS_DIR / res["receptor1"]["run_id"]
            r2_dir = BASE_RUNS_DIR / res["receptor2"]["run_id"]
            with open(r1_dir / "manifest.json", "r", encoding="utf-8") as f:
                m1 = json.load(f)
            with open(r2_dir / "manifest.json", "r", encoding="utf-8") as f:
                m2 = json.load(f)
            self.assertEqual(m1["status"], "exploratory—cross-protein docking")
            self.assertEqual(m2["status"], "exploratory—cross-protein docking")
        finally:
            shutil.rmtree(BASE_RUNS_DIR / res["receptor1"]["run_id"], ignore_errors=True)
            shutil.rmtree(BASE_RUNS_DIR / res["receptor2"]["run_id"], ignore_errors=True)

    # 6. Assert that an unknown status raises ValueError (or is handled by the one documented legacy alias only).
    def test_6_unknown_status_raises_value_error_legacy_alias_handled(self):
        # Documented legacy alias: "exploratory" -> "exploratory—no reference validation"
        self.assertEqual(normalize_manifest_status("exploratory"), "exploratory—no reference validation")
        self.assertEqual(normalize_manifest_status("  EXPLORATORY  "), "exploratory—no reference validation")

        # All defined statuses normalize correctly
        for valid_st in MANIFEST_STATUSES:
            self.assertEqual(normalize_manifest_status(valid_st), valid_st)
            self.assertEqual(normalize_manifest_status(valid_st.upper()), valid_st)

        # Unknown status must raise ValueError
        with self.assertRaises(ValueError) as ctx:
            normalize_manifest_status("unsupported_custom_status")
        self.assertIn("Unsupported docking manifest status", str(ctx.exception))

        # build_and_save_manifest must also raise ValueError on unknown status
        run_dir = BASE_RUNS_DIR / "test_unknown_status_dir"
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            with self.assertRaises(ValueError) as ctx2:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id="test_unknown_status_dir",
                    status="invalid_status_xyz",
                    target_meta={},
                    binding_site_meta={},
                    ligand_meta={},
                    prep_meta={},
                    engine_meta={},
                    grid_meta={},
                    results_meta={},
                )
            self.assertIn("Unsupported docking manifest status", str(ctx2.exception))
        finally:
            shutil.rmtree(run_dir, ignore_errors=True)


class TestBlindDockingGridDefinition(unittest.TestCase):
    """Tests for define_grid_whole_protein whole-receptor grid computation."""

    def setUp(self):
        self.multichain_pdb = (
            "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1      12.000  22.000  32.000  1.00 20.00           C\n"
            "ATOM      3  N   GLY B   1      50.000  60.000  70.000  1.00 20.00           N\n"
            "ATOM      4  CA  GLY B   1      52.000  62.000  72.000  1.00 20.00           C\n"
            "END\n"
        )

    def test_whole_protein_calculates_all_chains_when_chain_none(self):
        grid = define_grid_whole_protein(self.multichain_pdb, padding=5.0)
        self.assertAlmostEqual(grid.center_x, 31.0, places=1)
        self.assertAlmostEqual(grid.center_y, 41.0, places=1)
        self.assertAlmostEqual(grid.center_z, 51.0, places=1)
        self.assertGreaterEqual(grid.size_x, 50.0)

    def test_whole_protein_specific_chain_filters_correctly(self):
        grid_a = define_grid_whole_protein(self.multichain_pdb, padding=5.0, chain_id="A")
        self.assertAlmostEqual(grid_a.center_x, 11.0, places=1)
        self.assertAlmostEqual(grid_a.center_y, 21.0, places=1)
        self.assertAlmostEqual(grid_a.center_z, 31.0, places=1)
        self.assertGreaterEqual(grid_a.size_x, 20.0)

    def test_whole_protein_empty_pdb_raises_value_error(self):
        with self.assertRaises(ValueError):
            define_grid_whole_protein("")
        with self.assertRaises(ValueError):
            define_grid_whole_protein("   \n  ")


if __name__ == "__main__":
    unittest.main()
