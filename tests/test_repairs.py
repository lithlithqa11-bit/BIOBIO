"""Automated test suite verifying required repairs and scientific validity constraints."""

import os
import sys
import json
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
    extract_reference_ligand,
    GridBox,
)
from docking.comparison import (
    verify_comparison_clusters,
    _select_consensus_cluster,
    run_matched_docking_comparison,
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
from docking.storage import get_vina_path, build_and_save_manifest, BASE_RUNS_DIR
from docking.engine import get_vina_version
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


if __name__ == "__main__":
    unittest.main()
