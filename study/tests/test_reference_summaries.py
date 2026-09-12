"""Artificial-only tests, optionally including TSVs from the artificial R test."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import summarize_reference as s


def fixture():
    rows = []
    for unit in ["u1", "u2"]:
        for regime in s.REGIMES:
            for method in s.METHODS:
                for cell in ["A", "B"]:
                    for group in ["group0", "group1"]:
                        for sample in [1, 2]:
                            value = dict(bayesprism=.5, unico=.7, tca=.6, epic=.8, reference=.3, bulk=.2)[method]
                            if regime == s.REGIMES[1]:
                                value += .05 if method == "reference" else .1
                            rows.append(dict(task_id=unit, scenario="artificial", replicate_id=int(unit[-1]),
                                complete_four_both_regimes=unit == "u2", reference_regime=regime, method=method,
                                celltype=cell, celltype_design_changed=cell == "A", group=group,
                                sample_id=f"{group}_s{sample}", gene_stratum="ALL", genes=20, spearman=value))
    df = pd.DataFrame(rows)
    df.loc[(df.task_id == "u1") & (df.reference_regime == s.REGIMES[0]) & (df.method == "tca"), "spearman"] = np.nan
    return df


class EndpointSummaryTests(unittest.TestCase):
    def test_equal_group_then_equal_cell(self):
        x = fixture()
        tables = s.endpoint_summaries(x, ["spearman"])
        u = tables["UNIT"]
        self.assertEqual(len(u), 24)
        self.assertTrue(u.independent_unit.eq("original_simulation_unit").all())
        self.assertTrue(u.loc[(u.method == "tca") & (u.task_id == "u1") &
            (u.reference_regime == s.REGIMES[0]), "value"].isna().all())
        b = tables["PAIRED_BASELINE_UNIT"]
        row = b[(b.task_id == "u2") & (b.method == "unico") & (b.baseline == "reference") &
                (b.reference_regime == s.REGIMES[0])].iloc[0]
        self.assertAlmostEqual(row.value, .4)
        br = tables["PAIRED_BASELINE_REGIME_UNIT"]
        row = br[(br.task_id == "u2") & (br.method == "unico") & (br.baseline == "reference")].iloc[0]
        self.assertAlmostEqual(row.value, .05)
        d = tables["SCENARIO"]
        self.assertTrue(d.loc[d.scope == "ALL_120_REGISTERED_SLOTS", "planned_units"].eq(2).all())
        self.assertTrue(d.loc[d.scope == "COMPLETE_FOUR_METHODS_BOTH_REGIMES", "planned_units"].eq(1).all())

    def test_row_pairing_precedes_aggregation(self):
        x = fixture()
        # Disjoint arm availability within every group: no paired sample exists.
        mask = (x.task_id == "u2") & (x.method == "unico")
        x.loc[mask & (x.reference_regime == s.REGIMES[0]) & x.sample_id.str.endswith("s1"), "spearman"] = np.nan
        x.loc[mask & (x.reference_regime == s.REGIMES[1]) & x.sample_id.str.endswith("s2"), "spearman"] = np.nan
        tables = s.endpoint_summaries(x, ["spearman"])
        arms = tables["UNIT"]
        self.assertTrue(arms.loc[(arms.task_id == "u2") & (arms.method == "unico"), "value"].notna().all())
        paired = tables["PAIRED_REGIME_UNIT"]
        self.assertTrue(paired.loc[(paired.task_id == "u2") & (paired.method == "unico"), "value"].isna().all())

    def test_missing_one_group_is_not_single_group_mean(self):
        x = fixture()
        mask = (x.task_id == "u2") & (x.method == "unico") & (x.group == "group0")
        x.loc[mask, "spearman"] = np.nan
        tables = s.endpoint_summaries(x, ["spearman"])
        u = tables["UNIT"]
        self.assertTrue(u.loc[(u.task_id == "u2") & (u.method == "unico"), "value"].isna().all())

    def test_duplicate_or_missing_pair_rejected(self):
        x = fixture()
        with self.assertRaisesRegex(ValueError, "DUPLICATE"):
            s.paired_regime(pd.concat([x, x.iloc[[0]]]), "spearman")
        with self.assertRaisesRegex(ValueError, "REGIME_SLOT_MISSING"):
            s.paired_regime(x.iloc[1:], "spearman")
        with self.assertRaisesRegex(ValueError, "BASELINE_SLOT_MISSING"):
            s.paired_baseline(x[x.method != "bulk"], "spearman")

    def test_support_mismatch_rejected(self):
        x = fixture(); x.loc[x.index[0], "genes"] = 19
        with self.assertRaisesRegex(ValueError, "SUPPORT_MISMATCH"):
            s.paired_regime(x, "spearman")
        with self.assertRaisesRegex(ValueError, "SUPPORT_MISMATCH"):
            s.paired_baseline(x, "spearman")

    def test_existing_direction_class_semantics(self):
        x = pd.DataFrame(dict(truth_direction=s.CLASSES, estimated_direction=["zero"] * 3))
        y = s.legacy_direction_summary(x)
        self.assertAlmostEqual(y.balanced_accuracy_macro_recall, 1/3)
        self.assertEqual(y.nonzero_sign_accuracy, 0)
        self.assertEqual(y.zero_truth_recall, 1)
        x.loc[1, "estimated_direction"] = np.nan
        y = s.legacy_direction_summary(x)
        self.assertEqual(y.evaluable_records, 2)
        self.assertTrue(np.isnan(y.balanced_accuracy_macro_recall))
        self.assertTrue(np.isnan(y.zero_truth_recall))

    def test_both_120_slots_and_failures(self):
        rows = []
        for i in range(120):
            for regime in s.REGIMES:
                for method in s.METHODS:
                    missing = i < 2 and regime == s.REGIMES[0] and method == "tca"
                    rows.append(dict(task_id=f"u{i}", reference_regime=regime, method=method,
                        paired_common_genes=19, own_regime_common_genes=20, returned_genes=0 if missing else 20,
                        available=not missing, failure="PRESERVED_ORIGINAL_TCA_FAILURE" if missing else "",
                        complete_four_both_regimes=i >= 2))
        x = pd.DataFrame(rows); s.validate_availability(x)
        with self.assertRaisesRegex(ValueError, "SLOTS_DROPPED"):
            s.validate_availability(x.iloc[1:])

    def test_real_shape_slot_validator_artificial_data(self):
        base = fixture().drop_duplicates(["task_id", "reference_regime", "method"])
        rows = []
        for row in base.to_dict("records"):
            for cell in ["A", "B", "C", "D"]:
                for group in ["group0", "group1"]:
                    for sample in range(20):
                        for stratum in ["ALL", "DIRECT_SIGNAL", "NEGATIVE_CONTROL", "BACKGROUND", "DESIGN_CHANGED", "DESIGN_UNCHANGED"]:
                            rows.append({**row, "celltype": cell, "group": group,
                                         "sample_id": f"{group}_s{sample}", "gene_stratum": stratum})
        x = pd.DataFrame(rows)
        av = base[["task_id", "complete_four_both_regimes"]].drop_duplicates()
        s.validate_long_roster(x, av)
        with self.assertRaisesRegex(ValueError, "PLANNED_RECORDS_DROPPED"):
            s.validate_long_roster(x.iloc[1:], av)
        x.loc[0, "complete_four_both_regimes"] = not x.loc[0, "complete_four_both_regimes"]
        with self.assertRaisesRegex(ValueError, "COMPLETENESS_FLAG_DRIFT"):
            s.validate_long_roster(x, av)


def exercise_artificial_r(directory):
    d = Path(directory)
    av = s.read_tsv(d / "METHOD_AVAILABILITY.tsv"); s.validate_availability(av, expected_units=1)
    for filename, metrics, sample in [
        ("GENE_RANK_RECOVERY_LONG.tsv", ["spearman"], True),
        ("SAME_SCALE_ERROR_LONG.tsv", ["rmse", "sne"], True),
        ("SCALE_SPECIFIC_CONTRAST_LONG.tsv", ["truth_contrast_L2", "estimated_contrast_L2", "contrast_error_L2"], False)]:
        frame = s.read_tsv(d / filename)
        out = s.endpoint_summaries(frame, metrics, sample)
        assert out["UNIT"].task_id.nunique() == 1
        assert out["PAIRED_REGIME_UNIT"].reference_regime.eq("PAIRED_STATE0_MINUS_ORIGINAL").all()
    dr = s.read_tsv(d / "SCALE_SPECIFIC_DIRECTION_LONG.tsv")
    out = s.direction_summaries(dr)
    assert out["UNIT_CLASS_ACCURACY_COVERAGE"].task_id.nunique() == 1
    print("ARTIFICIAL_R_TO_PYTHON_TABLE_SCHEMA_PASS", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--artificial-r-directory")
    parser.add_argument("--evidence-output")
    args, remaining = parser.parse_known_args()
    result = unittest.main(argv=[__file__, *remaining], exit=False).result
    if not result.wasSuccessful():
        raise SystemExit(1)
    if args.artificial_r_directory:
        exercise_artificial_r(args.artificial_r_directory)
    if args.evidence_output:
        output = Path(args.evidence_output)
        assert not output.exists(), "TEST_EVIDENCE_EXISTS_NO_OVERWRITE"
        source_paths = [Path(__file__).resolve(), Path(s.__file__).resolve()]
        evidence = dict(status="PASS", python_tests=result.testsRun,
                        completed_utc=datetime.now(timezone.utc).isoformat(),
                        artificial_r_schema="PASS" if args.artificial_r_directory else "NOT_RUN",
                        real_model_fits=0, real_nnls_executions=0, real_endpoint_units=0,
                        sources=[dict(path=str(p), sha256=s.sha(p)) for p in source_paths],
                        versions=dict(numpy=np.__version__, pandas=pd.__version__))
        if args.artificial_r_directory:
            evidence["r_test_evidence"] = [dict(path=str(p.resolve()), sha256=s.sha(p)) for p in
                [Path(args.artificial_r_directory) / "ARTIFICIAL_R_ENDPOINT_TESTS.tsv",
                 Path(args.artificial_r_directory) / "ARTIFICIAL_TEST_SOURCE_HASHES.tsv"]]
        with output.open("x") as f:
            json.dump(evidence, f, indent=2); f.write("\n")
    print("ARTIFICIAL_REFERENCE_SUMMARY_PASS REAL_FITS_0_REAL_ENDPOINT_UNITS_0", flush=True)
