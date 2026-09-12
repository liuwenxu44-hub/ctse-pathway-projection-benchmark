"""Deterministic paired, unit-first summaries; no fits, p-values or winner scores.

The rowwise contrast is formed BEFORE averaging. In particular a difference of
unpaired arm means is never substituted for the mean of paired differences.
Only legacy endpoint quantities and their explicitly labelled contrasts occur.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REGIMES = ["ORIGINAL_POOLED_4800", "STATE0_ONLY_2400"]
METHODS = ["bayesprism", "unico", "tca", "epic", "reference", "bulk"]
CLASSES = ["negative", "zero", "positive"]
UNIT = ["task_id", "scenario", "replicate_id", "complete_four_both_regimes"]
DESIGN = ["celltype", "celltype_design_changed"]


def require(ok, code):
    if not bool(ok):
        raise ValueError(code)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=["NA"])


def unique(frame, keys, name):
    require(not frame.duplicated(keys).any(), f"DUPLICATE_{name}_KEY")


def row_keys(frame, sample):
    return UNIT + (["context"] if "context" in frame else []) + DESIGN + ["gene_stratum"] + (
        ["group", "sample_id"] if sample else [])


def dimensions(frame):
    return [x for x in ["reference_regime", "context", "method", "baseline", "gene_stratum",
                        "metric", "difference_definition"] if x in frame]


def aggregate(frame, value="value", sample=True):
    """The original equal-group then equal-celltype aggregation, with coverage."""
    keys = UNIT + dimensions(frame)
    if sample:
        bygroup = frame.groupby(keys + DESIGN + ["group"], dropna=False, sort=True)[value].agg(
            value="mean", evaluable_samples="count", planned_samples="size").reset_index()
        # A type with any evaluable data still needs both originally defined groups.
        bycell = bygroup.groupby(keys + DESIGN, dropna=False, sort=True).agg(
            value=("value", lambda x: x.mean() if x.notna().sum() == 2 else np.nan),
            evaluable_groups=("value", "count"), planned_groups=("value", "size"),
            evaluable_samples=("evaluable_samples", "sum"), planned_samples=("planned_samples", "sum")
        ).reset_index()
        require((bycell.planned_groups == 2).all(), "ORIGINAL_GROUP_SLOT_MISSING")
    else:
        bygroup = pd.DataFrame()
        unique(frame, keys + DESIGN, "NON_SAMPLE_RECORD")
        bycell = frame[keys + DESIGN + [value]].rename(columns={value: "value"}).copy()
    unit = bycell.groupby(keys, dropna=False, sort=True).agg(
        value=("value", "mean"), evaluable_celltypes=("value", "count"),
        planned_celltypes=("value", "size")).reset_index()
    if sample:
        counts = bycell.groupby(keys, dropna=False, sort=True)[["evaluable_samples", "planned_samples"]].sum().reset_index()
        unit = unit.merge(counts, on=keys, validate="one_to_one")
    unit["independent_unit"] = "original_simulation_unit"
    return bygroup, bycell, unit


def distribution(unit):
    keys = ["scenario"] + dimensions(unit) + (["celltype_design_changed"] if "celltype_design_changed" in unit else [])
    roster = unit[keys].drop_duplicates()
    out = []
    for scope, z in [("ALL_120_REGISTERED_SLOTS", unit),
                     ("COMPLETE_FOUR_METHODS_BOTH_REGIMES", unit[unit.complete_four_both_regimes])]:
        x = z.groupby(keys, dropna=False, sort=True).agg(
            planned_units=("value", "size"), evaluable_units=("value", "count"),
            mean=("value", "mean"), median=("value", "median"),
            q25=("value", lambda v: v.quantile(.25)), q75=("value", lambda v: v.quantile(.75))
        ).reset_index()
        x = roster.merge(x, on=keys, how="left", validate="one_to_one")
        for column in ["planned_units", "evaluable_units"]:
            x[column] = x[column].fillna(0).astype(int)
        x["scope"] = scope
        x["independent_unit"] = "original_simulation_unit"
        out.append(x)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def paired_baseline(frame, value, sample=True):
    keys = row_keys(frame, sample) + ["reference_regime"]
    unique(frame, keys + ["method"], "ENDPOINT")
    results = []
    for baseline in ["reference", "bulk"]:
        b = frame[frame.method == baseline][keys + [value, "genes"]].rename(
            columns={value: "baseline_value", "genes": "baseline_genes"})
        unique(b, keys, "BASELINE")
        z = frame[frame.method != baseline].merge(b, on=keys, how="left", validate="many_to_one", indicator=True)
        require((z._merge == "both").all(), "BASELINE_SLOT_MISSING")
        require((z.genes == z.baseline_genes).all(), "BASELINE_PAIRED_SUPPORT_MISMATCH")
        z["value"] = z[value] - z.baseline_value
        z["baseline"] = baseline
        z["difference_definition"] = "method_minus_baseline_" + value
        results.append(z.drop(columns=["_merge"]))
    return pd.concat(results, ignore_index=True)


def paired_regime(frame, value, sample=True, baseline=False):
    keys = row_keys(frame, sample) + ["method"] + (["baseline"] if baseline else [])
    unique(frame, keys + ["reference_regime"], "REGIME")
    old = frame[frame.reference_regime == REGIMES[0]][keys + [value, "genes"]].rename(
        columns={value: "original_value", "genes": "original_genes"})
    new = frame[frame.reference_regime == REGIMES[1]].copy()
    z = new.merge(old, on=keys, how="outer", validate="one_to_one", indicator=True)
    require((z._merge == "both").all(), "REGIME_SLOT_MISSING")
    require((z.genes == z.original_genes).all(), "REGIME_PAIRED_SUPPORT_MISMATCH")
    z["value"] = z[value] - z.original_value
    z["reference_regime"] = "PAIRED_STATE0_MINUS_ORIGINAL"
    z["difference_definition"] = ("state0_minus_original_of_method_minus_baseline_" if baseline else
                                    "state0_minus_original_") + (str(frame.metric.iloc[0]) if "metric" in frame else value)
    return z.drop(columns=["_merge"])


def endpoint_summaries(frame, metrics, sample=True):
    """Both arms, row-paired regime changes, baselines and paired baseline changes."""
    outputs = {}
    chunks = {name: [] for name in ["CELLTYPE_GROUP", "CELLTYPE_GROUPBALANCED", "UNIT", "SCENARIO",
                                   "UNIT_BY_DESIGN_ROLE", "SCENARIO_BY_DESIGN_ROLE",
                                   "PAIRED_BASELINE_UNIT", "PAIRED_BASELINE_SCENARIO",
                                   "PAIRED_REGIME_UNIT", "PAIRED_REGIME_SCENARIO",
                                   "PAIRED_BASELINE_REGIME_UNIT", "PAIRED_BASELINE_REGIME_SCENARIO"]}
    for metric in metrics:
        z = frame.copy(); z["metric"] = metric; z["value"] = z[metric]
        a, c, u = aggregate(z, sample=sample)
        if not a.empty:
            chunks["CELLTYPE_GROUP"].append(a)
        chunks["CELLTYPE_GROUPBALANCED"].append(c)
        chunks["UNIT"].append(u); chunks["SCENARIO"].append(distribution(u))
        role = c.groupby(UNIT + dimensions(c) + ["celltype_design_changed"], dropna=False, sort=True).agg(
            value=("value", "mean"), evaluable_celltypes=("value", "count"), planned_celltypes=("value", "size")
        ).reset_index()
        role["independent_unit"] = "original_simulation_unit"
        chunks["UNIT_BY_DESIGN_ROLE"].append(role)
        chunks["SCENARIO_BY_DESIGN_ROLE"].append(distribution(role))
        b = paired_baseline(z, metric, sample)
        _, _, bu = aggregate(b, sample=sample)
        chunks["PAIRED_BASELINE_UNIT"].append(bu)
        chunks["PAIRED_BASELINE_SCENARIO"].append(distribution(bu))
        r = paired_regime(z, metric, sample)
        _, _, ru = aggregate(r, sample=sample)
        chunks["PAIRED_REGIME_UNIT"].append(ru)
        chunks["PAIRED_REGIME_SCENARIO"].append(distribution(ru))
        br = paired_regime(b, "value", sample, baseline=True)
        _, _, bru = aggregate(br, sample=sample)
        chunks["PAIRED_BASELINE_REGIME_UNIT"].append(bru)
        chunks["PAIRED_BASELINE_REGIME_SCENARIO"].append(distribution(bru))
    for key, value in chunks.items():
        if value:
            outputs[key] = pd.concat(value, ignore_index=True)
    return outputs


def legacy_direction_summary(z):
    """Same original class_summary formulas, restricted to one independent unit.

    Invalid truth is kept in coverage but cannot count as a correct prediction.
    Absent true classes leave macro recall undefined, as in the original code.
    """
    truth = z.truth_direction; pred = z.estimated_direction
    valid_truth = truth.isin(CLASSES); ok = pred.isin(CLASSES) & valid_truth
    recalls = []
    for cl in CLASSES:
        take = (truth == cl) & ok
        recalls.append((pred[take] == cl).mean() if take.any() else np.nan)
    nz = truth.isin(["negative", "positive"]); zero = truth == "zero"
    return pd.Series(dict(
        planned_records=len(z), valid_truth_records=int(valid_truth.sum()), evaluable_records=int(ok.sum()),
        coverage=float(ok.mean()), overall_accuracy=(pred[ok] == truth[ok]).mean() if ok.any() else np.nan,
        balanced_accuracy_macro_recall=np.mean(recalls) if np.isfinite(recalls).all() else np.nan,
        nonzero_planned=int(nz.sum()), nonzero_evaluable=int((nz & ok).sum()),
        nonzero_sign_accuracy=(pred[nz & ok] == truth[nz & ok]).mean() if (nz & ok).any() else np.nan,
        zero_planned=int(zero.sum()), zero_evaluable=int((zero & ok).sum()),
        zero_truth_recall=(pred[zero & ok] == "zero").mean() if (zero & ok).any() else np.nan,
        recall_negative=recalls[0], recall_zero=recalls[1], recall_positive=recalls[2]))


def direction_summaries(frame):
    keys = UNIT + ["reference_regime", "context", "method", "pathway"]
    unique(frame, keys + DESIGN, "DIRECTION")
    units = []
    confusion = []
    for key, z in frame.groupby(keys, dropna=False, sort=True):
        info = dict(zip(keys, key)); units.append({**info, **legacy_direction_summary(z).to_dict()})
        for cl in CLASSES:
            for prediction in CLASSES + ["NO_OUTPUT"]:
                take = (z.truth_direction == cl) & ((~z.estimated_direction.isin(CLASSES)) if
                        prediction == "NO_OUTPUT" else (z.estimated_direction == prediction))
                confusion.append({**info, "truth_class": cl, "prediction_class": prediction,
                                  "records": int(take.sum()), "independent_unit": "original_simulation_unit"})
    u = pd.DataFrame(units)
    metriccols = ["coverage", "overall_accuracy", "balanced_accuracy_macro_recall", "nonzero_sign_accuracy",
                  "zero_truth_recall", "recall_negative", "recall_zero", "recall_positive"]
    long = u.melt(id_vars=keys, value_vars=metriccols, var_name="metric", value_name="value")
    # Pathways remain separate. No cell or pathway record is an independent n.
    long["gene_stratum"] = long.pop("pathway")
    scenario = distribution(long).rename(columns={"gene_stratum": "pathway"})
    counts = pd.DataFrame(confusion)
    pooled_counts = counts.groupby(["reference_regime", "context", "scenario", "method", "pathway",
                                   "truth_class", "prediction_class"], sort=True).agg(
        records=("records", "sum"), planned_units=("task_id", "nunique")).reset_index()
    pooled_counts["interpretation"] = "RECORD_COUNTS_ONLY_NOT_INDEPENDENT_REPLICATES"
    # Paired label accuracy on the same cell/path rows. This is existing overall
    # accuracy (0/1 correctness averaged within a unit), not a new composite score.
    z = frame.copy()
    z["overall_accuracy"] = np.where(z.truth_direction.isin(CLASSES) & z.estimated_direction.isin(CLASSES),
                                      (z.truth_direction == z.estimated_direction).astype(float), np.nan)
    z["gene_stratum"] = z.pop("pathway")
    pairs = endpoint_summaries(z, ["overall_accuracy"], sample=False)
    outputs = {"UNIT_CLASS_ACCURACY_COVERAGE": u, "SCENARIO_UNIT_DISTRIBUTION": scenario,
               "UNIT_CONFUSION_COUNTS": counts, "SCENARIO_RECORD_CONFUSION_COUNTS": pooled_counts}
    for name, val in pairs.items():
        if name.startswith("PAIRED_"):
            outputs[name] = val.rename(columns={"gene_stratum": "pathway"})
    return outputs


def validate_availability(av, expected_units=120):
    require(set(av.reference_regime) == set(REGIMES) and set(av.method) == set(METHODS), "AVAILABILITY_ROSTER")
    unique(av, ["task_id", "reference_regime", "method"], "AVAILABILITY")
    require(av.task_id.nunique() == expected_units and len(av) == expected_units * 12, "PLANNED_SLOTS_DROPPED")
    require((av.groupby(["reference_regime", "method"]).size() == expected_units).all(), "METHOD_SLOTS_DROPPED")
    require((av.paired_common_genes > 0).all() and
            (av.paired_common_genes <= av.own_regime_common_genes).all(), "SUPPORT_COUNTS_INVALID")
    require((av.groupby("task_id").paired_common_genes.nunique() == 1).all(), "SUPPORT_NOT_SHARED_ACROSS_REGIMES_METHODS")
    require((av.loc[~av.available, "returned_genes"] == 0).all(), "UNAVAILABLE_HAS_RETURNED_GENES")
    require(av.loc[~av.available, "failure"].fillna("").str.len().gt(0).all(), "UNAVAILABLE_REASON_MISSING")
    method = av[av.method.isin(METHODS[:4])]
    flags = method.groupby("task_id").available.all()
    require((av.complete_four_both_regimes == av.task_id.map(flags)).all(), "COMPLETE_FOUR_FLAG_WRONG")


def validate_long_roster(frame, av, sample=True, direction=False):
    """Fail closed on missing unit/method/group/type/stratum slots, even if NA."""
    require(set(frame.task_id) == set(av.task_id), "ENDPOINT_UNIT_ROSTER_INCOMPLETE")
    require(set(frame.reference_regime) == set(REGIMES), "ENDPOINT_REGIME_ROSTER_INCOMPLETE")
    stratum = "pathway" if direction else "gene_stratum"
    expected_strata = {"signal", "negative_control"} if direction else {
        "ALL", "DIRECT_SIGNAL", "NEGATIVE_CONTROL", "BACKGROUND", "DESIGN_CHANGED", "DESIGN_UNCHANGED"}
    require(set(frame[stratum]) == expected_strata, "ENDPOINT_STRATUM_ROSTER_CHANGED")
    contexts = {"LINEAR_CPM_SECONDARY": {"unico", "tca", "reference", "bulk"},
                "EPIC_NATIVE_LOG_SECONDARY": {"epic", "reference", "bulk"}}
    groups = frame.groupby("context", sort=True) if "context" in frame else [(None, frame)]
    if "context" in frame:
        require(set(frame.context) == set(contexts), "SCALE_CONTEXT_ROSTER_CHANGED")
    for ctx, z in groups:
        methods = set(METHODS) if ctx is None else contexts[ctx] | ({"always_zero"} if direction else set())
        require(set(z.method) == methods, "ILLEGAL_SCALE_METHOD_COMPARISON")
        keys = ["task_id", "reference_regime", "method", stratum]
        counts = z.groupby(keys, sort=True).size()
        require(len(counts) == av.task_id.nunique() * 2 * len(methods) * len(expected_strata) and
                (counts == (160 if sample else 4)).all(), "ENDPOINT_PLANNED_RECORDS_DROPPED")
        unique(z, keys + ["celltype"] + (["sample_id"] if sample else []), "LONG_RECORD")
        require((z.groupby(keys).celltype.nunique() == 4).all(), "ENDPOINT_CELLTYPE_SLOTS_CHANGED")
        if sample:
            require(set(z.group) == {"group0", "group1"} and
                    (z.groupby(keys + ["group"]).size() == 80).all() and
                    (z.groupby(keys).sample_id.nunique() == 40).all(), "ENDPOINT_GROUP_SAMPLE_SLOTS_CHANGED")
        flags = av.drop_duplicates("task_id").set_index("task_id").complete_four_both_regimes
        require((z.complete_four_both_regimes == z.task_id.map(flags)).all(), "ENDPOINT_COMPLETENESS_FLAG_DRIFT")


def main(directory):
    d = Path(directory).resolve(); require((d / "DERIVED_BUILD_COMPLETE").is_file(), "ENDPOINT_BUILD_INCOMPLETE")
    out = d / "SUMMARIES"; require(not out.exists(), "OUTPUT_EXISTS_NO_OVERWRITE")
    names = ["METHOD_AVAILABILITY.tsv", "GENE_RANK_RECOVERY_LONG.tsv", "SAME_SCALE_ERROR_LONG.tsv",
             "SCALE_SPECIFIC_CONTRAST_LONG.tsv", "SCALE_SPECIFIC_DIRECTION_LONG.tsv", "ACTUAL_ENDPOINT_READS.tsv"]
    # Paths are relative to the derived directory so independent A/B rebuilds
    # with identical scientific inputs also have byte-identical provenance.
    provenance = [{"path": name, "sha256": sha(d / name)} for name in names]
    av = read_tsv(d / names[0]); validate_availability(av)
    old_fail = av[(av.reference_regime == REGIMES[0]) & ~av.available]
    require(len(old_fail) == 2 and set(old_fail.method) == {"tca"} and
            set(old_fail.failure) == {"PRESERVED_ORIGINAL_TCA_FAILURE"}, "ORIGINAL_FAILURE_SLOTS_CHANGED")
    require(av.scenario.nunique() == 6 and
            (av.drop_duplicates("task_id").groupby("scenario").size() == 20).all(), "SCENARIO_DENOMINATORS_CHANGED")
    out.mkdir()
    def save(df, name):
        df.to_csv(out / name, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")
    coverage = av.groupby(["reference_regime", "scenario", "method"], sort=True).agg(
        planned_units=("task_id", "size"), available_units=("available", "sum"),
        min_returned=("returned_genes", "min"), max_returned=("returned_genes", "max"),
        min_own_regime_common=("own_regime_common_genes", "min"), max_own_regime_common=("own_regime_common_genes", "max"),
        min_paired_common=("paired_common_genes", "min"), max_paired_common=("paired_common_genes", "max"),
        complete_four_this_regime=("complete_four_this_regime", "sum"),
        complete_four_both_regimes=("complete_four_both_regimes", "sum")).reset_index()
    coverage["coverage"] = coverage.available_units / coverage.planned_units
    save(coverage, "METHOD_AVAILABILITY_SUMMARY.tsv")
    save(av[UNIT].drop_duplicates().sort_values("task_id"), "COMPLETE_FOUR_PAIRED_UNIT_ROSTER.tsv")
    for filename, prefix, metrics, sample in [
        ("GENE_RANK_RECOVERY_LONG.tsv", "RANK", ["spearman"], True),
        ("SAME_SCALE_ERROR_LONG.tsv", "ERROR", ["rmse", "sne"], True),
        ("SCALE_SPECIFIC_CONTRAST_LONG.tsv", "CONTRAST", ["truth_contrast_L2", "estimated_contrast_L2", "contrast_error_L2"], False)]:
        frame = read_tsv(d / filename)
        validate_long_roster(frame, av, sample)
        for name, table in endpoint_summaries(frame, metrics, sample).items():
            save(table, prefix + "_" + name + ".tsv")
        print(prefix, "UNIT_FIRST_PAIRED_SUMMARIES_COMPLETE", flush=True)
    dr = read_tsv(d / "SCALE_SPECIFIC_DIRECTION_LONG.tsv")
    validate_long_roster(dr, av, sample=False, direction=True)
    for name, table in direction_summaries(dr).items():
        save(table, "DIRECTION_" + name + ".tsv")
    for entry in provenance:
        require(sha(d / entry["path"]) == entry["sha256"], "ENDPOINT_CHANGED_DURING_SUMMARY")
    metadata = dict(status="DETERMINISTIC_DERIVED_ONLY", inputs=provenance,
                    source_path="study/summarize_reference.py", source_sha256=sha(__file__),
                    versions={"numpy": np.__version__, "pandas": pd.__version__},
                    independent_unit="original_simulation_unit", pair_order="form rowwise differences before unit aggregation",
                    scopes=["all 120 planned unit slots", "complete four methods in both regimes"],
                    limitations=["Reference state coverage and reference cell count both differ.",
                                 "Rank is gene-vector shape recovery; BayesPrism is not amplitude-comparable.",
                                 "Old-arm endpoints use the new paired support and do not replace frozen primary tables.",
                                 "No new p-values, no pooled gene/sample/pathway replicate count, no unified scale ranking."])
    (out / "SUMMARY_PROVENANCE.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (out / "SUMMARY_COMPLETE").write_text("DETERMINISTIC_PAIRED_UNIT_FIRST_NO_FITS_NO_NEW_TESTS\n")
    print("REFERENCE_SUMMARY_COMPLETE", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("derived_directory")
    main(parser.parse_args().derived_directory)
