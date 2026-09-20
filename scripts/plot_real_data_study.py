"""Render real-study figures and transparent metric denominators from saved results."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from spatial_collab.store import Project


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--style-helper", type=Path)
    args = parser.parse_args()
    helper = None
    if args.style_helper:
        spec = importlib.util.spec_from_file_location("scientific_figure_pro", args.style_helper)
        helper = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = helper
        spec.loader.exec_module(helper)
        helper.apply_publication_style(helper.FigureStyle(font_size=11, axes_linewidth=1.2))
    else:
        plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
    def load(name):
        return json.loads((args.study / name).read_text(encoding="utf-8"))
    receipt, plan = load("study-receipt.json"), load("frozen-plan.json")
    cells = Project(receipt["working_project"]).cells()
    runs = {scope: load(f"sensitivity-{scope}.json") for scope in ("whole_slice", "roi_induced")}
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    a, b, c, d = axes.flat
    colors = {"B_program": "#0072B2", "T_program": "#D55E00", "Mixed_program": "#8D5AA5", "Unresolved": "#C8CDD2"}
    for label in ("Unresolved", "B_program", "T_program", "Mixed_program"):
        items = [x for x in cells if x["label"] == label]
        a.scatter([x["x"] for x in items], [x["y"] for x in items], s=3, c=colors[label],
                  label=f"{label} (n={len(items):,})", linewidths=0, rasterized=True)
    roi = plan["analysis_roi_um"]
    for ax in (a, b):
        ax.add_patch(Rectangle((roi[0], roi[1]), roi[2]-roi[0], roi[3]-roi[1], fill=False, linewidth=1.4, edgecolor="black"))
        ax.set(xlim=(2700, 3300), ylim=(3300, 2700), xlabel="Native x (µm)", ylabel="Native y (µm)", aspect="equal")
    a.set_title("a  Declared marker programs, not cell-type truth", loc="left", fontsize=12)
    a.legend(markerscale=3, fontsize=9, loc="upper left", bbox_to_anchor=(0, -.2), ncol=2)
    b.scatter([x["x"] for x in cells], [x["y"] for x in cells], s=3, color="#D8DCDF", linewidths=0, rasterized=True)
    for nucleus, color, caption in [(0, "#CC79A7", "No assigned nucleus (n=98)"), (2, "#009E73", "More than one nucleus (n=492)")]:
        items = [x for x in cells if (x["attributes"]["nucleus_count"] == 0 if nucleus == 0 else x["attributes"]["nucleus_count"] > 1)]
        b.scatter([x["x"] for x in items], [x["y"] for x in items], s=6, color=color, linewidths=0, label=caption, rasterized=True)
    b.set_title("b  Source QC observations; no error verdict", loc="left", fontsize=12)
    b.legend(markerscale=2, fontsize=9, loc="upper left", bbox_to_anchor=(0, -.2))
    variants = {"non_single_nucleus_uncertain": ("Nucleus count ≠ 1: uncertain", "#CC79A7"),
                "non_single_nucleus_excluded": ("Nucleus count ≠ 1: excluded", "#009E73"),
                "require_three_markers": ("Require all 3 markers", "#0072B2")}
    for ax, scope, title in ((c, "whole_slice", "c  Halo retained; imported-window background"),
                            (d, "roi_induced", "d  ROI-only graph and ROI background")):
        entries = runs[scope]["results"]
        baseline = [x for x in entries if x["variant"] == "non_single_nucleus_uncertain"]
        ax.plot([x["radius_um"] for x in baseline], [100*x["before"]["excess_over_abundance"] for x in baseline],
                marker="o", color="black", linewidth=2, label="Working-label baseline")
        for variant, (label, color) in variants.items():
            points = [x for x in entries if x["variant"] == variant]
            ax.plot([x["radius_um"] for x in points], [100*x["after"]["excess_over_abundance"] for x in points],
                    marker="o", color=color, label=label)
        ax.axhline(0, color="#7A8188", linewidth=.8)
        ax.axhline(plan["min_effect"]*100, color="#B64342", linestyle="--", linewidth=1.1, label="Declared threshold: +10 pp")
        ax.set(xlabel="Declared radius (µm)", ylabel="Neighbor fraction − background (pp)",
               xticks=plan["radii_um"], ylim=(-6, 12))
        ax.set_title(title, loc="left", fontsize=12)
        ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -.2))
    fig.suptitle("Real Xenium reactive lymph node | 7,725 cells, 4,624 genes\nFrozen ROI: 3,383 cells + 100 µm halo | one donor, descriptive hypotheses only", fontsize=14)
    if helper:
        helper.finalize_figure(fig, args.study / "real-data-sensitivity", formats=["png", "pdf"], dpi=300)
    else:
        for extension in ("png", "pdf"):
            fig.savefig(args.study / f"real-data-sensitivity.{extension}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    rows = []
    for scope, run in runs.items():
        for entry in run["results"]:
            for state in ("before", "after"):
                metric = entry[state]
                rows.append({"scope": scope, "variant": entry["variant"], "radius_um": entry["radius_um"], "state": state,
                             **metric["graph"], "observed_fraction": metric["observed_fraction"], "background_fraction": metric["null_fraction"],
                             "excess": metric["excess_over_abundance"], "classification": metric["classification"]})
    with (args.study / "metric-denominators.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"png": str(args.study / "real-data-sensitivity.png"), "rows": len(rows)}))


if __name__ == "__main__":
    main()
