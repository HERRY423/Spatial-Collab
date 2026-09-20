"""Real server outputs through production JS renderers; DOM fixture, not browser acceptance."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from spatial_collab.analysis import compare
from spatial_collab.exploration import compare_regions
from spatial_collab.quality import inspect_quality
from spatial_collab.sensitivity import run_sensitivity
from spatial_collab.server import ToolService
from spatial_collab.store import Project


@pytest.mark.skipif(shutil.which("node") is None, reason="Node runtime unavailable for optional DOM unit smoke")
def test_workbench_qc_hypothesis_and_crop_renderers(tmp_path):
    project = Project.create(tmp_path / "project", [
        {"cell_id": "s", "x": 0, "y": 0, "label": "T", "counts": {"G": 1},
         "attributes": {"nucleus_area": None, "nucleus_count": 0, "segmentation_method": "source method"}},
        {"cell_id": "t", "x": 1, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "b", "x": 10, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "o", "x": 20, "y": 0, "label": "O", "counts": {}}],
        {"name": "Synthetic UI fixture", "source_kind": "synthetic", "slice_id": "declared-slice-id", "units": "micrometer",
         "source_files": [{"path": "fixture/cells<untrusted>.csv", "sha256": "a" * 64}],
         "working_annotation": {"kind": "automated_exploratory_marker_rule", "rules": {"markers": ["G"]},
                                "researcher_approval": "NOT_ESTABLISHED", "validated_cell_identity": False},
         "coordinate_system": "s-xy", "panel_genes": ["G"], "import_scope": {
             "kind": "spatial_window", "bounds": [-10, -10, 30, 30], "units": "micrometer",
             "sampling": "none", "selected_observation_count": 4,
             "selection_rule": "all cell centroids with xmin <= x <= xmax and ymin <= y <= ymax"}})
    revision = project.summary()["head_revision"]
    selected = project.set_selection(revision, cell_ids=["s"])
    variant = {"name": "hypothetical source exclusion", "cell_ids": ["s"],
               "changes": {"included": False}, "rationale": "Synthetic what-if; no approval"}
    hypothesis = run_sensitivity(project, revision, selected["selection_id"], [variant], [2.0], "T", "B")
    neighborhood = compare(project, revision, revision, selected["selection_id"], 2.0, "whole_slice", "T", "B")
    region = compare_regions(project, revision, revision, selected["selection_id"], genes=["G"])
    payload = {"project": ToolService(project.root).open_project(), "quality": inspect_quality(project, revision),
               "sensitivity": hypothesis, "neighborhood": neighborhood, "region": region}
    fixture = tmp_path / "payload.json"
    fixture.write_text(json.dumps(payload), encoding="utf-8")
    root = Path(__file__).parents[1]
    result = subprocess.run([shutil.which("node"), str(Path(__file__).with_name("workbench_render_smoke.cjs")),
                             str(root / "src/spatial_collab/static/workbench.js"),
                             str(root / "src/spatial_collab/static/workbench.html"), str(fixture)],
                            capture_output=True, text=True, encoding="utf-8", timeout=30)
    assert result.returncode == 0, result.stderr
    assert all(json.loads(result.stdout).values())
