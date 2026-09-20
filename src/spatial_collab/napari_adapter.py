"""Optional napari view/controller over the same revisioned Project as MCP.

The controller and reader have no Qt or napari dependency. Only the dock factory
imports the viewer libraries. Arrays follow napari (y, x); the store uses (x, y).
Viewer edits never write coordinates, labels or segmentations into source data.
"""
from __future__ import annotations

from copy import deepcopy
import json
from numbers import Integral
from pathlib import Path
from uuid import uuid4

import numpy as np

from .store import Project, SpatialError

MAX_LAYER_CELLS = 250_000
FEATURES = ("cell_id", "label", "included", "region")
KEY = "spatial_collab"


def _project_root(path):
    if isinstance(path, (list, tuple)):
        if len(path) != 1:
            raise SpatialError("Open one Spatial Collab project at a time.")
        path = path[0]
    root = Path(path).resolve()
    return root.parent if root.name == "project.sqlite3" else root


def _world_affine(layer):
    """Read the public 2D napari data-to-world mapping, including extra affine.

    The mapping is sampled at basis points and verified on additional probes.
    napari layer transforms are affine. Unknown/nonlinear or 3D mappings fail.
    """
    if layer.ndim != 2:
        raise SpatialError("Only two-dimensional calibrated layers are supported.")
    probes = np.array([[0., 0.], [1., 0.], [0., 1.], [2., 3.], [-2., 5.]])
    mapped = np.asarray([layer.data_to_world(point) for point in probes], dtype=float)
    if mapped.shape != probes.shape or not np.isfinite(mapped).all():
        raise SpatialError("Layer transform is not a finite two-dimensional affine.")
    linear = np.column_stack((mapped[1] - mapped[0], mapped[2] - mapped[0]))
    matrix = np.eye(3)
    matrix[:2, :2], matrix[:2, 2] = linear, mapped[0]
    if not np.isfinite(np.linalg.det(linear)) or abs(np.linalg.det(linear)) < 1e-12:
        raise SpatialError("Layer transform must be invertible.")
    if not np.allclose(probes @ linear.T + mapped[0], mapped, rtol=1e-10, atol=1e-10):
        raise SpatialError("Unsupported nonlinear layer transform.")
    return matrix


def capture_points(layer):
    """Copy GUI-owned arrays on the GUI thread, then validate off-thread."""
    try:
        features = {name: list(layer.features[name]) for name in FEATURES}
    except (KeyError, TypeError) as exc:
        raise SpatialError("Missing point identity features; refresh this project.") from exc
    return {"data": np.array(layer.data, dtype=float, copy=True),
            "features": features, "metadata": deepcopy(layer.metadata),
            "selected_indices": list(layer.selected_data), "transform": _world_affine(layer)}


def capture_polygon(layer):
    """Freeze exactly one selected rectangle/polygon in world (y, x)."""
    selected = list(layer.selected_data)
    if len(selected) != 1:
        raise SpatialError("Select exactly one polygon or rectangle in the ROI layer.")
    index = selected[0]
    if layer.shape_type[index] not in {"polygon", "rectangle"}:
        raise SpatialError("Only polygons and rectangles define an exact centroid ROI.")
    vertices = np.array(layer.data[index], dtype=float, copy=True)
    if vertices.ndim != 2 or vertices.shape[1] != 2 or not 3 <= len(vertices) <= 1001:
        raise SpatialError("ROI needs 3 to 1000 vertices in two dimensions.")
    transform = _world_affine(layer)
    world = vertices @ transform[:2, :2].T + transform[:2, 2]
    if not np.isfinite(world).all():
        raise SpatialError("ROI coordinates must be finite.")
    return {"metadata": deepcopy(layer.metadata), "polygon": world[:, ::-1].tolist()}


def capture_asset_layer(layer):
    """Freeze only a bounded materialized 2D image/labels layer on the GUI thread."""
    from .assets import validate_array_budget
    data, _, kind = layer.as_layer_data_tuple()
    if kind not in {"image", "labels"} or not isinstance(data, np.ndarray):
        raise SpatialError("Select a materialized NumPy image/labels plane; lazy/multiscale layers are not densified.")
    validate_array_budget(data.shape, data.dtype)
    return {"data": np.array(data, copy=True), "metadata": deepcopy(layer.metadata),
            "transform": _world_affine(layer), "kind": "segmentation" if kind == "labels" else "image"}


def add_asset_to_viewer(viewer, layer_data):
    """Use a detached rendering buffer; native Labels needs a writable buffer.

    Original assets and the core-loaded array remain read-only. GUI paint is
    disabled, while console/data replacement is detected at every inspection.
    """
    data, kwargs, kind = layer_data
    layer = getattr(viewer, "add_" + kind)(np.array(data, copy=True), **kwargs)
    if kind == "labels":
        layer.editable = False
    return layer


class NapariController:
    """Headless controller, serialized by the dock's background job queue.

    A snapshot pins source identity, revision, frame, point order and features.
    Refresh is explicit: a selection from another client is never silently used
    in an already prepared annotation proposal.
    """

    def __init__(self, project_path):
        self.project = Project(_project_root(project_path))
        self.preview = None
        self.refresh()

    def refresh(self):
        summary = self.project.summary()
        if summary["cell_count"] > MAX_LAYER_CELLS:
            raise SpatialError(f"napari point budget exceeded ({MAX_LAYER_CELLS}); no cells were sampled.")
        cells = self.project.cells(summary["head_revision"])
        self.summary = summary
        self.revision_id = summary["head_revision"]
        self.coordinates = np.asarray([[cell["y"], cell["x"]] for cell in cells], dtype=float)
        self.features = {field: [cell[field] for cell in cells] for field in FEATURES}
        meta = summary["metadata"]
        self.binding = {"schema": "spatial-collab.napari-frame.v1", "snapshot_id": uuid4().hex,
                        "project_path": str(self.project.root), "project_id": summary["project_id"],
                        "source_sha256": summary["source_sha256"], "revision_id": self.revision_id,
                        "slice_id": meta["slice_id"], "coordinate_system": meta["coordinate_system"],
                        "units": meta["units"], "axis_order": ["y", "x"], "complete": True}
        selection = summary["selection"]
        self.selection_id = selection["selection_id"] if selection and not selection["stale"] else None
        self.preview = None
        return summary

    def layer_data(self):
        return [(self.coordinates.copy(),
                 {"name": "Spatial Collab · observations", "features": deepcopy(self.features),
                  "metadata": {KEY: deepcopy(self.binding)}, "size": 6,
                  "face_color": "label", "face_color_cycle":
                  ["#4db6ac", "#f4b942", "#ef7d8b", "#8d9ff0", "#b0ca70", "#cb9bea"]}, "points"),
                ([], {"name": "Spatial Collab · ROI", "ndim": 2,
                      "metadata": {KEY: deepcopy(self.binding)}, "edge_color": "yellow",
                      "face_color": "transparent", "edge_width": 2}, "shapes")]

    def _current(self):
        summary = self.project.summary()
        if (summary["head_revision"] != self.revision_id
                or summary["source_sha256"] != self.binding["source_sha256"]):
            raise SpatialError("Stale napari view: another client changed the project. Refresh and review again.")
        return summary

    def _bound(self, metadata):
        if metadata.get(KEY) != self.binding:
            raise SpatialError("Layer identity, revision or coordinate frame changed. Refresh the project.")
        self._current()

    def validate_points(self, capture):
        self._bound(capture["metadata"])
        data = capture["data"]
        if data.shape != self.coordinates.shape or not np.array_equal(data, self.coordinates):
            raise SpatialError("Points were moved, added, deleted or reordered. Refresh; source coordinates are immutable.")
        if capture["features"] != self.features:
            raise SpatialError("Point identities or annotations were edited in the viewer. Refresh and use a revision preview.")
        if not np.allclose(capture["transform"], np.eye(3), rtol=0, atol=1e-12):
            raise SpatialError("The calibrated points layer was transformed. Refresh before sharing a selection.")

    def publish_points(self, capture, name="napari point selection"):
        self.validate_points(capture)
        indices = capture["selected_indices"]
        if any(isinstance(i, bool) or not isinstance(i, Integral) or i < 0
               or i >= len(self.coordinates) for i in indices) or len(set(indices)) != len(indices):
            raise SpatialError("Invalid point indices.")
        ids = [self.features["cell_id"][index] for index in sorted(indices)]
        selection = self.project.set_selection(self.revision_id, cell_ids=ids, name=name)
        self.selection_id, self.preview = selection["selection_id"], None
        return selection

    def publish_polygon(self, shape_capture, points_capture, name="napari polygon ROI"):
        # A moved points layer would visually imply a different membership.
        self.validate_points(points_capture)
        self._bound(shape_capture["metadata"])
        selection = self.project.set_selection(self.revision_id, polygon=shape_capture["polygon"], name=name)
        self.selection_id, self.preview = selection["selection_id"], None
        return selection

    def inspect(self, genes=None):
        self._current()
        if not self.selection_id:
            raise SpatialError("Share a selection or refresh to adopt the other client's selection first.")
        result = self.project.inspect_selection(genes)
        if result["selection"]["selection_id"] != self.selection_id:
            raise SpatialError("Another client changed the selection. Refresh to inspect its exact cells.")
        return result

    def propose(self, changes, rationale):
        self.preview = None
        self._current()
        if not self.selection_id:
            raise SpatialError("Share a selection first.")
        self.preview = self.project.propose_revision(self.revision_id, self.selection_id,
                                                      changes, rationale, actor="napari researcher")
        return deepcopy(self.preview)

    def load_proposal(self, proposal_id=None):
        """Load another host's immutable preview without changing shared selection.

        The caller must highlight the returned selection, not whichever selection
        happens to be active in another host. Loading does not grant confirmation.
        """
        self.preview = None
        self._current()
        if proposal_id is None:
            proposal_id = self.project.context()["latest_proposal_id"]
        if not proposal_id:
            raise SpatialError("There is no shared proposal to review.")
        proposal = self.project.get_proposal(proposal_id)
        if proposal["stale"] or proposal["base_revision"] != self.revision_id:
            raise SpatialError("This shared proposal is stale. Refresh and ask for a new preview.")
        if proposal["applied_revision"]:
            raise SpatialError("This shared proposal was already applied.")
        selection = self.project.get_selection(proposal["selection_id"])
        if selection is None or selection["stale"]:
            raise SpatialError("The shared proposal's exact selection is unavailable or stale.")
        if any(selection[key] != self.binding[key] for key in
               ("source_sha256", "slice_id", "coordinate_system", "units", "revision_id")):
            raise SpatialError("The shared proposal's selection does not match this project frame.")
        self._current()
        self.preview, self.selection_id = proposal, selection["selection_id"]
        return {"proposal": deepcopy(proposal), "selection": selection,
                "shared_selection_changed": False,
                "review_notice": "Highlight shows this proposal's frozen selection. "
                                 "Review the concrete delta and explicitly confirm before committing."}

    def commit(self, reviewer, confirmation):
        if self.preview is None:
            raise SpatialError("Create and inspect an annotation preview first.")
        self._current()
        result = self.project.apply_revision(self.preview["proposal_id"], self.revision_id, reviewer, confirmation)
        self.refresh()
        return result

    def restore(self, revision_id, reviewer, confirmation):
        self._current()
        result = self.project.revert_revision(revision_id, self.revision_id, reviewer, confirmation)
        self.refresh()
        return result

    def compare(self, base_revision, target_revision, selection_id=None, **parameters):
        if self.binding["units"] != "micrometer":
            raise SpatialError("Physical neighborhood analysis requires explicitly calibrated micrometer coordinates.")
        from .analysis import compare
        return compare(self.project, base_revision, target_revision, selection_id, **parameters)

    def list_assets(self):
        from .assets import list_assets
        return list_assets(self.project)

    def asset_layer(self, asset_id):
        from .assets import load_asset_layer
        self._current()
        return load_asset_layer(self.project, asset_id)

    def inspect_asset(self, asset_id, asset_capture, points_capture):
        from .assets import inspect_asset, validate_layer_capture
        self.validate_points(points_capture)
        validate_layer_capture(self.project, asset_id, asset_capture)
        indices = points_capture["selected_indices"]
        if len(indices) > 100:
            raise SpatialError("Select at most 100 exact observations for pixel/mask inspection.")
        if any(type(i) is bool or not isinstance(i, Integral) or i < 0 or i >= len(self.coordinates) for i in indices):
            raise SpatialError("Invalid point indices.")
        ids = [self.features["cell_id"][index] for index in sorted(indices)]
        return inspect_asset(self.project, asset_id, ids)

    def register_layer_snapshot(self, capture, *, source_sha256, slice_id, coordinate_system,
                                units, registration_note, label_to_cell_id=None, origin_sources=None, name=None,
                                sample_id=None):
        from .assets import register_array_snapshot
        self._current()
        swap = np.array([[0., 1, 0], [1, 0, 0], [0, 0, 1]])
        return register_array_snapshot(self.project, capture["data"], kind=capture["kind"],
                                       source_sha256=source_sha256, slice_id=slice_id,
                                       coordinate_system=coordinate_system, units=units,
                                       pixel_to_world=swap @ capture["transform"] @ swap,
                                       registration_note=registration_note, label_to_cell_id=label_to_cell_id,
                                       origin_sources=origin_sources or [], name=name, sample_id=sample_id)


def napari_get_reader(path):
    """Cheap probing only; never claim arbitrary SQLite files or raw formats."""
    try:
        root = _project_root(path)
        if (root / "project.sqlite3").is_file():
            return read_project_layers
    except (TypeError, ValueError, OSError):
        pass
    return None


def read_project_layers(path):
    return NapariController(path).layer_data()


class SpatialCollabWidget:
    """Lazy npe2 constructor with a host-recognized viewer injection signature.

    napari inspects the class's __init__ before calling it. Returning the actual
    QWidget from __new__ avoids importing Qt when merely discovering plugins.
    The explicit host viewer also keeps two simultaneously open viewers apart.
    """

    def __new__(cls, napari_viewer):
        return make_dock_widget(napari_viewer)

    def __init__(self, napari_viewer):
        pass


def make_dock_widget(napari_viewer=None):
    """npe2 widget factory; imports Qt only when a viewer requests the dock."""
    import napari
    from napari.qt.threading import thread_worker
    from qtpy.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout,
                               QLabel, QLineEdit, QPushButton, QTabWidget, QTextEdit, QVBoxLayout, QWidget)

    viewer = napari_viewer if napari_viewer is not None else napari.current_viewer()
    if viewer is None:
        raise SpatialError("Open a napari viewer before adding Spatial Collab.")

    class CollaborationDock(QWidget):
        def __init__(self):
            super().__init__()
            self.controller = self.points = self.shapes = self.worker = None
            self.restore_preview = None
            self.asset_layers = {}
            self.buttons = []
            layout = QVBoxLayout(self)
            self.status = QLabel("Open the same project folder used by your agent.")
            self.status.setWordWrap(True)
            layout.addWidget(self.status)
            form = QFormLayout()
            layout.addLayout(form)
            self.path = QLineEdit()
            form.addRow("Project folder", self.path)
            self.button(layout, "Open / refresh shared project", self.open_project)
            tabs = QTabWidget()
            layout.addWidget(tabs)
            explore, review, comparison = QWidget(), QWidget(), QWidget()
            for page, title in ((explore, "1 · Explore"), (review, "2 · Revise"),
                                (comparison, "3 · Compare")):
                tabs.addTab(page, title)
            explore_layout, review_layout, comparison_layout = (
                QVBoxLayout(explore), QVBoxLayout(review), QVBoxLayout(comparison))
            self.button(explore_layout, "Share selected point IDs", self.share_points)
            self.button(explore_layout, "Share selected polygon / rectangle", self.share_polygon)
            self.genes = QLineEdit()
            self.genes.setPlaceholderText("CD3D, LYZ (blank: first 20 panel genes)")
            explore_form = QFormLayout()
            explore_layout.addLayout(explore_form)
            explore_form.addRow("Markers", self.genes)
            self.button(explore_layout, "Inspect shared selection", self.inspect)
            self.asset_choice = QComboBox()
            explore_form.addRow("Registered image / mask", self.asset_choice)
            self.button(explore_layout, "Refresh registered asset list", self.refresh_assets)
            self.button(explore_layout, "Load registered image / mask", self.load_asset)
            self.button(explore_layout, "Inspect selected IDs in visible asset", self.inspect_visible_asset)
            explore_layout.addStretch()
            self.label = QLineEdit()
            self.included = QComboBox()
            self.included.addItems(["Keep inclusion", "Include", "Exclude"])
            self.rationale = QLineEdit()
            review_form = QFormLayout()
            review_layout.addLayout(review_form)
            self.shared_proposal = QLineEdit()
            self.shared_proposal.setPlaceholderText("Proposal ID from agent; blank = latest proposal")
            review_form.addRow("Shared proposal", self.shared_proposal)
            self.button(review_layout, "Load shared proposal for review", self.load_proposal)
            review_form.addRow("Proposed label", self.label)
            review_form.addRow("Inclusion", self.included)
            review_form.addRow("Reason", self.rationale)
            self.button(review_layout, "Preview annotation change", self.propose)
            self.reviewer = QLineEdit()
            review_form.addRow("Researcher name", self.reviewer)
            self.confirm = QCheckBox("I reviewed this preview / restore target and confirm")
            review_layout.addWidget(self.confirm)
            self.button(review_layout, "Commit reviewed preview", self.commit)
            self.restore_target = QComboBox()
            review_form.addRow("Restore target", self.restore_target)
            self.button(review_layout, "Preview restore target", self.preview_restore)
            self.button(review_layout, "Restore before as a new revision", self.restore)
            review_layout.addStretch()
            self.base = QComboBox()
            self.target = QComboBox()
            comparison_form = QFormLayout()
            comparison_layout.addLayout(comparison_form)
            comparison_form.addRow("Before / restore target", self.base)
            comparison_form.addRow("After", self.target)
            self.radius = QDoubleSpinBox()
            self.radius.setRange(0.001, 1e6)
            self.radius.setDecimals(3)
            self.radius.setValue(35)
            self.source_label = QLineEdit("T cell")
            self.target_label = QLineEdit("Myeloid")
            self.scope = QComboBox()
            self.scope.addItems(["roi_induced", "whole_slice"])
            self.threshold = QDoubleSpinBox()
            self.threshold.setRange(0.001, 1)
            self.threshold.setDecimals(3)
            self.threshold.setValue(0.1)
            self.comparison_roi = QLineEdit()
            self.comparison_roi.setPlaceholderText("Selection ID; blank = all imported observations")
            comparison_form.addRow("Frozen selection ID", self.comparison_roi)
            comparison_form.addRow("Radius (µm)", self.radius)
            comparison_form.addRow("Source label", self.source_label)
            comparison_form.addRow("Target label", self.target_label)
            comparison_form.addRow("Graph scope", self.scope)
            comparison_form.addRow("Descriptive effect threshold", self.threshold)
            self.compare_button = self.button(comparison_layout, "Compare before / after", self.compare)
            self.output = QTextEdit()
            self.output.setReadOnly(True)
            layout.addWidget(self.output)
            notice = QLabel("Image/mask registrations require explicit source, frame and cell-ID mapping. "
                            "Missing evidence stays unknown; overlays do not prove merged cells or true coexpression. "
                            "Named confirmation is attribution, not authenticated scientific approval.")
            notice.setWordWrap(True)
            layout.addWidget(notice)

        def button(self, layout, text, action):
            button = QPushButton(text)
            button.clicked.connect(lambda _=False: self.safe(action))
            layout.addWidget(button)
            self.buttons.append(button)
            return button

        def safe(self, action):
            try:
                action()
            except Exception as exc:
                self.failed(exc)

        def require(self):
            if self.controller is None or self.points not in viewer.layers or self.shapes not in viewer.layers:
                raise SpatialError("Open or refresh a project first.")

        def display(self, result):
            # Full preview available by scrolling; never hide proposal deltas.
            self.output.setPlainText(json.dumps(result, ensure_ascii=False, indent=2, default=str))
            if self.controller is not None:
                self.status.setText(f"Revision {self.controller.revision_id} · refresh after agent changes")

        def failed(self, exc):
            self.status.setText(str(exc))

        def run_job(self, operation, done):
            if self.worker is not None:
                raise SpatialError("Wait for the current operation to finish.")
            for button in self.buttons:
                button.setEnabled(False)
            self.status.setText("Working… the viewer remains interactive.")
            worker = thread_worker(operation)()
            self.worker = worker
            worker.returned.connect(lambda result: self.safe(lambda: done(result)))
            worker.errored.connect(self.failed)
            worker.finished.connect(self.job_finished)
            worker.start()

        def job_finished(self):
            self.worker = None
            for button in self.buttons:
                button.setEnabled(True)
            if self.controller is not None:
                self.compare_button.setEnabled(self.controller.binding["units"] == "micrometer")

        def open_project(self):
            path = self.path.text().strip()
            if not path:
                raise SpatialError("Enter the shared project folder.")
            self.run_job(lambda: NapariController(path), self.opened)

        def opened(self, controller):
            # Reuse this dock's layers only. Unrelated viewer layers stay intact.
            if self.controller is not None and self.controller.binding["source_sha256"] != controller.binding["source_sha256"]:
                for layer in self.asset_layers.values():
                    if layer in viewer.layers:
                        viewer.layers.remove(layer)
                self.asset_layers.clear()
            self.controller = controller
            self.render()
            self.display(controller.summary)

        def render(self):
            layers = self.controller.layer_data()
            for attr, (data, kwargs, kind) in zip(("points", "shapes"), layers):
                layer = getattr(self, attr)
                if layer is None or layer not in viewer.layers:
                    layer = getattr(viewer, "add_" + kind)(data, **kwargs)
                    setattr(self, attr, layer)
                else:
                    layer.data = data
                    layer.metadata = kwargs["metadata"]
                    layer.scale, layer.translate = (1, 1), (0, 0)
                    layer.rotate, layer.shear, layer.affine = 0, np.eye(2), np.eye(3)
                    if kind == "points":
                        layer.features = kwargs["features"]
                        layer.face_color = "label"
                if kind == "points":
                    layer.mode = "select"
            summary = self.controller.summary
            selection = summary["selection"]
            if selection and not selection["stale"]:
                selected = set(selection["cell_ids"])
                self.points.selected_data = {i for i, cid in enumerate(self.controller.features["cell_id"])
                                             if cid in selected}
            else:
                self.points.selected_data = set()
            before, after = self.base.currentText(), self.target.currentText()
            revisions = [item["revision_id"] for item in summary["revisions"]]
            for combo, previous, default in ((self.base, before, revisions[0]),
                                             (self.target, after, revisions[-1])):
                combo.clear()
                combo.addItems(revisions)
                combo.setCurrentText(previous if previous in revisions else default)
            self.target.setCurrentText(summary["head_revision"])
            self.restore_target.clear()
            self.restore_target.addItems(revisions)
            self.restore_preview = None
            self.confirm.setChecked(False)
            self.compare_button.setEnabled(self.controller.binding["units"] == "micrometer")
            self.status.setText(f"Revision {self.controller.revision_id} · {summary['cell_count']} observations · "
                                f"coordinates: {self.controller.binding['units']} · "
                                "refresh after agent changes")

        def refresh_assets(self):
            self.require()

            def listed(result):
                previous = self.asset_choice.currentData()
                self.asset_choice.clear()
                for asset in result["assets"]:
                    self.asset_choice.addItem(f"{asset['kind']} · {asset['name']} · {asset['mapping_status']}", asset["asset_id"])
                index = self.asset_choice.findData(previous)
                if index >= 0:
                    self.asset_choice.setCurrentIndex(index)
                self.display(result)

            self.run_job(self.controller.list_assets, listed)

        def load_asset(self):
            self.require()
            asset_id = self.asset_choice.currentData()
            if not asset_id:
                raise SpatialError("Refresh the asset list and select a registered image or mask first.")

            def loaded(layer_data):
                data, kwargs, kind = layer_data
                previous = self.asset_layers.get(asset_id)
                if previous is not None and previous in viewer.layers:
                    viewer.layers.remove(previous)
                layer = add_asset_to_viewer(viewer, layer_data)
                self.asset_layers[asset_id] = layer
                self.display(kwargs["metadata"])

            self.run_job(lambda: self.controller.asset_layer(asset_id), loaded)

        def inspect_visible_asset(self):
            self.require()
            asset_id = self.asset_choice.currentData()
            layer = self.asset_layers.get(asset_id)
            if layer is None or layer not in viewer.layers:
                raise SpatialError("Load this registered asset before inspecting its visible cell correspondence.")
            points, asset = capture_points(self.points), capture_asset_layer(layer)
            self.run_job(lambda: self.controller.inspect_asset(asset_id, asset, points), self.display)

        def selection_done(self, result):
            self.display(result)
            chosen = set(result["cell_ids"])
            self.points.selected_data = {i for i, cid in enumerate(self.controller.features["cell_id"])
                                         if cid in chosen}
            self.comparison_roi.setText(result["selection_id"])
            self.confirm.setChecked(False)
            self.status.setText(f"Shared {result['cell_count']} exact observations · {result['selection_id']}")

        def share_points(self):
            self.require()
            capture = capture_points(self.points)
            self.run_job(lambda: self.controller.publish_points(capture), self.selection_done)

        def share_polygon(self):
            self.require()
            points, shape = capture_points(self.points), capture_polygon(self.shapes)
            self.run_job(lambda: self.controller.publish_polygon(shape, points), self.selection_done)

        def inspect(self):
            self.require()
            genes = [value.strip() for value in self.genes.text().split(",") if value.strip()] or None
            self.run_job(lambda: self.controller.inspect(genes), self.display)

        def propose(self):
            self.require()
            changes = {}
            if self.label.text().strip():
                changes["label"] = self.label.text().strip()
            if self.included.currentIndex():
                changes["included"] = self.included.currentIndex() == 1
            rationale = self.rationale.text().strip()
            self.confirm.setChecked(False)
            self.restore_preview = None
            self.run_job(lambda: self.controller.propose(changes, rationale), self.display)

        def load_proposal(self):
            self.require()
            proposal_id = self.shared_proposal.text().strip() or None
            self.confirm.setChecked(False)
            self.controller.preview = None
            self.restore_preview = None
            points = capture_points(self.points)

            def loaded(result):
                self.display(result)
                chosen = set(result["selection"]["cell_ids"])
                self.points.selected_data = {i for i, cid in enumerate(self.controller.features["cell_id"])
                                             if cid in chosen}
                self.shared_proposal.setText(result["proposal"]["proposal_id"])
                self.comparison_roi.setText(result["selection"]["selection_id"])
                self.status.setText(f"Loaded exact proposal selection: {len(chosen)} observations. "
                                    "Review the delta; confirmation is required.")

            def load():
                self.controller.validate_points(points)
                return self.controller.load_proposal(proposal_id)

            self.run_job(load, loaded)

        def committed(self, result):
            self.render()
            self.display(result)

        def commit(self):
            self.require()
            reviewer, confirm = self.reviewer.text().strip(), self.confirm.isChecked()
            self.confirm.setChecked(False)
            points = capture_points(self.points)

            def commit():
                self.controller.validate_points(points)
                return self.controller.commit(reviewer, confirm)

            self.run_job(commit, self.committed)

        def preview_restore(self):
            self.require()
            target = self.restore_target.currentText()
            revision = self.controller.revision_id
            self.confirm.setChecked(False)
            self.controller.preview = None

            def ready(result):
                self.restore_preview = (target, revision)
                self.display(result)

            self.run_job(lambda: self.controller.project.get_revision(target), ready)

        def restore(self):
            self.require()
            revision = self.restore_target.currentText()
            if self.restore_preview != (revision, self.controller.revision_id):
                raise SpatialError("Preview this restore target before confirming it.")
            reviewer, confirm = self.reviewer.text().strip(), self.confirm.isChecked()
            self.confirm.setChecked(False)
            self.run_job(lambda: self.controller.restore(revision, reviewer, confirm), self.committed)

        def compare(self):
            self.require()
            base, target = self.base.currentText(), self.target.currentText()
            selection = self.comparison_roi.text().strip() or None
            parameters = dict(radius_um=self.radius.value(), graph_scope=self.scope.currentText(),
                              source_label=self.source_label.text().strip(),
                              target_label=self.target_label.text().strip(), min_effect=self.threshold.value())
            self.run_job(lambda: self.controller.compare(base, target, selection, **parameters), self.display)

    return CollaborationDock()
