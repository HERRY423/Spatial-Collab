# v0.2 implementation contract

Continue v0.1, preserve existing project compatibility. Version 0.2.0-alpha.1 / 0.2.0a1. No irreversible migrations or overwriting source datasets.

## Data semantics

Canonical Project API remains `cells()` / cell_id for backwards compatibility, but observations may be cells, spots or bins. Metadata extends with `platform`, `observation_unit` (`cell`, `spot`, `bin`), `label_semantics` (`cell_type`, `spot_annotation`, `bin_annotation`), and explicit coordinate frame/units. Legacy projects default cell+unknown platform in presentation, never rewrite original source/hash. Current computation remains within a single slice. Spot labels are regional/domain/working annotations, not pure cell identities or deconvolution.

Root owns store core enhancements, SpatialData adapter, server/UI integration, analysis semantic extension, main CLI, pyproject/manifests/docs integration. Import agent owns import_registry.py, format_readers.py, tests/test_formats.py and limited compatibility changes importers.py as notified. Napari agent owns napari_adapter.py, napari.yaml and tests. Benchmark agent owns benchmarks.py, benchmark task specs and associated tests/docs after source review.

## Import registry interface

Registry metadata plain JSON-friendly. `list_formats() -> list[dict]` capabilities with id, description, supported platform, observation units, required inputs/dependency/ceilings. `probe_source(path) -> dict` side-effect-free structured candidate readers and missing metadata; may inspect headers/filenames only; no full load, no confidence fake numbers, no automatic reader choice on ambiguity. `import_source(format_id, path, destination, **options) -> Project` calls explicit reader, lazy imports. Built-ins: anndata_h5ad, xenium (existing), tenx_mtx (generic matrix+positions+annotations), visium, cosmx_csv, merscope_csv, spatialdata_zarr (root implements). Expose register_reader declarative adapter so root can register SpatialData without altering importer-owned dispatch. Registration explicit trusted Python code; never arbitrary remote execution or inferred plugins from user files.

Formats supported must declare exact flavor and required user mappings. `platform` and `observation_unit` reflect actual reader. Xenium matrix fallback MTX requested. Generic 10x MTX positions CSV cell_id,x,y; labels CSV cell_id,label; mandatory units='micrometer', coordinate_system,slice_id. Visium supports legacy standard spot output, explicit microns_per_pixel (no inferred spot diameter), positions convention fullres x=pxl_col_in_fullres,y=pxl_row_in_fullres. Annotation sidecar exact alignment; explicit in_tissue subset handling, no undeclared filtering. CSV imaging profiles use global coordinates only or explicit mappings; FOV-local coordinate columns must reject pending transform, never merge overlapping FOVs. Preserve existing source hashes and pre/post source integrity checks; matrix shape, duplicate IDs, raw counts, sparse accumulation, resource limits.

New source import may allow explicitly unannotated labels (user option allow_unannotated) so discovery and marker exploration can start before full labels, but no fake annotation/comparison evidence. Existing import_h5ad/import_xenium behavior remains strict by default.

## Collaboration service additions (root)

Add lightweight revision/selection event cursor polling, cancellation/resumable task metadata if feasible, paged revision-pinned observation discovery/filtering, row/marker inspection, and capability manifest. Heavy import stays CLI, paths pinned for MCP. Keep old 11 tools working; add few semantic tools. Compare metrics must name observation_unit in evidence ceiling/limitations; never report spot neighborhoods as cell contacts. Further benchmark-oriented ROI composition/marker contrasts and radius sensitivity should be genuine computations with explicitly descriptive scope.

## Napari

Lazy optional npe2 reader/widget, model/controller separated from viewer. Same persistent selection and CAS proposal/commit. Points y,x conversion and raw cell_id feature mapping validated against project revision/frame. Event-bound layer metadata and refresh make stale viewer detectable. Shapes are ROI geometry, not segmentation. No source writes or viewer-selected-index guesses after transform/data mutation.
