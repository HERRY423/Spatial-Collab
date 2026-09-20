# Editable, frozen hypothesis plans

The v2 plan API keeps editable computational plans separate from source data and annotation revisions. A researcher and agent can inspect or revise a draft, freeze its exact inputs, and compare every declared alternative without creating an approved annotation or a fictitious scientific reviewer. The existing `run_sensitivity` v1 API and its exported results retain their original behavior.

## Lifecycle and persistence

`src/spatial_collab/hypotheses.py` exposes:

| Function | Purpose |
| --- | --- |
| `create_hypothesis_plan(project, spec)` | Create version 1, status `draft`. |
| `revise_hypothesis_plan(project, plan_id, expected_version, spec)` | Append a new draft from the current version, including after a previous freeze. |
| `freeze_hypothesis_plan(project, plan_id, expected_version)` | Append a frozen version with pinned source, revision, selections and resolved observation IDs. |
| `get_hypothesis_plan(project, plan_id, version=None)` | Read an exact version, or the current version. |
| `list_hypothesis_plans(project)` | List plan heads and retained version histories. |
| `hypothesis_context(project)` | Read a small state token and latest plan status without loading count payloads. |
| `run_hypothesis_plan(project, plan_id, version)` | Execute one explicitly frozen version across the complete declared grid. |
| `run_hypothesis_plan(project, frozen_plan=record)` | Recompute from embedded frozen parameters, without requiring a sidecar ledger. |

Each record has `schema=spatial-collab.hypothesis-plan.v2`, `plan_id`, integer `version`, `status`, `parent_version`, `parent_sha256`, `created_at`, `spec`, `resolved` and `plan_sha256`. The SHA-256 covers the canonical record without its own digest. A frozen version is never edited. Revising it creates a newer draft while the old version stays addressable and executable.

The project-local `hypotheses.sqlite3` sidecar stores immutable version rows. A transaction compares `expected_version` with the current head before appending; a stale editor must refresh. It also binds the ledger to project ID and source digest. Freeze and plan editing do not alter `Project.source`, scientific revision history or head. Runs use ordinary saved-run persistence; the freeze record is embedded in `parameters.frozen_plan` for bundle replay. Plan freezing records an explicit computational choice, not expert adjudication, human authentication or independently verified preregistration.

## Minimal specification

Replace the example revision, selection and cell IDs with actual values from the shared workspace. All shown top-level fields are required.

```json
{
  "name": "Nucleus-count sensitivity",
  "rationale": "Check whether a declared segmentation-related alternative changes the local descriptive metric.",
  "revision_id": "rev_actual",
  "selection_id": "sel_actual",
  "source_label": "T_program",
  "target_label": "B_program",
  "radii_um": [15, 35, 75],
  "graph_scopes": ["whole_slice", "roi_induced"],
  "min_effect": 0.1,
  "background": {"kind": "neighbor_universe"},
  "variants": [
    {
      "name": "non_single_nucleus_excluded",
      "rationale": "Computational exclusion assumption, not a confirmed error or a recommended QC rule.",
      "selector": {"predicate": {"all": [
        {"field": "attributes.nucleus_count", "op": "ne", "value": 1}
      ]}},
      "changes": {"included": false}
    }
  ]
}
```

A selector can instead be `{"cell_ids":["exact-observation-id"]}`. Predicate fields are `label`, `included`, `region`, `attributes.KEY`, or `counts.EXACT_FEATURE_ID`; operators are `eq`, `ne`, `gt`, `ge`, `lt`, `le`, and `in`. Predicates are bounded data declarations, not executable expressions. They support up to eight AND clauses. Type-incompatible comparisons, missing/null attributes and unmeasured feature IDs yield unknown. Integers and floats are comparable, but booleans only compare with booleans and strings only with strings; the string `"0"` is never coerced to numeric zero. For `in`, a compatible match determines membership; without a match, any incompatible candidate retains unknown. In an AND predicate a false clause determines nonmembership; otherwise unknown remains unknown. An absent sparse count is zero only when that exact feature occurs in the measured panel. Alias resolution is not guessed.

At freeze, predicates resolve against the pinned revision into sorted matched IDs and unknown IDs. Both lists are stored and hashed. A predicate with unknown observations makes the corresponding after-state unknown in every grid cell; known matches are not silently applied as though they were complete. A fully resolved empty predicate is a recorded zero-change alternative. Invalid explicit IDs prevent freezing. Run and replay recheck that the frozen resolution agrees with the pinned inputs.

At most eight variants, six distinct positive radii, and two graph scopes are supported. Every variant independently starts from the same baseline, and can only change `label` and/or `included`.

## Graph and background are separate choices

`whole_slice` means all included observations in the imported project, which may be a cropped window. `roi_induced` restricts candidate neighbor observations to the frozen ROI. Sources always remain included, selected observations bearing the declared source label.

Background is independently declared as:

- `{"kind":"neighbor_universe"}`: use the candidate neighbor universe for each graph scope.
- `{"kind":"imported_universe"}`: use all included imported observations, including when the graph is ROI-induced.
- `{"kind":"selection","selection_id":"sel_background"}`: use a second exact selection captured in the same source, coordinate frame and revision. Its IDs are frozen too.

After each hypothetical modification, label abundance and inclusion in the declared background are recomputed. The source is removed from its own background only if it belongs to that background. For a same-label comparison, its target membership is also removed. This makes overlap and changing denominators explicit; it does not make the background biologically representative.

A defined background with eligible observations but zero target labels has an observed zero background fraction, provided the target group exists in the neighbor universe. It is different from an empty background or a zero denominator after removing the source, which remains unknown.

## Two statistical estimands

For each nonisolated source observation *i*, let `d_i` be its number of eligible nonself radius neighbors, `t_i` its target-label neighbor count, and `q_i` the target fraction in its eligible nonself declared background.

| Method | Observed fraction | Background fraction | Interpretation |
| --- | --- | --- | --- |
| `edge_weighted` | `sum(t_i) / sum(d_i)` | `sum(d_i*q_i) / sum(d_i)` | One vote per directed neighbor edge; dense source neighborhoods have more weight. |
| `source_equal_weighted` | `mean(t_i/d_i)` | `mean(q_i)` | One vote per source with at least one eligible neighbor. |

Both report `excess_over_abundance = observed_fraction - null_fraction` and apply the declared positive excess threshold. Equal-source weighting excludes isolated sources from the estimand; their number is explicitly retained alongside total source count and number of sources with neighbors. It is a reference estimand, not an independent annotation truth or a method-accuracy ranking. Different answers can reveal sensitivity to source degree, density and background composition.

The reference test fixture deliberately has source degrees 3, 1 and 0. With one target neighbor at each nonisolated source, edge weighting gives `2/4 = 0.5`; equal-source weighting gives `(1/3 + 1)/2 = 2/3`. The isolated third source is counted but not silently assigned a zero fraction.

## Complete outputs and failure retention

Runs use `analysis_schema=spatial-collab.hypothesis-plan-run.v2`. Each `results` row has a stable `row_id`, variant, radius, graph/background scopes, matched/unknown/change counts, before/after calculations, per-method comparisons and reference differences. A calculation contains `status`, `reason`, graph denominators, selected composition, and `methods.edge_weighted` / `methods.source_equal_weighted`. A method reports observed/background fractions, excess, classification and its denominator.

Row status is `computed`, `unknown`, or `failed`. Missing source/target groups, unavailable background, unknown selector membership or absent physical calibration remain unknown. Physical units are never inferred from an uncalibrated frame; such a plan may freeze, but all physical graph cells are unknown without evaluating distances. A graph resource or numeric failure is retained in that cell while the remaining grid continues. Malformed frozen plans or mismatched immutable inputs are rejected before execution.

`difference_order` contains **every** row ID, sorted by maximum absolute before/after excess change across both methods, then stable row ID; rows without calculable differences stay at the end. There is no winner-only list. `status_counts` makes incomplete results visible. Method disagreement and threshold crossings are descriptive, not significance, causal interaction, correct identity or evidence of no biological relationship.

Replay consumes the embedded frozen plan using the same pinned revision and selections. It does not need the sidecar database or current plan head. The original v1 schemas remain separately replayable.
