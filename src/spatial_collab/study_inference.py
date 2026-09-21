"""Subject-aware continuous-outcome inference; explicit designs and failed fits."""
from collections import defaultdict
from importlib.metadata import version
import warnings

import numpy as np
from scipy import stats

from . import objects
from .spatial_statistics import bh
from .store import SpatialError, _hash, _text


def infer(project, study_id, records, plan):
    design = objects.get(project, study_id, "study")
    allowed = {"method", "control", "case", "metrics", "roi_class", "adjust_batch", "covariates", "missing_policy", "outcome_scale", "rationale"}
    if not isinstance(plan, dict) or set(plan)-allowed:
        raise SpatialError("Unknown study inference plan fields.")
    method = plan.get("method")
    if method not in {"welch", "paired_t", "mixedlm"}:
        raise SpatialError("Choose welch, paired_t or mixedlm.")
    for k in ("control", "case", "roi_class", "rationale"):
        _text(plan.get(k), k)
    if plan["control"] == plan["case"]:
        raise SpatialError("Contrast requires distinct conditions.")
    if plan.get("outcome_scale") != "continuous_section_summary":
        raise SpatialError("These Gaussian models require continuous section summaries; raw counts/proportions need an appropriate model or explicit upstream transformation.")
    metrics = plan.get("metrics")
    if not isinstance(metrics, list) or not 1 <= len(metrics) <= 100 or len(set(metrics)) != len(metrics):
        raise SpatialError("Declare 1..100 unique metrics as one complete testing family.")
    for m in metrics:
        _text(m, "metric")
    if plan.get("missing_policy", "reject") not in {"reject", "complete_subjects"}:
        raise SpatialError("Declare reject or complete_subjects missing policy.")
    adjust = plan.get("adjust_batch", True)
    if type(adjust) is not bool:
        raise SpatialError("adjust_batch must be boolean.")
    covariates = plan.get("covariates", [])
    if not isinstance(covariates, list) or len(set(covariates)) != len(covariates) or len(covariates) > 8:
        raise SpatialError("Choose up to eight explicit numeric covariates.")
    if method != "mixedlm" and covariates:
        raise SpatialError("Covariate adjustment requires mixedlm.")
    sections = {s["section_id"]: s for s in design["sections"] if s["condition"] in {plan["control"], plan["case"]} and s["roi_class"] == plan["roi_class"]}
    if not sections or {s["condition"] for s in sections.values()} != {plan["control"], plan["case"]}:
        raise SpatialError("Both conditions must exist in the selected ROI class.")
    sample_owners = defaultdict(set)
    for s in sections.values():
        sample_owners[s["sample_id"]].add((s["subject_id"], s["condition"]))
    if any(len(owners) != 1 for owners in sample_owners.values()):
        raise SpatialError("One sample ID cannot belong to different subjects or conditions.")
    if any(s["roi_origin"] == "posthoc_exploratory" for s in sections.values()):
        exploratory = True
    else:
        exploratory = False
    if not isinstance(records, list) or len(records) > 100000:
        raise SpatialError("Provide bounded per-section outcome records.")
    measured = {}
    for r in records:
        sid, metric = r.get("section_id"), r.get("metric")
        if sid not in sections or metric not in metrics or (sid, metric) in measured:
            raise SpatialError("Unknown, out-of-scope or duplicate section/metric record.")
        s = sections[sid]
        if r.get("source_sha256") != s["source_sha256"] or r.get("run_id") not in s["run_ids"]:
            raise SpatialError("Outcome does not match registered source and run identity.")
        for k in ("units", "method_version"):
            _text(r.get(k), k)
        value = r.get("value")
        if type(value) not in (int, float) or not np.isfinite(value):
            raise SpatialError("Finite continuous outcome required; missing values must remain missing records.")
        measured[(sid, metric)] = r
    output = []
    for metric in metrics:
        missing = [sid for sid in sections if (sid, metric) not in measured]
        excluded = {sections[s]["subject_id"] for s in missing}
        result = {"metric": metric, "status": "not_tested", "p_value": None, "q_value": None, "missing_sections": missing, "excluded_subjects": sorted(excluded), "effect_definition": f"{plan['case']} minus {plan['control']}"}
        output.append(result)
        if missing and plan.get("missing_policy", "reject") == "reject":
            result["reason"] = "missing_declared_sections"
            continue
        selected = [s for s in sections.values() if s["subject_id"] not in excluded]
        try:
            _fit(selected, measured, metric, plan, result, method, adjust, covariates)
        except (SpatialError, np.linalg.LinAlgError, ValueError) as exc:
            result.update(status="not_tested", reason=str(exc), p_value=None)
    # Failed hypotheses remain in the family as p=1; cannot improve other q's
    # by silently dropping difficult/nonconvergent outcomes.
    corrected = bh([r["p_value"] if r["p_value"] is not None else 1.0 for r in output])
    for row, q in zip(output, corrected):
        row["q_value"] = float(q) if row["p_value"] is not None else None
    def finite(value):
        if isinstance(value, dict):
            return {k: finite(v) for k, v in value.items()}
        if isinstance(value, list):
            return [finite(v) for v in value]
        return None if isinstance(value, float) and not np.isfinite(value) else value
    return objects.put(project, "studyresult", {"schema": "spatial-collab.study-inference.v1", "study_id": study_id, "study_sha256": design["object_sha256"], "records": records, "records_sha256": _hash(records), "plan": plan, "results": finite(output), "family_size": len(metrics), "multiple_testing": "BH over declared metrics; failed tests retained as p=1", "analysis_unit": "subject; sections are repeated measurements", "roi_inference": "exploratory_posthoc" if exploratory else "declared_prespecified_or_independent", "environment": {"numpy": version("numpy"), "scipy": version("scipy"), **({"statsmodels": version("statsmodels")} if method == "mixedlm" else {})}, "limitations": ["Researcher-declared subject identities and source/run bindings are not independent experimental verification.", "Observational contrasts do not establish causality; Gaussian assumptions must be reviewed.", "No small-sample Satterthwaite or Kenward-Roger correction for MixedLM; Wald inference is asymptotic."], "scientific_authorization": "NOT_ESTABLISHED"})


def _fit(sections, measured, metric, plan, result, method, adjust, covariates):
    grouped = defaultdict(list)
    units, versions = set(), set()
    for s in sections:
        r = measured[(s["section_id"], metric)]
        grouped[(s["subject_id"], s["condition"])].append(r["value"])
        units.add(r["units"])
        versions.add(r["method_version"])
    if len(units) != 1 or len(versions) != 1:
        raise SpatialError("Units and generating method version must match within a metric.")
    control, case = plan["control"], plan["case"]
    ids = {c: {s for s, condition in grouped if condition == c} for c in (control, case)}
    result.update(subjects_per_condition={c: len(ids[c]) for c in ids}, section_count=len(sections), units=next(iter(units)), method_version=next(iter(versions)))
    if min(map(len, ids.values())) < 3:
        raise SpatialError("At least three subjects per condition/pairs are required; more may be needed for reliable inference.")
    batches = sorted({s["batch"] for s in sections})
    x = np.column_stack([np.ones(len(sections)), [float(s["condition"] == case) for s in sections], *[[float(s["batch"] == b) for s in sections] for b in batches[1:]]])
    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise SpatialError("Condition and batch design is rank deficient; contrast not identifiable.")
    if method != "mixedlm" and len(batches) > 1:
        raise SpatialError("Use mixedlm for multiple batches; unadjusted t tests are not silently substituted.")
    if method == "welch":
        if ids[control] & ids[case]:
            raise SpatialError("Subjects occur in both conditions; use paired_t or mixedlm.")
        a = np.array([np.mean(grouped[(s, control)]) for s in sorted(ids[control])])
        b = np.array([np.mean(grouped[(s, case)]) for s in sorted(ids[case])])
        test = stats.ttest_ind(b, a, equal_var=False)
        effect = float(b.mean()-a.mean())
        ci = test.confidence_interval()
        result.update(effect=effect, statistic=float(test.statistic), degrees_of_freedom=float(test.df), confidence_interval=[float(ci.low), float(ci.high)], subject_values={control: a.tolist(), case: b.tolist()}, weighting="equal section mean within subject, equal subjects", p_value=float(test.pvalue))
    elif method == "paired_t":
        if ids[control] != ids[case]:
            raise SpatialError("Paired contrast requires the same subjects in both conditions; no implicit unmatched-subject exclusion.")
        ordered = sorted(ids[control])
        a = np.array([np.mean(grouped[(s, control)]) for s in ordered])
        b = np.array([np.mean(grouped[(s, case)]) for s in ordered])
        if np.std(b-a, ddof=1) <= 1e-14:
            raise SpatialError("Paired differences have zero residual variance.")
        test = stats.ttest_rel(b, a)
        ci = test.confidence_interval()
        result.update(effect=float((b-a).mean()), statistic=float(test.statistic), degrees_of_freedom=float(test.df), confidence_interval=[float(ci.low), float(ci.high)], subject_ids=ordered, subject_differences=(b-a).tolist(), weighting="equal section mean within subject/condition, paired by subject_id", p_value=float(test.pvalue))
    else:
        from statsmodels.regression.mixed_linear_model import MixedLM
        subject_ids = [s["subject_id"] for s in sections]
        if len(set(subject_ids)) < 6 or len(set(subject_ids)) == len(sections):
            raise SpatialError("Random intercept requires at least six subjects and repeated sections/conditions.")
        if not adjust and len(batches) > 1:
            raise SpatialError("Batch adjustment is required for a multibatch mixed model.")
        terms = ["intercept", "case_vs_control", *[f"batch:{b}" for b in batches[1:]]]
        for name in covariates:
            values = [s.get("covariates", {}).get(name) for s in sections]
            if any(type(v) not in (float, int) or not np.isfinite(v) for v in values):
                raise SpatialError(f"Missing/nonfinite numeric covariate: {name}")
            values = np.asarray(values, dtype=float)
            if np.std(values) == 0:
                raise SpatialError("Constant covariate is not identifiable.")
            x = np.column_stack([x, (values-values.mean())/values.std()])
            terms.append(name)
        if len(sections) <= x.shape[1] or np.linalg.matrix_rank(x) != x.shape[1]:
            raise SpatialError("Fixed effects are rank deficient or have no residual degrees of freedom.")
        y = np.array([measured[(s["section_id"], metric)]["value"] for s in sections])
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            fit = MixedLM(y, x, groups=subject_ids).fit(reml=True, method="lbfgs", maxiter=1000, disp=False)
        warning_text = [str(w.message) for w in captured]
        if len(set(subject_ids)) < 20:
            warning_text.append("Fewer than 20 subjects: asymptotic Wald intervals and p-values may be poorly calibrated; seek small-sample or design-specific validation.")
        var = float(fit.cov_re[0, 0])
        result.update(converged=bool(fit.converged), random_intercept_variance=var, residual_variance=float(fit.scale), warnings=warning_text, fixed_effect_terms=terms, fixed_effect_estimates=np.asarray(fit.fe_params).tolist(), subject_count=len(set(subject_ids)), inference="asymptotic_two_sided_normal_Wald_REML", design_rank=int(np.linalg.matrix_rank(x)))
        if not fit.converged or var <= 1e-8*max(float(fit.scale), 1e-12) or any("singular" in w.lower() or "positive definite" in w.lower() for w in warning_text):
            raise SpatialError("Mixed model nonconvergent or singular; effect diagnostics retained, inference withheld.")
        ci = fit.conf_int()[1]
        result.update(effect=float(fit.fe_params[1]), standard_error=float(fit.bse_fe[1]), confidence_interval=ci.tolist(), p_value=float(fit.pvalues[1]), weighting="section outcomes modeled with a subject random intercept, not independent cells")
    if result["p_value"] is None or not np.isfinite(result["p_value"]) or not np.isfinite(result["confidence_interval"]).all():
        raise SpatialError("Undefined test; zero variance or degenerate design.")
    result["status"] = "tested"
