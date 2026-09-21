"""Durable scientific workflow execution in a configured, separate Python process."""

import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from . import objects
from .integration_jobs import _alive
from .store import SpatialError, _hash, _id, _json, _now


def _packages():
    return dict(sorted((d.metadata["Name"], d.version) for d in metadata.distributions()))


def runtime(project, refresh=False):
    config = project.root / "analysis-runtime.json"
    record = (
        json.loads(config.read_text(encoding="utf-8"))
        if config.exists()
        else {"python": sys.executable, "packages": _packages(), "environment": "host_default"}
    )
    return _probe(record["python"]) if refresh else record


def _probe(python):
    executable = Path(python).resolve(strict=True)
    code = "import importlib.metadata as m,json,sys; print(json.dumps({'python':sys.executable,'packages':dict(sorted((d.metadata['Name'],d.version) for d in m.distributions()))}))"
    probe = subprocess.run(
        [str(executable), "-c", code],
        capture_output=True,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        text=True,
        check=True,
        timeout=60,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
    )
    record = json.loads(probe.stdout)
    record["configured_at"] = _now()
    return record


def configure(project, python):
    record = _probe(python)
    (project.root / "analysis-runtime.json").write_text(_json(record), encoding="utf-8")
    return record


def _schema(db):
    db.execute(
        "CREATE TABLE IF NOT EXISTS analysis_jobs (id TEXT PRIMARY KEY, recipe TEXT NOT NULL, cache_key TEXT NOT NULL, status TEXT NOT NULL, result_id TEXT, error TEXT, worker_pid INTEGER, created_at TEXT NOT NULL, updated_at TEXT NOT NULL)"
    )


def get(project, job_id):
    with project._db() as db:
        exists = db.execute("SELECT 1 FROM sqlite_master WHERE name='analysis_jobs'").fetchone()
        row = db.execute("SELECT * FROM analysis_jobs WHERE id=?", (job_id,)).fetchone() if exists else None
    if row is None:
        raise SpatialError("Unknown analysis job.")
    result = dict(row)
    result["recipe"] = json.loads(result["recipe"])
    return result


def list_jobs(project):
    with project._db() as db:
        exists = db.execute("SELECT 1 FROM sqlite_master WHERE name='analysis_jobs'").fetchone()
        ids = (
            db.execute("SELECT id FROM analysis_jobs ORDER BY rowid DESC LIMIT 50").fetchall()
            if exists
            else []
        )
    return {"jobs": [get(project, row[0]) for row in ids]}


def submit(project, spec, launch=True, force=False):
    from .workflow_methods import METHODS

    if (
        not isinstance(spec, dict)
        or spec.get("method") not in METHODS
        or not isinstance(spec.get("input_ids"), list)
    ):
        raise SpatialError("Choose an executable method and immutable input IDs.")
    recipe = {"seed": 0, "parameters": {}, **spec}
    if not isinstance(recipe["parameters"], dict):
        raise SpatialError("Method parameters must be an object.")
    rt = configure(project, runtime(project)["python"])
    key = _cache_key(project, recipe)
    return _enqueue(project, recipe, key, rt, launch, force)


def _cache_key(project, recipe):
    inputs = [objects.get(project, oid, "analysisinput")["object_sha256"] for oid in recipe["input_ids"]]
    rt = {k: v for k, v in runtime(project).items() if k != "configured_at"}
    code = {
        name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (
            "workflow_methods.py",
            "workflow_inputs.py",
            "spatial_statistics.py",
            "workflow_jobs.py",
            "workflow_validation.py",
            "paste_compat.py",
        )
    }
    assay = None
    if recipe["parameters"].get("assay_id"):
        from .proteomics import get_assay

        assay = get_assay(project, recipe["parameters"]["assay_id"])["assay_sha256"]
    return _hash({"recipe": recipe, "inputs": inputs, "assay": assay, "runtime": rt, "code": code})


def _enqueue(project, recipe, key, rt, launch, force):
    with project._db(write=True) as db:
        _schema(db)
        old = db.execute(
            "SELECT id FROM analysis_jobs WHERE cache_key=? AND status IN ('queued','running','succeeded') ORDER BY rowid DESC LIMIT 1",
            (key,),
        ).fetchone()
        if old and not force:
            return {**get(project, old[0]), "cache_hit": True}
        job_id = _id("analysisjob")
        db.execute(
            "INSERT INTO analysis_jobs VALUES(?,?,?,'queued',NULL,NULL,NULL,?,?)",
            (job_id, _json(recipe), key, _now(), _now()),
        )
    if launch:
        logs = project.root / "job_logs"
        if not logs.resolve().is_relative_to(project.root):
            raise SpatialError("Job logs must remain in the project.")
        logs.mkdir(exist_ok=True)
        env = {
            **os.environ,
            "PYTHONPATH": str(Path(__file__).resolve().parents[1]),
            "PYTHONIOENCODING": "utf-8",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
        try:
            with (logs / (job_id + ".log")).open("wb") as log:
                subprocess.Popen(
                    [rt["python"], "-m", "spatial_collab.workflow_worker", str(project.root), job_id],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    env=env,
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
                )
        except OSError as exc:
            _finish(project, job_id, "failed", error=str(exc))
    return {**get(project, job_id), "cache_hit": False}


def _finish(project, job_id, status, result_id=None, error=None):
    with project._db(write=True) as db:
        db.execute(
            "UPDATE analysis_jobs SET status=?,result_id=?,error=?,updated_at=? WHERE id=? AND status!='cancelled'",
            (status, result_id, error, _now(), job_id),
        )


def cancel(project, job_id):
    old = get(project, job_id)
    status = "cancelled" if old["status"] == "queued" or not _alive(old["worker_pid"]) else "cancelling"
    if old["status"] in {"queued", "running", "cancelling"}:
        _finish(project, job_id, status)
    return get(project, job_id)


def retry(project, job_id):
    old = get(project, job_id)
    if old["status"] not in {"failed", "cancelled"}:
        raise SpatialError("Only failed/cancelled jobs can be retried. Cancel an orphaned job first.")
    return submit(project, old["recipe"], force=True)


def execute(project, job_id):
    from .workflow_methods import run

    with project._db(write=True) as db:
        if db.execute("SELECT 1 FROM analysis_jobs WHERE status IN ('running','cancelling')").fetchone():
            return get(project, job_id)
        changed = db.execute(
            "UPDATE analysis_jobs SET status='running',worker_pid=?,updated_at=? WHERE id=? AND status='queued'",
            (os.getpid(), _now(), job_id),
        ).rowcount
    if not changed:
        return get(project, job_id)
    start = time.monotonic()

    def checkpoint():
        if get(project, job_id)["status"] in {"cancelled", "cancelling"}:
            raise SpatialError(
                "Cancelled at method checkpoint; training calls may finish their current phase first."
            )
        if time.monotonic() - start > 7200:
            raise SpatialError("Two-hour compute budget reached at method checkpoint.")

    try:
        job = get(project, job_id)
        if _cache_key(project, job["recipe"]) != job["cache_key"]:
            raise SpatialError(
                "Analysis code, inputs or configured environment changed after submission; retry to freeze a new task."
            )
        if Path(runtime(project)["python"]).resolve() != Path(sys.executable).resolve():
            raise SpatialError("Worker interpreter differs from the configured analysis environment.")
        if runtime(project).get("packages") != _packages():
            raise SpatialError("Installed analysis packages changed after submission; reconfigure and retry.")
        result = run(project, job["recipe"], checkpoint)
        checkpoint()
        if runtime(project).get("packages") != _packages():
            raise SpatialError(
                "Installed packages changed while running; retry with a freshly frozen environment."
            )
        if _cache_key(project, job["recipe"]) != job["cache_key"]:
            raise SpatialError(
                "Analysis dependencies changed while running; saved output is not accepted by this task."
            )
        _finish(project, job_id, "succeeded", result["object_id"])
    except Exception as exc:
        import traceback

        traceback.print_exc()
        cancelled = get(project, job_id)["status"] in {"cancelled", "cancelling"}
        _finish(project, job_id, "cancelled" if cancelled else "failed", error=f"{type(exc).__name__}: {exc}")
    return get(project, job_id)


def results(project):
    rows = []
    for oid in objects.catalog(project, "analysisresult"):
        r = objects.get(project, oid, "analysisresult")
        rows.append(
            {
                "result_id": oid,
                "method": r["method"],
                "task": r["task"],
                "input_ids": [x["input_id"] for x in r["inputs"]],
            }
        )
    return {"results": rows}


def inspect(project, result_id, offset=0, limit=100, field=None, input_index=0):
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 1000:
        raise SpatialError("Use a nonnegative offset and limit 1..1000.")
    r = objects.get(project, result_id, "analysisresult")
    inputs = r["inputs"]
    if type(input_index) is not int or not 0 <= input_index < len(inputs):
        raise SpatialError("Unknown input index.")
    source = objects.get(project, inputs[input_index]["input_id"], "analysisinput")
    out = r["output"]
    rows = out.get("rows")
    if rows is None:
        keys = (
            list(out.get("abundance", {}))
            if "abundance" in out
            else [
                key
                for key in (
                    "domains",
                    "embedding",
                    "rna_contribution_fraction",
                    "target_barycentric_coordinates",
                    "pathway_activity",
                    "pathway_pvalues",
                    "pathway_qvalues",
                )
                if key in out
            ]
        )
        field = field or (keys[0] if keys else None)
        if field not in keys:
            raise SpatialError("Choose a supported per-observation output field.")
        values = out["abundance"][field] if "abundance" in out else out[field]
        all_ids = (
            [cid for inp in inputs for cid in inp["observation_ids"]]
            if r["method"] == "harmony"
            else inputs[0]["observation_ids"]
        )
        if r["method"] == "harmony":
            start = sum(len(v["observation_ids"]) for v in inputs[:input_index])
            values = values[start : start + len(source["observation_ids"])]
            all_ids = source["observation_ids"]
        elif input_index != 0:
            raise SpatialError(
                "This output is defined on the first spatial input, not the reference/target axis."
            )
        rows = [
            {
                "cell_id": cid,
                "value": value,
                "xy": source["coordinates"][i] if source["coordinates"] else None,
            }
            for i, (cid, value) in enumerate(zip(all_ids, values))
        ]
    current = (
        source["provenance"].get("project_id") == project.context()["project_id"]
        and source["provenance"].get("revision_id") == project.context()["head_revision"]
    )
    metadata_out = {
        k: v
        for k, v in out.items()
        if k
        not in {
            "rows",
            "embedding",
            "uncorrected_embedding",
            "controls",
            "factor_loadings",
            "transport",
            "abundance",
            "rna_contribution_fraction",
            "target_barycentric_coordinates",
            "domains",
            "assignment_probabilities",
            "input_pca_embedding",
            "pathway_activity",
            "pathway_pvalues",
            "pathway_qvalues",
            "resource_network",
        }
    }
    return {
        "result_id": result_id,
        "method": r["method"],
        "task": r["task"],
        "field": field,
        "total": len(rows),
        "offset": offset,
        "rows": rows[offset : offset + limit],
        "next_offset": offset + limit if offset + limit < len(rows) else None,
        "metadata": metadata_out,
        "same_current_project_revision": current,
        "source_input_id": source["object_id"],
        "recipe": r["recipe"],
    }
