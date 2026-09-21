"""Durable jobs with isolated workers, cooperative cancellation and recipe caches."""
from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time

from . import objects
from .store import SpatialError, _hash, _id, _json, _now


def _schema(db):
    db.execute("CREATE TABLE IF NOT EXISTS integration_jobs (id TEXT PRIMARY KEY, recipe TEXT NOT NULL, "
               "cache_key TEXT NOT NULL, status TEXT NOT NULL, result_id TEXT, error TEXT, created_at TEXT NOT NULL, updated_at TEXT NOT NULL)")
    if "worker_pid" not in {r[1] for r in db.execute("PRAGMA table_info(integration_jobs)")}:
        db.execute("ALTER TABLE integration_jobs ADD COLUMN worker_pid INTEGER")


def submit(project, spec, launch=True):
    from .integration import BACKEND_VERSION, numerical_environment
    from .proteomics import get_assay
    allowed = {"revision_id", "assay_id", "rna_features", "protein_features", "observation_ids", "components", "clusters", "seed", "protein_transform", "backend", "feature_selection"}
    if not isinstance(spec, dict) or set(spec) - allowed or not {"revision_id", "assay_id"} <= set(spec):
        raise SpatialError("Job requires revision_id, assay_id and only supported baseline parameters.")
    defaults = {"components": 10, "clusters": 6, "seed": 0, "protein_transform": "log1p", "backend": "balanced_pca", "feature_selection": "variance"}
    spec = {**defaults, **spec}
    revision = project.get_revision(spec["revision_id"])
    assay = get_assay(project, spec["assay_id"])
    environment = numerical_environment()
    # Workers have a fixed thread policy, independent of the calling host.
    environment["thread_environment"] = dict.fromkeys(("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"), "1")
    code_names = ("integration.py", "smopca_adapter.py", "integration_jobs.py", "integration_runners.py")
    environment["implementation_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                            for name in code_names if Path(__file__).with_name(name).exists()}
    if spec["backend"] == "smopca":
        from .smopca_adapter import provenance
        environment["smopca"] = provenance()
    cache_key = _hash({"spec": spec, "source": project.context()["source_sha256"],
                       "revision": revision["revision_sha256"], "assay": assay["assay_sha256"],
                       "backend": BACKEND_VERSION, "environment": environment})
    with project._db(write=True) as db:
        _schema(db)
        old = db.execute("SELECT id FROM integration_jobs WHERE cache_key=? AND status IN ('queued','running','succeeded') ORDER BY rowid DESC LIMIT 1", (cache_key,)).fetchone()
        if old:
            identifier = old[0]
        else:
            identifier = _id("integrationjob")
            db.execute("INSERT INTO integration_jobs(id,recipe,cache_key,status,result_id,error,created_at,updated_at) VALUES(?,?,?,'queued',NULL,NULL,?,?)",
                       (identifier, _json(spec), cache_key, _now(), _now()))
    if old:
        result = get_job(project, identifier)
        if result["status"] == "succeeded":
            objects.get(project, result["result_id"], "integration")
        return {**result, "cache_hit": result["status"] == "succeeded"}
    if launch:
        logs = project.root / "job_logs"
        if not logs.resolve().is_relative_to(project.root):
            raise SpatialError("Job logs must remain inside the project.")
        logs.mkdir(exist_ok=True)
        env = os.environ.copy()
        env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        try:
            with (logs / (identifier + ".log")).open("wb") as stream:
                subprocess.Popen([sys.executable, "-m", "spatial_collab.integration_worker", str(project.root), identifier],
                                 stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=env,
                                 creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        except OSError as exc:
            _update(project, identifier, "failed", error=str(exc))
    return {**get_job(project, identifier), "cache_hit": False}


def get_job(project, job_id):
    with project._db() as db:
        exists = db.execute("SELECT 1 FROM sqlite_master WHERE name='integration_jobs'").fetchone()
        row = db.execute("SELECT * FROM integration_jobs WHERE id=?", (job_id,)).fetchone() if exists else None
    if row is None:
        raise SpatialError("Unknown integration job.")
    result = dict(row)
    result["recipe"] = json.loads(result["recipe"])
    result["cancellation"] = "Cooperative at phase boundaries; an active linear algebra call may finish first."
    result["recovery"] = "Cancel a stranded queued/running job, then retry; previous failure remains recorded."
    return result


def _update(project, identifier, status, result_id=None, error=None):
    with project._db(write=True) as db:
        db.execute("UPDATE integration_jobs SET status=?,result_id=?,error=?,updated_at=? WHERE id=? AND status != 'cancelled' AND (status != 'cancelling' OR ?='cancelled')",
                   (status, result_id, error, _now(), identifier, status))


def _alive(pid):
    if not pid:
        return False
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.windll.kernel32
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.GetExitCodeProcess.restype = wintypes.BOOL
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel.CloseHandle.restype = wintypes.BOOL
        handle = kernel.OpenProcess(0x1000, False, pid)
        if not handle:
            return False
        code = ctypes.c_ulong()
        try:
            return bool(kernel.GetExitCodeProcess(handle, ctypes.byref(code))) and code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False


def cancel(project, job_id):
    old = get_job(project, job_id)
    with project._db(write=True) as db:
        db.execute("UPDATE integration_jobs SET status=CASE WHEN status='queued' OR ? THEN 'cancelled' ELSE 'cancelling' END,updated_at=? WHERE id=? AND status IN ('queued','running','cancelling')", (not _alive(old.get("worker_pid")), _now(), job_id))
    return get_job(project, job_id)


def retry(project, job_id):
    record = get_job(project, job_id)
    if record["status"] not in {"failed", "cancelled"}:
        raise SpatialError("Only failed/cancelled jobs can be retried; cancel stranded jobs first.")
    return submit(project, record["recipe"])


def execute(project, job_id):
    from .integration import run_baseline
    with project._db(write=True) as db:
        if db.execute("SELECT 1 FROM integration_jobs WHERE status IN ('running','cancelling')").fetchone():
            return get_job(project, job_id)
        changed = db.execute("UPDATE integration_jobs SET status='running',updated_at=?,worker_pid=? WHERE id=? AND status='queued'", (_now(), os.getpid(), job_id)).rowcount
    if not changed:
        return get_job(project, job_id)
    def checkpoint():
        state = get_job(project, job_id)
        if state["status"] in {"cancelled", "cancelling"}:
            raise SpatialError("Job cancelled by user.")
        if time.monotonic() - start > 1800:
            raise SpatialError("Job exceeded the 30 minute compute budget.")
        _update(project, job_id, "running")
    start = time.monotonic()
    try:
        result = run_baseline(project, **get_job(project, job_id)["recipe"], checkpoint=checkpoint)
        checkpoint()
        _update(project, job_id, "succeeded", result["object_id"])
    except Exception as exc:
        status = "cancelled" if get_job(project, job_id)["status"] == "cancelling" else "failed"
        _update(project, job_id, status, error=f"{type(exc).__name__}: {exc}")
    return get_job(project, job_id)
