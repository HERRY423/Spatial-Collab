"""Internal subprocess entrypoint; writes only to a pinned project."""
import sys

from .integration_jobs import execute
from .store import Project

if __name__ == "__main__":
    project, job_id = Project(sys.argv[1]), sys.argv[2]
    while job_id:
        result = execute(project, job_id)
        print(result["status"], flush=True)
        if result["status"] in {"queued", "running", "cancelling"}:
            break  # Another worker owns the single per-project compute slot.
        with project._db() as db:
            row = db.execute("SELECT id FROM integration_jobs WHERE status='queued' ORDER BY rowid LIMIT 1").fetchone()
            job_id = row[0] if row else None
