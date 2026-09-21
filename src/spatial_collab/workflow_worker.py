"""Subprocess entrypoint, serial queue per project."""

import sys
from .store import Project
from .workflow_jobs import execute

if __name__ == "__main__":
    p, identifier = Project(sys.argv[1]), sys.argv[2]
    while identifier:
        result = execute(p, identifier)
        print(result["status"], flush=True)
        if result["status"] in {"queued", "running", "cancelling"}:
            break
        with p._db() as db:
            row = db.execute(
                "SELECT id FROM analysis_jobs WHERE status='queued' ORDER BY rowid LIMIT 1"
            ).fetchone()
            identifier = row[0] if row else None
