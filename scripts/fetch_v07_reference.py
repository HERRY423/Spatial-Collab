"""Fetch the public scvi-tools spleen/lymph single-cell reference with byte verification."""
import argparse
import hashlib
from pathlib import Path
import urllib.request

URL = "https://exampledata.scverse.org/scvi-tools/sln_111.h5ad"
SHA256 = "2098ca88739a9d2d789cb092481c326d7bcbfa73ca815ce3e80b7668cffbe452"


def fetch(destination):
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        temporary = destination.with_suffix(destination.suffix + ".download")
        with urllib.request.urlopen(URL, timeout=120) as response, temporary.open("xb") as stream:
            while block := response.read(1024 * 1024):
                stream.write(block)
                if stream.tell() > 100_000_000:
                    raise ValueError("Unexpected reference download size; temporary file retained for inspection")
        with temporary.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != SHA256:
            raise ValueError("Public reference differs from validation bytes; download retained, not accepted")
        temporary.rename(destination)
    with destination.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != SHA256:
            raise ValueError("Existing reference differs from validation bytes; not overwritten")
    print(destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    fetch(parser.parse_args().destination)
