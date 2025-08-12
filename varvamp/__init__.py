from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
import sys

source_location = Path(__file__).resolve()
pyproject = source_location.parent.parent / "pyproject.toml"

if sys.version_info >= (3, 11) and pyproject.exists():
    import tomllib
    with open(pyproject, "rb") as f:
        __version__ = tomllib.load(f).get("project", {}).get("version")
else:
    try:
        __version__ = version("varvamp")
    except PackageNotFoundError:
        __version__ = "unknown"
