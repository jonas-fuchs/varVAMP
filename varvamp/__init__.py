import importlib.metadata, pathlib, tomllib

# get __version__ from pyproject.toml
source_location = pathlib.Path("__file__").parent
if (source_location.parent / "pyproject.toml").exists():
    with open(source_location.parent / "pyproject.toml", "rb") as f:
        __version__ = tomllib.load(f)['project']['version']
else:
    __version__ = importlib.metadata.version("varvamp")
