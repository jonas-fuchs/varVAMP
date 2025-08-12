from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("varvamp")
except PackageNotFoundError:
    __version__ = "unknown"
