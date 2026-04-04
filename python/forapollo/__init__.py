"""forApollo -- Python ctypes bindings for libforapollo (forapollo_* C ABI).

Universal state estimation, navigation, and guidance.
The same Kalman filter that guided Apollo to the Moon.
"""

import ctypes
import os
import platform
from pathlib import Path

__version__ = "0.1.0"

_LIB_NAMES = {
    "Linux": "libforapollo.so",
    "Darwin": "libforapollo.dylib",
    "Windows": "forapollo.dll",
}


def _find_library() -> ctypes.CDLL:
    lib_name = _LIB_NAMES.get(platform.system(), "libforapollo.so")
    search_dirs: list[Path] = []

    env_path = os.environ.get("FORAPOLLO_LIB_PATH")
    if env_path:
        search_dirs.append(Path(env_path))

    pkg_dir = Path(__file__).resolve().parent
    search_dirs.append(pkg_dir)
    search_dirs.append(pkg_dir.parent)

    repo_root = pkg_dir.parent.parent
    search_dirs.append(repo_root / "zig-out" / "lib")
    search_dirs.append(repo_root / "prebuilt" / f"{platform.machine()}-linux")

    if os.path.exists("/proc/device-tree/model"):
        search_dirs.append(Path("/usr/local/lib"))
        search_dirs.append(Path("/usr/lib/aarch64-linux-gnu"))

    for d in search_dirs:
        candidate = d / lib_name
        if candidate.exists():
            return ctypes.CDLL(str(candidate))

    return ctypes.CDLL(lib_name)


lib = _find_library()
