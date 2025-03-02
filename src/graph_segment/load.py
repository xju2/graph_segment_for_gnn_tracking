import ctypes
from pathlib import Path


PREFERRED_LOAD_FLAG = ctypes.RTLD_LOCAL
this_file_path = Path(__file__).resolve().parent

def _load_wheel_installation(soname: str) -> ctypes.CDLL:
    lib = this_file_path / "lib64" / soname
    if lib.exists():
        return ctypes.CDLL(lib, PREFERRED_LOAD_FLAG)
    return None

def load_library():

    libs_to_return = []
    for soname in ["libgs_core.so"]:
        lib = _load_wheel_installation(soname)
        if lib:
            libs_to_return.append(lib)
    return libs_to_return
