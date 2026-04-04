"""forApollo orbital module -- coordinate transforms (ECI, ECEF, geodetic)."""

import ctypes
from ctypes import c_int, c_double, POINTER

from forapollo import lib as _lib

_f64p = POINTER(c_double)


def _make_info():
    info = c_int(0)
    return info, ctypes.byref(info)


# ---------------------------------------------------------------------------
# Coordinate transforms (forapollo_coords_*)
# ---------------------------------------------------------------------------

_lib.forapollo_coords_transform.argtypes = [c_int, c_int, c_int, _f64p, _f64p, c_double, _f64p, c_int, POINTER(c_int)]
_lib.forapollo_coords_transform.restype = None

_lib.forapollo_coords_eci_to_ecef.argtypes = [_f64p, _f64p, c_double, POINTER(c_int)]
_lib.forapollo_coords_eci_to_ecef.restype = None

_lib.forapollo_coords_ecef_to_eci.argtypes = [_f64p, _f64p, c_double, POINTER(c_int)]
_lib.forapollo_coords_ecef_to_eci.restype = None

_lib.forapollo_coords_ecef_to_geodetic.argtypes = [_f64p, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_double, c_double, POINTER(c_int)]
_lib.forapollo_coords_ecef_to_geodetic.restype = None

_lib.forapollo_coords_geodetic_to_ecef.argtypes = [c_double, c_double, c_double, _f64p, c_double, c_double, POINTER(c_int)]
_lib.forapollo_coords_geodetic_to_ecef.restype = None


def coords_transform(from_id: int, to_id: int, n: int, x_in, x_out, t: float, params, np_: int) -> int:
    info, ip = _make_info()
    _lib.forapollo_coords_transform(from_id, to_id, n, x_in, x_out, t, params, np_, ip)
    return info.value

def coords_eci_to_ecef(x_eci, x_ecef, gmst: float) -> int:
    info, ip = _make_info()
    _lib.forapollo_coords_eci_to_ecef(x_eci, x_ecef, gmst, ip)
    return info.value

def coords_ecef_to_eci(x_ecef, x_eci, gmst: float) -> int:
    info, ip = _make_info()
    _lib.forapollo_coords_ecef_to_eci(x_ecef, x_eci, gmst, ip)
    return info.value

def coords_ecef_to_geodetic(x_ecef, a_body: float = 6378137.0, f_body: float = 1.0/298.257223563) -> tuple[int, float, float, float]:
    info, ip = _make_info()
    lat, lon, alt = c_double(), c_double(), c_double()
    _lib.forapollo_coords_ecef_to_geodetic(x_ecef, ctypes.byref(lat), ctypes.byref(lon), ctypes.byref(alt), a_body, f_body, ip)
    return info.value, lat.value, lon.value, alt.value

def coords_geodetic_to_ecef(lat: float, lon: float, alt: float, x_ecef, a_body: float = 6378137.0, f_body: float = 1.0/298.257223563) -> int:
    info, ip = _make_info()
    _lib.forapollo_coords_geodetic_to_ecef(lat, lon, alt, x_ecef, a_body, f_body, ip)
    return info.value
