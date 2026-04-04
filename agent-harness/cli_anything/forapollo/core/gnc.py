"""Core domain types for guidance, navigation, and control operations."""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class StateVector:
    """Universal state vector -- the same pattern that guided Apollo."""
    x: np.ndarray          # state (position, velocity, etc.)
    P: np.ndarray          # covariance (flat n*n)
    n: int                 # state dimension

    def to_dict(self) -> dict:
        return {
            "state": self.x.tolist(),
            "covariance_diag": np.diag(self.P.reshape(self.n, self.n)).tolist(),
            "n": self.n,
        }


@dataclass
class GuidanceCommand:
    acceleration: np.ndarray  # commanded acceleration
    info: int                 # 0=ok, nonzero=error
    law: str                  # "zem", "zev", "pn_pure", "lambert", etc.

    def to_dict(self) -> dict:
        return {
            "acceleration": self.acceleration.tolist(),
            "info": self.info,
            "law": self.law,
        }


@dataclass
class OrbitalElements:
    a: float       # semi-major axis (m)
    e: float       # eccentricity
    i: float       # inclination (rad)
    omega: float   # argument of periapsis (rad)
    Omega: float   # RAAN (rad)
    nu: float      # true anomaly (rad)

    def to_dict(self) -> dict:
        return {
            "a": self.a, "e": self.e, "i": self.i,
            "omega": self.omega, "Omega": self.Omega, "nu": self.nu,
        }


@dataclass
class GeodeticPosition:
    lat: float   # latitude (rad)
    lon: float   # longitude (rad)
    alt: float   # altitude (m)

    def to_dict(self) -> dict:
        return {"lat_deg": np.degrees(self.lat), "lon_deg": np.degrees(self.lon),
                "alt_m": self.alt}


@dataclass
class LambertSolution:
    v1: np.ndarray  # departure velocity
    v2: np.ndarray  # arrival velocity
    tof: float      # time of flight
    info: int

    def to_dict(self) -> dict:
        return {
            "v1": self.v1.tolist(),
            "v2": self.v2.tolist(),
            "tof": self.tof,
            "delta_v1": float(np.linalg.norm(self.v1)),
            "delta_v2": float(np.linalg.norm(self.v2)),
            "info": self.info,
        }
