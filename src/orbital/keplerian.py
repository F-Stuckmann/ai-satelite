"""Keplerian orbital mechanics calculations.

Provides two-body orbital mechanics for satellite orbit analysis.
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

from ..core.constants import EARTH_RADIUS_KM, EARTH_MU


# Earth gravitational parameter
MU_EARTH = 398600.4418  # km³/s²


@dataclass
class OrbitalElements:
    """Classical Keplerian orbital elements.

    Attributes:
        a: Semi-major axis (km)
        e: Eccentricity (0 = circular, 0-1 = elliptical)
        i: Inclination (degrees)
        raan: Right Ascension of Ascending Node (degrees)
        aop: Argument of Periapsis (degrees)
        ta: True Anomaly (degrees)
    """
    a: float  # Semi-major axis (km)
    e: float = 0.0  # Eccentricity
    i: float = 0.0  # Inclination (deg)
    raan: float = 0.0  # RAAN (deg)
    aop: float = 0.0  # Argument of periapsis (deg)
    ta: float = 0.0  # True anomaly (deg)

    @classmethod
    def from_altitude(
        cls,
        altitude_km: float,
        inclination_deg: float = 0.0,
        eccentricity: float = 0.0,
        raan_deg: float = 0.0,
        aop_deg: float = 0.0,
        ta_deg: float = 0.0,
    ) -> "OrbitalElements":
        """Create orbital elements from altitude above Earth.

        Args:
            altitude_km: Altitude above Earth's surface (km)
            inclination_deg: Orbital inclination (degrees)
            eccentricity: Orbital eccentricity (0 for circular)
            raan_deg: Right Ascension of Ascending Node (degrees)
            aop_deg: Argument of periapsis (degrees)
            ta_deg: True anomaly (degrees)
        """
        a = EARTH_RADIUS_KM + altitude_km
        return cls(
            a=a,
            e=eccentricity,
            i=inclination_deg,
            raan=raan_deg,
            aop=aop_deg,
            ta=ta_deg,
        )

    @property
    def period_seconds(self) -> float:
        """Orbital period in seconds."""
        return 2 * np.pi * np.sqrt(self.a**3 / MU_EARTH)

    @property
    def period_minutes(self) -> float:
        """Orbital period in minutes."""
        return self.period_seconds / 60

    @property
    def mean_motion(self) -> float:
        """Mean motion in radians per second."""
        return np.sqrt(MU_EARTH / self.a**3)

    @property
    def velocity_circular(self) -> float:
        """Circular orbital velocity in km/s."""
        return np.sqrt(MU_EARTH / self.a)

    @property
    def velocity_periapsis(self) -> float:
        """Velocity at periapsis in km/s."""
        r_p = self.a * (1 - self.e)
        return np.sqrt(MU_EARTH * (2/r_p - 1/self.a))

    @property
    def velocity_apoapsis(self) -> float:
        """Velocity at apoapsis in km/s."""
        r_a = self.a * (1 + self.e)
        return np.sqrt(MU_EARTH * (2/r_a - 1/self.a))

    @property
    def altitude_periapsis(self) -> float:
        """Altitude at periapsis in km."""
        return self.a * (1 - self.e) - EARTH_RADIUS_KM

    @property
    def altitude_apoapsis(self) -> float:
        """Altitude at apoapsis in km."""
        return self.a * (1 + self.e) - EARTH_RADIUS_KM

    @property
    def specific_energy(self) -> float:
        """Specific orbital energy in km²/s²."""
        return -MU_EARTH / (2 * self.a)

    @property
    def specific_angular_momentum(self) -> float:
        """Specific angular momentum in km²/s."""
        return np.sqrt(MU_EARTH * self.a * (1 - self.e**2))


def orbital_velocity(altitude_km: float) -> float:
    """Calculate circular orbital velocity at given altitude.

    Args:
        altitude_km: Altitude above Earth's surface (km)

    Returns:
        Orbital velocity in km/s
    """
    r = EARTH_RADIUS_KM + altitude_km
    return np.sqrt(MU_EARTH / r)


def orbital_period(altitude_km: float) -> float:
    """Calculate orbital period at given altitude.

    Args:
        altitude_km: Altitude above Earth's surface (km)

    Returns:
        Orbital period in minutes
    """
    a = EARTH_RADIUS_KM + altitude_km
    t_seconds = 2 * np.pi * np.sqrt(a**3 / MU_EARTH)
    return t_seconds / 60


def orbits_per_day(altitude_km: float) -> float:
    """Calculate number of orbits per day.

    Args:
        altitude_km: Altitude above Earth's surface (km)

    Returns:
        Number of orbits per day
    """
    period_min = orbital_period(altitude_km)
    return 24 * 60 / period_min


def ground_track_shift(altitude_km: float) -> float:
    """Calculate westward ground track shift per orbit.

    Due to Earth's rotation, the ground track shifts westward
    each orbit.

    Args:
        altitude_km: Altitude above Earth's surface (km)

    Returns:
        Ground track shift in degrees per orbit
    """
    period_min = orbital_period(altitude_km)
    # Earth rotates 360° in 24 hours = 0.25°/min
    earth_rotation_rate = 360 / (24 * 60)  # deg/min
    return earth_rotation_rate * period_min


def eclipse_duration(altitude_km: float, beta_angle_deg: float = 0) -> float:
    """Calculate eclipse duration per orbit.

    Args:
        altitude_km: Altitude above Earth's surface (km)
        beta_angle_deg: Solar beta angle (degrees)

    Returns:
        Eclipse duration in minutes (0 if no eclipse)
    """
    r = EARTH_RADIUS_KM + altitude_km

    # Eclipse half-angle
    sin_rho = EARTH_RADIUS_KM / r
    rho = np.arcsin(sin_rho)

    # Check if eclipse occurs at this beta angle
    beta_rad = np.radians(abs(beta_angle_deg))
    if beta_rad > (np.pi/2 - rho):
        return 0.0  # No eclipse

    # Eclipse fraction (simplified)
    eclipse_fraction = rho / np.pi * np.cos(beta_rad)

    period_min = orbital_period(altitude_km)
    return eclipse_fraction * period_min


def delta_v_hohmann(r1_km: float, r2_km: float) -> Tuple[float, float]:
    """Calculate delta-v for Hohmann transfer between circular orbits.

    Args:
        r1_km: Initial orbital radius (km)
        r2_km: Final orbital radius (km)

    Returns:
        Tuple of (delta_v1, delta_v2) in km/s
    """
    # Transfer orbit semi-major axis
    a_t = (r1_km + r2_km) / 2

    # Velocities
    v1_circular = np.sqrt(MU_EARTH / r1_km)
    v2_circular = np.sqrt(MU_EARTH / r2_km)

    v1_transfer = np.sqrt(MU_EARTH * (2/r1_km - 1/a_t))
    v2_transfer = np.sqrt(MU_EARTH * (2/r2_km - 1/a_t))

    delta_v1 = abs(v1_transfer - v1_circular)
    delta_v2 = abs(v2_circular - v2_transfer)

    return delta_v1, delta_v2


def delta_v_plane_change(velocity_km_s: float, angle_deg: float) -> float:
    """Calculate delta-v for orbital plane change.

    Args:
        velocity_km_s: Orbital velocity (km/s)
        angle_deg: Plane change angle (degrees)

    Returns:
        Delta-v required in km/s
    """
    angle_rad = np.radians(angle_deg)
    return 2 * velocity_km_s * np.sin(angle_rad / 2)


def atmospheric_drag_lifetime(
    altitude_km: float,
    ballistic_coefficient: float,
    solar_activity: str = "moderate",
) -> float:
    """Estimate orbital lifetime due to atmospheric drag.

    Simplified model for LEO satellites.

    Args:
        altitude_km: Orbital altitude (km)
        ballistic_coefficient: m/(Cd*A) in kg/m²
        solar_activity: "low", "moderate", or "high"

    Returns:
        Estimated lifetime in years
    """
    # Atmospheric density model (very simplified)
    # Real calculations would use NRLMSISE-00 or similar

    density_factors = {
        "low": 0.5,
        "moderate": 1.0,
        "high": 3.0,
    }
    factor = density_factors.get(solar_activity, 1.0)

    # Base density at reference altitudes (kg/m³)
    if altitude_km < 200:
        return 0.01  # Days, not years
    elif altitude_km < 300:
        base_lifetime = 0.1 * (altitude_km - 200) / 100
    elif altitude_km < 400:
        base_lifetime = 0.5 + 2 * (altitude_km - 300) / 100
    elif altitude_km < 500:
        base_lifetime = 2.5 + 5 * (altitude_km - 400) / 100
    elif altitude_km < 600:
        base_lifetime = 7.5 + 10 * (altitude_km - 500) / 100
    else:
        base_lifetime = 17.5 + 20 * (altitude_km - 600) / 100

    # Adjust for ballistic coefficient and solar activity
    lifetime = base_lifetime * (ballistic_coefficient / 50) / factor

    return max(0.01, lifetime)


def summary(altitude_km: float, inclination_deg: float = 0) -> dict:
    """Generate orbital summary for given altitude.

    Args:
        altitude_km: Orbital altitude (km)
        inclination_deg: Orbital inclination (degrees)

    Returns:
        Dictionary of orbital parameters
    """
    orbit = OrbitalElements.from_altitude(altitude_km, inclination_deg)

    return {
        "altitude_km": altitude_km,
        "inclination_deg": inclination_deg,
        "semi_major_axis_km": orbit.a,
        "orbital_period_min": round(orbit.period_minutes, 2),
        "velocity_km_s": round(orbit.velocity_circular, 3),
        "velocity_km_h": round(orbit.velocity_circular * 3600, 0),
        "orbits_per_day": round(orbits_per_day(altitude_km), 2),
        "ground_track_shift_deg": round(ground_track_shift(altitude_km), 2),
        "eclipse_duration_min": round(eclipse_duration(altitude_km), 2),
    }
