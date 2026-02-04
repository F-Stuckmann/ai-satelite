"""View factor calculations for orbital thermal environment.

View factors determine how much thermal radiation a surface receives from
or emits to various sources (Earth, Sun, deep space).

Key thermal inputs for LEO satellites:
- Direct solar flux: ~1361 W/m² when illuminated
- Earth IR: ~240 W/m² (varies with altitude and position)
- Earth albedo: Reflected sunlight, ~30% of incident solar
- Deep space: 3K sink temperature
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional

from ..core.constants import SOLAR_FLUX, EARTH_IR_FLUX, EARTH_ALBEDO, EARTH_RADIUS_KM


@dataclass
class ViewFactors:
    """Calculate view factors and thermal environment for a satellite.

    Attributes:
        altitude_km: Orbital altitude in km
        beta_angle_deg: Solar beta angle (angle between orbital plane and sun vector)
    """

    altitude_km: float
    beta_angle_deg: float = 0.0

    def __post_init__(self):
        """Validate inputs."""
        if self.altitude_km < 100:
            raise ValueError("Altitude must be at least 100 km (LEO)")
        if self.altitude_km > 50000:
            raise ValueError("Altitude above 50000 km not supported")

    @property
    def orbital_radius_km(self) -> float:
        """Distance from Earth center."""
        return EARTH_RADIUS_KM + self.altitude_km

    @property
    def earth_view_factor(self) -> float:
        """View factor to Earth (fraction of hemisphere occupied by Earth).

        For a flat plate facing nadir:
        F = sin²(ρ) where ρ = arcsin(R_earth / r_orbit)
        """
        sin_rho = EARTH_RADIUS_KM / self.orbital_radius_km
        return sin_rho**2

    @property
    def space_view_factor(self) -> float:
        """View factor to deep space (complement of Earth view)."""
        return 1.0 - self.earth_view_factor

    @property
    def orbital_period_minutes(self) -> float:
        """Orbital period using Kepler's third law."""
        # T = 2π × sqrt(a³/μ)
        # μ_earth = 398600 km³/s²
        a = self.orbital_radius_km  # km
        mu = 398600.4418  # km³/s²
        t_seconds = 2 * np.pi * np.sqrt(a**3 / mu)
        return t_seconds / 60

    @property
    def eclipse_fraction(self) -> float:
        """Fraction of orbit spent in eclipse.

        Simplified model assuming circular orbit and beta=0.
        More accurate: depends on beta angle.
        """
        # Eclipse half-angle
        sin_rho = EARTH_RADIUS_KM / self.orbital_radius_km
        rho = np.arcsin(sin_rho)

        # For beta = 0, eclipse fraction is 2*rho / (2*pi)
        # For higher beta angles, eclipse fraction decreases
        beta_rad = np.radians(abs(self.beta_angle_deg))

        if beta_rad > (np.pi / 2 - rho):
            # No eclipse (high beta angle)
            return 0.0

        # Approximate eclipse fraction accounting for beta
        base_eclipse = rho / np.pi  # Fraction at beta=0
        beta_factor = np.cos(beta_rad)

        return base_eclipse * beta_factor

    @property
    def sunlight_fraction(self) -> float:
        """Fraction of orbit in sunlight."""
        return 1.0 - self.eclipse_fraction

    def solar_flux_on_surface(
        self,
        absorptivity: float = 0.3,
        sun_angle_deg: float = 0.0,
    ) -> float:
        """Calculate absorbed solar flux on a surface.

        Args:
            absorptivity: Solar absorptivity (0-1)
            sun_angle_deg: Angle between surface normal and sun vector

        Returns:
            Absorbed solar flux in W/m² (when illuminated)
        """
        cos_angle = np.cos(np.radians(sun_angle_deg))
        if cos_angle < 0:
            return 0.0  # Surface facing away from sun

        return SOLAR_FLUX * absorptivity * cos_angle

    def earth_ir_flux_on_surface(
        self,
        emissivity: float = 0.9,
        nadir_angle_deg: float = 0.0,
    ) -> float:
        """Calculate absorbed Earth IR flux on a surface.

        Args:
            emissivity: IR emissivity (= absorptivity for thermal equilibrium)
            nadir_angle_deg: Angle between surface normal and nadir

        Returns:
            Absorbed Earth IR in W/m²
        """
        cos_angle = np.cos(np.radians(nadir_angle_deg))
        if cos_angle < 0:
            return 0.0  # Surface facing away from Earth

        # IR flux decreases with Earth view factor
        return EARTH_IR_FLUX * self.earth_view_factor * emissivity * cos_angle

    def albedo_flux_on_surface(
        self,
        absorptivity: float = 0.3,
        nadir_angle_deg: float = 0.0,
    ) -> float:
        """Calculate absorbed Earth albedo flux on a surface.

        Albedo is reflected sunlight - only present when Earth-facing
        side is illuminated.

        Args:
            absorptivity: Solar absorptivity (0-1)
            nadir_angle_deg: Angle between surface normal and nadir

        Returns:
            Absorbed albedo in W/m² (average over illuminated arc)
        """
        cos_angle = np.cos(np.radians(nadir_angle_deg))
        if cos_angle < 0:
            return 0.0

        # Albedo flux (simplified - assume 50% of orbit has albedo contribution)
        albedo_flux = SOLAR_FLUX * EARTH_ALBEDO * self.earth_view_factor * 0.5
        return albedo_flux * absorptivity * cos_angle

    def total_absorbed_heat(
        self,
        area_m2: float,
        solar_absorptivity: float = 0.3,
        ir_emissivity: float = 0.9,
        in_sunlight: bool = True,
    ) -> float:
        """Calculate total absorbed environmental heat.

        Args:
            area_m2: Surface area in m²
            solar_absorptivity: Solar absorptivity (0-1)
            ir_emissivity: IR emissivity (0-1)
            in_sunlight: Whether currently in sunlight

        Returns:
            Total absorbed heat in Watts
        """
        # Earth IR (always present)
        q_ir = self.earth_ir_flux_on_surface(ir_emissivity) * area_m2

        if in_sunlight:
            # Direct solar
            q_solar = self.solar_flux_on_surface(solar_absorptivity) * area_m2
            # Albedo
            q_albedo = self.albedo_flux_on_surface(solar_absorptivity) * area_m2
        else:
            q_solar = 0.0
            q_albedo = 0.0

        return q_ir + q_solar + q_albedo

    def orbit_average_absorbed_heat(
        self,
        area_m2: float,
        solar_absorptivity: float = 0.3,
        ir_emissivity: float = 0.9,
    ) -> float:
        """Calculate orbit-average absorbed heat.

        Args:
            area_m2: Surface area in m²
            solar_absorptivity: Solar absorptivity (0-1)
            ir_emissivity: IR emissivity (0-1)

        Returns:
            Orbit-average absorbed heat in Watts
        """
        q_sunlight = self.total_absorbed_heat(
            area_m2, solar_absorptivity, ir_emissivity, in_sunlight=True
        )
        q_eclipse = self.total_absorbed_heat(
            area_m2, solar_absorptivity, ir_emissivity, in_sunlight=False
        )

        return (
            q_sunlight * self.sunlight_fraction
            + q_eclipse * self.eclipse_fraction
        )

    def summary(self) -> dict:
        """Return summary of orbital thermal environment."""
        return {
            "altitude_km": self.altitude_km,
            "beta_angle_deg": self.beta_angle_deg,
            "orbital_period_min": round(self.orbital_period_minutes, 1),
            "earth_view_factor": round(self.earth_view_factor, 3),
            "space_view_factor": round(self.space_view_factor, 3),
            "eclipse_fraction": round(self.eclipse_fraction, 3),
            "sunlight_fraction": round(self.sunlight_fraction, 3),
            "eclipse_duration_min": round(
                self.eclipse_fraction * self.orbital_period_minutes, 1
            ),
        }


def create_starlink_environment() -> ViewFactors:
    """Create view factors for typical Starlink orbit."""
    return ViewFactors(altitude_km=550, beta_angle_deg=0)
