"""Grid fin radiator model based on SpaceX patent US 11,834,205 B1.

The patent describes a cellular grid fin structure for maximizing heat exchange
with outer space. Key features:
- Variable density grid pattern (denser at center, sparser at edges)
- Direct component-to-chassis thermal mounting
- Integrated structural/thermal design
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Literal

from .stefan_boltzmann import radiative_power, equilibrium_temperature
from ..core.constants import STEFAN_BOLTZMANN, T_SPACE, MATERIALS


@dataclass
class GridFinRadiator:
    """Model of a cellular grid fin radiator panel.

    The grid fin structure increases effective radiating area through
    fin surface area enhancement while maintaining structural rigidity.

    Attributes:
        base_area_m2: Base panel area (footprint) in m²
        fin_height_mm: Height of fins in mm
        fin_spacing_mm: Spacing between fins in mm
        fin_thickness_mm: Thickness of fins in mm
        material: Material key from MATERIALS database
        emissivity: Surface emissivity (with coating)
        density_profile: 'uniform', 'center_dense', or custom array
    """

    base_area_m2: float
    fin_height_mm: float = 20.0
    fin_spacing_mm: float = 10.0
    fin_thickness_mm: float = 1.0
    material: str = "aluminum_6061"
    emissivity: float = 0.9
    density_profile: Literal["uniform", "center_dense"] = "center_dense"

    # Computed properties
    _area_enhancement: float = field(init=False, repr=False)
    _mass_kg: float = field(init=False, repr=False)

    def __post_init__(self):
        """Calculate derived properties."""
        self._calculate_area_enhancement()
        self._calculate_mass()

    def _calculate_area_enhancement(self) -> None:
        """Calculate effective area enhancement from fin geometry.

        For grid fins, the enhancement factor depends on:
        - Fin height to spacing ratio (H/S)
        - Fin efficiency (thermal gradient along fin)
        """
        h = self.fin_height_mm
        s = self.fin_spacing_mm
        t = self.fin_thickness_mm

        # Number of fin walls per unit length (both directions for grid)
        fins_per_m = 1000 / s  # fins/m

        # Fin surface area per base area
        # Each fin wall adds 2 × height × length area
        # Grid has fins in both x and y directions
        fin_area_ratio = 4 * h / s  # Simplified model

        # Fin efficiency (accounts for temperature drop along fin)
        # η = tanh(mL) / mL where m = sqrt(2h/(k*t))
        mat = MATERIALS[self.material]
        k = mat["k"]  # Thermal conductivity
        h_conv = 0  # No convection in space, radiation only

        # For radiation-dominated fins, efficiency is higher
        # Simplified: assume 85% efficiency for typical geometries
        fin_efficiency = 0.85

        # Apply density profile factor
        if self.density_profile == "center_dense":
            # Center has 1.5x density, edges have 0.7x
            # Average enhancement is slightly higher
            profile_factor = 1.15
        else:
            profile_factor = 1.0

        self._area_enhancement = 1 + fin_area_ratio * fin_efficiency * profile_factor

    def _calculate_mass(self) -> None:
        """Calculate radiator mass including fins."""
        mat = MATERIALS[self.material]
        rho = mat["rho"]  # kg/m³

        # Base plate (assume 2mm thick)
        base_thickness_m = 0.002
        base_volume = self.base_area_m2 * base_thickness_m

        # Fin volume (simplified)
        h = self.fin_height_mm / 1000  # Convert to m
        t = self.fin_thickness_mm / 1000
        s = self.fin_spacing_mm / 1000

        # Grid fin volume per m² of base
        fin_volume_per_m2 = 2 * h * t / s  # Two directions

        total_volume = base_volume + self.base_area_m2 * fin_volume_per_m2
        self._mass_kg = total_volume * rho

    @property
    def effective_area_m2(self) -> float:
        """Effective radiating area accounting for fin enhancement."""
        return self.base_area_m2 * self._area_enhancement

    @property
    def mass_kg(self) -> float:
        """Total radiator mass in kg."""
        return self._mass_kg

    @property
    def area_to_mass_ratio(self) -> float:
        """Effective area per unit mass (m²/kg)."""
        return self.effective_area_m2 / self._mass_kg

    def heat_rejection_capacity(self, temperature_k: float) -> float:
        """Calculate heat rejection at given temperature.

        Args:
            temperature_k: Radiator temperature in Kelvin

        Returns:
            Heat rejection capacity in Watts
        """
        return radiative_power(
            temperature_k,
            self.effective_area_m2,
            self.emissivity,
            T_SPACE,
        )

    def equilibrium_temp(self, heat_input_w: float) -> float:
        """Calculate equilibrium temperature for given heat input.

        Args:
            heat_input_w: Heat input in Watts

        Returns:
            Equilibrium temperature in Kelvin
        """
        return equilibrium_temperature(
            heat_input_w,
            self.effective_area_m2,
            self.emissivity,
            T_SPACE,
        )

    def w_per_kg_limit(self, max_temp_k: float) -> float:
        """Calculate W/kg limit at maximum temperature.

        Args:
            max_temp_k: Maximum allowable temperature in K

        Returns:
            W/kg thermal limit
        """
        max_power = self.heat_rejection_capacity(max_temp_k)
        return max_power / self._mass_kg

    def summary(self) -> dict:
        """Return summary of radiator properties."""
        return {
            "base_area_m2": self.base_area_m2,
            "effective_area_m2": round(self.effective_area_m2, 3),
            "area_enhancement": round(self._area_enhancement, 2),
            "mass_kg": round(self._mass_kg, 3),
            "area_to_mass_m2_per_kg": round(self.area_to_mass_ratio, 2),
            "material": self.material,
            "emissivity": self.emissivity,
            "heat_rejection_at_40C_W": round(self.heat_rejection_capacity(273.15 + 40), 1),
            "heat_rejection_at_60C_W": round(self.heat_rejection_capacity(273.15 + 60), 1),
        }


def create_starlink_radiator() -> GridFinRadiator:
    """Create a radiator approximating SpaceX Starlink design.

    Based on patent US 11,834,205 B1:
    - 3300W power dissipation
    - Grid fin structure
    - Aluminum chassis
    """
    # Estimate: ~2.5 m² base area with grid fin enhancement
    return GridFinRadiator(
        base_area_m2=2.5,
        fin_height_mm=25,
        fin_spacing_mm=12,
        fin_thickness_mm=1.5,
        material="aluminum_6061",
        emissivity=0.9,
        density_profile="center_dense",
    )
