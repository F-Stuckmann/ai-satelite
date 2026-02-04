"""Material properties and structural calculations for satellite design."""

import numpy as np
from dataclasses import dataclass
from typing import Optional

from ..core.constants import MATERIALS


@dataclass
class MaterialProperties:
    """Material properties for structural analysis.

    Attributes:
        name: Material name
        E: Young's modulus (Pa)
        nu: Poisson's ratio
        rho: Density (kg/m³)
        sigma_yield: Yield strength (Pa)
        sigma_ultimate: Ultimate strength (Pa)
        alpha: Thermal expansion coefficient (1/K)
        k: Thermal conductivity (W/(m·K))
        cp: Specific heat capacity (J/(kg·K))
    """
    name: str
    E: float  # Young's modulus (Pa)
    nu: float  # Poisson's ratio
    rho: float  # Density (kg/m³)
    sigma_yield: float  # Yield strength (Pa)
    sigma_ultimate: Optional[float] = None  # Ultimate strength (Pa)
    alpha: float = 0.0  # Thermal expansion coefficient (1/K)
    k: float = 0.0  # Thermal conductivity (W/(m·K))
    cp: float = 0.0  # Specific heat capacity (J/(kg·K))

    @classmethod
    def from_database(cls, material_key: str) -> "MaterialProperties":
        """Load material from constants database.

        Args:
            material_key: Key in MATERIALS dictionary

        Returns:
            MaterialProperties instance
        """
        if material_key not in MATERIALS:
            raise ValueError(f"Unknown material: {material_key}")

        mat = MATERIALS[material_key]
        return cls(
            name=mat["name"],
            E=mat["E"],
            nu=mat["nu"],
            rho=mat["rho"],
            sigma_yield=mat["sigma_yield"],
            alpha=mat.get("alpha", 0),
            k=mat.get("k", 0),
            cp=mat.get("cp", 0),
        )

    @property
    def G(self) -> float:
        """Shear modulus (Pa)."""
        return self.E / (2 * (1 + self.nu))

    @property
    def K(self) -> float:
        """Bulk modulus (Pa)."""
        return self.E / (3 * (1 - 2 * self.nu))


def thermal_stress(
    material: MaterialProperties,
    delta_T: float,
    constrained: bool = True,
) -> float:
    """Calculate thermal stress for temperature change.

    Args:
        material: Material properties
        delta_T: Temperature change (K or °C)
        constrained: If True, fully constrained (max stress)

    Returns:
        Thermal stress in Pa
    """
    if not constrained:
        return 0.0  # Free expansion, no stress

    # σ = E × α × ΔT
    return material.E * material.alpha * abs(delta_T)


def thermal_strain(material: MaterialProperties, delta_T: float) -> float:
    """Calculate thermal strain for temperature change.

    Args:
        material: Material properties
        delta_T: Temperature change (K or °C)

    Returns:
        Thermal strain (dimensionless)
    """
    return material.alpha * delta_T


def safety_factor(
    applied_stress: float,
    material: MaterialProperties,
    use_ultimate: bool = False,
) -> float:
    """Calculate safety factor.

    Args:
        applied_stress: Applied stress (Pa)
        material: Material properties
        use_ultimate: Use ultimate strength instead of yield

    Returns:
        Safety factor (>1 is safe)
    """
    if applied_stress <= 0:
        return float('inf')

    allowable = material.sigma_ultimate if use_ultimate and material.sigma_ultimate else material.sigma_yield
    return allowable / applied_stress


def mass_of_panel(
    length_m: float,
    width_m: float,
    thickness_mm: float,
    material: MaterialProperties,
) -> float:
    """Calculate mass of a rectangular panel.

    Args:
        length_m: Panel length (m)
        width_m: Panel width (m)
        thickness_mm: Panel thickness (mm)
        material: Material properties

    Returns:
        Mass in kg
    """
    volume = length_m * width_m * (thickness_mm / 1000)
    return volume * material.rho


def mass_of_hollow_cylinder(
    outer_radius_m: float,
    inner_radius_m: float,
    length_m: float,
    material: MaterialProperties,
) -> float:
    """Calculate mass of a hollow cylinder (tube).

    Args:
        outer_radius_m: Outer radius (m)
        inner_radius_m: Inner radius (m)
        length_m: Length (m)
        material: Material properties

    Returns:
        Mass in kg
    """
    area = np.pi * (outer_radius_m**2 - inner_radius_m**2)
    volume = area * length_m
    return volume * material.rho


def buckling_stress_panel(
    length_m: float,
    width_m: float,
    thickness_mm: float,
    material: MaterialProperties,
    boundary: str = "simply_supported",
) -> float:
    """Calculate critical buckling stress for a flat panel.

    Args:
        length_m: Panel length (m)
        width_m: Panel width (m)
        thickness_mm: Panel thickness (mm)
        material: Material properties
        boundary: "simply_supported" or "clamped"

    Returns:
        Critical buckling stress in Pa
    """
    t = thickness_mm / 1000  # Convert to m
    a = max(length_m, width_m)  # Longer dimension
    b = min(length_m, width_m)  # Shorter dimension

    # Buckling coefficient
    k = 4.0 if boundary == "simply_supported" else 6.97

    # Flexural rigidity
    D = material.E * t**3 / (12 * (1 - material.nu**2))

    # Critical stress
    sigma_cr = k * np.pi**2 * D / (b**2 * t)

    return sigma_cr


def natural_frequency_panel(
    length_m: float,
    width_m: float,
    thickness_mm: float,
    material: MaterialProperties,
    boundary: str = "simply_supported",
) -> float:
    """Calculate fundamental natural frequency of a flat panel.

    Args:
        length_m: Panel length (m)
        width_m: Panel width (m)
        thickness_mm: Panel thickness (mm)
        material: Material properties
        boundary: "simply_supported" or "clamped"

    Returns:
        Natural frequency in Hz
    """
    t = thickness_mm / 1000  # Convert to m
    a = length_m
    b = width_m

    # Frequency coefficient (depends on aspect ratio and boundary)
    if boundary == "simply_supported":
        lambda_sq = np.pi**2 * ((1/a**2) + (1/b**2))
    else:  # clamped
        lambda_sq = 22.4 / (a**2)  # Approximate for square plate

    # Flexural rigidity per unit width
    D = material.E * t**3 / (12 * (1 - material.nu**2))

    # Mass per unit area
    m = material.rho * t

    # Natural frequency
    omega = np.sqrt(D * lambda_sq / m)
    f = omega / (2 * np.pi)

    return f


def cte_mismatch_stress(
    material1: MaterialProperties,
    material2: MaterialProperties,
    delta_T: float,
    thickness_ratio: float = 1.0,
) -> tuple:
    """Calculate interfacial stress from CTE mismatch between bonded materials.

    Args:
        material1: First material
        material2: Second material
        delta_T: Temperature change (K)
        thickness_ratio: t2/t1 thickness ratio

    Returns:
        Tuple of (stress_in_material1, stress_in_material2) in Pa
    """
    alpha1 = material1.alpha
    alpha2 = material2.alpha
    E1 = material1.E
    E2 = material2.E

    delta_alpha = abs(alpha2 - alpha1)
    delta_epsilon = delta_alpha * abs(delta_T)

    # Simplified: assume equal thickness
    # More accurate would account for thickness ratio
    sigma1 = E1 * delta_epsilon / 2
    sigma2 = E2 * delta_epsilon / 2

    return sigma1, sigma2


def structural_summary(
    chassis_dims: tuple,
    wall_thickness_mm: float,
    material_key: str,
    temp_range: tuple = (-100, 80),
) -> dict:
    """Generate structural summary for a chassis.

    Args:
        chassis_dims: (length, width, height) in meters
        wall_thickness_mm: Wall thickness in mm
        material_key: Material key from database
        temp_range: (min_temp, max_temp) in Celsius

    Returns:
        Dictionary of structural properties
    """
    material = MaterialProperties.from_database(material_key)
    length, width, height = chassis_dims

    # Surface area and mass
    surface_area = 2 * (length*width + length*height + width*height)
    volume = surface_area * (wall_thickness_mm / 1000)
    mass = volume * material.rho

    # Thermal stress from temperature cycling
    delta_T = temp_range[1] - temp_range[0]
    thermal_stress_val = thermal_stress(material, delta_T)
    thermal_sf = safety_factor(thermal_stress_val, material)

    # Natural frequency (simplified - treat as panel)
    nat_freq = natural_frequency_panel(
        length, width, wall_thickness_mm, material, "simply_supported"
    )

    # Buckling stress
    buckling = buckling_stress_panel(
        length, width, wall_thickness_mm, material, "simply_supported"
    )

    return {
        "material": material.name,
        "mass_kg": round(mass, 3),
        "surface_area_m2": round(surface_area, 3),
        "thermal_stress_MPa": round(thermal_stress_val / 1e6, 2),
        "thermal_safety_factor": round(thermal_sf, 2),
        "natural_frequency_Hz": round(nat_freq, 1),
        "buckling_stress_MPa": round(buckling / 1e6, 1),
        "yield_strength_MPa": round(material.sigma_yield / 1e6, 1),
    }
