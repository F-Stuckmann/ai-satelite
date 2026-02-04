"""Stefan-Boltzmann radiative heat transfer calculations.

The Stefan-Boltzmann law describes thermal radiation:
    P = ε × σ × A × T⁴

Where:
    P = Power radiated (W)
    ε = Emissivity (0 to 1, 1 for blackbody)
    σ = Stefan-Boltzmann constant = 5.67×10⁻⁸ W/(m²·K⁴)
    A = Surface area (m²)
    T = Absolute temperature (K)
"""

import numpy as np
from ..core.constants import STEFAN_BOLTZMANN, T_SPACE


def radiative_power(
    temperature_k: float,
    area_m2: float,
    emissivity: float = 0.9,
    t_sink_k: float = T_SPACE,
) -> float:
    """Calculate net radiative heat transfer.

    Args:
        temperature_k: Surface temperature in Kelvin
        area_m2: Radiating surface area in m²
        emissivity: Surface emissivity (0-1), default 0.9 for black paint
        t_sink_k: Sink temperature (default: space at 3K)

    Returns:
        Net power radiated in Watts (positive = heat rejection)
    """
    if temperature_k <= 0:
        raise ValueError("Temperature must be positive (Kelvin)")
    if area_m2 <= 0:
        raise ValueError("Area must be positive")
    if not 0 < emissivity <= 1:
        raise ValueError("Emissivity must be between 0 and 1")

    q_out = emissivity * STEFAN_BOLTZMANN * area_m2 * temperature_k**4
    q_in = emissivity * STEFAN_BOLTZMANN * area_m2 * t_sink_k**4

    return q_out - q_in


def equilibrium_temperature(
    heat_input_w: float,
    area_m2: float,
    emissivity: float = 0.9,
    t_sink_k: float = T_SPACE,
) -> float:
    """Calculate equilibrium temperature for given heat input.

    At equilibrium: Q_in = Q_radiated
    Therefore: T = (Q_in / (ε × σ × A) + T_sink⁴)^(1/4)

    Args:
        heat_input_w: Total heat input in Watts
        area_m2: Radiating surface area in m²
        emissivity: Surface emissivity (0-1)
        t_sink_k: Sink temperature (default: space at 3K)

    Returns:
        Equilibrium temperature in Kelvin
    """
    if heat_input_w < 0:
        raise ValueError("Heat input cannot be negative")
    if area_m2 <= 0:
        raise ValueError("Area must be positive")
    if not 0 < emissivity <= 1:
        raise ValueError("Emissivity must be between 0 and 1")

    # T⁴ = Q / (ε × σ × A) + T_sink⁴
    t4 = heat_input_w / (emissivity * STEFAN_BOLTZMANN * area_m2) + t_sink_k**4

    return t4**0.25


def required_radiator_area(
    heat_rejection_w: float,
    temperature_k: float,
    emissivity: float = 0.9,
    t_sink_k: float = T_SPACE,
) -> float:
    """Calculate required radiator area for given heat rejection.

    Rearranging Stefan-Boltzmann:
    A = Q / (ε × σ × (T⁴ - T_sink⁴))

    Args:
        heat_rejection_w: Required heat rejection in Watts
        temperature_k: Radiator temperature in Kelvin
        emissivity: Surface emissivity (0-1)
        t_sink_k: Sink temperature (default: space at 3K)

    Returns:
        Required radiator area in m²
    """
    if heat_rejection_w < 0:
        raise ValueError("Heat rejection cannot be negative")
    if temperature_k <= t_sink_k:
        raise ValueError("Radiator temperature must be above sink temperature")
    if not 0 < emissivity <= 1:
        raise ValueError("Emissivity must be between 0 and 1")

    delta_t4 = temperature_k**4 - t_sink_k**4
    area = heat_rejection_w / (emissivity * STEFAN_BOLTZMANN * delta_t4)

    return area


def temperature_from_celsius(temp_c: float) -> float:
    """Convert Celsius to Kelvin."""
    return temp_c + 273.15


def temperature_to_celsius(temp_k: float) -> float:
    """Convert Kelvin to Celsius."""
    return temp_k - 273.15


def w_per_kg_thermal_limit(
    radiator_area_m2: float,
    radiator_mass_kg: float,
    max_temp_k: float,
    emissivity: float = 0.9,
) -> float:
    """Calculate W/kg limit based on thermal rejection capability.

    This represents the maximum power that can be dissipated per kg
    of radiator mass.

    Args:
        radiator_area_m2: Radiator area in m²
        radiator_mass_kg: Radiator mass in kg
        max_temp_k: Maximum allowable radiator temperature in K
        emissivity: Surface emissivity

    Returns:
        Maximum W/kg for thermal system
    """
    max_power = radiative_power(max_temp_k, radiator_area_m2, emissivity)
    return max_power / radiator_mass_kg


# Convenience: typical values
def verify_spacex_claim() -> dict:
    """Verify SpaceX patent claim of 12.6 W/kg.

    Patent US 11,834,205 B1 claims:
    - Mass: 261 kg
    - Power: 3300 W
    - W/kg: 12.6

    Returns:
        Dictionary with verification results
    """
    mass_kg = 261
    power_w = 3300
    claimed_w_per_kg = 12.6

    # Calculate actual W/kg
    actual_w_per_kg = power_w / mass_kg

    # Estimate required radiator area at typical operating temp (40°C)
    t_radiator = temperature_from_celsius(40)
    required_area = required_radiator_area(power_w, t_radiator)

    # At higher temp (60°C), less area needed
    t_radiator_hot = temperature_from_celsius(60)
    required_area_hot = required_radiator_area(power_w, t_radiator_hot)

    return {
        "patent": "US 11,834,205 B1",
        "claimed_w_per_kg": claimed_w_per_kg,
        "calculated_w_per_kg": round(actual_w_per_kg, 2),
        "verified": abs(actual_w_per_kg - claimed_w_per_kg) < 0.1,
        "required_radiator_area_40C": round(required_area, 2),
        "required_radiator_area_60C": round(required_area_hot, 2),
    }
