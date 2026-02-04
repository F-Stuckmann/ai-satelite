"""Physical constants and material/component databases."""

import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴) - Stefan-Boltzmann constant
SOLAR_FLUX = 1361.0  # W/m² - Solar constant at 1 AU
T_SPACE = 3.0  # K - Cosmic microwave background temperature
EARTH_IR_FLUX = 240.0  # W/m² - Average Earth IR emission
EARTH_ALBEDO = 0.3  # Average Earth albedo

# Earth parameters
EARTH_RADIUS_KM = 6371.0  # km
EARTH_MU = 398600.4418  # km³/s² - Gravitational parameter
J2 = 1.08263e-3  # Earth's J2 coefficient

# =============================================================================
# MATERIAL PROPERTIES DATABASE
# =============================================================================
# E = Young's modulus (Pa)
# nu = Poisson's ratio
# rho = Density (kg/m³)
# sigma_yield = Yield strength (Pa)
# alpha = Thermal expansion coefficient (1/K)
# k = Thermal conductivity (W/(m·K))
# cp = Specific heat capacity (J/(kg·K))

MATERIALS = {
    "aluminum_6061": {
        "name": "Aluminum 6061-T6",
        "E": 68.9e9,
        "nu": 0.33,
        "rho": 2700,
        "sigma_yield": 276e6,
        "alpha": 23.6e-6,
        "k": 167,
        "cp": 896,
    },
    "aluminum_7075": {
        "name": "Aluminum 7075-T6",
        "E": 71.7e9,
        "nu": 0.33,
        "rho": 2810,
        "sigma_yield": 503e6,
        "alpha": 23.4e-6,
        "k": 130,
        "cp": 960,
    },
    "titanium_6al4v": {
        "name": "Titanium 6Al-4V",
        "E": 113.8e9,
        "nu": 0.34,
        "rho": 4430,
        "sigma_yield": 880e6,
        "alpha": 8.6e-6,
        "k": 6.7,
        "cp": 526,
    },
    "cfrp_quasi_isotropic": {
        "name": "CFRP Quasi-Isotropic",
        "E": 70e9,
        "nu": 0.30,
        "rho": 1600,
        "sigma_yield": 600e6,
        "alpha": 0.5e-6,
        "k": 5,
        "cp": 1000,
    },
    "cfrp_unidirectional": {
        "name": "CFRP Unidirectional (fiber direction)",
        "E": 135e9,
        "nu": 0.30,
        "rho": 1600,
        "sigma_yield": 1500e6,
        "alpha": -0.5e-6,  # Negative CTE in fiber direction
        "k": 120,  # High in fiber direction
        "cp": 1000,
    },
    "stainless_304": {
        "name": "Stainless Steel 304",
        "E": 193e9,
        "nu": 0.29,
        "rho": 8000,
        "sigma_yield": 215e6,
        "alpha": 17.3e-6,
        "k": 16.2,
        "cp": 500,
    },
    "copper": {
        "name": "Copper (OFHC)",
        "E": 117e9,
        "nu": 0.34,
        "rho": 8940,
        "sigma_yield": 70e6,
        "alpha": 16.5e-6,
        "k": 401,
        "cp": 385,
    },
    "graphene_composite": {
        "name": "Graphene-Enhanced Composite (theoretical)",
        "E": 200e9,
        "nu": 0.25,
        "rho": 1800,
        "sigma_yield": 800e6,
        "alpha": 1e-6,
        "k": 1500,  # Enhanced thermal conductivity
        "cp": 800,
    },
}

# =============================================================================
# AI CHIP SPECIFICATIONS
# =============================================================================
# tdp = Thermal Design Power (W)
# mass_kg = Mass including heatsink mount (kg)
# tflops = Compute performance (TFLOPS)
# t_max = Maximum junction temperature (°C)

AI_CHIPS = {
    "nvidia_h100_sxm": {
        "name": "NVIDIA H100 SXM5",
        "tdp": 700,
        "mass_kg": 2.5,
        "tflops_fp16": 1979,
        "tflops_fp32": 67,
        "t_max": 83,
    },
    "nvidia_h100_pcie": {
        "name": "NVIDIA H100 PCIe",
        "tdp": 350,
        "mass_kg": 1.8,
        "tflops_fp16": 1513,
        "tflops_fp32": 51,
        "t_max": 83,
    },
    "nvidia_b200": {
        "name": "NVIDIA B200",
        "tdp": 1000,
        "mass_kg": 3.0,
        "tflops_fp16": 4500,
        "tflops_fp32": 140,
        "t_max": 85,
    },
    "google_tpu_v5e": {
        "name": "Google TPU v5e",
        "tdp": 200,
        "mass_kg": 1.5,
        "tflops_bf16": 393,
        "tflops_int8": 786,
        "t_max": 85,
    },
    "custom_space_asic": {
        "name": "Custom Space-Rated ASIC",
        "tdp": 50,
        "mass_kg": 0.3,
        "tflops_int8": 100,
        "t_max": 125,  # Rad-hard, higher temp tolerance
    },
}

# =============================================================================
# SOLAR CELL SPECIFICATIONS
# =============================================================================
# efficiency = BOL efficiency at 28°C
# degradation_rate = Annual degradation in LEO (%/year)
# mass_per_area = kg/m² including coverglass
# temp_coefficient = Efficiency change per °C (negative)

SOLAR_CELLS = {
    "triple_junction_gaas": {
        "name": "Triple-Junction GaAs",
        "efficiency": 0.30,
        "degradation_rate": 0.025,
        "mass_per_area": 0.84,
        "temp_coefficient": -0.0045,  # -0.45%/°C
        "t_ref": 28,  # Reference temperature (°C)
    },
    "quad_junction": {
        "name": "Quad-Junction IMM",
        "efficiency": 0.35,
        "degradation_rate": 0.02,
        "mass_per_area": 0.90,
        "temp_coefficient": -0.004,
        "t_ref": 28,
    },
    "thin_film_cigs": {
        "name": "Thin-Film CIGS",
        "efficiency": 0.20,
        "degradation_rate": 0.03,
        "mass_per_area": 0.50,
        "temp_coefficient": -0.003,
        "t_ref": 28,
    },
    "perovskite_tandem": {
        "name": "Perovskite Tandem (experimental)",
        "efficiency": 0.33,
        "degradation_rate": 0.05,  # Higher degradation currently
        "mass_per_area": 0.40,
        "temp_coefficient": -0.003,
        "t_ref": 28,
    },
}

# =============================================================================
# BATTERY SPECIFICATIONS
# =============================================================================
# specific_energy = Wh/kg
# cycle_life = Cycles at 80% DoD to 80% capacity
# operating_temp = (min, max) in °C

BATTERIES = {
    "li_ion_standard": {
        "name": "Li-Ion (Standard Space)",
        "specific_energy": 150,
        "cycle_life": 20000,
        "operating_temp": (-20, 45),
        "efficiency": 0.95,
    },
    "li_ion_high_energy": {
        "name": "Li-Ion (High Energy)",
        "specific_energy": 250,
        "cycle_life": 10000,
        "operating_temp": (-10, 40),
        "efficiency": 0.93,
    },
    "solid_state": {
        "name": "Solid-State Li (experimental)",
        "specific_energy": 400,
        "cycle_life": 5000,
        "operating_temp": (-30, 60),
        "efficiency": 0.96,
    },
}

# =============================================================================
# SPACEX STARLINK REFERENCE (US 11,834,205 B1)
# =============================================================================

STARLINK_REFERENCE = {
    "patent": "US 11,834,205 B1",
    "mass_kg": 261,
    "power_w": 3300,
    "w_per_kg": 12.6,
    "radiator_type": "cellular_grid_fin",
    "chassis_material": "aluminum_alloy",
    "altitude_km": 550,
    "design_life_years": 5,
}
