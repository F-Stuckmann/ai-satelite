"""Solar Panel on Heatsink Concept.

Innovative design that mounts solar cells directly on the back of radiator panels,
using a single structure for both power generation and thermal management.
"""

SOLAR_HEATSINK_CONCEPT = {
    "name": "Solar-on-Heatsink",
    "patent": None,
    "description": """
    Mount solar cells directly on radiator panel backside.

    Configuration:
    - Front (space-facing): High-emissivity radiator surface (ε=0.9)
    - Back (sun-facing): Solar cells bonded to panel

    Heat flow path:
    1. Solar cells generate power + waste heat
    2. Heat conducts through panel substrate
    3. Heat radiates to space from front surface

    Benefits:
    - Eliminates separate solar array mass/structure
    - Shorter thermal path (cells → radiator direct)
    - Potentially higher W/kg

    Challenges:
    - Solar cell efficiency drops ~0.5%/°C above 25°C
    - CTE mismatch between solar cells (3 ppm/K) and Al (23 ppm/K)
    - Need careful thermal design to limit cell temperature
    """,
    "parameters": {
        # Structure
        "chassis_length_m": 1.0,
        "chassis_width_m": 0.6,
        "chassis_height_m": 0.15,
        "wall_thickness_mm": 2.5,
        "material": "aluminum_7075",  # Stronger for thinner walls
        # Thermal (integrated solar-radiator panel)
        "radiator_type": "flat_panel",  # Flat for solar mounting
        "radiator_base_area_m2": 6.0,  # Larger area serves dual purpose
        "emissivity": 0.9,
        "panel_thickness_mm": 4.0,  # Thicker for heat spreading
        # Power (solar on radiator back)
        "solar_cell_type": "triple_junction_gaas",
        "solar_area_m2": 6.0,  # Same as radiator area
        "battery_capacity_wh": 400,
        # Payload
        "ai_chip_type": "nvidia_h100_pcie",
        "ai_chip_count": 1,
        "duty_cycle": 0.4,  # Limited by thermal
        # Orbit
        "altitude_km": 550,
        "beta_angle_deg": 0,
    },
    "analysis_required": [
        "Solar cell temperature vs efficiency curve",
        "Thermal gradient across integrated panel",
        "Structural stress from CTE mismatch",
        "Optimal panel thickness for heat spreading",
        "Attitude constraints (radiator must face space)",
    ],
    "potential_w_per_kg": "15-20 (estimated)",
    "thermal_coupling_model": {
        "cell_to_substrate_resistance": 0.1,  # K·m²/W
        "substrate_conductivity_w_per_mk": 167,  # Aluminum
        "min_panel_thickness_mm": 3.0,  # For adequate spreading
    },
    "cte_analysis": {
        "solar_cell_cte_ppm_k": 3,
        "aluminum_cte_ppm_k": 23,
        "temp_range_c": [-100, 80],  # Eclipse to full sun
        "strain_per_cycle": 0.002,  # 0.2% strain
        "mitigation": "Use compliant adhesive or Ti substrate",
    },
    "notes": [
        "Key insight: radiator and solar panel are both large flat surfaces",
        "Integration saves deployment mechanism mass",
        "Must manage solar cell temperature for efficiency",
        "Consider titanium substrate (CTE=8.6) to reduce mismatch",
    ],
}
