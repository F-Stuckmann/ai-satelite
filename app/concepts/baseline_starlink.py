"""SpaceX Starlink Baseline Design - US 11,834,205 B1.

Reference implementation based on the SpaceX spacecraft chassis patent.
"""

STARLINK_BASELINE = {
    "name": "SpaceX Starlink (Patent Baseline)",
    "patent": "US 11,834,205 B1",
    "description": """
    Reference design based on SpaceX patent claims.

    Key innovations:
    - Cellular grid fin radiator structure
    - Direct component-to-chassis thermal mounting
    - Central placement of high-heat components
    - Integrated structural/thermal design

    This represents the current state-of-the-art at 12.6 W/kg.
    """,
    "parameters": {
        # Structure
        "chassis_length_m": 1.2,
        "chassis_width_m": 0.8,
        "chassis_height_m": 0.2,
        "wall_thickness_mm": 3.0,
        "material": "aluminum_6061",
        # Thermal
        "radiator_type": "grid_fin",
        "radiator_base_area_m2": 2.5,
        "emissivity": 0.9,
        # Power
        "solar_cell_type": "triple_junction_gaas",
        "solar_area_m2": 8.0,
        "battery_capacity_wh": 500,
        # Payload
        "ai_chip_type": None,  # Starlink doesn't have AI compute
        "ai_chip_count": 0,
        # Orbit
        "altitude_km": 550,
        "beta_angle_deg": 0,
    },
    "claimed_metrics": {
        "mass_kg": 261,
        "power_w": 3300,
        "w_per_kg": 12.6,
    },
    "notes": [
        "Grid fin radiator provides ~2.5x area enhancement",
        "Aluminum chassis acts as heat spreader",
        "High-heat components mounted at center for best thermal coupling",
    ],
}
