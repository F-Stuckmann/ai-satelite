"""Body-Mounted Solar Concept.

Simplified design with solar cells on all chassis faces - no deployable arrays.
Trades power for simplicity and reliability.
"""

BODY_MOUNTED_CONCEPT = {
    "name": "Body-Mounted Solar",
    "patent": None,
    "description": """
    Solar cells mounted directly on all chassis faces.

    Configuration:
    - All 6 faces have solar cells
    - No deployable mechanisms
    - Attitude control keeps optimal face toward sun

    Benefits:
    - No deployment risk (major failure mode)
    - Lower mass (no hinges, motors, latches)
    - More robust to debris impacts
    - Simpler thermal design

    Challenges:
    - ~60% less power than tracking array (geometry)
    - Requires precise attitude control
    - Some faces always in shadow
    """,
    "parameters": {
        # Structure (larger chassis for more solar area)
        "chassis_length_m": 1.5,
        "chassis_width_m": 1.0,
        "chassis_height_m": 0.3,
        "wall_thickness_mm": 2.0,
        "material": "cfrp_quasi_isotropic",  # Light for larger chassis
        # Thermal
        "radiator_type": "flat_panel",
        "radiator_base_area_m2": 1.5,  # Dedicated radiator on one face
        "emissivity": 0.9,
        # Power (body-mounted)
        "solar_cell_type": "thin_film_cigs",  # Lighter, flexible
        "solar_area_m2": 4.6,  # Total of all faces minus radiator
        "battery_capacity_wh": 600,  # Larger for longer eclipse
        # Payload
        "ai_chip_type": "custom_space_asic",  # Lower power
        "ai_chip_count": 2,
        "duty_cycle": 0.6,
        # Orbit
        "altitude_km": 500,
        "beta_angle_deg": 0,
    },
    "power_geometry_analysis": {
        "sun_facing_fraction": 0.4,  # Average fraction of area facing sun
        "effective_area_factor": 0.4,  # vs deployed tracking array
        "attitude_modes": ["sun_tracking", "nadir_pointing", "inertial"],
    },
    "reliability_benefits": {
        "deployment_failure_rate_avoided": 0.001,  # Per mechanism
        "mechanisms_eliminated": ["hinge", "motor", "latch", "damper"],
        "estimated_mass_savings_kg": 5,
    },
    "potential_w_per_kg": "10-15 (estimated)",
    "notes": [
        "Best for missions prioritizing reliability over power",
        "Consider for LEO where drag limits lifetime anyway",
        "Works well with low-power AI ASICs",
    ],
}
