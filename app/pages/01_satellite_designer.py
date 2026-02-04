"""Satellite Designer - Interactive design interface with parameter sliders."""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

st.set_page_config(page_title="Satellite Designer", page_icon="🎛️", layout="wide")

# Import simulation modules
import sys
sys.path.insert(0, ".")

from src.core.constants import MATERIALS, AI_CHIPS, SOLAR_CELLS, BATTERIES, STARLINK_REFERENCE
from src.thermal.stefan_boltzmann import (
    required_radiator_area,
    equilibrium_temperature,
    radiative_power,
    temperature_from_celsius,
    temperature_to_celsius,
)
from src.thermal.grid_fin_radiator import GridFinRadiator
from src.thermal.view_factors import ViewFactors
from src.thermal.heat_sources import AIChipHeatSource


def create_satellite_3d(
    chassis_size: tuple,
    radiator_area: float,
    solar_area: float,
    temp_c: float,
) -> go.Figure:
    """Create 3D satellite visualization with thermal coloring."""

    # Chassis dimensions
    lx, ly, lz = chassis_size

    # Vertices of chassis box
    vertices = np.array([
        [0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0],
        [0, 0, lz], [lx, 0, lz], [lx, ly, lz], [0, ly, lz],
    ])

    # Center the chassis
    vertices -= np.array([lx/2, ly/2, lz/2])

    # Faces (indices into vertices)
    faces = [
        [0, 1, 5, 4],  # Front
        [2, 3, 7, 6],  # Back
        [0, 3, 7, 4],  # Left
        [1, 2, 6, 5],  # Right
        [0, 1, 2, 3],  # Bottom
        [4, 5, 6, 7],  # Top (radiator)
    ]

    # Temperature-based color (blue=cold, red=hot)
    temp_normalized = (temp_c - 20) / 80  # Normalize 20-100°C to 0-1
    temp_normalized = max(0, min(1, temp_normalized))

    # Color gradient: blue (20°C) -> yellow (60°C) -> red (100°C)
    if temp_normalized < 0.5:
        r = int(255 * temp_normalized * 2)
        g = int(255 * temp_normalized * 2)
        b = int(255 * (1 - temp_normalized * 2))
    else:
        r = 255
        g = int(255 * (1 - (temp_normalized - 0.5) * 2))
        b = 0

    chassis_color = f"rgb({r},{g},{b})"

    fig = go.Figure()

    # Add chassis as mesh
    x_chassis = []
    y_chassis = []
    z_chassis = []
    i_faces = []
    j_faces = []
    k_faces = []

    for idx, face in enumerate(faces):
        # Triangulate quad face
        v0, v1, v2, v3 = face
        i_faces.extend([v0, v0])
        j_faces.extend([v1, v2])
        k_faces.extend([v2, v3])

    fig.add_trace(go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=i_faces,
        j=j_faces,
        k=k_faces,
        color=chassis_color,
        opacity=0.9,
        name="Chassis",
        hoverinfo="name",
    ))

    # Add solar panels (deployed on sides)
    solar_panel_width = np.sqrt(solar_area / 2)  # Assume 2:1 aspect ratio
    solar_panel_length = solar_panel_width * 2

    # Left solar panel
    panel_offset = lx/2 + 0.1
    fig.add_trace(go.Mesh3d(
        x=[-panel_offset - solar_panel_width, -panel_offset, -panel_offset, -panel_offset - solar_panel_width],
        y=[-solar_panel_length/2, -solar_panel_length/2, solar_panel_length/2, solar_panel_length/2],
        z=[0, 0, 0, 0],
        i=[0, 0],
        j=[1, 2],
        k=[2, 3],
        color="darkblue",
        opacity=0.8,
        name="Solar Panel L",
    ))

    # Right solar panel
    fig.add_trace(go.Mesh3d(
        x=[panel_offset, panel_offset + solar_panel_width, panel_offset + solar_panel_width, panel_offset],
        y=[-solar_panel_length/2, -solar_panel_length/2, solar_panel_length/2, solar_panel_length/2],
        z=[0, 0, 0, 0],
        i=[0, 0],
        j=[1, 2],
        k=[2, 3],
        color="darkblue",
        opacity=0.8,
        name="Solar Panel R",
    ))

    # Add radiator fins (on top)
    fin_height = 0.05
    for i in range(5):
        y_pos = -ly/2 + ly * i / 4
        fig.add_trace(go.Scatter3d(
            x=[-lx/2, lx/2],
            y=[y_pos, y_pos],
            z=[lz/2 + fin_height, lz/2 + fin_height],
            mode="lines",
            line=dict(color="gray", width=3),
            showlegend=False,
        ))

    fig.update_layout(
        scene=dict(
            xaxis_title="X (m)",
            yaxis_title="Y (m)",
            zaxis_title="Z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.0)
            ),
        ),
        title=f"Satellite 3D View (Radiator Temp: {temp_c:.1f}°C)",
        height=500,
        margin=dict(l=0, r=0, t=40, b=0),
    )

    return fig


def main():
    st.title("🎛️ Satellite Designer")
    st.markdown("**Interactive satellite design with real-time W/kg calculation**")

    # Reference metric
    st.info(f"🎯 Target: SpaceX Patent claims **{STARLINK_REFERENCE['w_per_kg']} W/kg**")

    # Create columns for layout
    col_sliders, col_results = st.columns([1, 1])

    with col_sliders:
        st.subheader("Design Parameters")

        # === STRUCTURE ===
        with st.expander("🏗️ Structure", expanded=True):
            material = st.selectbox(
                "Chassis Material",
                options=list(MATERIALS.keys()),
                format_func=lambda x: MATERIALS[x]["name"],
                index=0,
            )

            col1, col2 = st.columns(2)
            with col1:
                chassis_length = st.slider("Length (m)", 0.5, 3.0, 1.2, 0.1)
                chassis_width = st.slider("Width (m)", 0.3, 2.0, 0.8, 0.1)
            with col2:
                chassis_height = st.slider("Height (m)", 0.1, 0.5, 0.2, 0.05)
                wall_thickness = st.slider("Wall (mm)", 1.0, 10.0, 3.0, 0.5)

            # Calculate chassis mass
            mat = MATERIALS[material]
            surface_area = 2 * (
                chassis_length * chassis_width +
                chassis_length * chassis_height +
                chassis_width * chassis_height
            )
            chassis_volume = surface_area * (wall_thickness / 1000)
            chassis_mass = chassis_volume * mat["rho"]

        # === THERMAL ===
        with st.expander("🌡️ Thermal System", expanded=True):
            radiator_type = st.selectbox(
                "Radiator Type",
                ["flat_panel", "grid_fin", "deployable"],
                format_func=lambda x: {"flat_panel": "Flat Panel", "grid_fin": "Grid Fin (SpaceX)", "deployable": "Deployable"}[x],
                index=1,
            )

            col1, col2 = st.columns(2)
            with col1:
                radiator_base_area = st.slider("Radiator Area (m²)", 0.5, 10.0, 2.5, 0.1)
            with col2:
                emissivity = st.slider("Emissivity", 0.7, 0.98, 0.9, 0.01)

            # Create radiator model
            if radiator_type == "grid_fin":
                radiator = GridFinRadiator(
                    base_area_m2=radiator_base_area,
                    material=material,
                    emissivity=emissivity,
                )
                radiator_effective_area = radiator.effective_area_m2
                radiator_mass = radiator.mass_kg
            else:
                # Simple flat panel
                radiator_effective_area = radiator_base_area
                panel_thickness = 0.003  # 3mm
                radiator_mass = radiator_base_area * panel_thickness * mat["rho"]

        # === POWER ===
        with st.expander("⚡ Power System", expanded=True):
            solar_cell_type = st.selectbox(
                "Solar Cell Type",
                options=list(SOLAR_CELLS.keys()),
                format_func=lambda x: SOLAR_CELLS[x]["name"],
                index=0,
            )

            col1, col2 = st.columns(2)
            with col1:
                solar_area = st.slider("Solar Array Area (m²)", 1.0, 30.0, 8.0, 0.5)
            with col2:
                battery_capacity = st.slider("Battery (Wh)", 100, 3000, 500, 50)

            solar_cell = SOLAR_CELLS[solar_cell_type]
            solar_power = solar_area * 1361 * solar_cell["efficiency"]
            solar_mass = solar_area * solar_cell["mass_per_area"]

            battery_type = "li_ion_standard"
            battery = BATTERIES[battery_type]
            battery_mass = battery_capacity / battery["specific_energy"]

        # === AI PAYLOAD ===
        with st.expander("🧠 AI Payload", expanded=True):
            chip_type = st.selectbox(
                "AI Chip Type",
                options=list(AI_CHIPS.keys()),
                format_func=lambda x: AI_CHIPS[x]["name"],
                index=0,
            )

            col1, col2 = st.columns(2)
            with col1:
                chip_count = st.slider("Number of Chips", 0, 8, 1, 1)
            with col2:
                duty_cycle = st.slider("Duty Cycle (%)", 10, 100, 50, 5) / 100

            chip = AI_CHIPS[chip_type]
            ai_power = chip["tdp"] * chip_count * duty_cycle
            ai_mass = chip["mass_kg"] * chip_count

        # === ORBIT ===
        with st.expander("🌍 Orbit", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                altitude = st.slider("Altitude (km)", 300, 1200, 550, 10)
            with col2:
                beta_angle = st.slider("Beta Angle (°)", 0, 90, 0, 5)

            orbit = ViewFactors(altitude_km=altitude, beta_angle_deg=beta_angle)

    # === CALCULATIONS ===
    # Total mass
    other_mass = 20  # Avionics, propulsion, etc.
    total_mass = chassis_mass + radiator_mass + solar_mass + battery_mass + ai_mass + other_mass

    # Total power generation (orbit average)
    orbit_avg_solar = solar_power * orbit.sunlight_fraction
    power_available = orbit_avg_solar * 0.9  # 90% efficiency to loads

    # Heat to dissipate
    heat_dissipation = ai_power + power_available * 0.05  # AI + 5% losses

    # Equilibrium temperature
    t_eq_k = equilibrium_temperature(heat_dissipation, radiator_effective_area, emissivity)
    t_eq_c = temperature_to_celsius(t_eq_k)

    # W/kg calculation
    w_per_kg = power_available / total_mass

    # === RESULTS ===
    with col_results:
        st.subheader("Results")

        # Primary metric: W/kg
        delta = w_per_kg - STARLINK_REFERENCE["w_per_kg"]
        st.metric(
            label="⚡ W/kg (Power-to-Mass Ratio)",
            value=f"{w_per_kg:.2f}",
            delta=f"{delta:+.2f} vs SpaceX",
            delta_color="normal" if delta >= 0 else "inverse",
        )

        # Status indicators
        col1, col2, col3 = st.columns(3)
        with col1:
            # Thermal status
            t_max = chip.get("t_max", 85)
            thermal_margin = t_max - t_eq_c
            if thermal_margin > 15:
                st.success(f"🌡️ Thermal: {t_eq_c:.1f}°C (+{thermal_margin:.0f}°C margin)")
            elif thermal_margin > 0:
                st.warning(f"🌡️ Thermal: {t_eq_c:.1f}°C (+{thermal_margin:.0f}°C margin)")
            else:
                st.error(f"🌡️ OVER TEMP: {t_eq_c:.1f}°C ({thermal_margin:.0f}°C)")

        with col2:
            # Power balance
            power_margin = power_available - ai_power
            if power_margin > 0:
                st.success(f"⚡ Power: +{power_margin:.0f}W margin")
            else:
                st.error(f"⚡ Power: {power_margin:.0f}W deficit")

        with col3:
            # Mass breakdown
            st.info(f"📦 Mass: {total_mass:.1f} kg")

        # 3D Visualization
        st.plotly_chart(
            create_satellite_3d(
                (chassis_length, chassis_width, chassis_height),
                radiator_effective_area,
                solar_area,
                t_eq_c,
            ),
            use_container_width=True,
        )

        # Mass breakdown pie chart
        st.subheader("Mass Breakdown")
        mass_data = {
            "Chassis": chassis_mass,
            "Radiator": radiator_mass,
            "Solar Array": solar_mass,
            "Battery": battery_mass,
            "AI Payload": ai_mass,
            "Other": other_mass,
        }

        fig_pie = go.Figure(data=[go.Pie(
            labels=list(mass_data.keys()),
            values=list(mass_data.values()),
            hole=0.4,
            textinfo="label+percent",
        )])
        fig_pie.update_layout(
            title=f"Total Mass: {total_mass:.1f} kg",
            height=300,
            margin=dict(l=20, r=20, t=40, b=20),
        )
        st.plotly_chart(fig_pie, use_container_width=True)

        # Detailed metrics
        with st.expander("📊 Detailed Metrics"):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Thermal**")
                st.write(f"- Radiator effective area: {radiator_effective_area:.2f} m²")
                st.write(f"- Heat dissipation: {heat_dissipation:.0f} W")
                st.write(f"- Equilibrium temp: {t_eq_c:.1f}°C")
                st.write(f"- Max allowable: {t_max}°C")

            with col2:
                st.markdown("**Power**")
                st.write(f"- Solar generation: {solar_power:.0f} W (peak)")
                st.write(f"- Orbit average: {orbit_avg_solar:.0f} W")
                st.write(f"- AI consumption: {ai_power:.0f} W")
                st.write(f"- Eclipse fraction: {orbit.eclipse_fraction*100:.1f}%")

            st.markdown("**Compute**")
            if chip_count > 0:
                tflops = chip.get("tflops_fp16", chip.get("tflops_bf16", 0)) * chip_count * duty_cycle
                st.write(f"- Effective TFLOPS: {tflops:.0f}")
                st.write(f"- TFLOPS/kg: {tflops/total_mass:.1f}")
                st.write(f"- TFLOPS/W: {tflops/ai_power:.2f}" if ai_power > 0 else "- TFLOPS/W: N/A")


if __name__ == "__main__":
    main()
