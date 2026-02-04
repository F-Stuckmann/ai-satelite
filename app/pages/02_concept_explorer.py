"""Concept Explorer - Browse and compare design concepts."""

import streamlit as st
import plotly.graph_objects as go

st.set_page_config(page_title="Concept Explorer", page_icon="💡", layout="wide")

import sys
sys.path.insert(0, ".")

from app.concepts import CONCEPTS
from src.core.constants import MATERIALS, AI_CHIPS, SOLAR_CELLS, STARLINK_REFERENCE
from src.thermal.stefan_boltzmann import equilibrium_temperature, temperature_to_celsius
from src.thermal.grid_fin_radiator import GridFinRadiator
from src.thermal.view_factors import ViewFactors


def calculate_concept_metrics(concept: dict) -> dict:
    """Calculate W/kg and other metrics for a concept."""
    params = concept["parameters"]

    # Get material properties
    material = params.get("material", "aluminum_6061")
    mat = MATERIALS.get(material, MATERIALS["aluminum_6061"])

    # Calculate chassis mass
    l, w, h = params["chassis_length_m"], params["chassis_width_m"], params["chassis_height_m"]
    wall_t = params["wall_thickness_mm"] / 1000
    surface_area = 2 * (l*w + l*h + w*h)
    chassis_mass = surface_area * wall_t * mat["rho"]

    # Calculate radiator mass and effective area
    radiator_area = params["radiator_base_area_m2"]
    if params["radiator_type"] == "grid_fin":
        radiator = GridFinRadiator(base_area_m2=radiator_area, material=material)
        radiator_effective = radiator.effective_area_m2
        radiator_mass = radiator.mass_kg
    else:
        radiator_effective = radiator_area
        radiator_mass = radiator_area * 0.003 * mat["rho"]

    # Calculate solar mass and power
    solar_type = params.get("solar_cell_type", "triple_junction_gaas")
    solar_cell = SOLAR_CELLS.get(solar_type, SOLAR_CELLS["triple_junction_gaas"])
    solar_area = params["solar_area_m2"]
    solar_mass = solar_area * solar_cell["mass_per_area"]
    solar_power_peak = solar_area * 1361 * solar_cell["efficiency"]

    # Orbit parameters
    altitude = params.get("altitude_km", 550)
    orbit = ViewFactors(altitude_km=altitude)
    solar_power_avg = solar_power_peak * orbit.sunlight_fraction * 0.9

    # Battery mass
    battery_wh = params.get("battery_capacity_wh", 500)
    battery_mass = battery_wh / 150  # Assume 150 Wh/kg

    # AI payload
    chip_type = params.get("ai_chip_type")
    chip_count = params.get("ai_chip_count", 0)
    duty_cycle = params.get("duty_cycle", 0.5)

    if chip_type and chip_count > 0 and chip_type in AI_CHIPS:
        chip = AI_CHIPS[chip_type]
        ai_mass = chip["mass_kg"] * chip_count
        ai_power = chip["tdp"] * chip_count * duty_cycle
    else:
        ai_mass = 0
        ai_power = 0

    # Other mass
    other_mass = 15

    # Totals
    total_mass = chassis_mass + radiator_mass + solar_mass + battery_mass + ai_mass + other_mass
    power_available = solar_power_avg

    # Thermal
    heat_load = ai_power + solar_power_avg * 0.05
    if heat_load > 0:
        t_eq_k = equilibrium_temperature(heat_load, radiator_effective, 0.9)
        t_eq_c = temperature_to_celsius(t_eq_k)
    else:
        t_eq_c = 20

    # W/kg
    w_per_kg = power_available / total_mass if total_mass > 0 else 0

    return {
        "total_mass_kg": total_mass,
        "power_available_w": power_available,
        "w_per_kg": w_per_kg,
        "temp_c": t_eq_c,
        "mass_breakdown": {
            "Chassis": chassis_mass,
            "Radiator": radiator_mass,
            "Solar": solar_mass,
            "Battery": battery_mass,
            "AI Payload": ai_mass,
            "Other": other_mass,
        },
    }


def main():
    st.title("💡 Concept Explorer")
    st.markdown("**Compare different satellite design concepts**")

    # Concept selection with radio buttons
    st.subheader("Select Design Concept")

    concept_names = {
        "starlink_baseline": "🛰️ Starlink Baseline (SpaceX Patent)",
        "solar_heatsink": "☀️ Solar-on-Heatsink (Your Idea!)",
        "body_mounted": "📦 Body-Mounted Solar",
        "graphene_radiator": "⚡ Graphene Radiator (Future Tech)",
    }

    selected_key = st.radio(
        "Choose a concept to explore:",
        options=list(CONCEPTS.keys()),
        format_func=lambda x: concept_names.get(x, x),
        horizontal=True,
    )

    concept = CONCEPTS[selected_key]
    metrics = calculate_concept_metrics(concept)

    # Layout
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader(concept["name"])

        # W/kg comparison
        w_per_kg = metrics["w_per_kg"]
        delta = w_per_kg - STARLINK_REFERENCE["w_per_kg"]

        st.metric(
            label="W/kg",
            value=f"{w_per_kg:.2f}",
            delta=f"{delta:+.2f} vs SpaceX baseline",
        )

        # Status indicators
        col_a, col_b = st.columns(2)
        with col_a:
            st.metric("Total Mass", f"{metrics['total_mass_kg']:.1f} kg")
        with col_b:
            st.metric("Power", f"{metrics['power_available_w']:.0f} W")

        st.metric("Operating Temp", f"{metrics['temp_c']:.1f}°C")

        # Description
        st.markdown("---")
        st.markdown("**Description:**")
        st.markdown(concept["description"])

        # Notes
        if "notes" in concept:
            st.markdown("**Key Points:**")
            for note in concept["notes"]:
                st.markdown(f"- {note}")

    with col2:
        # Mass breakdown chart
        st.subheader("Mass Breakdown")
        mass_data = metrics["mass_breakdown"]

        fig = go.Figure(data=[go.Pie(
            labels=list(mass_data.keys()),
            values=list(mass_data.values()),
            hole=0.4,
        )])
        fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20))
        st.plotly_chart(fig, use_container_width=True)

        # Parameters
        st.subheader("Parameters")
        params = concept["parameters"]

        with st.expander("Structure"):
            st.write(f"- Dimensions: {params['chassis_length_m']}×{params['chassis_width_m']}×{params['chassis_height_m']} m")
            st.write(f"- Material: {MATERIALS.get(params['material'], {}).get('name', params['material'])}")
            st.write(f"- Wall thickness: {params['wall_thickness_mm']} mm")

        with st.expander("Thermal"):
            st.write(f"- Radiator type: {params['radiator_type']}")
            st.write(f"- Radiator area: {params['radiator_base_area_m2']} m²")

        with st.expander("Power"):
            solar_type = params.get("solar_cell_type", "triple_junction_gaas")
            st.write(f"- Solar cells: {SOLAR_CELLS.get(solar_type, {}).get('name', solar_type)}")
            st.write(f"- Solar area: {params['solar_area_m2']} m²")
            st.write(f"- Battery: {params['battery_capacity_wh']} Wh")

        # Analysis required (for new concepts)
        if "analysis_required" in concept:
            st.subheader("Analysis Required")
            for item in concept["analysis_required"]:
                st.markdown(f"- ⚠️ {item}")

        # Potential W/kg
        if "potential_w_per_kg" in concept:
            st.info(f"**Estimated potential:** {concept['potential_w_per_kg']} W/kg")

    # Comparison chart
    st.markdown("---")
    st.subheader("W/kg Comparison")

    comparison_data = []
    for key, c in CONCEPTS.items():
        m = calculate_concept_metrics(c)
        comparison_data.append({
            "name": c["name"],
            "w_per_kg": m["w_per_kg"],
            "color": "green" if m["w_per_kg"] >= STARLINK_REFERENCE["w_per_kg"] else "orange",
        })

    fig_bar = go.Figure()
    fig_bar.add_trace(go.Bar(
        x=[d["name"] for d in comparison_data],
        y=[d["w_per_kg"] for d in comparison_data],
        marker_color=["green" if d["w_per_kg"] >= STARLINK_REFERENCE["w_per_kg"] else "orange" for d in comparison_data],
    ))

    # Add reference line
    fig_bar.add_hline(
        y=STARLINK_REFERENCE["w_per_kg"],
        line_dash="dash",
        line_color="red",
        annotation_text=f"SpaceX Patent: {STARLINK_REFERENCE['w_per_kg']} W/kg",
    )

    fig_bar.update_layout(
        yaxis_title="W/kg",
        height=400,
    )
    st.plotly_chart(fig_bar, use_container_width=True)


if __name__ == "__main__":
    main()
