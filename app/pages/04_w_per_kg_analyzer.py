"""W/kg Analyzer - Power-to-Mass Ratio Optimization Tool."""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.core.constants import (
    MATERIALS,
    SOLAR_CELLS,
    BATTERIES,
    AI_CHIPS,
    STARLINK_REFERENCE,
)
from src.thermal.stefan_boltzmann import (
    equilibrium_temperature,
    required_radiator_area,
)


st.set_page_config(page_title="W/kg Analyzer", page_icon="⚡", layout="wide")

st.title("⚡ W/kg Power-to-Mass Analyzer")
st.markdown("""
Explore and optimize the **power-to-mass ratio** (W/kg) for your AI satellite design.

**Reference**: SpaceX Starlink (US 11,834,205 B1) achieves **12.6 W/kg** (3,300W / 261kg)
""")

# Sidebar configuration
st.sidebar.header("Design Parameters")

# Power System
st.sidebar.subheader("Power System")
solar_type = st.sidebar.selectbox(
    "Solar Cell Type",
    list(SOLAR_CELLS.keys()),
    format_func=lambda x: SOLAR_CELLS[x]["name"],
)
solar_area = st.sidebar.slider("Solar Array Area (m²)", 1.0, 30.0, 8.5, 0.5)
solar_efficiency = SOLAR_CELLS[solar_type]["efficiency"]

battery_type = st.sidebar.selectbox(
    "Battery Type",
    list(BATTERIES.keys()),
    format_func=lambda x: BATTERIES[x]["name"],
)
battery_capacity = st.sidebar.slider("Battery Capacity (Wh)", 100, 5000, 1000, 100)

# Structure
st.sidebar.subheader("Structure")
material_type = st.sidebar.selectbox(
    "Chassis Material",
    list(MATERIALS.keys()),
    format_func=lambda x: MATERIALS[x]["name"],
)
chassis_volume = st.sidebar.slider("Chassis Volume (L)", 50, 500, 200, 10)
wall_thickness = st.sidebar.slider("Wall Thickness (mm)", 1.0, 10.0, 3.0, 0.5)

# Thermal System
st.sidebar.subheader("Thermal System")
radiator_area = st.sidebar.slider("Radiator Area (m²)", 1.0, 20.0, 6.0, 0.5)
emissivity = st.sidebar.slider("Radiator Emissivity", 0.7, 0.98, 0.9, 0.01)

# AI Payload
st.sidebar.subheader("AI Payload")
chip_type = st.sidebar.selectbox(
    "AI Chip Type",
    list(AI_CHIPS.keys()),
    format_func=lambda x: AI_CHIPS[x]["name"],
)
chip_count = st.sidebar.slider("Number of Chips", 0, 8, 2, 1)
duty_cycle = st.sidebar.slider("Compute Duty Cycle (%)", 10, 100, 50, 5)


# Calculate masses
def calculate_masses():
    """Calculate component masses."""
    # Solar array mass (kg/m²)
    solar_density = SOLAR_CELLS[solar_type]["mass_per_area"]
    solar_mass = solar_area * solar_density

    # Battery mass
    battery_specific_energy = BATTERIES[battery_type]["specific_energy"]  # Wh/kg
    battery_mass = battery_capacity / battery_specific_energy

    # Chassis mass (simplified box model)
    material = MATERIALS[material_type]
    # Approximate surface area from volume (assume cube-ish)
    side = (chassis_volume / 1000) ** (1/3)  # Convert L to m³
    surface_area = 6 * side ** 2
    chassis_mass = surface_area * (wall_thickness / 1000) * material["rho"]

    # Radiator mass (assume aluminum honeycomb, ~2 kg/m²)
    radiator_mass = radiator_area * 2.0

    # AI chip mass
    chip_mass = chip_count * AI_CHIPS[chip_type]["mass_kg"]

    # Power electronics, wiring, etc (~10% of other masses)
    other_mass = 0.1 * (solar_mass + battery_mass + chassis_mass)

    return {
        "solar_array": solar_mass,
        "battery": battery_mass,
        "chassis": chassis_mass,
        "radiator": radiator_mass,
        "ai_chips": chip_mass,
        "other": other_mass,
        "total": solar_mass + battery_mass + chassis_mass + radiator_mass + chip_mass + other_mass,
    }


def calculate_power():
    """Calculate power budget."""
    # Solar power generation
    solar_flux = 1361  # W/m²
    solar_power_generated = solar_area * solar_efficiency * solar_flux

    # AI chip power consumption
    chip_tdp = AI_CHIPS[chip_type]["tdp"]
    ai_power = chip_count * chip_tdp * (duty_cycle / 100)

    # System overhead (~15% of total)
    system_overhead = 0.15 * ai_power

    # Thermal management power (~5%)
    thermal_power = 0.05 * ai_power

    total_consumption = ai_power + system_overhead + thermal_power

    return {
        "generated": solar_power_generated,
        "ai_chips": ai_power,
        "system_overhead": system_overhead,
        "thermal_management": thermal_power,
        "total_consumption": total_consumption,
        "margin": solar_power_generated - total_consumption,
    }


def check_thermal_feasibility():
    """Check if thermal design can reject the heat."""
    power = calculate_power()
    heat_to_reject = power["total_consumption"]

    # Equilibrium temperature
    temp = equilibrium_temperature(heat_to_reject, radiator_area, emissivity)

    # Required area to stay below 80°C (353K)
    max_temp = 353  # K
    required_area = required_radiator_area(heat_to_reject, max_temp, emissivity)

    return {
        "equilibrium_temp_k": temp,
        "equilibrium_temp_c": temp - 273.15,
        "required_area_m2": required_area,
        "thermal_margin": radiator_area - required_area,
        "feasible": temp < 353,
    }


# Main calculations
masses = calculate_masses()
power = calculate_power()
thermal = check_thermal_feasibility()

# Calculate W/kg
available_power = min(power["generated"], power["total_consumption"] + power["margin"])
w_per_kg = available_power / masses["total"]

# Display main metrics
col1, col2, col3, col4 = st.columns(4)

with col1:
    delta_vs_starlink = w_per_kg - 12.6
    st.metric(
        "W/kg",
        f"{w_per_kg:.1f}",
        f"{delta_vs_starlink:+.1f} vs Starlink",
        delta_color="normal" if delta_vs_starlink >= 0 else "inverse",
    )

with col2:
    st.metric("Total Mass", f"{masses['total']:.1f} kg")

with col3:
    st.metric("Available Power", f"{available_power:.0f} W")

with col4:
    temp_status = "✓" if thermal["feasible"] else "✗"
    st.metric(
        "Thermal Status",
        f"{thermal['equilibrium_temp_c']:.0f}°C {temp_status}",
    )

# Comparison with references
st.markdown("---")
st.subheader("Comparison with Reference Designs")

comparison_data = {
    "Design": ["Your Design", "SpaceX Starlink", "Boeing 702", "NASA Dawn", "NASA DS1"],
    "W/kg": [w_per_kg, 12.6, 3.0, 8.0, 5.0],
    "Status": [
        "Current",
        "Reference (US 11,834,205)",
        "Traditional GEO",
        "Deep Space",
        "Ion Propulsion",
    ],
}

fig_comparison = go.Figure()
colors = ["#2E86AB", "#A23B72", "#E94F37", "#F18F01", "#C73E1D"]
for i, (design, wpkg) in enumerate(zip(comparison_data["Design"], comparison_data["W/kg"])):
    fig_comparison.add_trace(go.Bar(
        x=[design],
        y=[wpkg],
        name=design,
        marker_color=colors[i],
        text=[f"{wpkg:.1f}"],
        textposition="auto",
    ))

fig_comparison.update_layout(
    title="W/kg Comparison",
    yaxis_title="W/kg",
    showlegend=False,
    height=400,
)
st.plotly_chart(fig_comparison, use_container_width=True)

# Mass and Power Breakdown
col1, col2 = st.columns(2)

with col1:
    st.subheader("Mass Breakdown")
    mass_labels = list(masses.keys())[:-1]  # Exclude total
    mass_values = [masses[k] for k in mass_labels]

    fig_mass = go.Figure(data=[go.Pie(
        labels=[l.replace("_", " ").title() for l in mass_labels],
        values=mass_values,
        hole=0.4,
        marker_colors=["#2E86AB", "#A23B72", "#E94F37", "#F18F01", "#C73E1D", "#6B4226"],
    )])
    fig_mass.update_layout(height=350)
    st.plotly_chart(fig_mass, use_container_width=True)

with col2:
    st.subheader("Power Budget")
    power_labels = ["AI Chips", "System Overhead", "Thermal Mgmt", "Margin"]
    power_values = [
        power["ai_chips"],
        power["system_overhead"],
        power["thermal_management"],
        max(0, power["margin"]),
    ]

    fig_power = go.Figure(data=[go.Pie(
        labels=power_labels,
        values=power_values,
        hole=0.4,
        marker_colors=["#2E86AB", "#A23B72", "#E94F37", "#4CAF50"],
    )])
    fig_power.update_layout(height=350)
    st.plotly_chart(fig_power, use_container_width=True)

# Thermal Analysis
st.markdown("---")
st.subheader("Thermal Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    **Heat to Reject**: {power['total_consumption']:.0f} W

    **Radiator Area**: {radiator_area:.1f} m²

    **Required Area**: {thermal['required_area_m2']:.1f} m² (for 80°C limit)

    **Thermal Margin**: {thermal['thermal_margin']:.1f} m²
    """)

    if thermal["feasible"]:
        st.success("✓ Thermal design is feasible")
    else:
        st.error(f"✗ Need {thermal['required_area_m2']:.1f} m² radiator area")

with col2:
    # Temperature vs radiator area curve
    areas = np.linspace(1, 20, 100)
    temps = [equilibrium_temperature(power["total_consumption"], a, emissivity) - 273.15 for a in areas]

    fig_temp = go.Figure()
    fig_temp.add_trace(go.Scatter(
        x=areas,
        y=temps,
        mode="lines",
        name="Temperature",
        line=dict(color="#2E86AB", width=2),
    ))
    fig_temp.add_hline(y=80, line_dash="dash", line_color="red", annotation_text="80°C Limit")
    fig_temp.add_vline(x=radiator_area, line_dash="dash", line_color="green", annotation_text="Current")

    fig_temp.update_layout(
        title="Equilibrium Temperature vs Radiator Area",
        xaxis_title="Radiator Area (m²)",
        yaxis_title="Temperature (°C)",
        height=300,
    )
    st.plotly_chart(fig_temp, use_container_width=True)

# Optimization suggestions
st.markdown("---")
st.subheader("Optimization Suggestions")

suggestions = []

if w_per_kg < 12.6:
    suggestions.append("• Consider using higher efficiency solar cells to increase power")
    suggestions.append("• Reduce chassis mass with lighter materials (CFRP)")
    suggestions.append("• Optimize battery sizing - you may have excess capacity")

if not thermal["feasible"]:
    suggestions.append("• Increase radiator area or use grid fin design for better heat rejection")
    suggestions.append("• Reduce AI chip count or duty cycle")
    suggestions.append("• Consider higher emissivity coating")

if power["margin"] < 0:
    suggestions.append("• Power deficit - increase solar array or reduce consumption")

if masses["battery"] > masses["total"] * 0.2:
    suggestions.append("• Battery mass is >20% - consider higher specific energy cells")

if not suggestions:
    suggestions.append("✓ Design looks well optimized!")

for s in suggestions:
    st.markdown(s)

# Sensitivity Analysis
st.markdown("---")
st.subheader("Sensitivity Analysis")

st.markdown("How W/kg changes with key parameters:")

# Create sensitivity plots
fig_sens = make_subplots(rows=1, cols=3, subplot_titles=[
    "Solar Area Impact",
    "Chassis Material Impact",
    "AI Chip Count Impact",
])

# Solar area sensitivity
areas_sens = np.linspace(2, 30, 20)
wpkg_solar = []
for a in areas_sens:
    power_gen = a * solar_efficiency * 1361
    mass_add = a * SOLAR_CELLS[solar_type]["mass_per_area"]
    base_mass = masses["total"] - masses["solar_array"]
    wpkg_solar.append(power_gen / (base_mass + mass_add))

fig_sens.add_trace(
    go.Scatter(x=areas_sens, y=wpkg_solar, mode="lines", name="Solar Area"),
    row=1, col=1
)

# Material sensitivity
mat_names = []
mat_wpkg = []
for mat_key, mat_data in MATERIALS.items():
    side = (chassis_volume / 1000) ** (1/3)
    surface_area = 6 * side ** 2
    chassis_m = surface_area * (wall_thickness / 1000) * mat_data["rho"]
    base_mass = masses["total"] - masses["chassis"]
    mat_wpkg.append(available_power / (base_mass + chassis_m))
    mat_names.append(mat_data["name"][:10])

fig_sens.add_trace(
    go.Bar(x=mat_names, y=mat_wpkg, name="Material"),
    row=1, col=2
)

# Chip count sensitivity
chip_counts = range(0, 9)
wpkg_chips = []
for n in chip_counts:
    chip_m = n * AI_CHIPS[chip_type]["mass_kg"]
    base_mass = masses["total"] - masses["ai_chips"]
    wpkg_chips.append(available_power / (base_mass + chip_m))

fig_sens.add_trace(
    go.Scatter(x=list(chip_counts), y=wpkg_chips, mode="lines+markers", name="Chips"),
    row=1, col=3
)

fig_sens.update_layout(height=350, showlegend=False)
fig_sens.update_xaxes(title_text="Area (m²)", row=1, col=1)
fig_sens.update_xaxes(title_text="Material", row=1, col=2)
fig_sens.update_xaxes(title_text="Chip Count", row=1, col=3)
fig_sens.update_yaxes(title_text="W/kg", row=1, col=1)

st.plotly_chart(fig_sens, use_container_width=True)

# Export summary
st.markdown("---")
with st.expander("Export Design Summary"):
    summary = f"""
# AI Satellite Design Summary

## Key Metrics
- **W/kg**: {w_per_kg:.2f}
- **Total Mass**: {masses['total']:.1f} kg
- **Available Power**: {available_power:.0f} W
- **Equilibrium Temperature**: {thermal['equilibrium_temp_c']:.0f}°C

## Power System
- Solar Cell Type: {SOLAR_CELLS[solar_type]['name']}
- Solar Array Area: {solar_area} m²
- Solar Efficiency: {solar_efficiency*100:.0f}%
- Battery: {BATTERIES[battery_type]['name']} ({battery_capacity} Wh)

## Structure
- Material: {MATERIALS[material_type]['name']}
- Chassis Volume: {chassis_volume} L
- Wall Thickness: {wall_thickness} mm

## Thermal
- Radiator Area: {radiator_area} m²
- Emissivity: {emissivity}
- Thermal Margin: {thermal['thermal_margin']:.1f} m²

## AI Payload
- Chip Type: {AI_CHIPS[chip_type]['name']}
- Chip Count: {chip_count}
- Duty Cycle: {duty_cycle}%

## Mass Breakdown
{chr(10).join([f'- {k.replace("_", " ").title()}: {v:.1f} kg' for k, v in masses.items()])}
"""
    st.code(summary, language="markdown")
