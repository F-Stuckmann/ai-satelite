"""Structural Analysis - Thermal Stress and Material Properties."""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.core.constants import MATERIALS
from src.structural.materials import (
    MaterialProperties,
    thermal_stress,
    thermal_strain,
    safety_factor,
    buckling_stress_panel,
    natural_frequency_panel,
    cte_mismatch_stress,
    structural_summary,
)


st.set_page_config(page_title="Structural Analysis", page_icon="🔧", layout="wide")

st.title("🔧 Structural Analysis")
st.markdown("""
Analyze **thermal stress**, **buckling**, and **vibration** characteristics of satellite structures.

Key concerns for space structures:
- Thermal cycling stress (eclipse transitions: -100°C to +80°C)
- Launch vibration loads
- CTE mismatch between bonded materials
""")

# Sidebar configuration
st.sidebar.header("Structure Parameters")

# Material selection
material_key = st.sidebar.selectbox(
    "Primary Material",
    list(MATERIALS.keys()),
    format_func=lambda x: MATERIALS[x]["name"],
)

# Panel dimensions
st.sidebar.subheader("Panel Dimensions")
length_m = st.sidebar.slider("Length (m)", 0.1, 2.0, 0.5, 0.05)
width_m = st.sidebar.slider("Width (m)", 0.1, 2.0, 0.3, 0.05)
thickness_mm = st.sidebar.slider("Thickness (mm)", 0.5, 10.0, 2.0, 0.5)

# Thermal environment
st.sidebar.subheader("Thermal Environment")
temp_min = st.sidebar.slider("Min Temperature (°C)", -150, 0, -100, 10)
temp_max = st.sidebar.slider("Max Temperature (°C)", 0, 150, 80, 10)

# Boundary conditions
boundary = st.sidebar.selectbox(
    "Boundary Condition",
    ["simply_supported", "clamped"],
    format_func=lambda x: x.replace("_", " ").title(),
)

with st.sidebar.expander("ℹ️ Boundary Conditions"):
    st.markdown("""
**Simply Supported**: Edges can rotate freely but cannot translate.
- Lower buckling resistance
- Lower natural frequency
- Common for panels attached with hinges/pins

**Clamped (Fixed)**: Edges cannot rotate or translate.
- Higher buckling resistance (~1.7x)
- Higher natural frequency
- Common for welded/bonded edges

*Satellite panels are typically between these extremes.*
""")

# Load material
material = MaterialProperties.from_database(material_key)

# Calculate thermal properties
delta_T = temp_max - temp_min
thermal_stress_val = thermal_stress(material, delta_T)
thermal_strain_val = thermal_strain(material, delta_T)
thermal_sf = safety_factor(thermal_stress_val, material)

# Calculate structural properties
buckling = buckling_stress_panel(length_m, width_m, thickness_mm, material, boundary)
nat_freq = natural_frequency_panel(length_m, width_m, thickness_mm, material, boundary)
buckling_sf = safety_factor(thermal_stress_val, material, use_ultimate=False)

# Main metrics display
st.subheader("Key Results")

col1, col2, col3, col4 = st.columns(4)

with col1:
    sf_color = "normal" if thermal_sf > 2.0 else ("inverse" if thermal_sf < 1.5 else "off")
    st.metric(
        "Thermal Safety Factor",
        f"{thermal_sf:.2f}",
        "Safe" if thermal_sf > 1.5 else "Warning",
        delta_color=sf_color,
    )

with col2:
    st.metric(
        "Thermal Stress",
        f"{thermal_stress_val/1e6:.1f} MPa",
    )

with col3:
    st.metric(
        "Natural Frequency",
        f"{nat_freq:.1f} Hz",
    )

with col4:
    st.metric(
        "Buckling Stress",
        f"{buckling/1e6:.1f} MPa",
    )

# Material properties display
st.markdown("---")
col1, col2 = st.columns(2)

with col1:
    st.subheader(f"Material: {material.name}")

    mat_props = {
        "Young's Modulus (E)": f"{material.E/1e9:.0f} GPa",
        "Density (ρ)": f"{material.rho:.0f} kg/m³",
        "Yield Strength (σ_y)": f"{material.sigma_yield/1e6:.0f} MPa",
        "Poisson's Ratio (ν)": f"{material.nu:.2f}",
        "CTE (α)": f"{material.alpha*1e6:.1f} ppm/K",
        "Thermal Conductivity (k)": f"{material.k:.0f} W/(m·K)",
        "Shear Modulus (G)": f"{material.G/1e9:.1f} GPa",
        "Bulk Modulus (K)": f"{material.K/1e9:.1f} GPa",
    }

    for prop, value in mat_props.items():
        st.markdown(f"**{prop}**: {value}")

with col2:
    st.subheader("Panel Properties")

    panel_area = length_m * width_m
    panel_volume = panel_area * (thickness_mm / 1000)
    panel_mass = panel_volume * material.rho

    panel_props = {
        "Dimensions": f"{length_m*1000:.0f} × {width_m*1000:.0f} × {thickness_mm:.1f} mm",
        "Area": f"{panel_area*1e4:.0f} cm²",
        "Volume": f"{panel_volume*1e6:.1f} cm³",
        "Mass": f"{panel_mass*1000:.1f} g",
        "Boundary": boundary.replace("_", " ").title(),
    }

    for prop, value in panel_props.items():
        st.markdown(f"**{prop}**: {value}")

# Thermal stress analysis
st.markdown("---")
st.subheader("Thermal Stress Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    **Temperature Range**: {temp_min}°C to {temp_max}°C (ΔT = {delta_T}°C)

    **Thermal Strain**: {thermal_strain_val*1e6:.1f} με (microstrain)

    **Thermal Stress** (fully constrained): {thermal_stress_val/1e6:.2f} MPa

    **Yield Strength**: {material.sigma_yield/1e6:.0f} MPa

    **Safety Factor**: {thermal_sf:.2f}
    """)

    if thermal_sf > 2.0:
        st.success("✓ Excellent thermal margin")
    elif thermal_sf > 1.5:
        st.warning("⚠ Adequate but consider design margin")
    else:
        st.error("✗ Thermal stress exceeds safe limits")

with col2:
    # Stress vs temperature plot
    temps = np.linspace(-150, 150, 100)
    stresses = [thermal_stress(material, abs(t - 20)) / 1e6 for t in temps]  # From 20°C reference

    fig_stress = go.Figure()
    fig_stress.add_trace(go.Scatter(
        x=temps,
        y=stresses,
        mode="lines",
        name="Thermal Stress",
        line=dict(color="#2E86AB", width=2),
    ))
    fig_stress.add_hline(
        y=material.sigma_yield/1e6,
        line_dash="dash",
        line_color="red",
        annotation_text="Yield Strength",
    )
    fig_stress.add_vrect(
        x0=temp_min, x1=temp_max,
        fillcolor="lightblue", opacity=0.3,
        annotation_text="Operating Range",
    )

    fig_stress.update_layout(
        title="Thermal Stress vs Temperature",
        xaxis_title="Temperature (°C)",
        yaxis_title="Stress (MPa)",
        height=350,
    )
    st.plotly_chart(fig_stress, use_container_width=True)

# Buckling analysis
st.markdown("---")
st.subheader("Buckling Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    **Critical Buckling Stress**: {buckling/1e6:.1f} MPa

    **Applied Stress** (thermal): {thermal_stress_val/1e6:.2f} MPa

    **Buckling Safety Factor**: {buckling/thermal_stress_val:.2f}

    Buckling occurs when compressive stress exceeds the critical value.
    Thicker panels and clamped boundaries increase buckling resistance.
    """)

with col2:
    # Buckling vs thickness
    thicknesses = np.linspace(0.5, 10, 50)
    buckling_stresses = [
        buckling_stress_panel(length_m, width_m, t, material, boundary) / 1e6
        for t in thicknesses
    ]

    fig_buck = go.Figure()
    fig_buck.add_trace(go.Scatter(
        x=thicknesses,
        y=buckling_stresses,
        mode="lines",
        name="Critical Buckling Stress",
        line=dict(color="#E94F37", width=2),
    ))
    fig_buck.add_hline(
        y=thermal_stress_val/1e6,
        line_dash="dash",
        line_color="blue",
        annotation_text="Thermal Stress",
    )
    fig_buck.add_vline(
        x=thickness_mm,
        line_dash="dash",
        line_color="green",
        annotation_text="Current",
    )

    fig_buck.update_layout(
        title="Critical Buckling Stress vs Thickness",
        xaxis_title="Thickness (mm)",
        yaxis_title="Buckling Stress (MPa)",
        height=350,
    )
    st.plotly_chart(fig_buck, use_container_width=True)

# Vibration analysis
st.markdown("---")
st.subheader("Vibration Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    **Fundamental Frequency**: {nat_freq:.1f} Hz

    **Launch Environment Requirements**:
    - First mode > 35 Hz (typical launcher requirement)
    - Avoid resonance with launch vehicle frequencies

    **Current Status**: {"✓ Meets 35 Hz requirement" if nat_freq > 35 else "✗ Below 35 Hz - stiffen structure"}
    """)

    if nat_freq < 35:
        st.warning("Consider increasing thickness or using stiffer material")

with col2:
    # Frequency vs dimensions
    dims = np.linspace(0.1, 1.0, 50)
    freqs_length = [
        natural_frequency_panel(d, width_m, thickness_mm, material, boundary)
        for d in dims
    ]
    freqs_width = [
        natural_frequency_panel(length_m, d, thickness_mm, material, boundary)
        for d in dims
    ]

    fig_freq = go.Figure()
    fig_freq.add_trace(go.Scatter(
        x=dims,
        y=freqs_length,
        mode="lines",
        name="Varying Length",
        line=dict(color="#2E86AB", width=2),
    ))
    fig_freq.add_trace(go.Scatter(
        x=dims,
        y=freqs_width,
        mode="lines",
        name="Varying Width",
        line=dict(color="#A23B72", width=2),
    ))
    fig_freq.add_hline(
        y=35,
        line_dash="dash",
        line_color="red",
        annotation_text="35 Hz Requirement",
    )

    fig_freq.update_layout(
        title="Natural Frequency vs Panel Dimension",
        xaxis_title="Dimension (m)",
        yaxis_title="Frequency (Hz)",
        height=350,
    )
    st.plotly_chart(fig_freq, use_container_width=True)

# CTE Mismatch Analysis
st.markdown("---")
st.subheader("CTE Mismatch Analysis")

st.markdown("""
When bonding dissimilar materials (e.g., solar cells on aluminum radiator),
CTE mismatch causes interfacial stress during thermal cycling.
""")

col1, col2 = st.columns(2)

with col1:
    secondary_material = st.selectbox(
        "Secondary Material (bonded)",
        list(MATERIALS.keys()),
        index=1,  # Select second option by default
        format_func=lambda x: MATERIALS[x]["name"],
        key="secondary_mat",
    )

    mat2 = MaterialProperties.from_database(secondary_material)

    sigma1, sigma2 = cte_mismatch_stress(material, mat2, delta_T)

    st.markdown(f"""
    **Material 1**: {material.name} (CTE = {material.alpha*1e6:.1f} ppm/K)

    **Material 2**: {mat2.name} (CTE = {mat2.alpha*1e6:.1f} ppm/K)

    **CTE Difference**: {abs(material.alpha - mat2.alpha)*1e6:.1f} ppm/K

    **Interface Stresses**:
    - In {material.name}: {sigma1/1e6:.2f} MPa
    - In {mat2.name}: {sigma2/1e6:.2f} MPa
    """)

with col2:
    # CTE comparison chart
    materials_list = list(MATERIALS.keys())
    cte_values = [MATERIALS[m]["alpha"] * 1e6 for m in materials_list]
    mat_names = [MATERIALS[m]["name"] for m in materials_list]

    fig_cte = go.Figure()
    fig_cte.add_trace(go.Bar(
        x=mat_names,
        y=cte_values,
        marker_color=["#2E86AB", "#A23B72", "#E94F37", "#F18F01"],
        text=[f"{v:.1f}" for v in cte_values],
        textposition="auto",
    ))
    fig_cte.update_layout(
        title="CTE Comparison (ppm/K)",
        yaxis_title="CTE (ppm/K)",
        height=350,
    )
    st.plotly_chart(fig_cte, use_container_width=True)

# Material comparison
st.markdown("---")
st.subheader("Material Comparison")

# Create comparison table
comparison_data = []
for mat_key in MATERIALS.keys():
    mat = MaterialProperties.from_database(mat_key)
    stress = thermal_stress(mat, delta_T)
    sf = safety_factor(stress, mat)
    freq = natural_frequency_panel(length_m, width_m, thickness_mm, mat, boundary)
    buck = buckling_stress_panel(length_m, width_m, thickness_mm, mat, boundary)

    comparison_data.append({
        "Material": mat.name,
        "Density (kg/m³)": mat.rho,
        "E (GPa)": mat.E / 1e9,
        "Yield (MPa)": mat.sigma_yield / 1e6,
        "CTE (ppm/K)": mat.alpha * 1e6,
        "Thermal SF": sf,
        "Nat Freq (Hz)": freq,
        "Buckling (MPa)": buck / 1e6,
    })

# Display as table
import pandas as pd
df = pd.DataFrame(comparison_data)
df = df.round(2)
st.dataframe(df, use_container_width=True)

# Structural summary
st.markdown("---")
with st.expander("Structural Design Summary"):
    summary = structural_summary(
        (length_m, width_m, 0.1),  # Assuming 0.1m height
        thickness_mm,
        material_key,
        (temp_min, temp_max),
    )

    st.json(summary)
