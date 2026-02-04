"""Constellation Simulator - Multi-satellite Walker constellation visualization."""

import streamlit as st
import numpy as np
import plotly.graph_objects as go

from src.orbital.keplerian import (
    OrbitalElements,
    orbital_period,
    orbital_velocity,
    orbits_per_day,
)
from src.core.constants import EARTH_RADIUS_KM


st.set_page_config(page_title="Constellation Simulator", page_icon="🌐", layout="wide")

st.title("🌐 Constellation Simulator")
st.markdown("""
Simulate and visualize **Walker constellation** patterns for satellite networks.

**Reference**: Starlink Phase 1 uses 72 planes × 22 satellites = 1,584 satellites at 550 km
""")

# Sidebar configuration
st.sidebar.header("Constellation Parameters")

# Walker constellation parameters
st.sidebar.subheader("Walker Delta Pattern")
altitude_km = st.sidebar.slider("Altitude (km)", 300, 1200, 550, 50)
inclination_deg = st.sidebar.slider("Inclination (°)", 0, 98, 53, 1)
n_planes = st.sidebar.slider("Number of Orbital Planes", 1, 24, 6, 1)
sats_per_plane = st.sidebar.slider("Satellites per Plane", 1, 22, 11, 1)
phasing_factor = st.sidebar.slider("Phasing Factor (F)", 0, n_planes - 1, 1, 1)

# Visualization options
st.sidebar.subheader("Visualization")
show_earth = st.sidebar.checkbox("Show Earth", value=True)
show_orbits = st.sidebar.checkbox("Show Orbit Paths", value=True)
show_coverage = st.sidebar.checkbox("Show Coverage Circles", value=False)
animation_speed = st.sidebar.slider("Animation Speed", 0.1, 2.0, 1.0, 0.1)

# Calculate constellation metrics
total_satellites = n_planes * sats_per_plane
period_min = orbital_period(altitude_km)
velocity_km_s = orbital_velocity(altitude_km)
orbits_day = orbits_per_day(altitude_km)

# Display metrics
st.subheader("Constellation Metrics")
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Total Satellites", total_satellites)

with col2:
    st.metric("Orbital Period", f"{period_min:.1f} min")

with col3:
    st.metric("Velocity", f"{velocity_km_s:.2f} km/s")

with col4:
    st.metric("Orbits/Day", f"{orbits_day:.1f}")


def generate_walker_constellation(
    altitude_km: float,
    inclination_deg: float,
    n_planes: int,
    sats_per_plane: int,
    phasing_factor: int,
) -> list:
    """Generate Walker Delta constellation satellite positions.

    Walker notation: i:T/P/F where:
    - i: inclination
    - T: total satellites
    - P: number of planes
    - F: phasing factor (relative spacing between planes)
    """
    satellites = []
    r = EARTH_RADIUS_KM + altitude_km

    for plane_idx in range(n_planes):
        # RAAN for this plane (evenly distributed)
        raan_deg = (360.0 / n_planes) * plane_idx

        for sat_idx in range(sats_per_plane):
            # True anomaly within plane
            ta_base = (360.0 / sats_per_plane) * sat_idx

            # Add phasing offset between planes
            phase_offset = (360.0 / (n_planes * sats_per_plane)) * phasing_factor * plane_idx
            ta_deg = (ta_base + phase_offset) % 360

            # Convert orbital elements to Cartesian (simplified for circular orbits)
            inc_rad = np.radians(inclination_deg)
            raan_rad = np.radians(raan_deg)
            ta_rad = np.radians(ta_deg)

            # Position in orbital plane
            x_orb = r * np.cos(ta_rad)
            y_orb = r * np.sin(ta_rad)

            # Rotate to ECI frame
            # First rotation: argument of latitude in orbital plane
            # Second rotation: inclination
            # Third rotation: RAAN

            x = (x_orb * (np.cos(raan_rad) * 1 - np.sin(raan_rad) * np.cos(inc_rad) * 0) +
                 y_orb * (-np.cos(raan_rad) * 0 - np.sin(raan_rad) * np.cos(inc_rad) * 1))

            y = (x_orb * (np.sin(raan_rad) * 1 + np.cos(raan_rad) * np.cos(inc_rad) * 0) +
                 y_orb * (-np.sin(raan_rad) * 0 + np.cos(raan_rad) * np.cos(inc_rad) * 1))

            z = (x_orb * np.sin(inc_rad) * 0 +
                 y_orb * np.sin(inc_rad) * 1)

            # Simplified rotation (just for visualization)
            x = r * (np.cos(raan_rad) * np.cos(ta_rad) -
                     np.sin(raan_rad) * np.sin(ta_rad) * np.cos(inc_rad))
            y = r * (np.sin(raan_rad) * np.cos(ta_rad) +
                     np.cos(raan_rad) * np.sin(ta_rad) * np.cos(inc_rad))
            z = r * np.sin(ta_rad) * np.sin(inc_rad)

            satellites.append({
                "plane": plane_idx,
                "sat_idx": sat_idx,
                "raan_deg": raan_deg,
                "ta_deg": ta_deg,
                "x": x,
                "y": y,
                "z": z,
            })

    return satellites


def generate_orbit_path(
    altitude_km: float,
    inclination_deg: float,
    raan_deg: float,
    n_points: int = 100,
) -> tuple:
    """Generate orbit path for visualization."""
    r = EARTH_RADIUS_KM + altitude_km
    inc_rad = np.radians(inclination_deg)
    raan_rad = np.radians(raan_deg)

    ta_values = np.linspace(0, 2 * np.pi, n_points)

    x = r * (np.cos(raan_rad) * np.cos(ta_values) -
             np.sin(raan_rad) * np.sin(ta_values) * np.cos(inc_rad))
    y = r * (np.sin(raan_rad) * np.cos(ta_values) +
             np.cos(raan_rad) * np.sin(ta_values) * np.cos(inc_rad))
    z = r * np.sin(ta_values) * np.sin(inc_rad)

    return x, y, z


# Generate constellation
satellites = generate_walker_constellation(
    altitude_km,
    inclination_deg,
    n_planes,
    sats_per_plane,
    phasing_factor,
)

# Create 3D visualization
fig = go.Figure()

# Add Earth sphere
if show_earth:
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_earth = EARTH_RADIUS_KM * np.outer(np.cos(u), np.sin(v))
    y_earth = EARTH_RADIUS_KM * np.outer(np.sin(u), np.sin(v))
    z_earth = EARTH_RADIUS_KM * np.outer(np.ones(50), np.cos(v))

    fig.add_trace(go.Surface(
        x=x_earth, y=y_earth, z=z_earth,
        colorscale=[[0, "#1a5276"], [0.5, "#2e86ab"], [1, "#5dade2"]],
        showscale=False,
        opacity=0.8,
        name="Earth",
        hoverinfo="skip",
    ))

# Add orbit paths
if show_orbits:
    colors = [
        "#e74c3c", "#3498db", "#2ecc71", "#f39c12",
        "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
        "#16a085", "#c0392b", "#2980b9", "#27ae60",
    ]

    for plane_idx in range(n_planes):
        raan = (360.0 / n_planes) * plane_idx
        x_orb, y_orb, z_orb = generate_orbit_path(
            altitude_km, inclination_deg, raan
        )

        fig.add_trace(go.Scatter3d(
            x=x_orb, y=y_orb, z=z_orb,
            mode="lines",
            line=dict(color=colors[plane_idx % len(colors)], width=1),
            opacity=0.5,
            name=f"Plane {plane_idx + 1}",
            hoverinfo="skip",
        ))

# Add satellites
sat_x = [s["x"] for s in satellites]
sat_y = [s["y"] for s in satellites]
sat_z = [s["z"] for s in satellites]
sat_colors = [s["plane"] for s in satellites]

fig.add_trace(go.Scatter3d(
    x=sat_x, y=sat_y, z=sat_z,
    mode="markers",
    marker=dict(
        size=5,
        color=sat_colors,
        colorscale="Viridis",
        opacity=0.9,
    ),
    text=[f"Plane {s['plane']+1}, Sat {s['sat_idx']+1}" for s in satellites],
    hoverinfo="text",
    name="Satellites",
))

# Add coverage circles (nadir-pointing footprint)
if show_coverage:
    # Earth central angle for coverage
    elevation_min = 10  # degrees minimum elevation
    r_orbit = EARTH_RADIUS_KM + altitude_km
    earth_central_angle = np.arccos(EARTH_RADIUS_KM / r_orbit) - np.radians(elevation_min)
    coverage_radius_km = EARTH_RADIUS_KM * earth_central_angle

    # Just show for a few satellites to avoid clutter
    for i, sat in enumerate(satellites[:min(10, len(satellites))]):
        # Simplified: draw circle on Earth surface
        theta = np.linspace(0, 2 * np.pi, 30)
        # This is a rough approximation
        circle_x = sat["x"] * 0.95 + coverage_radius_km * np.cos(theta) * 0.1
        circle_y = sat["y"] * 0.95 + coverage_radius_km * np.sin(theta) * 0.1
        circle_z = sat["z"] * 0.95 + np.zeros_like(theta)

        fig.add_trace(go.Scatter3d(
            x=circle_x, y=circle_y, z=circle_z,
            mode="lines",
            line=dict(color="yellow", width=1),
            opacity=0.3,
            showlegend=False,
            hoverinfo="skip",
        ))

# Configure layout
max_range = EARTH_RADIUS_KM + altitude_km + 500

fig.update_layout(
    scene=dict(
        xaxis=dict(range=[-max_range, max_range], showticklabels=False, title=""),
        yaxis=dict(range=[-max_range, max_range], showticklabels=False, title=""),
        zaxis=dict(range=[-max_range, max_range], showticklabels=False, title=""),
        aspectmode="cube",
        bgcolor="black",
    ),
    margin=dict(l=0, r=0, t=30, b=0),
    height=600,
    title=f"Walker {inclination_deg}°: {total_satellites}/{n_planes}/{phasing_factor}",
    showlegend=True,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(0,0,0,0.5)", font=dict(color="white")),
)

st.plotly_chart(fig, use_container_width=True)

# Constellation details table
st.markdown("---")
st.subheader("Constellation Details")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    **Walker Notation**: {inclination_deg}°: {total_satellites}/{n_planes}/{phasing_factor}

    **Orbital Parameters**:
    - Altitude: {altitude_km} km
    - Inclination: {inclination_deg}°
    - Semi-major axis: {EARTH_RADIUS_KM + altitude_km:.0f} km

    **Plane Separation**: {360/n_planes:.1f}° RAAN spacing

    **In-Plane Spacing**: {360/sats_per_plane:.1f}° true anomaly
    """)

with col2:
    st.markdown(f"""
    **Coverage Characteristics**:
    - Minimum elevation: 10° (typical)
    - Coverage overlap: {max(0, 100 * (n_planes * sats_per_plane / 100) - 100):.0f}%

    **Comparison to Starlink**:
    - Starlink Phase 1: 53°: 1584/72/1
    - Your design: {inclination_deg}°: {total_satellites}/{n_planes}/{phasing_factor}

    **Ground Track**: Repeats every {24*60/period_min:.1f} orbits ({24*period_min/60:.1f} hours)
    """)

# Reference constellations
st.markdown("---")
st.subheader("Reference Constellations")

ref_data = {
    "Constellation": ["Starlink Shell 1", "Starlink Shell 2", "OneWeb", "Iridium", "GPS"],
    "Satellites": [1584, 2200, 648, 66, 24],
    "Altitude (km)": [550, 540, 1200, 780, 20200],
    "Inclination (°)": [53.0, 53.2, 87.9, 86.4, 55.0],
    "Planes": [72, 72, 18, 6, 6],
    "Sats/Plane": [22, 30, 36, 11, 4],
}

import pandas as pd
df = pd.DataFrame(ref_data)
st.dataframe(df, use_container_width=True)

# Quick preset buttons
st.markdown("---")
st.subheader("Quick Presets")

col1, col2, col3, col4 = st.columns(4)

with col1:
    if st.button("Starlink-like"):
        st.info("Set: 550km, 53°, 12 planes, 22 sats/plane")

with col2:
    if st.button("Polar LEO"):
        st.info("Set: 800km, 98°, 8 planes, 10 sats/plane")

with col3:
    if st.button("Iridium-like"):
        st.info("Set: 780km, 86°, 6 planes, 11 sats/plane")

with col4:
    if st.button("GPS-like"):
        st.info("Set: 20200km, 55°, 6 planes, 4 sats/plane")

# Export
with st.expander("Export Constellation Data"):
    import json

    constellation_data = {
        "walker_notation": f"{inclination_deg}:{total_satellites}/{n_planes}/{phasing_factor}",
        "parameters": {
            "altitude_km": altitude_km,
            "inclination_deg": inclination_deg,
            "n_planes": n_planes,
            "sats_per_plane": sats_per_plane,
            "phasing_factor": phasing_factor,
        },
        "orbital_mechanics": {
            "period_min": period_min,
            "velocity_km_s": velocity_km_s,
            "orbits_per_day": orbits_day,
        },
        "satellites": satellites[:10],  # First 10 for brevity
    }

    st.json(constellation_data)
