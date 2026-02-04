"""Orbital Viewer - 3D visualization of satellite orbiting Earth."""

import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.set_page_config(page_title="Orbital Viewer", page_icon="🌍", layout="wide")

import sys
sys.path.insert(0, ".")

from src.core.constants import EARTH_RADIUS_KM


def create_sphere(radius: float, resolution: int = 50) -> tuple:
    """Create sphere mesh vertices."""
    phi = np.linspace(0, np.pi, resolution)
    theta = np.linspace(0, 2 * np.pi, resolution)
    phi, theta = np.meshgrid(phi, theta)

    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    return x, y, z


def calculate_orbit_points(
    altitude_km: float,
    inclination_deg: float,
    raan_deg: float = 0,
    num_points: int = 360,
) -> tuple:
    """Calculate orbital path points.

    Args:
        altitude_km: Orbital altitude above Earth surface
        inclination_deg: Orbital inclination in degrees
        raan_deg: Right Ascension of Ascending Node in degrees
        num_points: Number of points in orbit

    Returns:
        (x, y, z) arrays of orbit coordinates in km
    """
    radius = EARTH_RADIUS_KM + altitude_km
    inc = np.radians(inclination_deg)
    raan = np.radians(raan_deg)

    # True anomaly (position in orbit)
    nu = np.linspace(0, 2 * np.pi, num_points)

    # Position in orbital plane
    x_orbital = radius * np.cos(nu)
    y_orbital = radius * np.sin(nu)
    z_orbital = np.zeros_like(nu)

    # Rotation matrices
    # 1. Rotate by inclination around x-axis
    x_inc = x_orbital
    y_inc = y_orbital * np.cos(inc) - z_orbital * np.sin(inc)
    z_inc = y_orbital * np.sin(inc) + z_orbital * np.cos(inc)

    # 2. Rotate by RAAN around z-axis
    x_final = x_inc * np.cos(raan) - y_inc * np.sin(raan)
    y_final = x_inc * np.sin(raan) + y_inc * np.cos(raan)
    z_final = z_inc

    return x_final, y_final, z_final


def get_satellite_position(
    altitude_km: float,
    inclination_deg: float,
    raan_deg: float,
    true_anomaly_deg: float,
) -> tuple:
    """Get satellite position at specific true anomaly."""
    radius = EARTH_RADIUS_KM + altitude_km
    inc = np.radians(inclination_deg)
    raan = np.radians(raan_deg)
    nu = np.radians(true_anomaly_deg)

    # Position in orbital plane
    x_orbital = radius * np.cos(nu)
    y_orbital = radius * np.sin(nu)

    # Apply rotations
    x_inc = x_orbital
    y_inc = y_orbital * np.cos(inc)
    z_inc = y_orbital * np.sin(inc)

    x_final = x_inc * np.cos(raan) - y_inc * np.sin(raan)
    y_final = x_inc * np.sin(raan) + y_inc * np.cos(raan)
    z_final = z_inc

    return x_final, y_final, z_final


def create_earth_texture_colors(x, y, z) -> np.ndarray:
    """Create simple Earth-like coloring based on latitude."""
    # Latitude from z coordinate
    lat = np.arcsin(z / EARTH_RADIUS_KM) * 180 / np.pi

    # Simple color scheme
    colors = np.zeros((*lat.shape, 3))

    # Ocean (blue) - default
    colors[:, :, 0] = 30  # R
    colors[:, :, 1] = 100  # G
    colors[:, :, 2] = 180  # B

    # Land masses (green) - simple approximation
    lon = np.arctan2(y, x) * 180 / np.pi

    # Add some "continents" based on longitude bands
    land_mask = (
        ((lon > -80) & (lon < -30) & (lat > -60) & (lat < 15)) |  # South America
        ((lon > -130) & (lon < -60) & (lat > 15) & (lat < 70)) |  # North America
        ((lon > -20) & (lon < 50) & (lat > -40) & (lat < 70)) |   # Europe/Africa
        ((lon > 60) & (lon < 150) & (lat > -50) & (lat < 70))     # Asia/Australia
    )

    colors[land_mask, 0] = 50   # R
    colors[land_mask, 1] = 150  # G
    colors[land_mask, 2] = 50   # B

    # Ice caps (white)
    ice_mask = (lat > 70) | (lat < -70)
    colors[ice_mask, 0] = 240
    colors[ice_mask, 1] = 240
    colors[ice_mask, 2] = 250

    return colors


def create_orbital_figure(
    altitude_km: float,
    inclination_deg: float,
    satellite_position_deg: float,
    show_orbit_path: bool = True,
    num_satellites: int = 1,
    satellites_per_plane: int = 1,
) -> go.Figure:
    """Create 3D figure with Earth and orbiting satellite(s)."""

    fig = go.Figure()

    # Create Earth sphere
    earth_x, earth_y, earth_z = create_sphere(EARTH_RADIUS_KM, resolution=60)
    colors = create_earth_texture_colors(earth_x, earth_y, earth_z)

    # Convert colors to plotly format
    surfacecolor = (
        colors[:, :, 0] * 65536 +
        colors[:, :, 1] * 256 +
        colors[:, :, 2]
    )

    fig.add_trace(go.Surface(
        x=earth_x,
        y=earth_y,
        z=earth_z,
        surfacecolor=surfacecolor,
        colorscale=[
            [0, 'rgb(30,100,180)'],    # Ocean
            [0.3, 'rgb(50,150,50)'],   # Land
            [0.7, 'rgb(50,150,50)'],   # Land
            [1, 'rgb(240,240,250)'],   # Ice
        ],
        showscale=False,
        name="Earth",
        hoverinfo="name",
        opacity=1.0,
    ))

    # Add atmosphere glow
    atm_x, atm_y, atm_z = create_sphere(EARTH_RADIUS_KM * 1.02, resolution=30)
    fig.add_trace(go.Surface(
        x=atm_x,
        y=atm_y,
        z=atm_z,
        colorscale=[[0, 'rgba(135,206,250,0.1)'], [1, 'rgba(135,206,250,0.1)']],
        showscale=False,
        name="Atmosphere",
        hoverinfo="skip",
        opacity=0.3,
    ))

    # Calculate orbital planes and satellites
    raan_spacing = 360 / max(1, num_satellites // satellites_per_plane) if num_satellites > satellites_per_plane else 0

    for plane_idx in range(max(1, num_satellites // satellites_per_plane)):
        raan = plane_idx * raan_spacing

        # Add orbit path
        if show_orbit_path:
            orbit_x, orbit_y, orbit_z = calculate_orbit_points(
                altitude_km, inclination_deg, raan
            )
            fig.add_trace(go.Scatter3d(
                x=orbit_x,
                y=orbit_y,
                z=orbit_z,
                mode="lines",
                line=dict(color="yellow", width=2),
                name=f"Orbit Plane {plane_idx + 1}" if num_satellites > satellites_per_plane else "Orbit",
                hoverinfo="name",
            ))

        # Add satellites in this plane
        for sat_idx in range(satellites_per_plane):
            if plane_idx * satellites_per_plane + sat_idx >= num_satellites:
                break

            # Spread satellites evenly in the plane
            sat_anomaly = satellite_position_deg + (sat_idx * 360 / satellites_per_plane)

            sat_x, sat_y, sat_z = get_satellite_position(
                altitude_km, inclination_deg, raan, sat_anomaly
            )

            # Satellite marker
            fig.add_trace(go.Scatter3d(
                x=[sat_x],
                y=[sat_y],
                z=[sat_z],
                mode="markers",
                marker=dict(
                    size=8,
                    color="red",
                    symbol="diamond",
                ),
                name=f"Satellite {plane_idx * satellites_per_plane + sat_idx + 1}",
                hovertemplate=(
                    f"Satellite {plane_idx * satellites_per_plane + sat_idx + 1}<br>"
                    f"Altitude: {altitude_km:.0f} km<br>"
                    f"<extra></extra>"
                ),
            ))

    # Camera and layout
    orbital_radius = EARTH_RADIUS_KM + altitude_km
    camera_distance = orbital_radius * 2.5

    fig.update_layout(
        scene=dict(
            xaxis=dict(
                visible=False,
                range=[-camera_distance, camera_distance],
            ),
            yaxis=dict(
                visible=False,
                range=[-camera_distance, camera_distance],
            ),
            zaxis=dict(
                visible=False,
                range=[-camera_distance, camera_distance],
            ),
            aspectmode="cube",
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=0.8),
            ),
            bgcolor="black",
        ),
        paper_bgcolor="black",
        title=dict(
            text=f"Satellite Orbit - {altitude_km:.0f} km altitude, {inclination_deg:.0f}° inclination",
            font=dict(color="white"),
        ),
        height=700,
        margin=dict(l=0, r=0, t=40, b=0),
        legend=dict(
            font=dict(color="white"),
            bgcolor="rgba(0,0,0,0.5)",
        ),
    )

    return fig


def main():
    st.title("🌍 Orbital Viewer")
    st.markdown("**Watch your satellite orbit Earth in 3D**")

    # Sidebar controls
    st.sidebar.header("Orbit Parameters")

    altitude_km = st.sidebar.slider(
        "Altitude (km)",
        min_value=200,
        max_value=2000,
        value=550,
        step=10,
        help="Height above Earth's surface",
    )

    inclination_deg = st.sidebar.slider(
        "Inclination (°)",
        min_value=0,
        max_value=98,
        value=53,
        step=1,
        help="Angle of orbit relative to equator",
    )

    st.sidebar.header("Visualization")

    show_orbit_path = st.sidebar.checkbox("Show orbit path", value=True)

    satellite_position = st.sidebar.slider(
        "Satellite Position (°)",
        min_value=0,
        max_value=359,
        value=0,
        step=5,
        help="Position along orbit (true anomaly)",
    )

    st.sidebar.header("Constellation")

    num_satellites = st.sidebar.slider(
        "Number of Satellites",
        min_value=1,
        max_value=24,
        value=1,
        step=1,
    )

    if num_satellites > 1:
        sats_per_plane = st.sidebar.slider(
            "Satellites per Plane",
            min_value=1,
            max_value=min(12, num_satellites),
            value=min(4, num_satellites),
            step=1,
        )
    else:
        sats_per_plane = 1

    # Create and display figure
    fig = create_orbital_figure(
        altitude_km=altitude_km,
        inclination_deg=inclination_deg,
        satellite_position_deg=satellite_position,
        show_orbit_path=show_orbit_path,
        num_satellites=num_satellites,
        satellites_per_plane=sats_per_plane,
    )

    st.plotly_chart(fig, use_container_width=True)

    # Orbit info
    col1, col2, col3 = st.columns(3)

    # Calculate orbital period
    orbital_radius = EARTH_RADIUS_KM + altitude_km
    mu = 398600.4418  # km³/s²
    period_seconds = 2 * np.pi * np.sqrt(orbital_radius**3 / mu)
    period_minutes = period_seconds / 60

    # Orbital velocity
    velocity = np.sqrt(mu / orbital_radius)  # km/s

    with col1:
        st.metric("Orbital Period", f"{period_minutes:.1f} min")

    with col2:
        st.metric("Orbital Velocity", f"{velocity:.2f} km/s")

    with col3:
        st.metric("Orbits per Day", f"{24 * 60 / period_minutes:.1f}")

    # Orbit type classification
    st.markdown("---")
    if altitude_km < 450:
        orbit_type = "Very Low Earth Orbit (VLEO)"
        note = "High drag, frequent reboosting needed"
    elif altitude_km < 1000:
        orbit_type = "Low Earth Orbit (LEO)"
        note = "Typical for Starlink, ISS"
    elif altitude_km < 2000:
        orbit_type = "LEO (upper)"
        note = "Lower drag, longer orbital lifetime"
    else:
        orbit_type = "Medium Earth Orbit (MEO)"
        note = "GPS satellites operate here"

    st.info(f"**{orbit_type}**: {note}")

    # Animation tip
    with st.expander("💡 Animation Tip"):
        st.markdown("""
        Use the **Satellite Position** slider to manually animate the satellite around its orbit.

        For auto-rotation of the 3D view:
        - Click and drag the Earth to rotate
        - Scroll to zoom in/out
        - Double-click to reset view
        """)


if __name__ == "__main__":
    main()
