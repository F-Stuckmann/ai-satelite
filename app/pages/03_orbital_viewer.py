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


def point_in_polygon(lat: float, lon: float, polygon: list) -> bool:
    """Check if a point is inside a polygon using ray casting."""
    n = len(polygon)
    inside = False
    j = n - 1
    for i in range(n):
        if ((polygon[i][0] > lat) != (polygon[j][0] > lat)) and \
           (lon < (polygon[j][1] - polygon[i][1]) * (lat - polygon[i][0]) /
            (polygon[j][0] - polygon[i][0]) + polygon[i][1]):
            inside = not inside
        j = i
    return inside


def is_land(lat: float, lon: float) -> bool:
    """Check if coordinates are on land using simplified continent polygons."""

    # Simplified continent boundaries (lat, lon pairs)
    # Format: list of (lat, lon) vertices for each continent

    continents = {
        # North America
        "north_america": [
            (49, -125), (60, -140), (70, -140), (72, -95), (70, -60),
            (50, -55), (45, -65), (30, -80), (25, -80), (20, -105),
            (30, -115), (35, -120), (49, -125)
        ],
        # South America
        "south_america": [
            (12, -70), (10, -75), (0, -80), (-5, -80), (-15, -75),
            (-23, -70), (-40, -65), (-55, -68), (-55, -63), (-40, -57),
            (-35, -55), (-25, -48), (-20, -40), (-5, -35), (5, -50),
            (10, -62), (12, -70)
        ],
        # Europe
        "europe": [
            (36, -10), (43, -10), (48, -5), (50, 2), (52, 5),
            (55, 8), (58, 10), (63, 10), (70, 20), (70, 30),
            (65, 40), (55, 40), (50, 40), (45, 35), (40, 25),
            (36, 25), (36, -10)
        ],
        # Africa
        "africa": [
            (37, -10), (35, 10), (30, 32), (22, 37), (12, 44),
            (0, 42), (-12, 40), (-25, 35), (-35, 20), (-35, 18),
            (-30, 27), (-20, 35), (-10, 40), (5, 10), (5, 0),
            (10, -15), (15, -17), (20, -17), (28, -13), (37, -10)
        ],
        # Asia (simplified)
        "asia": [
            (70, 30), (75, 70), (75, 100), (70, 140), (65, 170),
            (55, 165), (50, 140), (40, 130), (35, 130), (30, 120),
            (20, 110), (10, 105), (0, 100), (-8, 115), (5, 120),
            (20, 120), (25, 90), (20, 88), (22, 70), (25, 62),
            (30, 50), (40, 45), (45, 40), (55, 40), (65, 40), (70, 30)
        ],
        # Australia
        "australia": [
            (-12, 130), (-12, 142), (-18, 146), (-25, 153), (-35, 151),
            (-38, 147), (-38, 140), (-35, 135), (-32, 130), (-30, 115),
            (-23, 114), (-15, 125), (-12, 130)
        ],
        # Greenland
        "greenland": [
            (60, -45), (65, -55), (70, -55), (78, -70), (83, -40),
            (83, -30), (78, -20), (70, -22), (65, -40), (60, -45)
        ],
        # Antarctica (simplified)
        "antarctica": [
            (-65, -60), (-70, -60), (-75, -70), (-80, -90), (-85, -120),
            (-80, -150), (-75, 170), (-70, 160), (-68, 140), (-66, 100),
            (-68, 70), (-70, 30), (-75, 0), (-70, -30), (-65, -60)
        ],
        # Madagascar
        "madagascar": [
            (-12, 49), (-16, 50), (-24, 47), (-26, 45), (-22, 43),
            (-16, 45), (-12, 49)
        ],
        # UK/Ireland
        "uk": [
            (50, -6), (52, -10), (55, -8), (58, -6), (59, -3),
            (56, -2), (54, 0), (52, 2), (51, 1), (50, -6)
        ],
        # Japan
        "japan": [
            (31, 130), (33, 130), (35, 135), (38, 140), (42, 145),
            (45, 145), (44, 142), (40, 140), (35, 137), (33, 133),
            (31, 130)
        ],
        # New Zealand
        "new_zealand": [
            (-35, 173), (-38, 178), (-42, 177), (-47, 167), (-45, 167),
            (-41, 173), (-37, 175), (-35, 173)
        ],
        # Indonesia (simplified)
        "indonesia": [
            (-5, 95), (-6, 105), (-7, 110), (-8, 115), (-8, 120),
            (-5, 120), (-2, 115), (2, 110), (2, 105), (0, 100),
            (-5, 95)
        ],
        # Central America
        "central_america": [
            (20, -105), (18, -95), (15, -90), (10, -85), (8, -80),
            (10, -77), (15, -83), (18, -88), (22, -98), (20, -105)
        ],
    }

    for name, polygon in continents.items():
        if point_in_polygon(lat, lon, polygon):
            return True
    return False


def create_earth_texture_colors(x, y, z) -> np.ndarray:
    """Create Earth-like coloring with realistic continents."""
    # Calculate lat/lon from coordinates
    lat = np.arcsin(np.clip(z / EARTH_RADIUS_KM, -1, 1)) * 180 / np.pi
    lon = np.arctan2(y, x) * 180 / np.pi

    # Initialize colors array
    colors = np.zeros((*lat.shape, 3), dtype=np.uint8)

    # Vectorized ocean coloring with depth variation
    depth_factor = 0.8 + 0.2 * np.abs(lat) / 90  # Lighter near poles
    colors[:, :, 0] = (20 * depth_factor).astype(np.uint8)   # R
    colors[:, :, 1] = (80 * depth_factor).astype(np.uint8)   # G
    colors[:, :, 2] = (160 * depth_factor).astype(np.uint8)  # B

    # Check each point for land (vectorized where possible)
    for i in range(lat.shape[0]):
        for j in range(lat.shape[1]):
            if is_land(lat[i, j], lon[i, j]):
                # Vary green based on latitude (tropical vs temperate)
                lat_abs = abs(lat[i, j])
                if lat_abs < 25:
                    # Tropical - darker green
                    colors[i, j] = [34, 120, 34]
                elif lat_abs < 45:
                    # Temperate - medium green
                    colors[i, j] = [50, 140, 50]
                elif lat_abs < 60:
                    # Northern - lighter/browner
                    colors[i, j] = [80, 130, 60]
                else:
                    # Tundra - grayish
                    colors[i, j] = [140, 150, 130]

                # Add some desert coloring near tropics
                if 15 < lat_abs < 35:
                    # Check for desert regions (simplified)
                    if (20 < lon[i, j] < 50 and 15 < lat[i, j] < 35) or \
                       (-120 < lon[i, j] < -100 and 25 < lat[i, j] < 40) or \
                       (40 < lon[i, j] < 80 and 20 < lat[i, j] < 40) or \
                       (120 < lon[i, j] < 140 and -30 < lat[i, j] < -20):
                        colors[i, j] = [194, 178, 128]  # Sandy color

    # Ice caps and snow (overwrite land colors near poles)
    ice_mask = (lat > 75) | (lat < -65)
    colors[ice_mask] = [250, 250, 255]

    # Partial ice/snow
    partial_ice = ((lat > 65) & (lat <= 75)) | ((lat < -55) & (lat >= -65))
    colors[partial_ice, 0] = np.clip(colors[partial_ice, 0] + 80, 0, 255)
    colors[partial_ice, 1] = np.clip(colors[partial_ice, 1] + 80, 0, 255)
    colors[partial_ice, 2] = np.clip(colors[partial_ice, 2] + 80, 0, 255)

    return colors


def create_orbital_figure(
    altitude_km: float,
    inclination_deg: float,
    satellite_position_deg: float,
    show_orbit_path: bool = True,
    num_satellites: int = 1,
    satellites_per_plane: int = 1,
    animate: bool = False,
    animation_speed: int = 5,
) -> go.Figure:
    """Create 3D figure with Earth and orbiting satellite(s)."""

    fig = go.Figure()

    # Create Earth sphere
    earth_x, earth_y, earth_z = create_sphere(EARTH_RADIUS_KM, resolution=60)
    colors = create_earth_texture_colors(earth_x, earth_y, earth_z)

    # Convert colors to plotly format (convert to int32 to avoid overflow)
    surfacecolor = (
        colors[:, :, 0].astype(np.int32) * 65536 +
        colors[:, :, 1].astype(np.int32) * 256 +
        colors[:, :, 2].astype(np.int32)
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

    # Track trace indices for animation
    trace_count = 2  # Earth + Atmosphere
    satellite_trace_indices = []

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
            trace_count += 1

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
            satellite_trace_indices.append(trace_count)
            trace_count += 1

    # Camera and layout
    orbital_radius = EARTH_RADIUS_KM + altitude_km
    camera_distance = orbital_radius * 2.5

    # Create animation frames if animation is enabled
    if animate and num_satellites <= 4 and satellite_trace_indices:
        frames = []
        n_frames = 72  # 72 frames for full orbit (5° per frame)

        for frame_idx in range(n_frames):
            frame_position = (satellite_position_deg + frame_idx * 5) % 360
            frame_data = []

            # Recalculate satellite positions for this frame
            raan_spacing_anim = 360 / max(1, num_satellites // satellites_per_plane) if num_satellites > satellites_per_plane else 0

            for plane_idx in range(max(1, num_satellites // satellites_per_plane)):
                raan = plane_idx * raan_spacing_anim

                for sat_idx in range(satellites_per_plane):
                    if plane_idx * satellites_per_plane + sat_idx >= num_satellites:
                        break

                    sat_anomaly = frame_position + (sat_idx * 360 / satellites_per_plane)
                    sat_x, sat_y, sat_z = get_satellite_position(
                        altitude_km, inclination_deg, raan, sat_anomaly
                    )

                    frame_data.append(go.Scatter3d(
                        x=[sat_x],
                        y=[sat_y],
                        z=[sat_z],
                        mode="markers",
                        marker=dict(size=8, color="red", symbol="diamond"),
                    ))

            # Only update satellite traces, keep Earth and orbit paths
            frames.append(go.Frame(
                data=frame_data,
                name=str(frame_idx),
                traces=satellite_trace_indices,  # Only animate satellite traces
            ))

        fig.frames = frames

        # Add play/pause buttons
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    showactive=False,
                    y=0.1,
                    x=0.05,
                    xanchor="left",
                    buttons=[
                        dict(
                            label="▶ Play",
                            method="animate",
                            args=[
                                None,
                                dict(
                                    frame=dict(duration=80, redraw=True),
                                    fromcurrent=True,
                                    mode="immediate",
                                )
                            ],
                        ),
                        dict(
                            label="⏸ Pause",
                            method="animate",
                            args=[
                                [None],
                                dict(
                                    frame=dict(duration=0, redraw=True),
                                    mode="immediate",
                                )
                            ],
                        ),
                    ],
                    font=dict(color="white"),
                    bgcolor="rgba(50,50,50,0.8)",
                )
            ],
            sliders=[
                dict(
                    active=0,
                    yanchor="top",
                    xanchor="left",
                    currentvalue=dict(
                        font=dict(size=12, color="white"),
                        prefix="Position: ",
                        suffix="°",
                        visible=True,
                        xanchor="center",
                    ),
                    transition=dict(duration=0),
                    pad=dict(b=10, t=50),
                    len=0.9,
                    x=0.05,
                    y=0,
                    steps=[
                        dict(
                            args=[
                                [str(i)],
                                dict(
                                    frame=dict(duration=0, redraw=True),
                                    mode="immediate",
                                )
                            ],
                            label=str(i * 5),
                            method="animate",
                        )
                        for i in range(n_frames)
                    ],
                    font=dict(color="white"),
                    bgcolor="rgba(50,50,50,0.5)",
                )
            ],
        )

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

    # Animation toggle
    enable_animation = st.sidebar.checkbox(
        "🎬 Enable Animation",
        value=True,
        help="Enable Play/Pause controls in the 3D view"
    )

    satellite_position = st.sidebar.slider(
        "Initial Position (°)",
        min_value=0,
        max_value=359,
        value=0,
        step=5,
        help="Starting position along orbit",
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
        animate=enable_animation,
    )

    st.plotly_chart(fig, use_container_width=True)

    # Orbit info
    col1, col2, col3, col4 = st.columns(4)

    # Calculate orbital period
    orbital_radius = EARTH_RADIUS_KM + altitude_km
    mu = 398600.4418  # km³/s²
    period_seconds = 2 * np.pi * np.sqrt(orbital_radius**3 / mu)
    period_minutes = period_seconds / 60

    # Orbital velocity
    velocity_km_s = np.sqrt(mu / orbital_radius)  # km/s
    velocity_km_h = velocity_km_s * 3600  # km/h
    velocity_mach = velocity_km_s * 1000 / 343  # Mach number (vs sea level sound)

    with col1:
        st.metric("Orbital Period", f"{period_minutes:.1f} min")

    with col2:
        st.metric(
            "Velocity",
            f"{velocity_km_s:.2f} km/s",
            f"{velocity_km_h:,.0f} km/h",
        )

    with col3:
        st.metric("Speed (Mach)", f"{velocity_mach:.0f}x")

    with col4:
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

    # Ground track visualization
    st.markdown("---")
    st.subheader("🗺️ Ground Track")

    # Calculate ground track (subsatellite points)
    # Earth rotates 360° in ~24 hours = 0.25°/min
    earth_rotation_rate = 360 / (24 * 60)  # deg/min

    ground_track_lons = []
    ground_track_lats = []

    for minute in range(int(period_minutes * 3)):  # 3 orbits
        # Satellite position in orbit
        orbit_angle = (minute / period_minutes) * 360  # degrees
        orbit_angle_rad = np.radians(orbit_angle)

        # Calculate latitude (depends on inclination and position in orbit)
        lat = np.degrees(np.arcsin(np.sin(np.radians(inclination_deg)) * np.sin(orbit_angle_rad)))

        # Calculate longitude (RAAN + in-orbit position - Earth rotation)
        # Simplified: assume RAAN=0
        lon_in_orbit = np.degrees(np.arctan2(
            np.cos(np.radians(inclination_deg)) * np.sin(orbit_angle_rad),
            np.cos(orbit_angle_rad)
        ))
        earth_rotation = earth_rotation_rate * minute
        lon = (lon_in_orbit - earth_rotation + 180) % 360 - 180

        ground_track_lats.append(lat)
        ground_track_lons.append(lon)

    # Create ground track plot
    fig_ground = go.Figure()

    # Add world map outline (simplified coastlines)
    fig_ground.add_trace(go.Scattergeo(
        lon=ground_track_lons,
        lat=ground_track_lats,
        mode="lines",
        line=dict(color="red", width=2),
        name="Ground Track",
    ))

    # Current satellite position on ground track
    current_idx = int((satellite_position / 360) * len(ground_track_lons) / 3) % len(ground_track_lons)
    fig_ground.add_trace(go.Scattergeo(
        lon=[ground_track_lons[current_idx]],
        lat=[ground_track_lats[current_idx]],
        mode="markers",
        marker=dict(size=12, color="yellow", symbol="star"),
        name="Current Position",
    ))

    fig_ground.update_geos(
        showland=True,
        landcolor="rgb(50, 150, 50)",
        showocean=True,
        oceancolor="rgb(30, 100, 180)",
        showcountries=True,
        countrycolor="rgb(100, 100, 100)",
        showcoastlines=True,
        coastlinecolor="rgb(80, 80, 80)",
        projection_type="natural earth",
    )

    fig_ground.update_layout(
        title=f"Ground Track (3 orbits, {3 * period_minutes:.0f} min)",
        height=400,
        margin=dict(l=0, r=0, t=40, b=0),
    )

    st.plotly_chart(fig_ground, use_container_width=True)

    # Tips
    with st.expander("💡 Tips"):
        st.markdown("""
        **Animation Controls (in 3D view):**
        - Click **▶ Play** button below the Earth to start animation
        - Click **⏸ Pause** to stop
        - Use the slider to scrub through the orbit manually

        **3D View Controls:**
        - Click and drag to rotate the view
        - Scroll to zoom in/out
        - Double-click to reset view

        **Constellation Mode:**
        - Increase "Number of Satellites" to simulate a constellation
        - Satellites are distributed across orbital planes
        - Animation works best with 1-4 satellites
        """)


if __name__ == "__main__":
    main()
