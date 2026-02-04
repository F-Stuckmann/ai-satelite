"""AI Satellite Designer - Main Streamlit Application.

Run with: streamlit run app/main.py
"""

import streamlit as st

# Page config must be first Streamlit command
st.set_page_config(
    page_title="AI Satellite Designer",
    page_icon="🛰️",
    layout="wide",
    initial_sidebar_state="expanded",
)


def main():
    """Main application entry point."""
    st.title("🛰️ AI Satellite Constellation Simulator")

    st.markdown("""
    **Verify SpaceX Patent US 11,834,205 B1 and push W/kg limits**

    This tool lets you interactively design AI satellites and explore
    the theoretical limits of power-to-mass ratio (W/kg).

    ---

    ### Quick Reference: W/kg Benchmarks

    | Spacecraft | W/kg |
    |------------|------|
    | Boeing 702 | 3.0 |
    | NASA DS1 | 5.0 |
    | NASA Dawn | 8.0 |
    | **SpaceX Starlink (Patent)** | **12.6** |
    | Your Design | **???** |

    ---

    ### Getting Started

    Use the sidebar to navigate between pages:

    1. **🎛️ Satellite Designer** - Main design interface with parameter sliders
    2. **💡 Concept Explorer** - Pre-built design concepts to explore
    3. **🌍 Constellation Sim** - Multi-satellite simulation
    4. **📈 W/kg Optimizer** - Automated optimization

    ---
    """)

    # Quick demo calculation
    st.subheader("Quick Thermal Calculator")

    col1, col2 = st.columns(2)

    with col1:
        power_w = st.slider("Heat to Dissipate (W)", 100, 5000, 700, 50)
        temp_c = st.slider("Radiator Temperature (°C)", 20, 100, 50, 5)

    with col2:
        emissivity = st.slider("Emissivity", 0.5, 0.98, 0.9, 0.01)

    # Calculate required radiator area
    from src.thermal.stefan_boltzmann import required_radiator_area, temperature_from_celsius

    temp_k = temperature_from_celsius(temp_c)
    area_m2 = required_radiator_area(power_w, temp_k, emissivity)

    st.metric(
        label="Required Radiator Area",
        value=f"{area_m2:.2f} m²",
        delta=f"at {temp_c}°C",
    )

    # Show Stefan-Boltzmann equation
    with st.expander("📐 Stefan-Boltzmann Equation"):
        st.latex(r"P = \varepsilon \cdot \sigma \cdot A \cdot T^4")
        st.markdown("""
        Where:
        - P = Power radiated (W)
        - ε = Emissivity (0-1)
        - σ = 5.67×10⁻⁸ W/(m²·K⁴)
        - A = Area (m²)
        - T = Temperature (K)
        """)


if __name__ == "__main__":
    main()
