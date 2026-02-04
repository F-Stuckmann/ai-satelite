# AI Satellite Constellation Simulator

Full-fidelity simulation platform to verify and push beyond SpaceX Patent **US 11,834,205 B1** claims, exploring theoretical limits of W/kg for AI satellite constellations.

## Overview

This project provides an interactive web-based tool to:
- **Verify** SpaceX's claimed 12.6 W/kg power-to-mass ratio
- **Explore** design concepts to push W/kg even higher
- **Simulate** thermal, structural, and orbital dynamics
- **Visualize** satellites in 3D with real-time parameter adjustment

## Patent Reference

**US 11,834,205 B1** - SpaceX Spacecraft Chassis
- Claimed: **12.6 W/kg** (3,300W / 261kg)
- Key innovation: Cellular grid fin radiators with integrated thermal/structural design

| Spacecraft | W/kg |
|------------|------|
| Boeing 702 | 3.0 |
| NASA DS1 | 5.0 |
| NASA Dawn | 8.0 |
| **SpaceX Starlink** | **12.6** |
| **Your design?** | **???** |

## Features

### Interactive Design Tool
- Parameter sliders for chassis, thermal, power, payload, and orbit
- Real-time W/kg calculation
- 3D rotating satellite visualization
- Thermal and stress overlays

### Design Concepts Library
- **Starlink Baseline**: Reference SpaceX patent design
- **Solar-on-Heatsink**: Mount solar cells directly on radiator
- **Body-Mounted Solar**: No deployable arrays
- **Graphene Radiator**: Advanced materials
- **PCM Thermal Buffer**: Phase change for burst compute

### Physics Simulation
- Stefan-Boltzmann radiative cooling
- Grid fin radiator modeling
- View factors (Earth, Sun, space)
- Eclipse thermal cycling
- FEA structural analysis

### Orbital Mechanics
- Keplerian + SGP4 propagation
- J2 and higher-order perturbations
- Constellation topology (Walker patterns)

## Installation

```bash
# Clone the repo
git clone git@github.com:F-Stuckmann/ai-satelite.git
cd ai-satelite

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .
```

## Quick Start (Local)

```bash
# Run the interactive web app
streamlit run app/main.py

# Or run example scripts
python examples/verify_spacex_patent.py
```

## Remote Access (Run on Server, Access via IP)

If you want to run the app on a server and access it from another computer:

```bash
# 1. SSH into your server
ssh user@your-server-ip

# 2. Setup (first time only)
cd ai-satelite
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .

# 3. Run accessible from network
streamlit run app/main.py --server.address 0.0.0.0 --server.port 8501
```

Then open in your browser:
```
http://<your-server-ip>:8501
```

### Run in Background

```bash
# Option 1: nohup (persists after logout)
nohup streamlit run app/main.py --server.address 0.0.0.0 > streamlit.log 2>&1 &

# Option 2: screen (can reattach later)
screen -S satellite
streamlit run app/main.py --server.address 0.0.0.0
# Detach: Ctrl+A, D
# Reattach: screen -r satellite

# Option 3: tmux
tmux new -s satellite
streamlit run app/main.py --server.address 0.0.0.0
# Detach: Ctrl+B, D
# Reattach: tmux attach -t satellite
```

### Firewall Note

Make sure port 8501 is open:
```bash
# Check if listening
lsof -i :8501

# On Linux with ufw
sudo ufw allow 8501
```

## Project Structure

```
ai-satelite/
├── app/                    # Streamlit web application
│   ├── main.py
│   ├── pages/              # Multi-page app
│   ├── components/         # UI components
│   └── concepts/           # Design concept definitions
├── src/                    # Core simulation modules
│   ├── thermal/            # Thermal physics
│   ├── structural/         # Structural analysis
│   ├── orbital/            # Orbital mechanics
│   ├── constellation/      # Multi-satellite simulation
│   ├── power/              # Power systems
│   ├── ai_payload/         # AI chip modeling
│   ├── optimization/       # W/kg optimization
│   └── visualization/      # 3D rendering
├── tests/                  # Unit tests
├── examples/               # Example scripts
└── docs/                   # Documentation
```

## Key Equations

### Radiative Cooling (Stefan-Boltzmann)
```
P = ε × σ × A × T⁴

Where:
  P = Power radiated (W)
  ε = Emissivity (0.9 for black radiator)
  σ = 5.67×10⁻⁸ W/(m²·K⁴)
  A = Radiator area (m²)
  T = Temperature (K)
```

### W/kg Calculation
```
W/kg = Total_Power / Total_Mass

Where:
  Total_Power = Solar_Power × Duty_Cycle - Losses
  Total_Mass = Structure + Thermal + Power + Payload
```

## Contributing

Contributions welcome! Please read the contributing guidelines and submit PRs.

## License

MIT License - see [LICENSE](LICENSE)

## References

- [US 11,834,205 B1](https://patents.google.com/patent/US11834205B1/en) - SpaceX Spacecraft Chassis Patent
- [Poliastro](https://docs.poliastro.space/) - Python Astrodynamics Library
- [SGP4](https://pypi.org/project/sgp4/) - Satellite Propagation
