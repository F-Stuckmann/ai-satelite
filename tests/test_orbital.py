"""Tests for orbital mechanics module."""

import pytest
import numpy as np

from src.orbital.keplerian import (
    OrbitalElements,
    orbital_velocity,
    orbital_period,
    orbits_per_day,
    ground_track_shift,
    eclipse_duration,
    delta_v_hohmann,
    delta_v_plane_change,
    atmospheric_drag_lifetime,
    summary,
)


class TestOrbitalElements:
    """Tests for OrbitalElements dataclass."""

    def test_from_altitude_default(self):
        """Test creating orbit from altitude with defaults."""
        orbit = OrbitalElements.from_altitude(400)

        # Semi-major axis should be Earth radius + altitude
        assert orbit.a == pytest.approx(6371 + 400, rel=0.01)
        assert orbit.e == 0.0
        assert orbit.i == 0.0

    def test_from_altitude_with_inclination(self):
        """Test creating orbit with inclination."""
        orbit = OrbitalElements.from_altitude(550, inclination_deg=53.0)

        assert orbit.a == pytest.approx(6371 + 550, rel=0.01)
        assert orbit.i == 53.0

    def test_period_iss_orbit(self):
        """Test orbital period for ISS-like orbit (~400 km)."""
        orbit = OrbitalElements.from_altitude(400)

        # ISS period is ~92.6 minutes
        assert orbit.period_minutes == pytest.approx(92.6, rel=0.02)

    def test_velocity_circular(self):
        """Test circular orbital velocity."""
        orbit = OrbitalElements.from_altitude(400)

        # LEO velocity is ~7.66 km/s
        assert orbit.velocity_circular == pytest.approx(7.66, rel=0.02)

    def test_specific_energy_negative(self):
        """Test that bound orbits have negative specific energy."""
        orbit = OrbitalElements.from_altitude(400)

        assert orbit.specific_energy < 0

    def test_elliptical_orbit_velocities(self):
        """Test periapsis > apoapsis velocity for elliptical orbit."""
        orbit = OrbitalElements(a=10000, e=0.2)

        assert orbit.velocity_periapsis > orbit.velocity_apoapsis

    def test_altitude_periapsis_apoapsis(self):
        """Test periapsis < apoapsis altitude for elliptical orbit."""
        orbit = OrbitalElements(a=10000, e=0.2)

        assert orbit.altitude_periapsis < orbit.altitude_apoapsis


class TestOrbitalFunctions:
    """Tests for standalone orbital functions."""

    def test_orbital_velocity_decreases_with_altitude(self):
        """Test that orbital velocity decreases with altitude."""
        v_low = orbital_velocity(400)
        v_high = orbital_velocity(800)

        assert v_low > v_high

    def test_orbital_period_increases_with_altitude(self):
        """Test that orbital period increases with altitude."""
        p_low = orbital_period(400)
        p_high = orbital_period(800)

        assert p_low < p_high

    def test_orbits_per_day_leo(self):
        """Test orbits per day for LEO satellite."""
        # ISS makes ~15.5 orbits per day
        opd = orbits_per_day(400)

        assert opd == pytest.approx(15.5, rel=0.05)

    def test_ground_track_shift_positive(self):
        """Test that ground track shifts westward (positive degrees)."""
        shift = ground_track_shift(400)

        assert shift > 0
        # Should be roughly 22-23 degrees for LEO
        assert 20 < shift < 30


class TestEclipse:
    """Tests for eclipse calculations."""

    def test_eclipse_duration_positive(self):
        """Test eclipse duration is positive for LEO."""
        duration = eclipse_duration(400)

        # LEO eclipse is typically 30-40 minutes
        assert duration > 0
        assert duration < 60

    def test_eclipse_duration_increases_with_altitude(self):
        """Test eclipse duration generally increases with altitude (up to a point)."""
        dur_low = eclipse_duration(400)
        dur_mid = eclipse_duration(600)

        # This relationship holds for low-altitude range
        assert dur_mid > dur_low * 0.9  # Allow some variation

    def test_no_eclipse_high_beta(self):
        """Test no eclipse at high beta angles."""
        # High beta angle (near summer/winter solstice for sun-sync)
        duration = eclipse_duration(400, beta_angle_deg=70)

        assert duration == 0.0


class TestDeltaV:
    """Tests for delta-v calculations."""

    def test_hohmann_transfer_positive(self):
        """Test Hohmann transfer delta-v values are positive."""
        r1 = 6771  # 400 km altitude
        r2 = 7171  # 800 km altitude

        dv1, dv2 = delta_v_hohmann(r1, r2)

        assert dv1 > 0
        assert dv2 > 0

    def test_hohmann_leo_to_geo(self):
        """Test Hohmann transfer from LEO to GEO."""
        r1 = 6771  # 400 km LEO
        r2 = 42164  # GEO

        dv1, dv2 = delta_v_hohmann(r1, r2)

        # Total delta-v to GEO is ~3.9 km/s
        total_dv = dv1 + dv2
        assert total_dv == pytest.approx(3.9, rel=0.1)

    def test_plane_change_zero_angle(self):
        """Test zero plane change requires zero delta-v."""
        dv = delta_v_plane_change(7.0, 0.0)

        assert dv == pytest.approx(0.0, abs=1e-10)

    def test_plane_change_90_deg(self):
        """Test 90 degree plane change requires ~sqrt(2) * velocity."""
        v = 7.0
        dv = delta_v_plane_change(v, 90.0)

        # For 90 degree change: dv = 2 * v * sin(45) = sqrt(2) * v
        expected = np.sqrt(2) * v
        assert dv == pytest.approx(expected, rel=0.01)


class TestAtmosphericDrag:
    """Tests for atmospheric drag lifetime estimation."""

    def test_lifetime_increases_with_altitude(self):
        """Test orbital lifetime increases with altitude."""
        bc = 50  # ballistic coefficient

        lt_low = atmospheric_drag_lifetime(300, bc)
        lt_mid = atmospheric_drag_lifetime(500, bc)
        lt_high = atmospheric_drag_lifetime(700, bc)

        assert lt_low < lt_mid < lt_high

    def test_lifetime_increases_with_bc(self):
        """Test lifetime increases with ballistic coefficient."""
        alt = 400

        lt_low = atmospheric_drag_lifetime(alt, 25)
        lt_high = atmospheric_drag_lifetime(alt, 100)

        assert lt_low < lt_high

    def test_lifetime_decreases_high_solar_activity(self):
        """Test lifetime decreases with high solar activity."""
        alt = 400
        bc = 50

        lt_low = atmospheric_drag_lifetime(alt, bc, "low")
        lt_high = atmospheric_drag_lifetime(alt, bc, "high")

        assert lt_high < lt_low


class TestSummary:
    """Tests for orbital summary function."""

    def test_summary_contains_required_keys(self):
        """Test summary returns all required keys."""
        result = summary(400)

        required_keys = [
            "altitude_km",
            "inclination_deg",
            "semi_major_axis_km",
            "orbital_period_min",
            "velocity_km_s",
            "velocity_km_h",
            "orbits_per_day",
            "ground_track_shift_deg",
            "eclipse_duration_min",
        ]

        for key in required_keys:
            assert key in result

    def test_summary_altitude_matches(self):
        """Test summary altitude matches input."""
        result = summary(550)

        assert result["altitude_km"] == 550

    def test_summary_inclination_matches(self):
        """Test summary inclination matches input."""
        result = summary(400, inclination_deg=53.0)

        assert result["inclination_deg"] == 53.0
