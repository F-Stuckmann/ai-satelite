"""Tests for thermal physics module."""

import pytest
import sys
sys.path.insert(0, ".")

from src.thermal.stefan_boltzmann import (
    radiative_power,
    equilibrium_temperature,
    required_radiator_area,
    temperature_from_celsius,
    verify_spacex_claim,
)
from src.thermal.grid_fin_radiator import GridFinRadiator, create_starlink_radiator
from src.thermal.view_factors import ViewFactors
from src.core.constants import STEFAN_BOLTZMANN


class TestStefanBoltzmann:
    """Tests for Stefan-Boltzmann calculations."""

    def test_radiative_power_basic(self):
        """Test basic radiative power calculation."""
        # At 300K, 1 m², ε=1.0
        power = radiative_power(300, 1.0, 1.0, 0)
        expected = STEFAN_BOLTZMANN * 300**4
        assert abs(power - expected) < 1

    def test_radiative_power_with_sink(self):
        """Test radiative power with non-zero sink temperature."""
        power = radiative_power(350, 2.0, 0.9, 3)
        assert power > 0

    def test_equilibrium_temperature(self):
        """Test equilibrium temperature calculation."""
        # 1000W into 1m² should give certain temperature
        t_eq = equilibrium_temperature(1000, 1.0, 0.9)
        assert 350 < t_eq < 400  # Roughly 77-127°C

    def test_required_area_inverse_of_power(self):
        """Test that required_area is inverse of radiative_power."""
        power = 500
        temp = 350
        area = required_radiator_area(power, temp, 0.9)

        # Verify by calculating power with that area
        power_check = radiative_power(temp, area, 0.9)
        assert abs(power - power_check) < 1

    def test_temperature_conversion(self):
        """Test Celsius to Kelvin conversion."""
        assert temperature_from_celsius(0) == 273.15
        assert temperature_from_celsius(100) == 373.15

    def test_verify_spacex_claim(self):
        """Verify SpaceX patent claim calculation."""
        result = verify_spacex_claim()
        assert result["verified"] is True
        assert abs(result["calculated_w_per_kg"] - 12.64) < 0.1


class TestGridFinRadiator:
    """Tests for grid fin radiator model."""

    def test_area_enhancement(self):
        """Grid fins should enhance effective area."""
        radiator = GridFinRadiator(base_area_m2=1.0)
        assert radiator.effective_area_m2 > radiator.base_area_m2

    def test_mass_calculation(self):
        """Radiator should have positive mass."""
        radiator = GridFinRadiator(base_area_m2=2.0)
        assert radiator.mass_kg > 0

    def test_heat_rejection(self):
        """Test heat rejection capacity."""
        radiator = GridFinRadiator(base_area_m2=2.0)
        capacity = radiator.heat_rejection_capacity(temperature_from_celsius(50))
        assert capacity > 0

    def test_starlink_radiator(self):
        """Test Starlink reference radiator."""
        radiator = create_starlink_radiator()
        summary = radiator.summary()
        assert summary["base_area_m2"] == 2.5
        assert summary["area_enhancement"] > 2.0


class TestViewFactors:
    """Tests for orbital view factor calculations."""

    def test_earth_view_factor_decreases_with_altitude(self):
        """Higher altitude means smaller Earth view factor."""
        vf_low = ViewFactors(altitude_km=400)
        vf_high = ViewFactors(altitude_km=800)
        assert vf_low.earth_view_factor > vf_high.earth_view_factor

    def test_view_factors_sum_to_one(self):
        """Earth + space view factors should sum to 1."""
        vf = ViewFactors(altitude_km=550)
        assert abs(vf.earth_view_factor + vf.space_view_factor - 1.0) < 0.001

    def test_eclipse_fraction(self):
        """Eclipse fraction should be between 0 and 1."""
        vf = ViewFactors(altitude_km=550)
        assert 0 <= vf.eclipse_fraction <= 1
        assert 0 <= vf.sunlight_fraction <= 1
        assert abs(vf.eclipse_fraction + vf.sunlight_fraction - 1.0) < 0.001

    def test_orbital_period(self):
        """Test orbital period calculation."""
        vf = ViewFactors(altitude_km=550)
        # LEO period should be ~90-100 minutes
        assert 85 < vf.orbital_period_minutes < 100


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
