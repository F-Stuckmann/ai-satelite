"""Tests for structural analysis module."""

import pytest
import numpy as np

from src.structural.materials import (
    MaterialProperties,
    thermal_stress,
    thermal_strain,
    safety_factor,
    mass_of_panel,
    mass_of_hollow_cylinder,
    buckling_stress_panel,
    natural_frequency_panel,
    cte_mismatch_stress,
    structural_summary,
)


class TestMaterialProperties:
    """Tests for MaterialProperties dataclass."""

    def test_from_database_aluminum(self):
        """Test loading aluminum from database."""
        mat = MaterialProperties.from_database("aluminum_6061")

        assert mat.name == "Aluminum 6061-T6"
        assert mat.E > 0
        assert mat.rho > 0
        assert mat.sigma_yield > 0
        assert 0 < mat.nu < 0.5

    def test_from_database_invalid_raises(self):
        """Test loading invalid material raises error."""
        with pytest.raises(ValueError, match="Unknown material"):
            MaterialProperties.from_database("unobtainium")

    def test_shear_modulus_positive(self):
        """Test shear modulus is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        assert mat.G > 0

    def test_bulk_modulus_positive(self):
        """Test bulk modulus is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        assert mat.K > 0

    def test_shear_modulus_relation(self):
        """Test G = E / (2 * (1 + nu)) relationship."""
        mat = MaterialProperties.from_database("aluminum_6061")

        expected_G = mat.E / (2 * (1 + mat.nu))
        assert mat.G == pytest.approx(expected_G, rel=1e-10)

    def test_bulk_modulus_relation(self):
        """Test K = E / (3 * (1 - 2*nu)) relationship."""
        mat = MaterialProperties.from_database("titanium_6al4v")

        expected_K = mat.E / (3 * (1 - 2 * mat.nu))
        assert mat.K == pytest.approx(expected_K, rel=1e-10)


class TestThermalStress:
    """Tests for thermal stress calculations."""

    def test_thermal_stress_positive(self):
        """Test thermal stress is positive for positive delta-T."""
        mat = MaterialProperties.from_database("aluminum_6061")

        stress = thermal_stress(mat, 100)  # 100K change

        assert stress > 0

    def test_thermal_stress_symmetric(self):
        """Test thermal stress is same for heating and cooling."""
        mat = MaterialProperties.from_database("aluminum_6061")

        stress_heat = thermal_stress(mat, 100)
        stress_cool = thermal_stress(mat, -100)

        assert stress_heat == pytest.approx(stress_cool, rel=1e-10)

    def test_thermal_stress_unconstrained_zero(self):
        """Test unconstrained thermal expansion has zero stress."""
        mat = MaterialProperties.from_database("aluminum_6061")

        stress = thermal_stress(mat, 100, constrained=False)

        assert stress == 0.0

    def test_thermal_stress_formula(self):
        """Test thermal stress follows σ = E * α * ΔT."""
        mat = MaterialProperties.from_database("aluminum_6061")
        delta_T = 100

        stress = thermal_stress(mat, delta_T)
        expected = mat.E * mat.alpha * abs(delta_T)

        assert stress == pytest.approx(expected, rel=1e-10)


class TestThermalStrain:
    """Tests for thermal strain calculations."""

    def test_thermal_strain_positive_for_heating(self):
        """Test positive strain for heating (expansion)."""
        mat = MaterialProperties.from_database("aluminum_6061")

        strain = thermal_strain(mat, 100)

        assert strain > 0

    def test_thermal_strain_negative_for_cooling(self):
        """Test negative strain for cooling (contraction)."""
        mat = MaterialProperties.from_database("aluminum_6061")

        strain = thermal_strain(mat, -100)

        assert strain < 0

    def test_thermal_strain_formula(self):
        """Test thermal strain follows ε = α * ΔT."""
        mat = MaterialProperties.from_database("titanium_6al4v")
        delta_T = 50

        strain = thermal_strain(mat, delta_T)
        expected = mat.alpha * delta_T

        assert strain == pytest.approx(expected, rel=1e-10)


class TestSafetyFactor:
    """Tests for safety factor calculations."""

    def test_safety_factor_positive(self):
        """Test safety factor is positive for positive stress."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sf = safety_factor(100e6, mat)  # 100 MPa

        assert sf > 0

    def test_safety_factor_infinite_for_zero_stress(self):
        """Test safety factor is infinite for zero stress."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sf = safety_factor(0, mat)

        assert sf == float('inf')

    def test_safety_factor_decreases_with_stress(self):
        """Test safety factor decreases as stress increases."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sf_low = safety_factor(50e6, mat)
        sf_high = safety_factor(200e6, mat)

        assert sf_low > sf_high

    def test_safety_factor_yield_vs_ultimate(self):
        """Test safety factor using ultimate strength is higher."""
        mat = MaterialProperties(
            name="Test",
            E=70e9,
            nu=0.33,
            rho=2700,
            sigma_yield=200e6,
            sigma_ultimate=300e6,
        )

        sf_yield = safety_factor(100e6, mat, use_ultimate=False)
        sf_ultimate = safety_factor(100e6, mat, use_ultimate=True)

        assert sf_ultimate > sf_yield


class TestMassCalculations:
    """Tests for mass calculation functions."""

    def test_panel_mass_positive(self):
        """Test panel mass is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        mass = mass_of_panel(0.5, 0.3, 2.0, mat)

        assert mass > 0

    def test_panel_mass_proportional_to_dimensions(self):
        """Test panel mass scales with dimensions."""
        mat = MaterialProperties.from_database("aluminum_6061")

        mass_1x = mass_of_panel(0.5, 0.3, 2.0, mat)
        mass_2x = mass_of_panel(1.0, 0.3, 2.0, mat)  # Double length

        assert mass_2x == pytest.approx(2 * mass_1x, rel=1e-10)

    def test_hollow_cylinder_mass_positive(self):
        """Test hollow cylinder mass is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        mass = mass_of_hollow_cylinder(0.1, 0.095, 1.0, mat)

        assert mass > 0

    def test_hollow_cylinder_inner_less_than_outer(self):
        """Test hollow cylinder with inner > outer gives negative mass."""
        mat = MaterialProperties.from_database("aluminum_6061")

        # This would be physically impossible but tests the formula
        mass = mass_of_hollow_cylinder(0.095, 0.1, 1.0, mat)

        # The function doesn't validate, so mass will be negative
        assert mass < 0


class TestBucklingStress:
    """Tests for buckling stress calculations."""

    def test_buckling_stress_positive(self):
        """Test critical buckling stress is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sigma_cr = buckling_stress_panel(0.5, 0.3, 2.0, mat)

        assert sigma_cr > 0

    def test_buckling_stress_increases_with_thickness(self):
        """Test buckling stress increases with thickness (^2 relationship)."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sigma_thin = buckling_stress_panel(0.5, 0.3, 1.0, mat)
        sigma_thick = buckling_stress_panel(0.5, 0.3, 2.0, mat)

        # Buckling stress: σ_cr = k*π²*E*t² / (12*(1-ν²)*b²)
        # Doubling thickness gives 4x stress (t² relationship)
        ratio = sigma_thick / sigma_thin
        assert ratio == pytest.approx(4.0, rel=0.01)

    def test_clamped_higher_than_simply_supported(self):
        """Test clamped boundary has higher buckling stress."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sigma_ss = buckling_stress_panel(0.5, 0.3, 2.0, mat, "simply_supported")
        sigma_cl = buckling_stress_panel(0.5, 0.3, 2.0, mat, "clamped")

        assert sigma_cl > sigma_ss


class TestNaturalFrequency:
    """Tests for natural frequency calculations."""

    def test_natural_frequency_positive(self):
        """Test natural frequency is positive."""
        mat = MaterialProperties.from_database("aluminum_6061")

        freq = natural_frequency_panel(0.5, 0.3, 2.0, mat)

        assert freq > 0

    def test_frequency_increases_with_thickness(self):
        """Test frequency increases with thickness."""
        mat = MaterialProperties.from_database("aluminum_6061")

        freq_thin = natural_frequency_panel(0.5, 0.3, 1.0, mat)
        freq_thick = natural_frequency_panel(0.5, 0.3, 2.0, mat)

        assert freq_thick > freq_thin

    def test_frequency_decreases_with_panel_size(self):
        """Test frequency decreases with larger panel."""
        mat = MaterialProperties.from_database("aluminum_6061")

        freq_small = natural_frequency_panel(0.3, 0.2, 2.0, mat)
        freq_large = natural_frequency_panel(0.6, 0.4, 2.0, mat)

        assert freq_small > freq_large

    def test_clamped_higher_frequency(self):
        """Test clamped boundary has higher frequency."""
        mat = MaterialProperties.from_database("aluminum_6061")

        freq_ss = natural_frequency_panel(0.3, 0.3, 2.0, mat, "simply_supported")
        freq_cl = natural_frequency_panel(0.3, 0.3, 2.0, mat, "clamped")

        assert freq_cl > freq_ss


class TestCTEMismatch:
    """Tests for CTE mismatch stress calculations."""

    def test_cte_mismatch_produces_stress(self):
        """Test CTE mismatch between materials produces stress."""
        mat1 = MaterialProperties.from_database("aluminum_6061")
        mat2 = MaterialProperties.from_database("cfrp_quasi_isotropic")  # Very different CTE

        sigma1, sigma2 = cte_mismatch_stress(mat1, mat2, 100)

        assert sigma1 > 0
        assert sigma2 > 0

    def test_same_material_no_mismatch_stress(self):
        """Test same material has zero mismatch stress."""
        mat = MaterialProperties.from_database("aluminum_6061")

        sigma1, sigma2 = cte_mismatch_stress(mat, mat, 100)

        assert sigma1 == pytest.approx(0.0, abs=1e-10)
        assert sigma2 == pytest.approx(0.0, abs=1e-10)

    def test_mismatch_stress_increases_with_delta_t(self):
        """Test mismatch stress increases with temperature change."""
        mat1 = MaterialProperties.from_database("aluminum_6061")
        mat2 = MaterialProperties.from_database("titanium_6al4v")

        sigma1_small, _ = cte_mismatch_stress(mat1, mat2, 50)
        sigma1_large, _ = cte_mismatch_stress(mat1, mat2, 100)

        assert sigma1_large > sigma1_small


class TestStructuralSummary:
    """Tests for structural summary function."""

    def test_summary_contains_required_keys(self):
        """Test summary returns all required keys."""
        result = structural_summary(
            (0.5, 0.3, 0.1),  # Chassis dimensions
            2.0,  # Wall thickness
            "aluminum_6061",
        )

        required_keys = [
            "material",
            "mass_kg",
            "surface_area_m2",
            "thermal_stress_MPa",
            "thermal_safety_factor",
            "natural_frequency_Hz",
            "buckling_stress_MPa",
            "yield_strength_MPa",
        ]

        for key in required_keys:
            assert key in result

    def test_summary_mass_positive(self):
        """Test summary mass is positive."""
        result = structural_summary(
            (0.5, 0.3, 0.1),
            2.0,
            "aluminum_6061",
        )

        assert result["mass_kg"] > 0

    def test_summary_respects_temp_range(self):
        """Test summary uses provided temperature range."""
        result_small = structural_summary(
            (0.5, 0.3, 0.1),
            2.0,
            "aluminum_6061",
            temp_range=(-50, 50),  # 100K range
        )
        result_large = structural_summary(
            (0.5, 0.3, 0.1),
            2.0,
            "aluminum_6061",
            temp_range=(-100, 100),  # 200K range
        )

        # Larger temp range = larger thermal stress
        assert result_large["thermal_stress_MPa"] > result_small["thermal_stress_MPa"]
