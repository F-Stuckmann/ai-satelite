"""Heat source models for satellite components.

Models thermal power dissipation from various satellite components,
particularly AI chips and power electronics.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np

from ..core.constants import AI_CHIPS


@dataclass
class HeatSource:
    """Generic heat source model.

    Attributes:
        name: Identifier for the heat source
        power_w: Nominal power dissipation in Watts
        mass_kg: Mass of the component in kg
        position: (x, y, z) position on chassis in meters
        t_max_c: Maximum allowable temperature in Celsius
    """

    name: str
    power_w: float
    mass_kg: float
    position: tuple[float, float, float] = (0.0, 0.0, 0.0)
    t_max_c: float = 85.0

    @property
    def power_density_w_per_kg(self) -> float:
        """Power density in W/kg."""
        return self.power_w / self.mass_kg if self.mass_kg > 0 else 0.0

    def power_at_duty_cycle(self, duty_cycle: float) -> float:
        """Calculate average power at given duty cycle.

        Args:
            duty_cycle: Fraction of time at full power (0-1)

        Returns:
            Average power in Watts
        """
        return self.power_w * max(0.0, min(1.0, duty_cycle))


@dataclass
class AIChipHeatSource(HeatSource):
    """AI accelerator chip heat source with compute metrics.

    Extends HeatSource with AI-specific performance metrics.

    Attributes:
        chip_type: Key from AI_CHIPS database
        tflops: Compute performance in TFLOPS
        duty_cycle: Fraction of time at TDP (0-1)
    """

    chip_type: str = "nvidia_h100_sxm"
    tflops: float = 0.0
    duty_cycle: float = 1.0

    def __post_init__(self):
        """Load chip specifications from database if available."""
        if self.chip_type in AI_CHIPS:
            chip = AI_CHIPS[self.chip_type]
            if self.power_w == 0:
                self.power_w = chip["tdp"]
            if self.mass_kg == 0:
                self.mass_kg = chip["mass_kg"]
            if self.tflops == 0:
                # Use FP16 TFLOPS as primary metric
                self.tflops = chip.get("tflops_fp16", chip.get("tflops_bf16", 0))
            self.t_max_c = chip.get("t_max", 85)

    @property
    def tflops_per_watt(self) -> float:
        """Compute efficiency in TFLOPS/W."""
        return self.tflops / self.power_w if self.power_w > 0 else 0.0

    @property
    def tflops_per_kg(self) -> float:
        """Compute density in TFLOPS/kg."""
        return self.tflops / self.mass_kg if self.mass_kg > 0 else 0.0

    @property
    def effective_tflops(self) -> float:
        """Effective TFLOPS accounting for duty cycle."""
        return self.tflops * self.duty_cycle

    @property
    def average_power_w(self) -> float:
        """Average power considering duty cycle."""
        return self.power_w * self.duty_cycle

    def summary(self) -> dict:
        """Return summary of chip properties."""
        return {
            "name": self.name,
            "chip_type": self.chip_type,
            "tdp_w": self.power_w,
            "mass_kg": self.mass_kg,
            "tflops": self.tflops,
            "duty_cycle": self.duty_cycle,
            "average_power_w": round(self.average_power_w, 1),
            "effective_tflops": round(self.effective_tflops, 1),
            "tflops_per_watt": round(self.tflops_per_watt, 2),
            "tflops_per_kg": round(self.tflops_per_kg, 1),
            "t_max_c": self.t_max_c,
        }


class PowerElectronicsHeatSource(HeatSource):
    """Power electronics (DC-DC converters, etc.) heat source."""

    def __init__(
        self,
        name: str,
        power_throughput_w: float,
        efficiency: float = 0.95,
        mass_kg: float = 0.0,
        position: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ):
        """Initialize power electronics heat source.

        Args:
            name: Component identifier
            power_throughput_w: Power handled by the electronics
            efficiency: Conversion efficiency (0-1)
            mass_kg: Component mass
            position: Position on chassis
        """
        # Heat dissipation is the inefficiency
        power_loss = power_throughput_w * (1 - efficiency)

        # Estimate mass if not provided (typical: 1 kg per 500W throughput)
        if mass_kg == 0:
            mass_kg = power_throughput_w / 500

        super().__init__(
            name=name,
            power_w=power_loss,
            mass_kg=mass_kg,
            position=position,
            t_max_c=105,  # Power electronics typically tolerate higher temps
        )

        self.power_throughput_w = power_throughput_w
        self.efficiency = efficiency


def calculate_total_heat_load(sources: list[HeatSource]) -> dict:
    """Calculate total heat load from multiple sources.

    Args:
        sources: List of HeatSource objects

    Returns:
        Dictionary with total power, mass, and breakdown
    """
    total_power = sum(s.power_w for s in sources)
    total_mass = sum(s.mass_kg for s in sources)

    return {
        "total_power_w": total_power,
        "total_mass_kg": total_mass,
        "power_density_w_per_kg": total_power / total_mass if total_mass > 0 else 0,
        "source_count": len(sources),
        "breakdown": [
            {"name": s.name, "power_w": s.power_w, "mass_kg": s.mass_kg}
            for s in sources
        ],
    }


def create_ai_compute_payload(
    chip_type: str = "nvidia_h100_sxm",
    chip_count: int = 1,
    duty_cycle: float = 0.5,
) -> list[AIChipHeatSource]:
    """Create AI compute payload with multiple chips.

    Args:
        chip_type: Type of AI chip from database
        chip_count: Number of chips
        duty_cycle: Operating duty cycle (0-1)

    Returns:
        List of AIChipHeatSource objects
    """
    chips = []
    for i in range(chip_count):
        chip = AIChipHeatSource(
            name=f"AI_Chip_{i+1}",
            power_w=0,  # Will be loaded from database
            mass_kg=0,
            chip_type=chip_type,
            duty_cycle=duty_cycle,
        )
        chips.append(chip)

    return chips
