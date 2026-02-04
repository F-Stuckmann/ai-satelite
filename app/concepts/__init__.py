"""Design concepts library for satellite configurations."""

from .baseline_starlink import STARLINK_BASELINE
from .solar_heatsink import SOLAR_HEATSINK_CONCEPT
from .body_mounted_solar import BODY_MOUNTED_CONCEPT
from .graphene_radiator import GRAPHENE_RADIATOR_CONCEPT

CONCEPTS = {
    "starlink_baseline": STARLINK_BASELINE,
    "solar_heatsink": SOLAR_HEATSINK_CONCEPT,
    "body_mounted": BODY_MOUNTED_CONCEPT,
    "graphene_radiator": GRAPHENE_RADIATOR_CONCEPT,
}
