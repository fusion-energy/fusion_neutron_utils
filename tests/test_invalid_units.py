from fusion_neutron_utils import (
    neutron_energy_mean_and_std_dev,
    reactivity,
    relative_reaction_rates,
)
from pytest import raises


def test_reactivity_invalid_temperature_units():
    with raises(ValueError, match="Invalid temperature units"):
        reactivity(10, temperature_units="K")


def test_reactivity_invalid_reactivity_units():
    with raises(ValueError, match="Invalid reaction rate units"):
        reactivity(10, temperature_units="keV", reactivity_units="m3")


def test_reactivity_invalid_equation():
    with raises(ValueError, match="equations are supported"):
        reactivity(10, temperature_units="keV", equation="unknown")


def test_neutron_energy_invalid_temperature_units():
    with raises(ValueError, match="Invalid temperature units"):
        neutron_energy_mean_and_std_dev(10, temperature_units="K", reaction="D+T=n+a")


def test_neutron_energy_invalid_energy_units():
    with raises(ValueError, match="Invalid energy units"):
        neutron_energy_mean_and_std_dev(
            10, temperature_units="keV", neutron_energy_units="J", reaction="D+T=n+a"
        )


def test_relative_reaction_rates_invalid_temperature_units():
    with raises(ValueError, match="Invalid temperature units"):
        relative_reaction_rates(10, temperature_units="K")
