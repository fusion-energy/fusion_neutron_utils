from collections.abc import Sequence

__all__ = [
    "fusion_energy_per_reaction",
    "mean_neutron_energy",
    "neutron_energy_mean_and_std_dev",
    "neutron_energy_spectrum",
    "neutron_rate_from_fusion_power",
    "neutrons_per_reaction",
    "reactivity",
    "relative_neutron_rates",
    "relative_reaction_rates",
]


def neutron_energy_mean_and_std_dev(
    ion_temperature: float,
    temperature_units: str | None = ...,
    neutron_energy_units: str | None = ...,
    reaction: str | None = ...,
) -> tuple[float, float]: ...


def reactivity(
    ion_temperature: float,
    temperature_units: str | None = ...,
    reactivity_units: str | None = ...,
    reaction: str | None = ...,
    equation: str | None = ...,
) -> float: ...


def relative_reaction_rates(
    ion_temperature: float,
    temperature_units: str | None = ...,
    deuterium_fraction: float | None = ...,
    tritium_fraction: float | None = ...,
    equation: str | None = ...,
) -> list[float]: ...


def relative_neutron_rates(
    ion_temperature: float,
    temperature_units: str | None = ...,
    deuterium_fraction: float | None = ...,
    tritium_fraction: float | None = ...,
) -> list[float]: ...


def neutron_energy_spectrum(
    energies: Sequence[float],
    ion_temperature: float,
    temperature_units: str | None = ...,
    neutron_energy_units: str | None = ...,
    reaction: str | None = ...,
) -> list[float]: ...


def fusion_energy_per_reaction(
    reaction: str,
    energy_units: str | None = ...,
) -> float: ...


def neutrons_per_reaction(reaction: str) -> int: ...


def neutron_rate_from_fusion_power(
    fusion_power: float,
    ion_temperature: float,
    temperature_units: str | None = ...,
    deuterium_fraction: float | None = ...,
    tritium_fraction: float | None = ...,
) -> float: ...


def mean_neutron_energy(
    ion_temperature: float,
    temperature_units: str | None = ...,
    neutron_energy_units: str | None = ...,
    deuterium_fraction: float | None = ...,
    tritium_fraction: float | None = ...,
) -> float: ...
