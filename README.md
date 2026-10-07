[![CI Python testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml) [![CI Rust testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml)

A package for calculating neutron properties from DT, DD and TT fusion reactions.

This package accurately calculates the neutron energies and distributions by accounting for the plasma temperature.

- mean neutron energy
- neutron energy standard deviation
- thermal reactivity
- TT neutron energy spectrum

A Python :snake: package with a Rust :crab: backend.

The package makes use of the [Sadler–Van Belle formula](https://doi.org/10.1016/j.fusengdes.2012.02.025) and [Bosch-Hale](https://doi.org/10.1088/0029-5515%2F32%2F4%2FI07) parametrization for reactivity. For energy distributions [Ballabio](https://doi.org/10.1088/0029-5515/38/11/310) is used. The TT reactivity (Hale) and TT spectrum data are taken from [NeSST](https://github.com/aidancrilly/NeSST) (MIT license).

This package has been inspired by the [NeSST package](https://github.com/aidancrilly/NeSST).

## Install

```bash
pip install fusion_neutron_utils
```

## Usage

To find the average neutron energy and the standard deviation of that energy.

These results could then be input into [openmc.stats.Normal](https://docs.openmc.org/en/v0.12.1/pythonapi/generated/openmc.stats.Normal.html) distribution for a neutron source term in OpenMC

```python
from fusion_neutron_utils import neutron_energy_mean_and_std_dev
neutron_energy_mean_and_std_dev(
    reaction='D+T=n+a',
    ion_temperature=30e3,
    temperature_units='eV',
    neutron_energy_units='eV'
)
>>>>(14092196.942384735, 413861.375751198)
```


The relative reaction rates can be found for the different reactions, this can be useful for setting the relative ```Source.strength``` in OpenMC. Relative reaction rates returns the fractions of all reactions that are DT, DD (n+He3), DD (p+T) and TT, in that order, for a fuel with the given fractions of deuterium and tritium ions. Note the  DD (p+T)  does not emit a neutron but is there for completeness, and each TT reaction emits two neutrons.
```python
from fusion_neutron_utils import relative_reaction_rates
relative_reaction_rates(
    ion_temperature=30e3,
    temperature_units='eV',
    deuterium_fraction=0.5,
    tritium_fraction=0.5,
)
>>>[0.9896963934218174, 0.0039039462139493744, 0.0035022352914826, 0.002897425072750721]
```

The reactivity can also be found, this can be useful for finding the relative reaction rate for D+D, D+T or T+T at a specific temperature

```python
from fusion_neutron_utils import reactivity
reactivity(
    ion_temperature=30e3,
    temperature_units='eV',
    reactivity_units='m^3/s',
    reaction='D+D=p+T',
    equation='Bosch-Hale'
)
>>>4.728252714517674e-24
```

The TT reaction has three bodies in the final state, so its neutrons have a continuous energy spectrum rather than a peak. The spectrum dN/dE, normalised to integrate to 1, can be found at a list of energies. This can be used with [openmc.stats.Tabular](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.stats.Tabular.html).

```python
from fusion_neutron_utils import neutron_energy_spectrum
neutron_energy_spectrum(
    energies=[2.0, 4.0, 6.0],
    ion_temperature=10,
    temperature_units='keV',
    neutron_energy_units='MeV',
    reaction='T+T=2n+a',
)
>>>[0.12770950136868156, 0.13181225603982472, 0.12147036033626678]
```
