[![CI Python testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml) [![CI Rust testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml)

A package for calculating neutron properties from DT and DD fusion reactions.

This package accurately calculates the neutron energies and distributions by accounting for the plasma temperature.

- mean neutron energy
- neutron energy standard deviation
- thermal reactivity 

A Python :snake: package with a Rust :crab: backend.

The package makes use of the [Sadler–Van Belle formula](https://doi.org/10.1016/j.fusengdes.2012.02.025) and [Bosch-Hale](https://doi.org/10.1088/0029-5515%2F32%2F4%2FI07) parametrization for reactivity. For energy distributions [Ballabio](https://doi.org/10.1088/0029-5515/38/11/310) is used.

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


The relative reaction rates can be found for the different reactions, this can be useful for setting the relative ```Source.strength``` in OpenMC. Relative reaction rates returns the DT, DD (n+He3), and DD (p+T) reaction rates in that order. Note the  DD (p+T)  does not emit a neutron but is there for completeness.
```python
from fusion_neutron_utils import relative_reaction_rates
relative_reaction_rates(
    ion_temperature=30e3,
    temperature_units='eV',
    dt_fraction=0.1,
    dd_fraction=0.9,  
)
>>>[0.9989900631769298, 0.0005595866500573114, 0.00045035017301275736] 
```

The reactivity can also be found, this can be useful for finding the relative reaction rate for D+D or D+T at a specific temperature

```python
from fusion_neutron_utils import reactivity
reactivity(
    ion_temperature=30e3,
    temperature_units='eV',
    reactivity_units='m^3/s',
    reaction='D+D=p+T',
    equation='Bosch-Hale'
)
>>>2.2356373732343755e-53
```
