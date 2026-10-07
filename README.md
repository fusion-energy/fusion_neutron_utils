[![CI Python testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-python.yml) [![CI Rust testing](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml/badge.svg)](https://github.com/fusion-energy/fusion_neutron_utils/actions/workflows/ci-rust.yml)

A package for calculating neutron properties from DT, DD and TT fusion reactions.

This package accurately calculates the neutron energies and distributions by accounting for the plasma temperature.

- mean neutron energy and standard deviation for D+T and D+D
- T+T neutron energy spectrum
- thermal reactivity for D+T, both D+D branches and T+T
- relative reaction rates and relative neutron rates for a D/T fuel mix
- energy released and neutrons emitted per reaction (Q-values)
- neutron emission rate for a fusion power, used to normalise OpenMC tallies
- mean energy of all emitted neutrons for a D/T fuel mix

A Python :snake: package with a Rust :crab: backend.

The package makes use of the [Sadler-Van Belle formula](https://doi.org/10.1016/j.fusengdes.2012.02.025) and [Bosch-Hale](https://doi.org/10.1088/0029-5515%2F32%2F4%2FI07) parametrization for reactivity. For energy distributions [Ballabio](https://doi.org/10.1088/0029-5515/38/11/310) is used. The TT reactivity (Hale) and TT spectrum data are taken from [NeSST](https://github.com/aidancrilly/NeSST) (MIT license).

This package has been inspired by the [NeSST package](https://github.com/aidancrilly/NeSST).

## Install

```bash
pip install fusion_neutron_utils
```

## Usage

All functions take a single ion temperature and return Python floats or lists. Temperatures default to eV and energies to eV, use the `temperature_units` and `neutron_energy_units` arguments to change these. See the [Reference](#reference) section for all accepted values.

### Neutron energy for D+T and D+D

The mean neutron energy and its standard deviation, from the Ballabio fits. These can be used in an [openmc.stats.Normal](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.stats.Normal.html) energy distribution.

```python
from fusion_neutron_utils import neutron_energy_mean_and_std_dev
neutron_energy_mean_and_std_dev(
    reaction='D+T=n+a',
    ion_temperature=30e3,
    temperature_units='eV',
    neutron_energy_units='eV'
)
>>>(14092196.942891473, 413861.375751198)
```

### Neutron energy spectrum for T+T

The T+T reaction has three bodies in the final state, so its neutrons have a continuous energy spectrum rather than a peak. The spectrum dN/dE, normalised to integrate to 1, can be found at a list (or numpy array) of energies. This can be used in an [openmc.stats.Tabular](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.stats.Tabular.html) energy distribution.

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

### Reactivity

The thermal reactivity <sigma v> for a single reaction.

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

For D+T the older Sadler-Van Belle fit can be used instead of Bosch-Hale. It is within a few percent of Bosch-Hale.

```python
from fusion_neutron_utils import reactivity
reactivity(
    ion_temperature=30,
    temperature_units='keV',
    reactivity_units='m^3/s',
    reaction='D+T=n+a',
    equation='Sadler-Van Belle'
)
>>>6.5695456075579885e-22
```

T+T uses a tabulated reactivity from Hale, so it takes `equation='Hale'` or no equation.

```python
from fusion_neutron_utils import reactivity
reactivity(
    ion_temperature=30,
    temperature_units='keV',
    reactivity_units='cm^3/s',
    reaction='T+T=2n+a',
)
>>>3.911718324197931e-18
```

### Relative reaction rates and relative neutron rates

For a fuel with the given fractions of deuterium and tritium ions, these return a list for the DT, DD (n+He3), DD (p+T) and TT reactions, in that order.

`relative_reaction_rates` gives the fraction of all **reactions**:

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

`relative_neutron_rates` gives the fraction of all **neutrons**. D+D=p+T emits no neutron and each T+T reaction emits two, so these differ from the reaction fractions. Use these as the relative strengths of each reaction in an OpenMC source.

```python
from fusion_neutron_utils import relative_neutron_rates
relative_neutron_rates(
    ion_temperature=10,
    temperature_units='keV',
    deuterium_fraction=0.5,
    tritium_fraction=0.5,
)
>>>[0.9938422010961021, 0.0026341092281159653, 0.0, 0.003523689675781927]
```

### Neutron emission rate from fusion power

The number of neutrons per second for a fusion power in W, including all four reactions for the given ion temperature and fuel composition:

```python
from fusion_neutron_utils import neutron_rate_from_fusion_power
neutron_rate_from_fusion_power(
    fusion_power=500e6,
    ion_temperature=10,
    temperature_units='keV',
    deuterium_fraction=0.5,
    tritium_fraction=0.5,
)
>>>1.781276864853604e+20
```

The energy released per reaction (Q-value, from the AME2020 atomic masses) and the neutrons emitted per reaction are also available:

```python
from fusion_neutron_utils import fusion_energy_per_reaction, neutrons_per_reaction
fusion_energy_per_reaction('D+T=n+a', energy_units='MeV')
>>>17.589299865
neutrons_per_reaction('T+T=2n+a')
>>>2
```

### Mean energy of all emitted neutrons

Averaged over the neutrons from all reactions, using the Ballabio means for D+T and D+D and the mean of the T+T spectrum:

```python
from fusion_neutron_utils import mean_neutron_energy
mean_neutron_energy(
    ion_temperature=10,
    temperature_units='keV',
    neutron_energy_units='MeV',
    deuterium_fraction=0.5,
    tritium_fraction=0.5,
)
>>>13.999459130682638
```

### Comparing fuel mixes and temperatures

The functions take a single temperature and fuel composition, so loop to compare several. The share of neutrons from D+D grows quickly as the fuel becomes deuterium rich:

```python
from fusion_neutron_utils import relative_neutron_rates

for deuterium_fraction in [0.5, 0.7, 0.9, 1.0]:
    dt, dd, _, tt = relative_neutron_rates(
        ion_temperature=10,
        temperature_units='keV',
        deuterium_fraction=deuterium_fraction,
        tritium_fraction=1 - deuterium_fraction,
    )
    print(f"D fraction {deuterium_fraction:.1f}: DT {dt:.4f} DD {dd:.4f} TT {tt:.4f}")
```

```
D fraction 0.5: DT 0.9938 DD 0.0026 TT 0.0035
D fraction 0.7: DT 0.9924 DD 0.0061 TT 0.0015
D fraction 0.9: DT 0.9763 DD 0.0233 TT 0.0004
D fraction 1.0: DT 0.0000 DD 1.0000 TT 0.0000
```

and the neutrons per second per MW and the mean neutron energy change with temperature:

```python
from fusion_neutron_utils import mean_neutron_energy, neutron_rate_from_fusion_power

for ion_temperature in [5, 10, 20, 40]:
    rate = neutron_rate_from_fusion_power(fusion_power=1e6, ion_temperature=ion_temperature, temperature_units='keV')
    energy = mean_neutron_energy(ion_temperature=ion_temperature, temperature_units='keV', neutron_energy_units='MeV')
    print(f"{ion_temperature} keV: {rate:.4e} neutrons/s per MW, mean energy {energy:.3f} MeV")
```

```
5 keV: 3.5641e+17 neutrons/s per MW, mean energy 13.978 MeV
10 keV: 3.5626e+17 neutrons/s per MW, mean energy 13.999 MeV
20 keV: 3.5654e+17 neutrons/s per MW, mean energy 14.004 MeV
40 keV: 3.5778e+17 neutrons/s per MW, mean energy 13.973 MeV
```

## Making an OpenMC source

The energy distribution of a D/T plasma neutron source combines a Normal distribution for D+T and D+D with the tabulated T+T spectrum, weighted by the share of neutrons from each reaction. Use `relative_neutron_rates` for the weights, not `relative_reaction_rates`, as the source samples neutrons. This example is a point source; the same energy distribution can be used with any spatial distribution.

```python
import numpy as np
import openmc
from fusion_neutron_utils import (
    neutron_energy_mean_and_std_dev,
    neutron_energy_spectrum,
    neutron_rate_from_fusion_power,
    relative_neutron_rates,
)

ion_temperature = 10  # keV
deuterium_fraction = 0.5
tritium_fraction = 0.5

dt_mean, dt_std_dev = neutron_energy_mean_and_std_dev(
    ion_temperature, temperature_units='keV', neutron_energy_units='eV', reaction='D+T=n+a'
)
dd_mean, dd_std_dev = neutron_energy_mean_and_std_dev(
    ion_temperature, temperature_units='keV', neutron_energy_units='eV', reaction='D+D=n+He3'
)
tt_energies = np.linspace(1e4, 12e6, 1200)  # eV
tt_spectrum = neutron_energy_spectrum(
    tt_energies, ion_temperature, temperature_units='keV', neutron_energy_units='eV'
)

dt_share, dd_share, _, tt_share = relative_neutron_rates(
    ion_temperature,
    temperature_units='keV',
    deuterium_fraction=deuterium_fraction,
    tritium_fraction=tritium_fraction,
)

source = openmc.IndependentSource()
source.space = openmc.stats.Point((0, 0, 0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Mixture(
    [dt_share, dd_share, tt_share],
    [
        openmc.stats.Normal(dt_mean, dt_std_dev),
        openmc.stats.Normal(dd_mean, dd_std_dev),
        openmc.stats.Tabular(tt_energies, tt_spectrum),
    ],
)
# neutrons per second, useful when combining several sources
source.strength = neutron_rate_from_fusion_power(
    fusion_power=500e6,
    ion_temperature=ion_temperature,
    temperature_units='keV',
    deuterium_fraction=deuterium_fraction,
    tritium_fraction=tritium_fraction,
)
```

For a pure deuterium or pure tritium fuel some shares are zero, those entries can be left out of the mixture.

## Normalising OpenMC tallies

OpenMC tally results are per source neutron, whatever the source strength. Multiply them by the neutron emission rate to get rates. For example a flux tally (in particle-cm per source neutron) divided by the cell volume (cm^3) and multiplied by the neutrons per second gives the flux in neutrons per cm^2 per second:

```python
import openmc
from fusion_neutron_utils import neutron_rate_from_fusion_power

neutrons_per_second = neutron_rate_from_fusion_power(
    fusion_power=500e6,
    ion_temperature=10,
    temperature_units='keV',
    deuterium_fraction=0.5,
    tritium_fraction=0.5,
)

with openmc.StatePoint('statepoint.10.h5') as statepoint:
    flux_tally = statepoint.get_tally(name='flux')
    flux_per_source_neutron = flux_tally.mean.flatten()  # particle-cm per source neutron

cell_volume = 1e6  # cm^3, the volume of the tallied cell
flux = flux_per_source_neutron * neutrons_per_second / cell_volume  # neutrons/cm^2/s
```

The neutron rate depends on the ion temperature and fuel composition through the mix of reactions. For a plasma with temperature and density profiles, use a representative temperature or a source model that integrates over the profiles, such as [openmc-plasma-source](https://github.com/fusion-energy/openmc-plasma-source).

## Reference

### Reactions

| `reaction` | Energy released (MeV) | Neutrons per reaction | Neutron energy | Reactivity `equation` |
|---|---|---|---|---|
| `'D+T=n+a'` | 17.589 | 1 | `neutron_energy_mean_and_std_dev` | `'Bosch-Hale'` (default) or `'Sadler-Van Belle'` |
| `'D+D=n+He3'` | 3.269 | 1 | `neutron_energy_mean_and_std_dev` | `'Bosch-Hale'` |
| `'D+D=p+T'` | 4.033 | 0 | none (no neutron) | `'Bosch-Hale'` |
| `'T+T=2n+a'` | 11.332 | 2 | `neutron_energy_spectrum` | `'Hale'` |

### Units

| Argument | Accepted values | Default |
|---|---|---|
| `temperature_units` | `'eV'`, `'keV'`, `'MeV'`, `'GeV'` | `'eV'` |
| `neutron_energy_units`, `energy_units` | `'eV'`, `'keV'`, `'MeV'`, `'GeV'` | `'eV'` |
| `reactivity_units` | `'m^3/s'`, `'cm^3/s'`, `'mm^3/s'` | `'m^3/s'` |
| `fusion_power` | W | |

### Defaults

- `reaction` defaults to `'D+T=n+a'` in `reactivity` and `neutron_energy_mean_and_std_dev`, and to `'T+T=2n+a'` in `neutron_energy_spectrum`.
- `deuterium_fraction` and `tritium_fraction` default to 0.5 each and must sum to 1.

### Validity

- The Bosch-Hale parametrization is fitted for ion temperatures of 0.2 keV to 100 keV.
- The T+T reactivity is tabulated from 0.1 keV to 1000 keV. Outside this range `reactivity` for T+T raises a `ValueError`, and so do the functions that combine all four reactions (`relative_reaction_rates`, `relative_neutron_rates`, `neutron_rate_from_fusion_power` and `mean_neutron_energy`).
- Invalid reactions, units, equations, fractions and non-positive temperatures raise a `ValueError`.
