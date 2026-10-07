import numpy as np
from fusion_neutron_utils import neutron_energy_spectrum, reactivity, relative_reaction_rates
from pytest import approx, raises

# Hale T+T reactivity in m^3/s at tabulated temperatures, as given by NeSST reac_TT
hale_table = {1: 1.8436e-29, 10: 4.0283e-25, 100: 1.9340e-23}


def test_tt_reactivity_matches_hale_table():
    for ion_temperature_kev, expected in hale_table.items():
        result = reactivity(ion_temperature_kev, "keV", "m^3/s", "T+T=2n+a")
        assert result == approx(expected, rel=1e-4, abs=0)


def test_tt_reactivity_interpolates_between_table_points():
    low = reactivity(10, "keV", "m^3/s", "T+T=2n+a")
    high = reactivity(100, "keV", "m^3/s", "T+T=2n+a")
    middle = reactivity(30, "keV", "m^3/s", "T+T=2n+a")
    assert low < middle < high


def test_tt_reactivity_outside_table_raises():
    with raises(ValueError, match="tabulated"):
        reactivity(0.05, "keV", "m^3/s", "T+T=2n+a")
    with raises(ValueError, match="tabulated"):
        reactivity(2000, "keV", "m^3/s", "T+T=2n+a")


def test_tt_reactivity_only_hale_equation():
    assert reactivity(10, "keV", "m^3/s", "T+T=2n+a", "Hale") == reactivity(10, "keV", "m^3/s", "T+T=2n+a")
    with raises(ValueError, match="Hale"):
        reactivity(10, "keV", "m^3/s", "T+T=2n+a", "Bosch-Hale")


def test_tt_spectrum_matches_nesst():
    # NeSST dNdE_TT at 10 keV, converted to per MeV
    nesst = {2: 0.12771561460611056, 4: 0.1318181805159505, 6: 0.12147579209254876, 8: 0.08710691182086254}
    result = neutron_energy_spectrum(list(nesst.keys()), 10, "keV", "MeV")
    assert result == approx(list(nesst.values()), rel=1e-3, abs=0)


def test_tt_spectrum_normalised():
    energies_ev = np.linspace(0, 12e6, 4000)
    for ion_temperature_kev in (1, 10, 50):
        spectrum = neutron_energy_spectrum(list(energies_ev), ion_temperature_kev, "keV", "eV")
        assert np.trapezoid(spectrum, energies_ev) == approx(1.0, rel=1e-3)


def test_tt_spectrum_units():
    in_ev = neutron_energy_spectrum([4e6], 10e3, "eV", "eV")[0]
    in_mev = neutron_energy_spectrum([4.0], 10, "keV", "MeV")[0]
    assert in_mev == approx(in_ev * 1e6, rel=1e-12, abs=0)


def test_tt_spectrum_invalid_inputs():
    with raises(ValueError, match="T\\+T=2n\\+a"):
        neutron_energy_spectrum([14e6], 10e3, reaction="D+T=n+a")
    with raises(ValueError, match="non-negative"):
        neutron_energy_spectrum([-1.0], 10e3)
    with raises(ValueError, match="positive"):
        neutron_energy_spectrum([1e6], 0.0)


def test_tt_relative_reaction_rate():
    # rate_TT / rate_DT = (f_T^2 / 2 * <sigma v>_TT) / (f_D * f_T * <sigma v>_DT)
    dt, _, _, tt = relative_reaction_rates(10, "keV", deuterium_fraction=0.5, tritium_fraction=0.5)
    expected = 0.5 * reactivity(10, "keV", "m^3/s", "T+T=2n+a") / reactivity(10, "keV", "m^3/s", "D+T=n+a")
    assert tt / dt == approx(expected, rel=1e-9, abs=0)
