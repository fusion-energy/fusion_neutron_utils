from fusion_neutron_utils import reactivity, relative_reaction_rates
from pytest import approx, mark

# Bosch and Hale, Nucl. Fusion 32 (1992) 611, Table VIII, in cm^3/s.
# The same values are given by NeSST (reac_DT and reac_DD).
bosch_hale_table = {
    "D+T=n+a": {1: 6.857e-21, 2: 2.977e-19, 5: 1.366e-17, 10: 1.136e-16, 20: 4.330e-16, 50: 8.649e-16, 100: 8.448e-16},
    "D+D=n+He3": {1: 9.933e-23, 2: 3.110e-21, 5: 9.128e-20, 10: 6.023e-19, 20: 2.603e-18, 50: 1.133e-17, 100: 2.682e-17},
    "D+D=p+T": {1: 1.017e-22, 2: 3.150e-21, 5: 9.024e-20, 10: 5.781e-19, 20: 2.399e-18, 50: 9.838e-18, 100: 2.244e-17},
}


@mark.parametrize("reaction", bosch_hale_table.keys())
def test_reactivity_matches_bosch_hale_table(reaction):
    for ion_temperature_kev, expected in bosch_hale_table[reaction].items():
        result = reactivity(
            ion_temperature=ion_temperature_kev,
            temperature_units="keV",
            reactivity_units="cm^3/s",
            reaction=reaction,
        )
        assert result == approx(expected, rel=1e-3, abs=0)


def test_reactivity_units():
    cm3 = reactivity(10, "keV", "cm^3/s", "D+T=n+a")
    assert reactivity(10, "keV", "m^3/s", "D+T=n+a") == approx(cm3 * 1e-6, abs=0)
    assert reactivity(10, "keV", "mm^3/s", "D+T=n+a") == approx(cm3 * 1e3, abs=0)
    assert reactivity(10e3, "eV", "cm^3/s", "D+T=n+a") == approx(cm3, abs=0)


def test_reactivity_default_units_are_m3_per_s():
    assert reactivity(10e3) == approx(1.136e-22, rel=1e-3, abs=0)


def test_sadler_van_belle_close_to_bosch_hale():
    for ion_temperature_kev in (5, 10, 20, 50):
        svb = reactivity(ion_temperature_kev, "keV", "m^3/s", "D+T=n+a", "Sadler-Van Belle")
        bh = reactivity(ion_temperature_kev, "keV", "m^3/s", "D+T=n+a", "Bosch-Hale")
        assert svb == approx(bh, rel=0.05, abs=0)


def test_relative_reaction_rates_ratio():
    # rate_ij = f_i * f_j * <sigma v>_ij / (1 + delta_ij), so for a 50:50 fuel
    # each DD branch over DT is 0.5 * <sigma v>_DD / <sigma v>_DT
    dt, dd_n, dd_p = relative_reaction_rates(
        10, "keV", deuterium_fraction=0.5, tritium_fraction=0.5
    )
    assert dt + dd_n + dd_p == approx(1.0)
    assert dd_n / dt == approx(0.5 * 6.023e-19 / 1.136e-16, rel=1e-3, abs=0)
    assert dd_p / dt == approx(0.5 * 5.781e-19 / 1.136e-16, rel=1e-3, abs=0)
