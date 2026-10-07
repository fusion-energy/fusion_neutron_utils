import numpy as np
from fusion_neutron_utils import (
    fusion_energy_per_reaction,
    mean_neutron_energy,
    neutron_energy_mean_and_std_dev,
    neutron_energy_spectrum,
    neutron_rate_from_fusion_power,
    neutrons_per_reaction,
    relative_neutron_rates,
    relative_reaction_rates,
)
from pytest import approx, raises

REACTIONS = ["D+T=n+a", "D+D=n+He3", "D+D=p+T", "T+T=2n+a"]
MEV_IN_JOULES = 1.602176634e-13


def test_fusion_energy_per_reaction():
    expected_mev = {"D+T=n+a": 17.589, "D+D=n+He3": 3.269, "D+D=p+T": 4.033, "T+T=2n+a": 11.332}
    for reaction, q in expected_mev.items():
        assert fusion_energy_per_reaction(reaction, "MeV") == approx(q, abs=1e-3)
    assert fusion_energy_per_reaction("D+T=n+a") == approx(fusion_energy_per_reaction("D+T=n+a", "MeV") * 1e6, abs=0)


def test_neutrons_per_reaction():
    assert [neutrons_per_reaction(r) for r in REACTIONS] == [1, 1, 0, 2]


def test_unknown_reaction_raises():
    with raises(ValueError, match="reaction must be one of"):
        fusion_energy_per_reaction("p+B11=3a")
    with raises(ValueError, match="reaction must be one of"):
        neutrons_per_reaction("p+B11=3a")


def test_relative_neutron_rates_weights_reactions_by_neutrons():
    reactions = relative_reaction_rates(10, "keV")
    neutrons = relative_neutron_rates(10, "keV")
    weighted = [f * neutrons_per_reaction(r) for f, r in zip(reactions, REACTIONS)]
    assert neutrons == approx([w / sum(weighted) for w in weighted], rel=1e-12, abs=0)
    assert neutrons[2] == 0.0
    assert sum(neutrons) == approx(1.0)


def test_neutron_rate_pure_dt_limit():
    # with only D+T reactions 1 MW gives 1e6 / Q_DT neutrons per second
    rate = neutron_rate_from_fusion_power(1e6, 10, "keV")
    pure_dt = 1e6 / (fusion_energy_per_reaction("D+T=n+a", "MeV") * MEV_IN_JOULES)
    assert rate == approx(pure_dt, rel=0.01)


def test_neutron_rate_from_fusion_power_mixture():
    for d, t in ((0.5, 0.5), (0.9, 0.1), (1.0, 0.0), (0.0, 1.0)):
        fractions = relative_reaction_rates(10, "keV", d, t)
        joules = sum(f * fusion_energy_per_reaction(r, "MeV") * MEV_IN_JOULES for f, r in zip(fractions, REACTIONS))
        neutrons = sum(f * neutrons_per_reaction(r) for f, r in zip(fractions, REACTIONS))
        expected = 2e9 / joules * neutrons
        assert neutron_rate_from_fusion_power(2e9, 10, "keV", d, t) == approx(expected, rel=1e-12, abs=0)


def test_neutron_rate_negative_power_raises():
    with raises(ValueError, match="fusion_power"):
        neutron_rate_from_fusion_power(-1.0, 10, "keV")


def test_mean_neutron_energy_pure_deuterium_is_dd_mean():
    dd_mean, _ = neutron_energy_mean_and_std_dev(10, "keV", "MeV", "D+D=n+He3")
    assert mean_neutron_energy(10, "keV", "MeV", 1.0, 0.0) == approx(dd_mean, rel=1e-12, abs=0)


def test_mean_neutron_energy_pure_tritium_is_tt_spectrum_mean():
    energies = np.linspace(0, 12, 12001)
    spectrum = np.array(neutron_energy_spectrum(list(energies), 10, "keV", "MeV"))
    expected = np.trapezoid(energies * spectrum, energies) / np.trapezoid(spectrum, energies)
    assert mean_neutron_energy(10, "keV", "MeV", 0.0, 1.0) == approx(expected, rel=1e-4)


def test_mean_neutron_energy_mixture():
    fractions = relative_neutron_rates(10, "keV")
    dt_mean, _ = neutron_energy_mean_and_std_dev(10, "keV", "MeV", "D+T=n+a")
    dd_mean, _ = neutron_energy_mean_and_std_dev(10, "keV", "MeV", "D+D=n+He3")
    tt_mean = mean_neutron_energy(10, "keV", "MeV", 0.0, 1.0)
    expected = fractions[0] * dt_mean + fractions[1] * dd_mean + fractions[3] * tt_mean
    assert mean_neutron_energy(10, "keV", "MeV") == approx(expected, rel=1e-12, abs=0)
