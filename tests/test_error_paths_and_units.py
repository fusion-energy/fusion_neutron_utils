import math

from pytest import approx, mark, raises

from fusion_neutron_utils import (
    fusion_energy_per_reaction,
    mean_neutron_energy,
    neutron_energy_mean_and_std_dev,
    neutron_energy_spectrum,
    neutron_rate_from_fusion_power,
    reactivity,
    relative_neutron_rates,
    relative_reaction_rates,
)

# factor to convert a value in keV into each energy unit
ENERGY_UNIT_FROM_KEV = {"eV": 1e3, "keV": 1.0, "MeV": 1e-3, "GeV": 1e-6}


@mark.parametrize("reaction", ["D+D=p+T", "T+T=2n+a", "p+B11=3a", ""])
def test_neutron_energy_invalid_reaction(reaction):
    with raises(ValueError, match="reaction must be either"):
        neutron_energy_mean_and_std_dev(10, "keV", reaction=reaction)


@mark.parametrize("reaction", ["D+T=n+a", "D+D=n+He3"])
@mark.parametrize("units", ENERGY_UNIT_FROM_KEV.keys())
def test_neutron_energy_mean_and_std_dev_units(reaction, units):
    mean_kev, std_kev = neutron_energy_mean_and_std_dev(10, "keV", "keV", reaction)
    mean, std = neutron_energy_mean_and_std_dev(10, "keV", units, reaction)
    factor = ENERGY_UNIT_FROM_KEV[units]
    assert mean == approx(mean_kev * factor, rel=1e-12, abs=0)
    assert std == approx(std_kev * factor, rel=1e-12, abs=0)


@mark.parametrize("reaction", ["D+T=n+a", "D+D=n+He3"])
@mark.parametrize("temperature, units", [(10e3, "eV"), (0.01, "MeV"), (1e-5, "GeV")])
def test_neutron_energy_temperature_units(reaction, temperature, units):
    expected = neutron_energy_mean_and_std_dev(10, "keV", "eV", reaction)
    result = neutron_energy_mean_and_std_dev(temperature, units, "eV", reaction)
    assert result == approx(expected, rel=1e-12, abs=0)


@mark.parametrize("reaction", ["D+T=n+a", "D+D=n+He3"])
def test_std_dev_increases_with_ion_temperature(reaction):
    stds = [
        neutron_energy_mean_and_std_dev(t, "keV", "eV", reaction)[1]
        for t in (1, 5, 10, 20, 50)
    ]
    assert stds == sorted(stds)
    assert all(s > 0 for s in stds)


def test_neutron_energy_negative_temperature_raises():
    with raises(ValueError, match="positive"):
        neutron_energy_mean_and_std_dev(-1, "keV")


@mark.parametrize("ion_temperature", [0.0, -5.0])
def test_reactivity_non_positive_temperature_raises(ion_temperature):
    with raises(ValueError, match="positive"):
        reactivity(ion_temperature, "keV")


def test_reactivity_invalid_reaction():
    with raises(ValueError, match="reactions are supported"):
        reactivity(10, "keV", reaction="D+He3=p+a")


def test_reactivity_sadler_van_belle_only_supports_dt():
    for reaction in ("D+D=n+He3", "D+D=p+T"):
        with raises(ValueError, match="Sadler-Van Belle"):
            reactivity(10, "keV", reaction=reaction, equation="Sadler-Van Belle")


def test_reactivity_hale_equation_not_valid_for_dt():
    with raises(ValueError, match="equations are supported"):
        reactivity(10, "keV", reaction="D+T=n+a", equation="Hale")


@mark.parametrize("temperature, units", [(10e3, "eV"), (0.01, "MeV"), (1e-5, "GeV")])
def test_reactivity_temperature_units(temperature, units):
    expected = reactivity(10, "keV", "m^3/s", "D+T=n+a")
    assert reactivity(temperature, units, "m^3/s", "D+T=n+a") == approx(expected, rel=1e-12, abs=0)


def test_tt_reactivity_units():
    m3 = reactivity(10, "keV", "m^3/s", "T+T=2n+a")
    cm3 = reactivity(10, "keV", "cm^3/s", "T+T=2n+a")
    mm3 = reactivity(10, "keV", "mm^3/s", "T+T=2n+a")
    assert cm3 == approx(m3 * 1e6, rel=1e-12, abs=0)
    assert mm3 == approx(cm3 * 1e3, rel=1e-12, abs=0)
    assert reactivity(0.01, "MeV", "m^3/s", "T+T=2n+a") == approx(m3, rel=1e-12, abs=0)


@mark.parametrize("ion_temperature_kev", [0.1, 1000.0])
def test_tt_reactivity_table_endpoints_are_valid(ion_temperature_kev):
    result = reactivity(ion_temperature_kev, "keV", "m^3/s", "T+T=2n+a")
    assert math.isfinite(result)
    assert result > 0


def test_relative_reaction_rates_invalid_equation():
    with raises(ValueError, match="equations are supported"):
        relative_reaction_rates(10, "keV", equation="unknown")


def test_relative_reaction_rates_only_one_fraction_given_raises():
    # tritium_fraction defaults to 0.5, so 0.7 + 0.5 does not sum to 1
    with raises(ValueError, match="sum to 1.0"):
        relative_reaction_rates(10, "keV", deuterium_fraction=0.7)
    with raises(ValueError, match="sum to 1.0"):
        relative_reaction_rates(10, "keV", tritium_fraction=0.7)


def test_relative_reaction_rates_sum_tolerance_boundary():
    # the check allows a deviation of 1e-6 from 1.0
    relative_reaction_rates(10, "keV", deuterium_fraction=0.5, tritium_fraction=0.5 + 5e-7)
    with raises(ValueError, match="sum to 1.0"):
        relative_reaction_rates(10, "keV", deuterium_fraction=0.5, tritium_fraction=0.5 + 2e-6)


def test_relative_reaction_rates_sadler_van_belle_differs_from_bosch_hale():
    bh = relative_reaction_rates(10, "keV", equation="Bosch-Hale")
    svb = relative_reaction_rates(10, "keV", equation="Sadler-Van Belle")
    assert svb[0] != bh[0]
    assert svb[0] == approx(bh[0], rel=0.05, abs=0)


def test_relative_neutron_rates_pure_deuterium():
    # D+D=p+T emits no neutrons and there is no DT or TT
    result = relative_neutron_rates(10, "keV", deuterium_fraction=1.0, tritium_fraction=0.0)
    assert result == approx([0.0, 1.0, 0.0, 0.0], abs=1e-12)


def test_relative_neutron_rates_pure_tritium():
    result = relative_neutron_rates(10, "keV", deuterium_fraction=0.0, tritium_fraction=1.0)
    assert result == approx([0.0, 0.0, 0.0, 1.0], abs=1e-12)


def test_relative_neutron_rates_custom_fuel_weights_by_neutrons():
    reactions = relative_reaction_rates(10, "keV", 0.9, 0.1)
    neutrons = relative_neutron_rates(10, "keV", 0.9, 0.1)
    weights = [reactions[0], reactions[1], 0.0, 2 * reactions[3]]
    assert neutrons == approx([w / sum(weights) for w in weights], rel=1e-12, abs=0)
    assert sum(neutrons) == approx(1.0)


def test_relative_neutron_rates_temperature_units():
    expected = relative_neutron_rates(10, "keV")
    for temperature, units in ((10e3, "eV"), (0.01, "MeV")):
        assert relative_neutron_rates(temperature, units) == approx(expected, rel=1e-12, abs=0)


def test_relative_neutron_rates_sadler_van_belle():
    result = relative_neutron_rates(10, "keV", equation="Sadler-Van Belle")
    assert sum(result) == approx(1.0)


def test_relative_neutron_rates_error_propagation():
    with raises(ValueError, match="between 0 and 1"):
        relative_neutron_rates(10, "keV", 1.5, -0.5)
    with raises(ValueError, match="sum to 1.0"):
        relative_neutron_rates(10, "keV", 0.5, 0.4)
    with raises(ValueError, match="Invalid temperature units"):
        relative_neutron_rates(10, "K")
    with raises(ValueError, match="equations are supported"):
        relative_neutron_rates(10, "keV", equation="unknown")


def test_fusion_energy_per_reaction_units():
    assert fusion_energy_per_reaction("D+T=n+a", "keV") == approx(17589.299865, rel=1e-12, abs=0)
    assert fusion_energy_per_reaction("D+T=n+a", "GeV") == approx(17.589299865e-3, rel=1e-12, abs=0)
    assert fusion_energy_per_reaction("D+T=n+a", "eV") == approx(17589299.865, rel=1e-12, abs=0)


def test_fusion_energy_per_reaction_invalid_units():
    with raises(ValueError, match="Invalid energy units"):
        fusion_energy_per_reaction("D+T=n+a", "J")


def test_neutron_rate_zero_power_is_zero():
    assert neutron_rate_from_fusion_power(0.0, 10, "keV") == 0.0


def test_neutron_rate_nan_power_raises():
    with raises(ValueError, match="fusion_power"):
        neutron_rate_from_fusion_power(float("nan"), 10, "keV")


def test_neutron_rate_hard_coded_reference_for_1_gw_dt():
    # 1 GW / (17.589 MeV * 1.602e-13 J/MeV) = 3.548e20 reactions per second,
    # one neutron each. DD contributions at 10 keV are below 1 %.
    assert neutron_rate_from_fusion_power(1e9, 10, "keV") == approx(3.548e20, rel=0.01)


def test_neutron_rate_scales_linearly_with_power():
    base = neutron_rate_from_fusion_power(1e6, 10, "keV")
    assert neutron_rate_from_fusion_power(5e6, 10, "keV") == approx(5 * base, rel=1e-12, abs=0)


def test_neutron_rate_temperature_units_and_equation():
    expected = neutron_rate_from_fusion_power(1e9, 10, "keV")
    for temperature, units in ((10e3, "eV"), (0.01, "MeV")):
        assert neutron_rate_from_fusion_power(1e9, temperature, units) == approx(expected, rel=1e-12, abs=0)
    svb = neutron_rate_from_fusion_power(1e9, 10, "keV", equation="Sadler-Van Belle")
    assert svb == approx(expected, rel=0.05, abs=0)


@mark.parametrize("units", ENERGY_UNIT_FROM_KEV.keys())
def test_mean_neutron_energy_units(units):
    mean_kev = mean_neutron_energy(10, "keV", "keV")
    mean = mean_neutron_energy(10, "keV", units)
    assert mean == approx(mean_kev * ENERGY_UNIT_FROM_KEV[units], rel=1e-12, abs=0)


def test_mean_neutron_energy_invalid_units():
    with raises(ValueError, match="Invalid energy units"):
        mean_neutron_energy(10, "keV", "J")


def test_mean_neutron_energy_invalid_temperature_units():
    with raises(ValueError, match="Invalid temperature units"):
        mean_neutron_energy(10, "K")


def test_mean_neutron_energy_temperature_units():
    expected = mean_neutron_energy(10, "keV", "eV")
    for temperature, units in ((10e3, "eV"), (0.01, "MeV")):
        assert mean_neutron_energy(temperature, units, "eV") == approx(expected, rel=1e-12, abs=0)


def test_mean_neutron_energy_error_propagation():
    with raises(ValueError, match="sum to 1.0"):
        mean_neutron_energy(10, "keV", "MeV", 0.5, 0.4)
    with raises(ValueError, match="equations are supported"):
        mean_neutron_energy(10, "keV", "MeV", equation="unknown")


def test_tt_spectrum_invalid_units():
    with raises(ValueError, match="Invalid temperature units"):
        neutron_energy_spectrum([1e6], 10, temperature_units="K")
    with raises(ValueError, match="Invalid energy units"):
        neutron_energy_spectrum([1e6], 10, "keV", "J")


def test_tt_spectrum_negative_temperature_raises():
    with raises(ValueError, match="positive"):
        neutron_energy_spectrum([1e6], -10, "keV")


def test_tt_spectrum_nan_energy_raises():
    with raises(ValueError, match="non-negative"):
        neutron_energy_spectrum([1.0, float("nan")], 10, "keV", "MeV")


def test_tt_spectrum_empty_energies():
    assert neutron_energy_spectrum([], 10, "keV", "MeV") == []


def test_tt_spectrum_preserves_input_order():
    energies = [1.0, 6.0, 3.0, 8.0, 2.0]
    result = neutron_energy_spectrum(energies, 10, "keV", "MeV")
    for energy, value in zip(energies, result):
        assert neutron_energy_spectrum([energy], 10, "keV", "MeV")[0] == approx(value, rel=1e-12, abs=0)


def test_tt_spectrum_is_non_negative_and_zero_above_kinematic_limit():
    energies = [0.5 * i for i in range(21)]  # 0 to 10 MeV
    assert all(s >= 0.0 for s in neutron_energy_spectrum(energies, 10, "keV", "MeV"))
    # the center of mass spectrum ends at 9.4 MeV, far above that at 10 keV
    assert neutron_energy_spectrum([20.0], 10, "keV", "MeV") == [0.0]


def test_tt_spectrum_energy_units():
    in_mev = neutron_energy_spectrum([4.0], 10, "keV", "MeV")[0]
    in_kev = neutron_energy_spectrum([4000.0], 10, "keV", "keV")[0]
    assert in_kev == approx(in_mev / 1e3, rel=1e-12, abs=0)


# NOTE AdrGos: Known issues, written as the behavior we want. strict=True makes the
# suite flag them as soon as the code is fixed so the marker can be removed.

@mark.xfail(strict=True, reason="ion_temperature <= 0.0 is False for NaN and inf, so they pass validation")
@mark.parametrize("bad_temperature", [float("nan"), float("inf")])
@mark.parametrize(
    "function",
    [
        lambda t: reactivity(t, "keV"),
        lambda t: neutron_energy_mean_and_std_dev(t, "keV"),
    ],
)
def test_non_finite_temperature_raises(function, bad_temperature):
    with raises(ValueError):
        function(bad_temperature)


@mark.xfail(
    strict=True,
    raises=ValueError,
    reason="reaction_fractions always evaluates the TT reactivity, which raises below 0.1 keV even with no tritium",
)
def test_pure_deuterium_below_tt_table_range():
    result = relative_reaction_rates(0.05, "keV", deuterium_fraction=1.0, tritium_fraction=0.0)
    assert sum(result) == approx(1.0)