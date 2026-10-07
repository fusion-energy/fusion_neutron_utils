import pytest
from fusion_neutron_utils import relative_reaction_rates

def test_relative_reaction_rates_default_fractions():
    result = relative_reaction_rates(ion_temperature=10e3)
    assert isinstance(result, list)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_custom_fractions():
    result = relative_reaction_rates(ion_temperature=10e3, deuterium_fraction=0.7, tritium_fraction=0.3)
    assert isinstance(result, list)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_temperature_units():
    result = relative_reaction_rates(ion_temperature=10.0, temperature_units="keV")
    assert isinstance(result, list)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_custom_equation():
    result = relative_reaction_rates(ion_temperature=10e3, equation="Bosch-Hale")
    assert isinstance(result, list)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_invalid_fractions():
    with pytest.raises(TypeError):
        relative_reaction_rates(ion_temperature=10.0, deuterium_fraction="invalid", tritium_fraction="invalid")

def test_relative_reaction_rates_zero_fractions():
    with pytest.raises(ValueError):
        relative_reaction_rates(ion_temperature=10.0, deuterium_fraction=0.0, tritium_fraction=0.0)

def test_relative_reaction_rates_non_sum_1_fractions():
    with pytest.raises(ValueError):
        relative_reaction_rates(ion_temperature=10.0, deuterium_fraction=0.5, tritium_fraction=0.4)

if __name__ == "__main__":
    pytest.main()

def test_relative_reaction_rates_fuel_composition():
    # DT scales with f_D * f_T and DD with f_D^2 / 2
    dt_a, dd_a, _, _ = relative_reaction_rates(10.0, "keV", deuterium_fraction=0.5, tritium_fraction=0.5)
    dt_b, dd_b, _, _ = relative_reaction_rates(10.0, "keV", deuterium_fraction=0.9, tritium_fraction=0.1)
    ratio_a = dd_a / dt_a
    ratio_b = dd_b / dt_b
    expected = (0.9**2 / (0.9 * 0.1)) / (0.5**2 / (0.5 * 0.5))
    assert ratio_b / ratio_a == pytest.approx(expected, rel=1e-9)


def test_relative_reaction_rates_pure_deuterium():
    dt, dd_n, dd_p, tt = relative_reaction_rates(10.0, "keV", deuterium_fraction=1.0, tritium_fraction=0.0)
    assert dt == 0.0
    assert tt == 0.0
    assert dd_n + dd_p == pytest.approx(1.0)


def test_relative_reaction_rates_pure_tritium():
    assert relative_reaction_rates(10.0, "keV", deuterium_fraction=0.0, tritium_fraction=1.0) == [0.0, 0.0, 0.0, 1.0]


def test_relative_reaction_rates_fraction_out_of_range():
    with pytest.raises(ValueError, match="between 0 and 1"):
        relative_reaction_rates(10.0, "keV", deuterium_fraction=1.5, tritium_fraction=-0.5)


def test_relative_reaction_rates_sadler_van_belle():
    result = relative_reaction_rates(10.0, "keV", equation="Sadler-Van Belle")
    assert sum(result) == pytest.approx(1.0)
