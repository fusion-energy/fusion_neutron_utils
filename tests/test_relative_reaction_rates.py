import pytest
from fusion_neutron_utils import relative_reaction_rates

def test_relative_reaction_rates_default_fractions():
    result = relative_reaction_rates(ion_temperature=10.0)
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_custom_fractions():
    result = relative_reaction_rates(ion_temperature=10.0, dt_fraction=0.3, dd_fraction=0.7)
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_temperature_units():
    result = relative_reaction_rates(ion_temperature=10.0, temperature_units="keV")
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_custom_equation():
    result = relative_reaction_rates(ion_temperature=10.0, equation="Bosch-Hale")
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(x, float) for x in result)

def test_relative_reaction_rates_invalid_fractions():
    with pytest.raises(TypeError):
        relative_reaction_rates(ion_temperature=10.0, dt_fraction="invalid", dd_fraction="invalid")

def test_relative_reaction_rates_zero_fractions():
    with pytest.raises(ValueError):
        relative_reaction_rates(ion_temperature=10.0, dt_fraction=0.0, dd_fraction=0.0)

def test_relative_reaction_rates_non_sum_1_fractions():
    with pytest.raises(ValueError):
        relative_reaction_rates(ion_temperature=10.0, dt_fraction=0.5, dd_fraction=0.4)

if __name__ == "__main__":
    pytest.main()