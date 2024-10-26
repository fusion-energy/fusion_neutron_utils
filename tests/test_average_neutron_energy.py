import NeSST as nst
from fusion_neutron_utils import neutron_energy_mean_and_std_dev
from pytest import approx

def test_mean_energy_with_nesst():

    ion_temperature = 19e3

    nesst_dt_mean, nest_dt_std_dev, nest_dt_var = nst.DTprimspecmoments(ion_temperature)
    nesst_dd_mean, nest_dd_std_dev, nest_dd_var = nst.DDprimspecmoments(ion_temperature)

    fnu_dd_mean, fnu_dd_std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+D=n+He3',
        ion_temperature=ion_temperature,
        temperature_units='eV',
        neutron_energy_units='eV'
    )

    fnu_dt_mean, fnu_dt_std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+T=n+a',
        ion_temperature=ion_temperature,
        temperature_units='eV',
        neutron_energy_units='eV'
    )
    # assert nesst_dt_mean == fnu_dt_mean
    assert nesst_dd_mean == fnu_dd_mean
    assert nest_dt_std_dev == fnu_dt_std_dev
    assert nest_dd_std_dev == fnu_dd_std_dev

def test_mean_energy():
    mean, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+D=n+He3',
        ion_temperature=20e3,
        temperature_units='eV',
        neutron_energy_units='eV'
    )
    assert mean == approx(2.5e6, abs = 0.2e6)

    mean, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+D=n+He3',
        ion_temperature=20,
        temperature_units='keV',
        neutron_energy_units='MeV'
    )
    assert mean == approx(2.5, abs = 0.2)

    mean, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+T=n+a',
        ion_temperature=20e3,
        temperature_units='eV',
    )
    assert mean == approx(14.06e6, abs = 0.04e6)

def test_mean_energy_increases_with_ion_temperature():
    mean_cold, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+D=n+He3',
        ion_temperature=20e3,
        temperature_units='eV',
    )

    mean_hot, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+D=n+He3',
        ion_temperature=40e3,
        temperature_units='eV',
    )
    assert mean_cold < mean_hot

def test_mean_energy_units():
    mean_kev, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+T=n+a',
        ion_temperature=30,
        temperature_units='keV',
        neutron_energy_units='keV'
    )

    mean_ev, std_dev = neutron_energy_mean_and_std_dev(
        reaction='D+T=n+a',
        ion_temperature=30e3,
        temperature_units='eV',
        neutron_energy_units='eV'
    )
    assert mean_kev == mean_ev/1e3

# def test_