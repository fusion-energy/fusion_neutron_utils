use pyo3::prelude::*;

mod tt_data;

#[pyfunction(signature = (ion_temperature, temperature_units=None, neutron_energy_units=None, reaction=None))]
/// Calculate the average neutron energy for a given ion temperature and reaction.
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', neutron_energy_units='eV', reaction='D+T=n+a')")]
fn neutron_energy_mean_and_std_dev(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    neutron_energy_units: Option<&str>,
    reaction: Option<&str>,
) -> PyResult<(f64, f64)> {
    let reaction = reaction.unwrap_or("D+T=n+a");

    // values from Ballabio paper
    let (a_1, a_2, a_3, a_4, mean) = match reaction {
        "D+D=n+He3" => (4.69515, -0.040729, 0.47, 0.81844, 2.4486858678216934e6),
        "D+T=n+a" => (5.30509, 0.0024736, 1.84, 1.3818, 14028394.744466662),
        _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("reaction must be either 'D+D=n+He3' or 'D+T=n+a'")),
    };

    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units)?; // Ballabio equation accepts KeV units
    if ion_temperature_kev <= 0.0 {
        return Err(value_error("Ion temperature must be positive and non-zero"));
    }

    // units of mean_delta are in put into ev with the 1000 multiplication
    let mean_delta = 1000.0 *( a_1 * ion_temperature_kev.powf(2.0 / 3.0) / (1.0 + a_2 * ion_temperature_kev.powf(a_3)) + a_4 * ion_temperature_kev);

    let mean_adjusted = mean + mean_delta;  

    let mean_scaled =  scale_energy_in_kev_to_requested_units(mean_adjusted/1e3, neutron_energy_units)?;

    let (w_0, a_1, a_2, a_3, a_4) = match reaction {
        "D+D=n+He3" => (82.542, 1.7013e-3, 0.16888, 0.49, 7.9460e-4),
        "D+T=n+a" => (177.259, 5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5),
        _ => unreachable!(), // This case is already handled above
    };

    let delta = a_1 * ion_temperature_kev.powf(2.0 / 3.0) / (1.0 + a_2 * ion_temperature_kev.powf(a_3)) + a_4 * ion_temperature_kev;

    // 2.3548200450309493 on the line below comes from equation 2* math.sqrt(math.log(2)*2)
    let variance = ((w_0 * (1.0 + delta)).powi(2) * ion_temperature_kev) / 2.3548200450309493_f64.powi(2);
    // let variance = variance * 1e6; // converting keV^2 back to eV^2
    let std_dev = variance.sqrt();
    let std_dev = scale_energy_in_kev_to_requested_units(std_dev, neutron_energy_units)?;

    Ok((mean_scaled, std_dev))
}



#[pyfunction(signature = (ion_temperature, temperature_units=None, reactivity_units=None, reaction=None, equation=None))]
/// Thermal reactivity <sigma v> assuming a Maxwellian ion temperature distribution.
///
/// D+T=n+a, D+D=n+He3 and D+D=p+T use the Bosch-Hale parametrization (valid for
/// 0.2 keV to 100 keV), D+T=n+a can also use Sadler-Van Belle. T+T=2n+a uses
/// the tabulated reactivity from Hale (0.1 keV to 1000 keV, interpolated
/// log-log), taken from NeSST.
///
/// Args:
///     ion_temperature_kev (float): Ion temperature.
///     temperature_units
///
/// Returns:
///     float: The thermal reactivity in m^3/s.
///
/// Source:
/// Sec. 5.2, Eqn. (12) - (14), Table VII in
/// H.-s. Bosch and G. M. Dale,
/// "Improved Formulas for Fusion Cross-Sections and Thermal Reactivities",
/// Nucl. Fusion 32, 611 (1992)
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', reactivity_units='m^3/s', reaction='D+T=n+a', equation='Bosch-Hale')")]
fn reactivity(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    reactivity_units: Option<&str>,
    reaction: Option<&str>,
    equation: Option<&str>,
) -> PyResult<f64> {

    if ion_temperature <= 0.0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Ion temperature must be positive and non-zero"));
    }

    let reaction = reaction.unwrap_or("D+T=n+a");
    let equation_str = equation.unwrap_or("Bosch-Hale");

    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units)?;

    if reaction == "T+T=2n+a" {
        if !matches!(equation, None | Some("Hale")) {
            return Err(value_error("Only the 'Hale' equation is supported for the 'T+T=2n+a' reaction"));
        }
        return scale_reactivity_units(tt_reactivity(ion_temperature_kev)?, reactivity_units);
    }

    let sigma_thermal_reactivity_scaled = if equation_str == "Bosch-Hale"{
        let (c1, c2, c3, c4, c5, c6, c7, gamov, mrc2) = if reaction == "D+T=n+a" {
            (
                1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5,
                34.3827, 1_124_656.0,
            )
        } else if reaction == "D+D=p+T" {
            (
                5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0,
                31.3970, 937_814.0,
            )
        } else if reaction == "D+D=n+He3" {
            (
                5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0,
                31.3970, 937_814.0,
            )
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Only 'D+T=n+a', 'D+D=p+T', 'D+D=n+He3' and 'T+T=2n+a' reactions are supported"));
        };

        let sigma_thermal_reactivity =  bosch_and_hale_equations(c1, c2, c3, c4, c5, c6, c7, gamov, mrc2, ion_temperature_kev)?;
        scale_reactivity_units(sigma_thermal_reactivity, reactivity_units)
    }else if equation_str == "Sadler-Van Belle"{
        if reaction != "D+T=n+a" {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Only 'D+T=n+a' reaction is supported for 'Sadler-Van Belle' equation"));
        }
        let sigma_thermal_reactivity = sadler_van_belle(ion_temperature_kev)?;
        scale_reactivity_units(sigma_thermal_reactivity, reactivity_units)
    }else{
        Err(value_error("Only 'Bosch-Hale' and 'Sadler-Van Belle' equations are supported"))
    };
    sigma_thermal_reactivity_scaled
    
}


#[pyfunction(signature = (ion_temperature, temperature_units=None, deuterium_fraction=None, tritium_fraction=None, equation=None))]
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', deuterium_fraction=0.5, tritium_fraction=0.5, equation='Bosch-Hale')")]
/// Calculate the relative reaction rates for given ion temperature and fuel composition.
///
/// The reaction rate per unit volume between ion species i and j is
/// n_i * n_j * <sigma v>_ij / (1 + delta_ij), so for fuel ion fractions f_D
/// and f_T the DT rate is proportional to f_D * f_T * <sigma v>_DT, each
/// DD branch to f_D^2 / 2 * <sigma v>_DD and TT to f_T^2 / 2 * <sigma v>_TT.
///
/// Parameters
/// ----------
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// deuterium_fraction : float, optional
///     The fraction of fuel ions that are deuterium. Default is 0.5.
/// tritium_fraction : float, optional
///     The fraction of fuel ions that are tritium. Default is 0.5.
/// equation : str, optional
///     The equation used for the D+T reactivity, 'Bosch-Hale' or
///     'Sadler-Van Belle'. The D+D reactivities always use Bosch-Hale.
///     Default is 'Bosch-Hale'.
///
/// Returns
/// -------
/// List[float]
///     The fractions of all reactions that are DT, DD (n+He3), DD (p+T) and TT, in that order.
///
/// Examples
/// --------
/// >>> relative_reaction_rates(10e3)
/// [dt_rate, dd_n_rate, dd_p_rate, tt_rate]
///
/// >>> relative_reaction_rates(10.0, temperature_units='keV', deuterium_fraction=0.9, tritium_fraction=0.1)
/// [dt_rate, dd_n_rate, dd_p_rate, tt_rate]
fn relative_reaction_rates(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    deuterium_fraction: Option<f64>,
    tritium_fraction: Option<f64>,
    equation: Option<&str>,
) -> Result<Vec<f64>, PyErr> {
    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units)?;
    Ok(reaction_fractions(ion_temperature_kev, deuterium_fraction, tritium_fraction, equation)?.to_vec())
}

/// Fractions of all reactions that are DT, DD (n+He3), DD (p+T) and TT, in the
/// order of REACTIONS.
fn reaction_fractions(
    ion_temperature_kev: f64,
    deuterium_fraction: Option<f64>,
    tritium_fraction: Option<f64>,
    equation: Option<&str>,
) -> PyResult<[f64; 4]> {
    let deuterium_fraction = deuterium_fraction.unwrap_or(0.5);
    let tritium_fraction = tritium_fraction.unwrap_or(0.5);

    let equation_str = equation.unwrap_or("Bosch-Hale");

    if !(0.0..=1.0).contains(&deuterium_fraction) || !(0.0..=1.0).contains(&tritium_fraction) {
        return Err(value_error("deuterium_fraction and tritium_fraction must each be between 0 and 1"));
    }
    let total_fraction = deuterium_fraction + tritium_fraction;
    let tol : f64 = 0.000001;
    if !(total_fraction > 1. - tol && total_fraction < 1. + tol) {
        return Err(value_error("The deuterium_fraction + tritium_fraction do not sum to 1.0 and are not within a small tolerance (+/-0.000001)"));
    }

    let dt_reactivity = reactivity(ion_temperature_kev, Some("keV"), Some("m^3/s"), Some("D+T=n+a"), Some(equation_str))?;
    let dd_reactivity_1 = reactivity(ion_temperature_kev, Some("keV"), Some("m^3/s"), Some("D+D=n+He3"), Some("Bosch-Hale"))?;
    let dd_reactivity_2 = reactivity(ion_temperature_kev, Some("keV"), Some("m^3/s"), Some("D+D=p+T"), Some("Bosch-Hale"))?;
    let tt_reactivity = reactivity(ion_temperature_kev, Some("keV"), Some("m^3/s"), Some("T+T=2n+a"), None)?;

    let dt_rate = deuterium_fraction * tritium_fraction * dt_reactivity;
    let dd_rate_1 = 0.5 * deuterium_fraction.powi(2) * dd_reactivity_1;
    let dd_rate_2 = 0.5 * deuterium_fraction.powi(2) * dd_reactivity_2;
    let tt_rate = 0.5 * tritium_fraction.powi(2) * tt_reactivity;

    let total_rate = dt_rate + dd_rate_1 + dd_rate_2 + tt_rate;

    Ok([dt_rate / total_rate, dd_rate_1 / total_rate, dd_rate_2 / total_rate, tt_rate / total_rate])
}


fn sadler_van_belle(ion_temperature: f64) -> Result<f64, PyErr> {
    let c = [
        2.5663271e-18,
        19.983026,
        2.5077133e-2,
        2.5773408e-3,
        6.1880463e-5,
        6.6024089e-2,
        8.1215505e-3,
    ];

    let u = 1.0 - ion_temperature * (c[2] + ion_temperature * (c[3] - c[4] * ion_temperature))
        / (1.0 + ion_temperature * (c[5] + c[6] * ion_temperature));

    let val = c[0]
        * ((-c[1] * (u / ion_temperature).powf(1.0 / 3.0)).exp())
        / (u.powf(5.0 / 6.0) * ion_temperature.powf(2.0 / 3.0));

    // the fit gives m^3/s, converting to cm^3/s
    Ok(val * 1.0e6)
}

/// T+T thermal reactivity in cm^3/s, log-log interpolation of the Hale table.
fn tt_reactivity(ion_temperature_kev: f64) -> PyResult<f64> {
    let temperatures = &tt_data::TT_REACTIVITY_TEMPERATURE_KEV;
    let values = &tt_data::TT_REACTIVITY_CM3_PER_S;
    let (t_min, t_max) = (temperatures[0], temperatures[temperatures.len() - 1]);
    if !(t_min..=t_max).contains(&ion_temperature_kev) {
        return Err(value_error(&format!(
            "The T+T reactivity is tabulated from {t_min} keV to {t_max} keV, not {ion_temperature_kev} keV"
        )));
    }
    let i = temperatures.partition_point(|&t| t <= ion_temperature_kev).clamp(1, temperatures.len() - 1);
    let (t0, t1) = (temperatures[i - 1], temperatures[i]);
    let (v0, v1) = (values[i - 1], values[i]);
    let fraction = (ion_temperature_kev / t0).ln() / (t1 / t0).ln();
    Ok(v0 * (v1 / v0).powf(fraction))
}

fn bosch_and_hale_equations(c1: f64, c2: f64, c3: f64, c4: f64, c5: f64, c6: f64, c7: f64, gamov: f64, mrc2: f64, ion_temperature_kev: f64) -> Result<f64, PyErr> {
    // Equation 13
    let theta: f64 = ion_temperature_kev * (1.0 - (ion_temperature_kev * (c2 + ion_temperature_kev * (c4 + ion_temperature_kev * c6))) / (1.0 + ion_temperature_kev * (c3 + ion_temperature_kev * (c5 + ion_temperature_kev * c7)))) .powi(-1);

    // Equation 14
    let xi: f64 = (gamov.powi(2) / (4.0 * theta)).powf(1.0 / 3.0);

    // Equation 12, gives cm^3/s
    let sigma_thermal_reactivity: f64 = c1 * theta * (xi / (mrc2 * ion_temperature_kev.powi(3))).sqrt() * (-3.0 * xi).exp();

    Ok(sigma_thermal_reactivity)
}

// Triton to neutron mass ratio, CODATA 2018 mass energy equivalents
const TRITON_NEUTRON_MASS_RATIO: f64 = 2808.92113668 / 939.56542194;

#[pyfunction(signature = (energies, ion_temperature, temperature_units=None, neutron_energy_units=None, reaction=None))]
#[pyo3(text_signature = "(energies, ion_temperature, temperature_units='eV', neutron_energy_units='eV', reaction='T+T=2n+a')")]
/// Neutron energy spectrum dN/dE for the T+T=2n+a reaction.
///
/// The T+T reaction has three bodies in the final state, so the neutron
/// spectrum is a continuum rather than a peak. The centre of mass spectrum
/// (Brune fit, from NeSST) is broadened for the ion temperature following
/// Appelbe et al., High Energy Density Physics 2016, as done in NeSST's
/// dNdE_TT. For D+T and D+D use neutron_energy_mean_and_std_dev.
///
/// Parameters
/// ----------
/// energies : list of float
///     Neutron energies at which to evaluate the spectrum.
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// neutron_energy_units : str, optional
///     The units of the energies. Default is 'eV'.
/// reaction : str, optional
///     Only 'T+T=2n+a' is supported. Default is 'T+T=2n+a'.
///
/// Returns
/// -------
/// List[float]
///     dN/dE at each energy, normalised to integrate to 1 over energy in
///     neutron_energy_units.
fn neutron_energy_spectrum(
    energies: Vec<f64>,
    ion_temperature: f64,
    temperature_units: Option<&str>,
    neutron_energy_units: Option<&str>,
    reaction: Option<&str>,
) -> PyResult<Vec<f64>> {
    if reaction.unwrap_or("T+T=2n+a") != "T+T=2n+a" {
        return Err(value_error("Only the 'T+T=2n+a' reaction is supported, use neutron_energy_mean_and_std_dev for 'D+T=n+a' and 'D+D=n+He3'"));
    }
    let ion_temperature_ev = scale_temperature_units_to_kev(ion_temperature, temperature_units)? * 1e3;
    if ion_temperature_ev <= 0.0 {
        return Err(value_error("Ion temperature must be positive and non-zero"));
    }
    // energy of one neutron_energy_unit in eV
    let unit_in_ev = scale_energy_in_kev_to_requested_units(1.0, Some("eV"))? / scale_energy_in_kev_to_requested_units(1.0, neutron_energy_units)?;
    if energies.iter().any(|&e| !(e >= 0.0)) {
        return Err(value_error("Neutron energies must be non-negative"));
    }
    let energies_ev: Vec<f64> = energies.iter().map(|e| e * unit_in_ev).collect();
    Ok(tt_spectrum_per_ev(&energies_ev, ion_temperature_ev).iter().map(|s| s * unit_in_ev).collect())
}

/// TT neutron spectrum dN/dE in 1/eV at energies in eV, normalised to 1.
fn tt_spectrum_per_ev(energies_ev: &[f64], ion_temperature_ev: f64) -> Vec<f64> {
    // Appelbe et al. HEDP 2016, as implemented in NeSST TT_spectrum_model.spec.
    // The kernel exp(-a (sqrt(E) - sqrt(E'))^2) is a Gaussian in sqrt(E) with
    // standard deviation 1 / sqrt(2 a).
    let a = 2.0 * TRITON_NEUTRON_MASS_RATIO / ion_temperature_ev;

    // Grid points per kernel width, the minimum and maximum refinement of the
    // data grid, and the number of kernel standard deviations summed over.
    // Chosen so the spectrum is within 6e-5 of the peak of a fully converged
    // calculation from 0.1 keV to 300 keV.
    const REFINE_PER_KERNEL_WIDTH: f64 = 4.0;
    const MIN_REFINE: usize = 8;
    const MAX_REFINE: usize = 64;
    const KERNEL_CUTOFF: f64 = 6.0;

    // The tabulated centre of mass spectrum is linearly interpolated onto a
    // finer grid. The kernel width in energy at 1 MeV, 2 sqrt(1 MeV) / sqrt(2 a),
    // sets the refinement, which is highest at low ion temperatures where the
    // kernel is narrow compared with the 37 keV data spacing. A minimum
    // refinement resolves the 1 / sqrt(E') factor at low E' when the kernel
    // is wide.
    let data_spacing = (tt_data::TT_COM_SPECTRUM_ENERGY_MEV[1] - tt_data::TT_COM_SPECTRUM_ENERGY_MEV[0]) * 1e6;
    let kernel_width_at_1_mev = 2.0 * 1e3 / (2.0 * a).sqrt();
    let refine = ((REFINE_PER_KERNEL_WIDTH * data_spacing / kernel_width_at_1_mev).ceil() as usize).clamp(MIN_REFINE, MAX_REFINE);

    // Integrals use the trapezium rule. The E = 0 point is skipped as the
    // kernel has a 1 / sqrt(E) factor and the spectrum is zero there.
    let com_data = &tt_data::TT_COM_SPECTRUM;
    let d_energy = data_spacing / refine as f64;
    let n_fine = (com_data.len() - 1) * refine + 1;
    let mut com_sqrt_energies = Vec::with_capacity(n_fine - 1);
    let mut com_weights = Vec::with_capacity(n_fine - 1);
    for k in 1..n_fine {
        let (i, step) = (k / refine, (k % refine) as f64 / refine as f64);
        let value = if step == 0.0 { com_data[i] } else { com_data[i] + step * (com_data[i + 1] - com_data[i]) };
        let trapezium = if k == n_fine - 1 { 0.5 } else { 1.0 };
        com_sqrt_energies.push((k as f64 * d_energy).sqrt());
        com_weights.push(value * trapezium * d_energy / (k as f64 * d_energy).sqrt());
    }
    let integral: f64 = com_weights.iter().zip(com_sqrt_energies.iter()).map(|(w, s)| w * s).sum();
    let norm = 0.5 * (a / std::f64::consts::PI).sqrt() / integral;

    // Only centre of mass energies within KERNEL_CUTOFF standard deviations
    // contribute, the kernel is below exp(-KERNEL_CUTOFF^2 / 2) beyond that.
    let window = KERNEL_CUTOFF / (2.0 * a).sqrt();

    energies_ev
        .iter()
        .map(|&e| {
            let sqrt_e = e.sqrt();
            let lo = com_sqrt_energies.partition_point(|&s| s < sqrt_e - window);
            let hi = com_sqrt_energies.partition_point(|&s| s <= sqrt_e + window);
            let per_ev: f64 = com_sqrt_energies[lo..hi]
                .iter()
                .zip(com_weights[lo..hi].iter())
                .map(|(&sqrt_e_com, &w)| (-a * (sqrt_e - sqrt_e_com).powi(2)).exp() * w)
                .sum();
            norm * per_ev
        })
        .collect()
}

/// Mean energy in eV of the TT neutron spectrum, by the trapezium rule on a
/// 10 keV grid up to 12 MeV (the centre of mass spectrum ends at 9.4 MeV).
fn tt_mean_energy_ev(ion_temperature_ev: f64) -> f64 {
    let energies_ev: Vec<f64> = (0..=1200).map(|i| i as f64 * 1e4).collect();
    let spectrum = tt_spectrum_per_ev(&energies_ev, ion_temperature_ev);
    let (mut first_moment, mut zeroth_moment) = (0.0, 0.0);
    for i in 0..energies_ev.len() - 1 {
        let width = energies_ev[i + 1] - energies_ev[i];
        first_moment += 0.5 * (energies_ev[i] * spectrum[i] + energies_ev[i + 1] * spectrum[i + 1]) * width;
        zeroth_moment += 0.5 * (spectrum[i] + spectrum[i + 1]) * width;
    }
    first_moment / zeroth_moment
}

// The reactions in the order returned by relative_reaction_rates, with the
// energy released per reaction (keV) and neutrons emitted per reaction. The
// energies are Q-values from the AME2020 atomic mass excesses,
// M. Wang et al 2021 Chinese Phys. C 45 030003.
const REACTIONS: [(&str, f64, u32); 4] = [
    ("D+T=n+a", 17589.299865, 1),
    ("D+D=n+He3", 3268.90885, 1),
    ("D+D=p+T", 4032.663826, 0),
    ("T+T=2n+a", 11332.06981, 2),
];

const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

fn reaction_index(reaction: &str) -> PyResult<usize> {
    REACTIONS
        .iter()
        .position(|(name, _, _)| *name == reaction)
        .ok_or_else(|| value_error("reaction must be one of 'D+T=n+a', 'D+D=n+He3', 'D+D=p+T' or 'T+T=2n+a'"))
}

#[pyfunction(signature = (reaction, energy_units=None))]
#[pyo3(text_signature = "(reaction, energy_units='eV')")]
/// Energy released per fusion reaction (the Q-value), from the AME2020 atomic
/// mass excesses.
///
/// Parameters
/// ----------
/// reaction : str
///     One of 'D+T=n+a', 'D+D=n+He3', 'D+D=p+T' or 'T+T=2n+a'.
/// energy_units : str, optional
///     The units of the returned energy, 'eV', 'keV', 'MeV' or 'GeV'. Default is 'eV'.
///
/// Returns
/// -------
/// float
///     The energy released per reaction.
fn fusion_energy_per_reaction(reaction: &str, energy_units: Option<&str>) -> PyResult<f64> {
    scale_energy_in_kev_to_requested_units(REACTIONS[reaction_index(reaction)?].1, energy_units)
}

#[pyfunction]
#[pyo3(text_signature = "(reaction)")]
/// Number of neutrons emitted per fusion reaction.
///
/// Parameters
/// ----------
/// reaction : str
///     One of 'D+T=n+a', 'D+D=n+He3', 'D+D=p+T' or 'T+T=2n+a'.
///
/// Returns
/// -------
/// int
///     1 for D+T=n+a and D+D=n+He3, 0 for D+D=p+T and 2 for T+T=2n+a.
fn neutrons_per_reaction(reaction: &str) -> PyResult<u32> {
    Ok(REACTIONS[reaction_index(reaction)?].2)
}

/// Fractions of all neutrons that come from DT, DD (n+He3), DD (p+T) and TT.
fn neutron_fractions(reaction_fractions: &[f64; 4]) -> [f64; 4] {
    let per_reaction: Vec<f64> = reaction_fractions
        .iter()
        .zip(REACTIONS.iter())
        .map(|(fraction, (_, _, neutrons))| fraction * *neutrons as f64)
        .collect();
    let total: f64 = per_reaction.iter().sum();
    [per_reaction[0] / total, per_reaction[1] / total, per_reaction[2] / total, per_reaction[3] / total]
}

#[pyfunction(signature = (ion_temperature, temperature_units=None, deuterium_fraction=None, tritium_fraction=None, equation=None))]
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', deuterium_fraction=0.5, tritium_fraction=0.5, equation='Bosch-Hale')")]
/// Fractions of all neutrons emitted by each reaction, for a given ion
/// temperature and fuel composition.
///
/// Unlike relative_reaction_rates this counts neutrons, so D+D=p+T contributes
/// none and each T+T=2n+a reaction contributes two. These are the relative
/// source strengths to use for each reaction in OpenMC.
///
/// Parameters
/// ----------
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// deuterium_fraction : float, optional
///     The fraction of fuel ions that are deuterium. Default is 0.5.
/// tritium_fraction : float, optional
///     The fraction of fuel ions that are tritium. Default is 0.5.
/// equation : str, optional
///     The equation used for the D+T reactivity. Default is 'Bosch-Hale'.
///
/// Returns
/// -------
/// List[float]
///     The fractions of all neutrons from DT, DD (n+He3), DD (p+T) and TT, in that order.
fn relative_neutron_rates(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    deuterium_fraction: Option<f64>,
    tritium_fraction: Option<f64>,
    equation: Option<&str>,
) -> PyResult<Vec<f64>> {
    let ion_temperature_kev = scale_temperature_units_to_kev(ion_temperature, temperature_units)?;
    let fractions = reaction_fractions(ion_temperature_kev, deuterium_fraction, tritium_fraction, equation)?;
    Ok(neutron_fractions(&fractions).to_vec())
}

#[pyfunction(signature = (fusion_power, ion_temperature, temperature_units=None, deuterium_fraction=None, tritium_fraction=None, equation=None))]
#[pyo3(text_signature = "(fusion_power, ion_temperature, temperature_units='eV', deuterium_fraction=0.5, tritium_fraction=0.5, equation='Bosch-Hale')")]
/// Neutron emission rate for a given fusion power, ion temperature and fuel
/// composition.
///
/// The fusion power is divided by the mean energy released per reaction, for
/// the mix of reactions at this temperature and fuel composition, and
/// multiplied by the mean number of neutrons per reaction. All four reactions
/// are included. This is the number to multiply OpenMC tally results (per
/// source neutron) by to get rates.
///
/// Parameters
/// ----------
/// fusion_power : float
///     The fusion power in W.
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// deuterium_fraction : float, optional
///     The fraction of fuel ions that are deuterium. Default is 0.5.
/// tritium_fraction : float, optional
///     The fraction of fuel ions that are tritium. Default is 0.5.
/// equation : str, optional
///     The equation used for the D+T reactivity. Default is 'Bosch-Hale'.
///
/// Returns
/// -------
/// float
///     Neutrons emitted per second.
fn neutron_rate_from_fusion_power(
    fusion_power: f64,
    ion_temperature: f64,
    temperature_units: Option<&str>,
    deuterium_fraction: Option<f64>,
    tritium_fraction: Option<f64>,
    equation: Option<&str>,
) -> PyResult<f64> {
    if !(fusion_power >= 0.0) {
        return Err(value_error("fusion_power must be non-negative"));
    }
    let ion_temperature_kev = scale_temperature_units_to_kev(ion_temperature, temperature_units)?;
    let fractions = reaction_fractions(ion_temperature_kev, deuterium_fraction, tritium_fraction, equation)?;
    let mut energy_per_reaction_joules = 0.0;
    let mut neutrons_per_reaction = 0.0;
    for (fraction, (_, q_kev, neutrons)) in fractions.iter().zip(REACTIONS.iter()) {
        energy_per_reaction_joules += fraction * q_kev * 1e3 * ELEMENTARY_CHARGE;
        neutrons_per_reaction += fraction * *neutrons as f64;
    }
    Ok(fusion_power / energy_per_reaction_joules * neutrons_per_reaction)
}

#[pyfunction(signature = (ion_temperature, temperature_units=None, neutron_energy_units=None, deuterium_fraction=None, tritium_fraction=None, equation=None))]
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', neutron_energy_units='eV', deuterium_fraction=0.5, tritium_fraction=0.5, equation='Bosch-Hale')")]
/// Mean energy of the emitted neutrons for a given ion temperature and fuel
/// composition, averaged over all neutrons from all reactions.
///
/// The D+T and D+D means use the Ballabio fits (as in
/// neutron_energy_mean_and_std_dev) and the T+T mean is the mean of the
/// spectrum from neutron_energy_spectrum.
///
/// Parameters
/// ----------
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// neutron_energy_units : str, optional
///     The units of the returned energy. Default is 'eV'.
/// deuterium_fraction : float, optional
///     The fraction of fuel ions that are deuterium. Default is 0.5.
/// tritium_fraction : float, optional
///     The fraction of fuel ions that are tritium. Default is 0.5.
/// equation : str, optional
///     The equation used for the D+T reactivity. Default is 'Bosch-Hale'.
///
/// Returns
/// -------
/// float
///     The mean neutron energy.
fn mean_neutron_energy(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    neutron_energy_units: Option<&str>,
    deuterium_fraction: Option<f64>,
    tritium_fraction: Option<f64>,
    equation: Option<&str>,
) -> PyResult<f64> {
    let ion_temperature_kev = scale_temperature_units_to_kev(ion_temperature, temperature_units)?;
    let fractions = neutron_fractions(&reaction_fractions(ion_temperature_kev, deuterium_fraction, tritium_fraction, equation)?);
    let (dt_mean_ev, _) = neutron_energy_mean_and_std_dev(ion_temperature_kev, Some("keV"), Some("eV"), Some("D+T=n+a"))?;
    let (dd_mean_ev, _) = neutron_energy_mean_and_std_dev(ion_temperature_kev, Some("keV"), Some("eV"), Some("D+D=n+He3"))?;
    let tt_mean_ev = if fractions[3] > 0.0 { tt_mean_energy_ev(ion_temperature_kev * 1e3) } else { 0.0 };
    let mean_ev = fractions[0] * dt_mean_ev + fractions[1] * dd_mean_ev + fractions[3] * tt_mean_ev;
    scale_energy_in_kev_to_requested_units(mean_ev / 1e3, neutron_energy_units)
}

fn value_error(message: &str) -> PyErr {
    PyErr::new::<pyo3::exceptions::PyValueError, _>(message.to_string())
}

fn scale_reactivity_units(sigma_thermal_reactivity: f64, reactivity_units: Option<&str>) -> PyResult<f64> {
    match reactivity_units.unwrap_or("m^3/s") {
        "m^3/s" => Ok(sigma_thermal_reactivity * 1.0e-6),
        "cm^3/s" => Ok(sigma_thermal_reactivity),
        "mm^3/s" => Ok(sigma_thermal_reactivity * 1.0e3),
        _ => Err(value_error("Invalid reaction rate units, accepted values are 'm^3/s', 'cm^3/s', 'mm^3/s'")),
    }
}

fn scale_temperature_units_to_kev(ion_temperature: f64, temperature_units: Option<&str>) -> PyResult<f64> {
    match temperature_units.unwrap_or("eV") {
        "keV" => Ok(ion_temperature),
        "eV" => Ok(ion_temperature * 1e-3),
        "MeV" => Ok(ion_temperature * 1e3),
        "GeV" => Ok(ion_temperature * 1e6),
        _ => Err(value_error("Invalid temperature units, accepted values are 'eV', 'keV', 'MeV' or 'GeV'")),
    }
}

fn scale_energy_in_kev_to_requested_units(energy_in_kev: f64, energy_units: Option<&str>) -> PyResult<f64> {
    match energy_units.unwrap_or("eV") {
        "eV" => Ok(energy_in_kev * 1e3), // converting keV to eV
        "keV" => Ok(energy_in_kev),
        "MeV" => Ok(energy_in_kev / 1e3), // converting keV to MeV
        "GeV" => Ok(energy_in_kev / 1e6), // converting keV to GeV
        _ => Err(value_error("Invalid energy units, accepted values are 'eV', 'keV', 'MeV' or 'GeV'")),
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn fusion_neutron_utils(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reactivity, m)?)?;
    m.add_function(wrap_pyfunction!(relative_reaction_rates, m)?)?;
    m.add_function(wrap_pyfunction!(neutron_energy_mean_and_std_dev, m)?)?;
    m.add_function(wrap_pyfunction!(neutron_energy_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(fusion_energy_per_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(neutrons_per_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(relative_neutron_rates, m)?)?;
    m.add_function(wrap_pyfunction!(neutron_rate_from_fusion_power, m)?)?;
    m.add_function(wrap_pyfunction!(mean_neutron_energy, m)?)?;
    Ok(())
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scale_temperature_units_to_kev_eV() {
        let ion_temperature = 1000.0; // 1000 eV
        let temperature_units = Some("eV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units).unwrap();
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_keV() {
        let ion_temperature = 1.0; // 1 keV
        let temperature_units = Some("keV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units).unwrap();
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_MeV() {
        let ion_temperature = 0.001; // 1 MeV
        let temperature_units = Some("MeV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units).unwrap();
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_GeV() {
        let ion_temperature = 0.000001; // 1 GeV
        let temperature_units = Some("GeV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units).unwrap();
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_invalid_units() {
        let ion_temperature = 1000.0;
        let temperature_units = Some("K");
        assert!(scale_temperature_units_to_kev(ion_temperature, temperature_units).is_err());
    }
}


#[cfg(test)]
mod tests2 {
    use super::*;

    #[test]
    fn test_scale_to_ev() {
        let energy_in_kev = 1.0;
        let target_unit = "eV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit)).unwrap();
        assert_eq!(result, 1000.0); // 1 keV = 1000 eV
    }

    #[test]
    fn test_scale_to_mev() {
        let energy_in_kev = 1000.0;
        let target_unit = "MeV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit)).unwrap();
        assert_eq!(result, 1.0); // 1000 keV = 1 MeV
    }

    #[test]
    fn test_scale_to_gev() {
        let energy_in_kev = 1_000_000.0;
        let target_unit = "GeV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit)).unwrap();
        assert_eq!(result, 1.0); // 1,000,000 keV = 1 GeV
    }

    #[test]
    fn test_invalid_unit() {
        let energy_in_kev = 1.0;
        let target_unit = "invalid";
        assert!(scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit)).is_err());
    }
}

#[cfg(test)]
mod tests3 {
    use super::*;

    #[test]
    fn test_scale_to_m3_per_s() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("m^3/s");
        let result = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units).unwrap();
        assert_eq!(result, 1.0 * 1.0e-6); // 1.0 * 1.0e-6 = 1.0e-6
    }

    #[test]
    fn test_scale_to_cm3_per_s() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("cm^3/s");
        let result = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units).unwrap();
        assert_eq!(result, 1.0); // 1.0 cm^3/s = 1.0 cm^3/s
    }


    #[test]
    fn test_invalid_unit() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("invalid");
        assert!(scale_reactivity_units(sigma_thermal_reactivity, reactivity_units).is_err());
    }
}