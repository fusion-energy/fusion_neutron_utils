use pyo3::prelude::*;

#[pyfunction(signature = (ion_temperature, temperature_units=None, neutron_energy_units=None, reaction=None))]
/// Calculate the average neutron energy for a given ion temperature and reaction.
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', neutron_energy_units='eV', reaction='D+T=n+a')")]
fn neutron_energy_mean_and_std_dev(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    neutron_energy_units: Option<&str>,
    reaction: Option<&str>,
) -> PyResult<(f64, f64)> {
    // values from Ballabio paper
    let (a_1, a_2, a_3, a_4, mean) = match reaction {
        Some("D+D=n+He3") => (4.69515, -0.040729, 0.47, 0.81844, 2.4486858678216934e6),
        Some("D+T=n+a") => (5.30509, 0.0024736, 1.84, 1.3818, 14028394.744466662),
        _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("reaction must be either 'D+D=n+He3' or 'D+T=n+a'")),
    };

    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units); // Ballabio equation accepts KeV units

    // units of mean_delta are in put into ev with the 1000 multiplication
    let mean_delta = 1000.0 *( a_1 * ion_temperature_kev.powf(0.66666666) / (1.0 + a_2 * ion_temperature_kev.powf(a_3)) + a_4 * ion_temperature_kev);

    let mean_adjusted = mean + mean_delta;  

    let mean_scaled =  scale_energy_in_kev_to_requested_units(mean_adjusted/1e3, neutron_energy_units);

    let (w_0, a_1, a_2, a_3, a_4) = match reaction {
        Some("D+D=n+He3") => (82.542, 1.7013e-3, 0.16888, 0.49, 7.9460e-4),
        Some("D+T=n+a") => (177.259, 5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5),
        _ => unreachable!(), // This case is already handled above
    };

    let delta = a_1 * ion_temperature_kev.powf(2.0 / 3.0) / (1.0 + a_2 * ion_temperature_kev.powf(a_3)) + a_4 * ion_temperature_kev;

    // 2.3548200450309493 on the line below comes from equation 2* math.sqrt(math.log(2)*2)
    let variance = ((w_0 * (1.0 + delta)).powi(2) * ion_temperature_kev) / 2.3548200450309493_f64.powi(2);
    // let variance = variance * 1e6; // converting keV^2 back to eV^2
    let std_dev = variance.sqrt();
    let std_dev = scale_energy_in_kev_to_requested_units(std_dev, neutron_energy_units);

    Ok((mean_scaled, std_dev))
}



#[pyfunction(signature = (ion_temperature, temperature_units=None, reactivity_units=None, reaction=None, equation=None))]
/// Bosch-Hale parametrization of D+T thermal reactivity, assuming a Maxwellian ion
/// temperature distribution. If ion_temperature_kev is given in keV, the returned
/// <sigma v> is in m^3/s. This is valid for 0.2 keV <= ion_temperature_kev <= 100 keV.
///
/// D + T -> T(3.56 MeV) + n(14.03 MeV)
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

    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units);

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
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Only 'D+T=n+a', 'D+D=p+T', and 'D+D=n+He3' reactions are supported"));
        };

        let sigma_thermal_reactivity =  bosch_and_hale_equations(c1, c2, c3, c4, c5, c6, c7, gamov, mrc2, ion_temperature_kev)?;
        let scaling_factor = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units);
        Ok(sigma_thermal_reactivity*scaling_factor)
    }else if equation_str == "Sadler-Van Belle"{
        if reaction != "D+T=n+a" {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Only 'D+T=n+a' reaction is supported for 'Sadler-Van Belle' equation"));
        }
        let sigma_thermal_reactivity = sadler_van_belle(ion_temperature_kev)?;
        let scaling_factor = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units);
        Ok(sigma_thermal_reactivity*scaling_factor)
    }else{
        panic!("Only 'Bosch-Hale' and 'Sadler-Van Belle' equations are supported");
    };
    sigma_thermal_reactivity_scaled
    
}


#[pyfunction(signature = (ion_temperature, temperature_units=None, dt_fraction=None, dd_fraction=None, equation=None))]
#[pyo3(text_signature = "(ion_temperature, temperature_units='eV', dt_fraction=0.5, dd_fraction=0.5 equation='Bosch-Hale')")]
/// Calculate the relative reaction rates for given ion temperature and isotope fractions.
///
/// Parameters
/// ----------
/// ion_temperature : float
///     The ion temperature.
/// temperature_units : str, optional
///     The units of the ion temperature. Default is 'eV'.
/// dt_fraction : float, optional
///     The fraction of DT reactions. Default is 0.5.
/// dd_fraction : float, optional
///     The fraction of DD reactions. Default is 0.5.
/// equation : str, optional
///     The equation to use for reactivity calculations. Default is 'Bosch-Hale'.
///
/// Returns
/// -------
/// List[float]
///     A list containing the relative reaction rates for DT, DD (n+He3), and DD (p+T) reactions.
///
/// Examples
/// --------
/// >>> relative_reaction_rates(10.0)
/// [dt_reactivity, dd_reactivity_1, dd_reactivity_2]
///
/// >>> relative_reaction_rates(10.0, temperature_units='K', dt_fraction=0.3, dd_fraction=0.7, equation='Custom-Equation')
/// [dt_reactivity, dd_reactivity_1, dd_reactivity_2]
fn relative_reaction_rates(
    ion_temperature: f64,
    temperature_units: Option<&str>,
    dt_fraction: Option<f64>,
    dd_fraction: Option<f64>,
    equation: Option<&str>,
) -> Result<Vec<f64>, PyErr> {

    let ion_temperature_kev: f64 = scale_temperature_units_to_kev(ion_temperature, temperature_units);

    let dt_fraction = dt_fraction.unwrap_or(0.5);
    let dd_fraction = dd_fraction.unwrap_or(0.5);

    let equation_str = equation.unwrap_or("Bosch-Hale");
    
    let total_fraction = dt_fraction + dd_fraction;
    let tol : f64 = 0.000001;
    if !(total_fraction > 1. - tol && total_fraction < 1. + tol) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("The dt_fraction + dd_fraction do not sum to 1.0 and are not within a small tolerance (+/-0.000001)"));
    }

    let dt_reactivity = reactivity(
        ion_temperature_kev,
        Some("keV"),
        Some("m^3/s"),
        Some("D+T=n+a"),
        Some(equation_str)
        // Some("Bosch-Hale")
    )?;

    let dd_reactivity_1 = reactivity(
        ion_temperature_kev,
        Some("keV"),
        Some("m^3/s"),
        Some("D+D=n+He3"),
        Some(equation_str)
        // Some("Bosch-Hale")
    )?;

    let dd_reactivity_2 = reactivity(
        ion_temperature_kev,
        Some("keV"),
        Some("m^3/s"),
        Some("D+D=p+T"),
        Some(equation_str)
        // Some("Bosch-Hale")
    )?;

    let total_reactivity = dt_reactivity * dt_fraction + dd_reactivity_1 * dd_fraction + dd_reactivity_2 * dd_fraction;
    let dt_normalized = dt_reactivity * dt_fraction / total_reactivity;
    let dd1_normalized = dd_reactivity_1 * dd_fraction / total_reactivity;
    let dd2_normalized = dd_reactivity_2 * dd_fraction / total_reactivity;

    Ok(vec![dt_normalized, dd1_normalized, dd2_normalized])

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

    Ok(val)
}

fn bosch_and_hale_equations(c1: f64, c2: f64, c3: f64, c4: f64, c5: f64, c6: f64, c7: f64, gamov: f64, mrc2: f64, ion_temperature_kev: f64) -> Result<f64, PyErr> {
    // Equation 13
    let theta: f64 = ion_temperature_kev * (1.0 - (ion_temperature_kev * (c2 + ion_temperature_kev * (c4 + ion_temperature_kev * c6))) / (1.0 + ion_temperature_kev * (c3 + ion_temperature_kev * (c5 + ion_temperature_kev * c7)))) .powi(-1);

    // Equation 14
    let xi: f64 = (gamov.powi(2) / (4.0 * theta)).powf(1.0 / 3.0);

    // Equation 12
    let sigma_thermal_reactivity: f64 = c1 * theta * (xi / (mrc2 * ion_temperature_kev.powi(3))).sqrt() * (-3.0 * xi).exp();

    Ok(sigma_thermal_reactivity * 1.0e-6)
}

fn scale_reactivity_units(sigma_thermal_reactivity: f64, reactivity_units: Option<&str>) -> f64 {
    match reactivity_units.unwrap_or("m^3/s") {
        "m^3/s" => sigma_thermal_reactivity * 1.0e-6,
        "cm^3/s" => sigma_thermal_reactivity,
        "mm^3/s" => sigma_thermal_reactivity * 1.0e3,
        _ => panic!("Invalid reaction rate units, accepted values are 'm^3/s', 'cm^3/s', 'mm^3/s'")
    }
}

fn scale_temperature_units_to_kev(ion_temperature: f64, temperature_units: Option<&str>) -> f64 {
    match temperature_units.unwrap_or("eV") {
        "keV" => ion_temperature,
        "eV" => ion_temperature * 1e-3,
        "MeV" => ion_temperature * 1e3,
        "GeV" => ion_temperature * 1e6,
        _ => panic!("Invalid temperature units, accepted values are 'eV', 'keV', 'MeV' or 'GeV'")
    } 
}

fn scale_energy_in_kev_to_requested_units(energy_in_kev: f64, temperature_units: Option<&str>) -> f64 {
    match temperature_units.unwrap_or("eV") {
        "eV" => energy_in_kev * 1e3, // converting keV to eV
        "keV" => energy_in_kev, // converting keV to MeV
        "MeV" => energy_in_kev / 1e3, // converting keV to MeV
        "GeV" => energy_in_kev / 1e6, // converting keV to GeV
        _ => panic!("Unsupported temperature units, accepted values are 'eV', 'MeV', 'GeV'"),
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn fusion_neutron_utils(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reactivity, m)?)?;
    m.add_function(wrap_pyfunction!(relative_reaction_rates, m)?)?;
    m.add_function(wrap_pyfunction!(neutron_energy_mean_and_std_dev, m)?)?;
    Ok(())
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scale_temperature_units_to_kev_eV() {
        let ion_temperature = 1000.0; // 1000 eV
        let temperature_units = Some("eV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units);
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_keV() {
        let ion_temperature = 1.0; // 1 keV
        let temperature_units = Some("keV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units);
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_MeV() {
        let ion_temperature = 0.001; // 1 MeV
        let temperature_units = Some("MeV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units);
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    fn test_scale_temperature_units_to_kev_GeV() {
        let ion_temperature = 0.000001; // 1 GeV
        let temperature_units = Some("GeV");
        let result = scale_temperature_units_to_kev(ion_temperature, temperature_units);
        assert_eq!(result, 1.0); // 1 keV
    }

    #[test]
    #[should_panic(expected = "Invalid temperature units")]
    fn test_scale_temperature_units_to_kev_invalid_units() {
        let ion_temperature = 1000.0;
        let temperature_units = Some("K");
        scale_temperature_units_to_kev(ion_temperature, temperature_units);
    }
}


#[cfg(test)]
mod tests2 {
    use super::*;

    #[test]
    fn test_scale_to_ev() {
        let energy_in_kev = 1.0;
        let target_unit = "eV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit));
        assert_eq!(result, 1000.0); // 1 keV = 1000 eV
    }

    #[test]
    fn test_scale_to_mev() {
        let energy_in_kev = 1000.0;
        let target_unit = "MeV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit));
        assert_eq!(result, 1.0); // 1000 keV = 1 MeV
    }

    #[test]
    fn test_scale_to_gev() {
        let energy_in_kev = 1_000_000.0;
        let target_unit = "GeV";
        let result = scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit));
        assert_eq!(result, 1.0); // 1,000,000 keV = 1 GeV
    }

    #[test]
    #[should_panic(expected = "Unsupported temperature units, accepted values are 'eV', 'MeV', 'GeV'")]
    fn test_invalid_unit() {
        let energy_in_kev = 1.0;
        let target_unit = "invalid";
        scale_energy_in_kev_to_requested_units(energy_in_kev, Some(target_unit));
    }
}

#[cfg(test)]
mod tests3 {
    use super::*;

    #[test]
    fn test_scale_to_m3_per_s() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("m^3/s");
        let result = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units);
        assert_eq!(result, 1.0 * 1.0e-6); // 1.0 * 1.0e-6 = 1.0e-6
    }

    #[test]
    fn test_scale_to_cm3_per_s() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("cm^3/s");
        let result = scale_reactivity_units(sigma_thermal_reactivity, reactivity_units);
        assert_eq!(result, 1.0); // 1.0 cm^3/s = 1.0 cm^3/s
    }


    #[test]
    #[should_panic(expected = "Invalid reaction rate units, accepted values are 'm^3/s', 'cm^3/s', 'mm^3/s'")]
    fn test_invalid_unit() {
        let sigma_thermal_reactivity = 1.0;
        let reactivity_units = Some("invalid");
        scale_reactivity_units(sigma_thermal_reactivity, reactivity_units);
    }
}