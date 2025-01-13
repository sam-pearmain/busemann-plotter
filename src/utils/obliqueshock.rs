#![allow(dead_code)]

use std::f64::consts::PI; 
use super::numerics::bisection;
use super::error::FlowError;


pub enum Input {
    UpstreamMach(f64),
    NormalUpstreamMach(f64),
    DeflectionAngle(f64), // potentially need to change to allow for strong or weak oblique shocks, currently only does weak shocks
    ShockAngle(f64),
}

pub enum Output {
    DownstreamMach,
    DeflectionAngle,
    ShockAngle,
    PressureRatio,
    DensityRatio,
    TemperatureRatio,
    StagnationPressureRatio,
    NormalUpstreamMach,
    NormalDownstreamMach,
}

#[derive(Debug)]
pub struct ObliqueShock {
    upstream_mach: f64,             // M1
    downstream_mach: f64,           // M2
    deflection_angle: f64,          // θ
    shock_angle: f64,               // β
    pressure_ratio: f64,            // p2 / p1
    density_ratio: f64,             // ρ2 / ρ1
    temperature_ratio: f64,         // T2 / T1
    stagnation_pressure_ratio: f64, // p02 / p01
    normal_upstream_mach: f64,      // M1n
    normal_downstream_mach: f64,    // M2n
}

impl ObliqueShock {
    pub fn from_mach_and_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, FlowError> {
        let deflection_angle: f64 = calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio)?;
        let downstream_mach: f64 = calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)?;
        let pressure_ratio: f64 = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let density_ratio: f64 = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let temperature_ratio: f64 = calc_temperature_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let stagnation_pressure_ratio: f64 = calc_stagnation_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let normal_upstream_mach: f64 = calc_normal_upstream_mach(upstream_mach, shock_angle)?;
        let normal_downstream_mach: f64 = calc_normal_downstream_mach(downstream_mach, shock_angle, deflection_angle)?;

        Ok(ObliqueShock {
            upstream_mach,
            downstream_mach,
            deflection_angle,
            shock_angle,
            pressure_ratio,
            density_ratio,
            temperature_ratio,
            stagnation_pressure_ratio,
            normal_upstream_mach,
            normal_downstream_mach,
        })
    }
    
    pub fn from_mach_and_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, FlowError> {
        let shock_angle = calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio)?;
        ObliqueShock::from_mach_and_shock_angle(upstream_mach, shock_angle, specific_heat_ratio)
    }

    pub fn from_mach_and_normal_mach(upstream_mach: f64, normal_upstream_mach: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, FlowError> {
        if upstream_mach <= 1.0 {
            return Err(FlowError::InvalidMachNumber);
        }
        if normal_upstream_mach <= 1.0 || normal_upstream_mach > upstream_mach {
            return Err(FlowError::InvalidMachNumber);
        }

        let shock_angle = (normal_upstream_mach / upstream_mach).asin();
        ObliqueShock::from_mach_and_shock_angle(upstream_mach, shock_angle, specific_heat_ratio)
    }
}

pub fn calculate(input: Vec<Input>, output: Output, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let oblique_shock = match input.as_slice() {
        [Input::UpstreamMach(upstream_mach), Input::ShockAngle(shock_angle)] => {
            ObliqueShock::from_mach_and_shock_angle(*upstream_mach, *shock_angle, specific_heat_ratio)?
        }
        [Input::UpstreamMach(upstream_mach), Input::DeflectionAngle(deflection_angle)] => {
            ObliqueShock::from_mach_and_deflection_angle(*upstream_mach, *deflection_angle, specific_heat_ratio)?
        }
        [Input::UpstreamMach(upstream_mach), Input::DeflectionAngle(_), Input::ShockAngle(shock_angle)] => {
            ObliqueShock::from_mach_and_shock_angle(*upstream_mach, *shock_angle, specific_heat_ratio)?
        }
        [Input::UpstreamMach(upstream_mach), Input::NormalUpstreamMach(normal_upstream_mach)] => {
            ObliqueShock::from_mach_and_normal_mach(*upstream_mach, *normal_upstream_mach, specific_heat_ratio)?
        }
        _ => {
            return Err(FlowError::MathError);
        }
    };

    match output {
        Output::DownstreamMach => Ok(oblique_shock.downstream_mach),
        Output::DeflectionAngle => Ok(oblique_shock.deflection_angle),
        Output::ShockAngle => Ok(oblique_shock.shock_angle),
        Output::PressureRatio => Ok(oblique_shock.pressure_ratio),
        Output::DensityRatio => Ok(oblique_shock.density_ratio),
        Output::TemperatureRatio => Ok(oblique_shock.temperature_ratio),
        Output::StagnationPressureRatio => Ok(oblique_shock.stagnation_pressure_ratio),
        Output::NormalUpstreamMach => Ok(oblique_shock.normal_upstream_mach),
        Output::NormalDownstreamMach => Ok(oblique_shock.normal_downstream_mach),
    }
}

fn calc_downstream_mach(upstream_mach: f64, shock_angle: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let normal_upstream_mach: f64 = calc_normal_upstream_mach(upstream_mach, shock_angle)?;

    // this is wrong
    let normal_downstream_mach: f64 = (
        (1.0 + (specific_heat_ratio - 1.0) * normal_upstream_mach.powi(2) / 2.0)
        / (specific_heat_ratio * normal_upstream_mach.powi(2) - (specific_heat_ratio - 1.0) / 2.0)
    ).sqrt();

    let downstream_mach: f64 = normal_downstream_mach / (shock_angle - deflection_angle).sin();
    Ok(downstream_mach)
}

pub fn calc_downstream_mach_from_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let deflection_angle = calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio)?;
    calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)
}

pub fn calc_downstream_mach_from_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let shock_angle = calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio)?;
    calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)
}

pub fn calc_deflection_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let tan_deflection_angle = 
        2.0 / shock_angle.tan() * 
        (upstream_mach.powi(2) * shock_angle.sin().powi(2) - 1.0) / 
        (upstream_mach.powi(2) * (specific_heat_ratio + (2.0 * shock_angle).cos()) + 2.0);
    Ok(tan_deflection_angle.atan())
}

pub fn calc_pressure_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let pressure_ratio: f64 = 
        (2.0 * specific_heat_ratio * upstream_mach.powi(2) *
        shock_angle.sin().powi(2) - (specific_heat_ratio - 1.0)) / 
        (specific_heat_ratio + 1.0);
    Ok(pressure_ratio)
}

pub fn calc_density_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let density_ratio: f64 = 
        (specific_heat_ratio + 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) /
        ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) + 2.0);
    Ok(density_ratio)
    
}

pub fn calc_temperature_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let pressure_ratio: f64 = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let density_ratio: f64 = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let temperature_ratio: f64 = pressure_ratio * (1.0 / density_ratio);
    Ok(temperature_ratio)
}

pub fn calc_stagnation_pressure_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    let stagnation_pressure_ratio: f64 =
        calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?.powf(specific_heat_ratio / (specific_heat_ratio - 1.0)) *
        (1.0 / calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?).powf(1.0 / (specific_heat_ratio - 1.0));
    Ok(stagnation_pressure_ratio)
}

pub fn calc_shock_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    if upstream_mach <= 1.0 {
        return Err(FlowError::InvalidMachNumber);
    }
    let f = |shock_angle: f64| {
        let calculated_deflection_angle = match calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio) {
            Ok(value) => value,
            Err(_) => panic!("erm what"),
        };
        return calculated_deflection_angle - deflection_angle
    };

    let lower_bound: f64 = deflection_angle;
    let upper_bound: f64 = PI / 2.0;

    let shock_angle: f64 = bisection(&f, lower_bound, upper_bound, None, None);

    if shock_angle.is_nan() {
        return Err(FlowError::MathError);
    }

    Ok(shock_angle)
}

pub fn calc_max_shock_angle(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, FlowError> {
    if upstream_mach <= 1.0 {
        return Err(FlowError::InvalidMachNumber);
    }
    
    let sin_max_shock_angle: f64 = 
        ((1.0 / (specific_heat_ratio * upstream_mach.powi(2))) * 
        (1.0 +
            ((specific_heat_ratio + 1.0) * (
                (specific_heat_ratio + 1.0) * upstream_mach.powi(4) / 16.0 +
                (specific_heat_ratio - 1.0) * upstream_mach.powi(2) / 2.0 +
                1.0
            ).sqrt())
        )).sqrt();

    if sin_max_shock_angle > 1.0 || sin_max_shock_angle < 0.0 {
        return Err(FlowError::MathError);
    }

    let shock_angle: f64 = sin_max_shock_angle.asin();
    Ok(shock_angle)
}

pub fn calc_normal_upstream_mach(upstream_mach: f64, shock_angle: f64) -> Result<f64, FlowError> {
    if upstream_mach <= 1.0 {
        return Err(FlowError::InvalidMachNumber);
    }
    Ok(upstream_mach * shock_angle.sin())
}

pub fn calc_normal_downstream_mach(downstream_mach: f64, shock_angle: f64, deflection_angle: f64) -> Result<f64, FlowError> {
    if downstream_mach <= 1.0 {
        return Err(FlowError::InvalidMachNumber);
    }
    Ok(downstream_mach * (shock_angle - deflection_angle).sin())
}