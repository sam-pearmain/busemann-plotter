#[allow(dead_code)]

use crate::utils::error::FlowError;

pub fn calc_deflection_angle(
    upstream_mach: f64, 
    shock_angle: f64,
    gamma: f64
) -> Result<f64, FlowError> {
    let tan_deflection_angle: f64 =
        2.0 / shock_angle.tan() * 
        (upstream_mach.powi(2) * shock_angle.sin().powi(2) - 1.0) / 
        (upstream_mach.powi(2) * (gamma + (2.0 * shock_angle).cos()) + 2.0);
    Ok(tan_deflection_angle.atan())
}

pub fn calc_pressure_ratio(
    k_squared: f64,
    gamma: f64,
) -> f64 {
    ((gamma + 2.0) * k_squared / ((gamma - 1.0) * k_squared + 2.0)).powf(gamma / (gamma - 1.0)) *
    ((gamma + 1.0) / ((2.0 * gamma * k_squared) - gamma + 1.0)).powf(1.0 / (gamma - 1.0))
}

pub fn calc_downstream_mach(
    k_squared: f64,
    upstream_mach: f64,
    gamma: f64,
) -> f64 {
    let numerator: f64 = 
        (gamma + 1.0).powi(2) * upstream_mach.powi(2) * k_squared -
        (4.0 * (k_squared - 1.0) * (gamma * k_squared + 1.0));
    let denominator: f64 = 
        (2.0 * gamma * k_squared - (gamma - 1.0)) *
        ((gamma - 1.0) * k_squared + 2.0);
    numerator / denominator
}