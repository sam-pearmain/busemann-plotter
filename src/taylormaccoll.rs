#[allow(dead_code)]

use crate::utils::error::FlowError;

pub struct VelocityVector {
    radial_component: f64,      // u
    tangential_component: f64,  // v
}

pub struct VelocityVectorDerivative {
    radial_derivative: f64,     // du / dθ
    tangential_derivative: f64, // dv / dθ
}

pub fn taylor_maccoll(
    mach_vector: VelocityVector,
    theta: f64,
    gamma: f64,
) -> Result<VelocityVectorDerivative, FlowError>{
    // get radial and tangential velocity components
    let u: f64 = mach_vector.radial_component;
    let v: f64 = mach_vector.tangential_component;
    
    let mach: f64 = (u.powi(2) + v.powi(2)).sqrt();
    if mach < 1.0 {
        return Err(FlowError::InvalidMachNumber)
    }

    let u_derivative: f64 = 
        v + ((gamma - 1.0) / 2.0 * u * v) *
            (u + (v * 1.0 / theta.tan())) / 
            (v.powi(2) - 1.0);
    let v_derivative: f64 = 
        - u + (1.0 + ((gamma - 1.0) / 2.0) * v.powi(2)) *
            (u + (v * 1.0 / theta.tan())) / 
            (v.powi(2) - 1.0);

    Ok(VelocityVectorDerivative{
        radial_derivative: u_derivative,
        tangential_derivative: v_derivative,
    })
}