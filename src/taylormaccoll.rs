#[allow(dead_code)]

use std::collections::HashMap;
use crate::utils::error::FlowError;

pub struct VelocityVector {
    radial_component: f64,      // u
    tangential_component: f64,  // v
}

impl VelocityVector {
    pub fn valid_mach_number(&self) -> bool {
        self.get_mach_number() >= 1.0
    }

    pub fn get_mach_number(&self) -> f64 {
        (self.radial_component.powi(2) + self.tangential_component.powi(2)).sqrt()
    }
}

pub struct VelocityVectorDerivative {
    radial_derivative: f64,     // du / dθ
    tangential_derivative: f64, // dv / dθ
}

pub fn taylor_maccoll(
    velocity_vector: VelocityVector,
    theta: f64,
    gamma: f64,
) -> Result<VelocityVectorDerivative, FlowError>{
    // get radial and tangential velocity components
    let u: f64 = velocity_vector.radial_component;
    let v: f64 = velocity_vector.tangential_component;
    
    if !velocity_vector.valid_mach_number() {
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

pub fn solve_taylor_maccoll(
    initial_velocity_vector: VelocityVector,
    initial_theta: f64,
    final_theta: f64,

) -> Result<HashMap<f64, VelocityVector, f64>, FlowError> {
    if !initial_velocity_vector.valid_mach_number() {
        return Err(FlowError::InvalidMachNumber);
    }
    
}