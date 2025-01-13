#[allow(dead_code)]

use crate::utils::error::FlowError;

#[derive(Debug, Clone)]
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

#[derive(Debug, Clone)]
pub struct VelocityVectorDerivative {
    radial_derivative: f64,     // du / dθ
    tangential_derivative: f64, // dv / dθ
}

pub struct TaylorMaccollResult {
    // a struct to organise the results from integrating taylor maccoll equations
    velocity_vector: VelocityVector,
    radial_distance: f64,
    theta: f64,
}

pub fn streamline(
    velocity_vector: &VelocityVector,
    r: f64 // the radial distance 
) -> Result<f64, FlowError> {
    if !velocity_vector.valid_mach_number() {
        return Err(FlowError::InvalidMachNumber);
    }

    let r_derivative = 
        r * velocity_vector.radial_component / velocity_vector.tangential_component;

    Ok(r_derivative)
}

pub fn taylor_maccoll(
    velocity_vector: &VelocityVector,
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
    initial_r: f64,
    gamma: f64,
    steps: u32,
) -> Result<Vec<TaylorMaccollResult>, FlowError> {
    // 4th order runge kutta integration of taylor maccoll equations
    if !initial_velocity_vector.valid_mach_number() {
        return Err(FlowError::InvalidMachNumber);
    }

    // set step size
    let h: f64 = (final_theta - initial_theta) / steps as f64;
    
    // vector to store results
    let mut results: Vec<TaylorMaccollResult> = Vec::new();

    // push initial values to results
    results.push(TaylorMaccollResult{
        velocity_vector: initial_velocity_vector.clone(),
        radial_distance: initial_r,
        theta: initial_theta,
    });

    // starting conditions for integration
    let mut current_radial_velocity: f64 = initial_velocity_vector.radial_component;
    let mut current_tangential_velocity: f64 = initial_velocity_vector.tangential_component;
    let mut current_radial_distance: f64 = initial_r;
    let mut current_theta: f64 = initial_theta;

    for _ in 0..steps {
        // first runge-kutta constant
        let k1_velocity_vector: VelocityVector = 
            VelocityVector {
                radial_component: current_radial_velocity,
                tangential_component: current_tangential_velocity,
            };
        let k1: VelocityVectorDerivative = 
            taylor_maccoll(
                &k1_velocity_vector,
                current_theta,
                gamma,
            )?;
        let k1_radial: f64 = h * k1.radial_derivative;
        let k1_tangential: f64 = h * k1.tangential_derivative;
        let k1_contour: f64 = h * streamline(&k1_velocity_vector, current_radial_distance)?;

        // second runge-kutta constant
        let k2_velocity_vector: VelocityVector = 
            VelocityVector {
                radial_component: current_radial_velocity + (0.5 * k1_radial),
                tangential_component: current_tangential_velocity + (0.5 * k1_radial),
            };
        let k2: VelocityVectorDerivative = 
            taylor_maccoll(
                &k2_velocity_vector,
                current_theta + (0.5 * h),
                gamma,
            )?;
    }

    Ok(results)
}