#[allow(dead_code)]

use crate::utils::error::FlowError;

#[derive(Debug, Clone)]
pub struct VelocityVector {
    pub radial_component: f64,      // u
    pub tangential_component: f64,  // v
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
    pub velocity_vector: VelocityVector,
    pub radial_distance: f64,
    pub theta: f64,
}

pub fn streamline(
    velocity_vector: &VelocityVector,
    r: f64 // the radial distance 
) -> Result<f64, FlowError> {
    if !velocity_vector.valid_mach_number() {
        return Err(FlowError::InvalidMachNumber);
    }

    let r_derivative: f64 = 
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
        let k1: VelocityVectorDerivative = taylor_maccoll(&k1_velocity_vector, current_theta, gamma)?;
        let k1_radial: f64 = h * k1.radial_derivative;
        let k1_tangential: f64 = h * k1.tangential_derivative;
        let k1_contour: f64 = h * streamline(&k1_velocity_vector, current_radial_distance)?;

        // second runge-kutta constant
        let k2_velocity_vector: VelocityVector = 
            VelocityVector {
                radial_component: current_radial_velocity + (0.5 * k1_radial),
                tangential_component: current_tangential_velocity + (0.5 * k1_tangential),
            };
        let k2: VelocityVectorDerivative = taylor_maccoll(&k2_velocity_vector, current_theta + (0.5 * h), gamma)?;
        let k2_radial: f64 = h * k2.radial_derivative;
        let k2_tangential: f64 = h * k2.tangential_derivative;
        let k2_contour: f64 = h * streamline(&k2_velocity_vector, current_radial_distance + (0.5 * k1_contour))?;

        // third runge-kutta constant
        let k3_velocity_vector: VelocityVector = 
            VelocityVector {
                radial_component: current_radial_velocity + (0.5 * k2_radial),
                tangential_component: current_tangential_velocity + (0.5 * k2_tangential),
            };
        let k3: VelocityVectorDerivative = taylor_maccoll(&k3_velocity_vector, current_theta + (0.5 * h), gamma)?;
        let k3_radial: f64 = h * k3.radial_derivative;
        let k3_tangential: f64 = h * k3.tangential_derivative;
        let k3_conour: f64 = h * streamline(&k3_velocity_vector, current_radial_distance + (0.5 * k2_contour))?;

        // fourth runge-kutta constant
        let k4_velocity_vector: VelocityVector =
            VelocityVector {
                radial_component: current_radial_velocity + k3_radial,
                tangential_component: current_tangential_velocity + k3_tangential,
            };
        let k4: VelocityVectorDerivative = taylor_maccoll(&k4_velocity_vector, current_theta + h, gamma)?;
        let k4_radial: f64 = h * k4.radial_derivative;
        let k4_tangential: f64 = h * k4.tangential_derivative;
        let k4_contour: f64 = h * streamline(&k4_velocity_vector, current_radial_distance + k3_conour)?;

        // calculate subsequent velocity components, radial distance and theta
        let next_radial_velocity: f64 = 
            current_radial_velocity + (1.0 / 6.0) *
            (k1_radial + 2.0 * k2_radial + 2.0 * k3_radial + k4_radial);
        let next_tangential_velocity: f64 = 
            current_tangential_velocity + (1.0 / 6.0) *
            (k1_tangential + 2.0 * k2_tangential + 2.0 * k3_tangential + k4_tangential);
        let next_radial_distance: f64 = 
            current_radial_distance + (1.0 / 6.0) *
            (k1_contour + 2.0 * k2_contour + 2.0 * k3_conour + k4_contour);
        let next_theta: f64 = current_theta + h;

        // break clause
        let cross_stream_mach: f64 = 
            next_radial_velocity * next_theta.sin() +
            next_tangential_velocity * next_theta.cos();
        if cross_stream_mach >= 0.0 {
            // freestream condition reached
            break;
        }

        // append results to results vec
        results.push(
            TaylorMaccollResult {
                velocity_vector: VelocityVector {
                    radial_component: next_radial_velocity,
                    tangential_component: next_tangential_velocity,
                },
                radial_distance: next_radial_distance,
                theta: next_theta,
            }
        );

        // update current values with their subsequent value and loop
        current_radial_velocity = next_radial_velocity;
        current_tangential_velocity = next_tangential_velocity;
        current_radial_distance = next_radial_distance;
        current_theta = next_theta;
    }

    Ok(results)
}