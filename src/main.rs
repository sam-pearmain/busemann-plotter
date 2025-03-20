#[allow(dead_code)]

use std::f64::consts::PI;
use csv::Writer;

mod utils;
mod taylormaccoll;

fn main() {
    let gamma: f64 = 1.0;
    let design_efficiency: f64 = 0.99;
    let design_chamber_mach: f64 = 2.0;

    // determine k^2 and M2 from design pressure efficiency
    let k_squared = solve_for_k_squared(design_efficiency, gamma, 1.0, 5.0, None, None);
    let upstream_mach = solve_for_upstream_mach(design_chamber_mach, k_squared, gamma, 0.0, 5.0, None, None);
    
    // get the shock angle for the conical shock at the throat
    let shock_angle: f64 = (k_squared.sqrt() / upstream_mach).asin();
    let initial_deflection_angle: f64 = utils::equations::calc_deflection_angle(upstream_mach, shock_angle, gamma)
        .expect("couldn't get deflection angle");
    let initial_theta: f64 = shock_angle - initial_deflection_angle;
    let initial_radial_velocity = upstream_mach * shock_angle.cos();
    let initial_tangential_velocity: f64 = -upstream_mach * shock_angle.sin();

    // solve taylor maccoll
    match taylormaccoll::solve_taylor_maccoll(
        taylormaccoll::VelocityVector{
            radial_component: initial_radial_velocity,
            tangential_component: initial_tangential_velocity, 
        }, 
        initial_theta, 
        PI, 
        1.0, 
        gamma, 
        8000,
    ) {
        Ok(results) => {
            let mut wtr = Writer::from_path("results.csv")
            .expect("failed to create csv");

            // header
            wtr.write_record(&["theta (rad)", "r (-)", "radial mach", "tangential mach", "mach number", "x", "y"])
                .expect("Failed to write header");

            // write line to csv
            let mut x_data: Vec<f64> = Vec::new();
            let mut y_data: Vec<f64> = Vec::new();
            for result in results {
                let mach_number = (
                    result.velocity_vector.radial_component.powi(2) +
                    result.velocity_vector.tangential_component.powi(2)
                ).sqrt();

                let x_coord: f64 = -result.radial_distance * result.theta.cos();
                let y_coord: f64 = result.radial_distance * result.theta.sin();
                x_data.push(x_coord);
                y_data.push(y_coord);
                
                wtr.write_record(&[
                    format!("{:.6}", result.theta),
                    format!("{:.6}", result.radial_distance),
                    format!("{:.6}", result.velocity_vector.radial_component),
                    format!("{:.6}", result.velocity_vector.tangential_component),
                    format!("{:.6}", mach_number),
                    format!("{:.6}", x_coord),
                    format!("{:.6}", y_coord),
                ])
                .expect("failed to write line");
            }
            
            // fit a polynomial
            let poly = utils::polyfit::polyfit(&x_data, &y_data, 20);
            println!("fitted polynomial: {}", poly);
            poly.plot("inlet.png", (x_data[0], *x_data.last().expect("erm")))
                .expect("idk what happened");

            // flush the csv writer
            wtr.flush().expect("flush failed");
        },
        Err(e) => {
            panic!("there was an error: {:?}", e);
        }
    }
}

fn solve_for_k_squared(
    target_pressure_ratio: f64,
    gamma: f64,
    x1: f64, // lower bound for k^2
    x2: f64, // upper bound for k^2
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    let func = |k_squared: f64| {
        utils::equations::calc_pressure_ratio(k_squared, gamma) - target_pressure_ratio
    };

    utils::numerics::bisection(&func, x1, x2, tolerance, max_iters)
}

fn solve_for_upstream_mach(
    target_mach: f64,
    k_squared: f64,
    gamma: f64,
    x1: f64,
    x2: f64,
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    let func = |upstream_mach: f64| {
        utils::equations::calc_downstream_mach(k_squared, upstream_mach, gamma) - target_mach
    };

    utils::numerics::bisection(&func, x1, x2, tolerance, max_iters)
}

