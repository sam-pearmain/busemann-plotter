#[allow(dead_code)]

mod utils;
mod taylormaccoll;

fn main() {
    let gamma: f64 = 1.4;
    let design_efficiency: f64 = 0.9;

    // determine k^2 from design pressure efficiency
    let k_squared = solve_for_k_squared(
        design_efficiency, 
        gamma, 
        x1, // erm
        x2, 
        None, 
        None,
    );
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
        utils::equations::calc_downstream_mach(k_squared, upstream_mach, gamma)
    };
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
