pub fn bisection (
    f: &impl Fn(f64) -> f64,
    x1: f64, // 1st solution bound
    x2: f64, // 2nd solution bound
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    // default tolerance 1e-9 unless otherwise given
    let tolerance = tolerance.unwrap_or(1e-9);
    let max_iters = max_iters.unwrap_or(200);

    // initialise lower and upper bound according to given bounds
    let (mut lowerbound, mut upperbound) = if x1 < x2 { (x1, x2) } else { (x2, x1) };

    // iterate
    for _ in 0..max_iters {
        let midpoint = (upperbound + lowerbound) / 2.0;
        
        // check convergence
        if f(midpoint).abs() < tolerance || (upperbound - lowerbound) / 2.0 < tolerance {
            return midpoint;
        }

        // update bounds
        if (f(midpoint) * f(lowerbound)) > 0.0 {
            lowerbound = midpoint;
        } else {
            upperbound = midpoint;
        }
    }

    panic!("solution not converged");
}
