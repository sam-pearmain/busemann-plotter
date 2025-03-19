#![allow(dead_code)]

#[derive(Debug)]
pub struct Polynomial {
    order: usize,
    coeffs: Vec<f64>,
}

impl Polynomial {
    pub fn new(order: usize) -> Self {
        Polynomial { order, coeffs: Vec::with_capacity(order) }
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "polynomial order: {}, coeffs: ", self.order)?;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if i == 0 {
                write!(f, "{}", coeff)?;
            } else {
                write!(f, " {}x^{}", coeff, i)?;
            }
        }
        Ok(())
    }
}

pub fn polyfit(x_data: &[f64], y_data: &[f64], order: usize) -> Polynomial {
    use nalgebra::{DMatrix, DVector};
    
    let a_order = order + 1;
    let n = y_data.len();

    let mut a_vec = Vec::with_capacity(n * a_order);
    for &x in x_data.iter().take(n) {
        let mut power = 1.0;
        for _ in 0..a_order {
            a_vec.push(power);
            power *= x;
        }
    }

    let a_matrix = DMatrix::from_row_slice(n, a_order, &a_vec);
    let b_vector = DVector::from_row_slice(&y_data[..n]);

    let ata = a_matrix.transpose() * &a_matrix;
    let atb = a_matrix.transpose() * b_vector;
    let coeffs_vec = ata.lu().solve(&atb)
        .expect("failed to solve the normal equations");

    let mut poly = Polynomial::new(order);
    poly.coeffs = coeffs_vec.iter().cloned().collect();
    poly
}