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

    pub fn eval(&self, x: f64) -> f64 {
        self.coeffs
            .iter()
            .enumerate()
            .fold(0.0, |acc, (i, &c)| acc + c * x.powi(i as i32))
    }

    pub fn plot(&self, filename: &str, x_range: (f64, f64)) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;
        let root = BitMapBackend::new(filename, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let steps = 100;
        let mut y_min = f64::INFINITY;
        let mut y_max = f64::NEG_INFINITY;
        for i in 0..=steps {
            let x = x_range.0 + (x_range.1 - x_range.0) * i as f64 / steps as f64;
            let y = self.eval(x);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
        let mut chart = ChartBuilder::on(&root)
            .caption("Polynomial Plot", ("sans-serif", 20).into_font())
            .margin(30)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(x_range.0..x_range.1, y_min..y_max)?;
        chart.configure_mesh().draw()?;
        chart.draw_series(LineSeries::new(
            (0..=steps).map(|i| {
                let x = x_range.0 + (x_range.1 - x_range.0) * i as f64 / steps as f64;
                (x, self.eval(x))
            }),
            &RED,
        ))?
        .label("Polynomial")
        .legend(|(x, y)| Path::new(vec![(x, y), (x + 20, y)], &RED));
        chart.configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
        root.present()?;
        Ok(())
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
    let coeffs_vec = ata.lu().solve(&atb).expect("failed to solve the normal equations");
    let mut poly = Polynomial::new(order);
    poly.coeffs = coeffs_vec.iter().cloned().collect();
    poly
}