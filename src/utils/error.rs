#[derive(Debug)]
pub enum FlowError {
    InvalidSpecificHeatRatio,
    InvalidMachNumber,
    MissingValues,
    MathError,
}
