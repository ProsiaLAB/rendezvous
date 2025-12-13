use prosia_extensions::types::{RMatrix, RVector};

pub trait OdeSystem {
    /// Number of components / dimensions in the ODE system
    fn dimension(&self) -> usize;

    fn derivatives(&mut self);
    fn scale(&self);
    fn pre_time_step(&mut self, dt: f64);
    fn post_time_step(&mut self, dt: f64);
    fn needs_nbody(&self) -> bool;
}

pub struct OdeState {
    pub y: RVector,
}

pub struct OdeInternalState {
    scale: RVector,
    c: RVector,
    d: RMatrix,
    y1: RVector,
    y0_dot: RVector,
    y_dot: RVector,
    y_temp: RVector,
}
