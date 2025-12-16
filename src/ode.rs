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
    pub scale: RVector,
    pub c: RVector,
    pub d: RMatrix,
    pub y1: RVector,
    pub y0_dot: RVector,
    pub y_dot: RVector,
    pub y_temp: RVector,
}
