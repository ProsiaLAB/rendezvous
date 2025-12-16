//! The Gragg-Bulirsch-Stoer integrator (short BS for Bulirsch-Stoer)
//! is an adaptive integrator which uses Richardson extrapolation and
//! the modified midpoint method to obtain solutions to ordinary
//! differential equations.

use crate::integrator::{ForceSplitIntegrator, StepContext};

pub struct Gbs {
    /// Allowed absolute scalar error
    pub eps_abs: f64,
    /// Allowed relative scalar error
    pub eps_rel: f64,
    pub min_dt: f64,
    pub max_dt: f64,

    pub sequence: Vec<usize>,
    pub cost_per_step: Vec<usize>,
    pub cost_per_time_unit: Vec<f64>,
    pub optimal_step: Vec<f64>,
    pub coeff: Vec<f64>,
    pub dt_proposed: f64,
    pub first_or_last_step: bool,
    pub previous_rejected: bool,
    pub target_iter: usize,
    pub user_ode_needs_nbody: bool,
}

impl ForceSplitIntegrator for Gbs {
    fn pre_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }
}
