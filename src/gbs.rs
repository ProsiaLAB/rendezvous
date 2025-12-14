//! The Gragg-Bulirsch-Stoer integrator (short BS for Bulirsch-Stoer)
//! is an adaptive integrator which uses Richardson extrapolation and
//! the modified midpoint method to obtain solutions to ordinary
//! differential equations.

use crate::{collision::CollisionResolver, integrator::Integrator, rendezvous::Simulation};

pub struct Gbs {
    /// Allowed absolute scalar error
    pub eps_abs: f64,
    /// Allowed relative scalar error
    pub eps_rel: f64,
    pub min_dt: f64,
    pub max_dt: f64,

    pub(crate) sequence: Vec<usize>,
    pub(crate) cost_per_step: Vec<usize>,
    pub(crate) cost_per_time_unit: Vec<f64>,
    pub(crate) optimal_step: Vec<f64>,
    pub(crate) coeff: Vec<f64>,
    pub(crate) dt_proposed: f64,
    pub(crate) first_or_last_step: bool,
    pub(crate) previous_rejected: bool,
    pub(crate) target_iter: usize,
    pub(crate) user_ode_needs_nbody: bool,
}
