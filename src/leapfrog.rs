//! This is the standard leap frog integrator. It is second order
//! and symplectic. No configuration is available
//! (the timestep is set in the simulation structure).

use crate::integrator::{ForceSplitIntegrator, StepContext};

pub struct LeapFrog;

impl ForceSplitIntegrator for LeapFrog {
    fn pre_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }
}
