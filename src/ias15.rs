use crate::integrator::{ForceSplitIntegrator, StepContext};

pub struct Ias15 {
    pub epsilon: f64,
    pub min_dt: f64,
    pub adaptive_mode: AdaptiveMode,
    pub iterations_max_exceeded: usize,
    pub n_allocated: usize,
}

impl Ias15 {
    pub fn reset(&mut self) {
        // Reset internal state
        todo!()
    }
}

impl ForceSplitIntegrator for Ias15 {
    fn pre_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }
}

pub enum AdaptiveMode {
    Individual,
    Global,
    Prs23,
    Aarseth85,
}
