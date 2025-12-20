use crate::integrator::ForceSplit;
use crate::integrator::Reset;
use crate::integrator::StepContext;

pub struct Ias15 {
    pub epsilon: f64,
    pub min_dt: f64,
    pub adaptive_mode: AdaptiveMode,
    pub iterations_max_exceeded: usize,
    pub n_allocated: usize,
}

impl Reset for Ias15 {
    fn reset(&mut self) {
        // Reset internal state
        todo!()
    }
}

impl ForceSplit for Ias15 {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}

pub enum AdaptiveMode {
    Individual,
    Global,
    Prs23,
    Aarseth85,
}
