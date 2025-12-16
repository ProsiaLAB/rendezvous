use crate::integrator::{ForceSplitIntegrator, StepContext};

pub struct Sei;

impl ForceSplitIntegrator for Sei {
    fn pre_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }
}
