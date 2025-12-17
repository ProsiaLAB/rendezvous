use crate::integrator::{ForceSplitIntegrator, StepContext, SyncContext, Synchronizable};

pub struct Eos;

impl Synchronizable for Eos {
    fn synchronize(&mut self, _ctx: SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplitIntegrator for Eos {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
