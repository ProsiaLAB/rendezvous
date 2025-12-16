use crate::integrator::{ForceSplitIntegrator, StepContext, SyncContext, Synchronizable};

pub struct Saba;

impl Synchronizable for Saba {
    fn synchronize(&mut self, _ctx: SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplitIntegrator for Saba {
    fn pre_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: StepContext<'_>) {
        todo!()
    }
}
