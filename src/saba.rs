use crate::integrator::{ForceSplitIntegrator, Synchronizable};
use crate::integrator::{StepContext, SyncContext};

pub struct Saba;

impl Synchronizable for Saba {
    fn synchronize(&mut self, _ctx: SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplitIntegrator for Saba {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
