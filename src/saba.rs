use crate::integrator::{ForceSplit, Synchronize};
use crate::integrator::{StepContext, SyncContext};

pub struct Saba;

impl Synchronize for Saba {
    fn synchronize(&mut self, _ctx: &mut SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplit for Saba {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
