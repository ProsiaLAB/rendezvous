use crate::integrator::{ForceSplit, Synchronize};
use crate::integrator::{StepContext, SyncContext};

pub struct Janus;

impl Synchronize for Janus {
    fn synchronize(&mut self, _ctx: &mut SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplit for Janus {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
