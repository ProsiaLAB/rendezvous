use crate::integrator::{ForceSplit, Synchronize};
use crate::integrator::{StepContext, SyncContext};

pub struct Eos;

impl Synchronize for Eos {
    fn synchronize(&mut self, _ctx: &mut SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplit for Eos {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
