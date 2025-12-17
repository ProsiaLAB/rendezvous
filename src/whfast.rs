use crate::integrator::{ForceSplitIntegrator, StepContext, SyncContext, Synchronizable};

pub struct WHFast {
    pub recalculate_coordinates_this_time_step: bool,
}

impl Synchronizable for WHFast {
    fn synchronize(&mut self, _ctx: SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplitIntegrator for WHFast {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
