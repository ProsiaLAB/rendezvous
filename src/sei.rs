use crate::integrator::ForceSplit;
use crate::integrator::StepContext;

pub struct Sei {
    pub omega: f64,
    pub omega_z: f64,
}

impl ForceSplit for Sei {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}
