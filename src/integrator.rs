use crate::gravity::IgnoreGravityTerms;
use crate::particle::Particle;

use crate::eos::Eos;
use crate::gbs::Gbs;
use crate::ias15::Ias15;
use crate::janus::Janus;
use crate::leapfrog::LeapFrog;
use crate::mercurius::Mercurius;
use crate::saba::Saba;
use crate::sei::Sei;
use crate::trace::Trace;
use crate::whfast::WHFast;

pub enum Integrator {
    Ias15(Ias15),
    WHFast(WHFast),
    LeapFrog(LeapFrog),
    Gbs(Gbs),
    Sei(Sei),
    Janus(Janus),
    Mercurius(Mercurius),
    Saba(Saba),
    Eos(Eos),
    Trace(Trace),
    None,
}

impl ForceSplitIntegrator for Integrator {
    fn pre_force(&mut self, ctx: &mut StepContext<'_>) {
        match self {
            Integrator::WHFast(i) => i.pre_force(ctx),
            Integrator::Saba(i) => i.pre_force(ctx),
            Integrator::Mercurius(i) => i.pre_force(ctx),
            Integrator::Janus(i) => i.pre_force(ctx),
            Integrator::Eos(i) => i.pre_force(ctx),
            Integrator::Ias15(i) => i.pre_force(ctx),
            Integrator::LeapFrog(i) => i.pre_force(ctx),
            Integrator::Sei(i) => i.pre_force(ctx),
            Integrator::Gbs(i) => i.pre_force(ctx),
            Integrator::Trace(i) => i.pre_force(ctx),
            Integrator::None => {}
        }
    }

    fn post_force(&mut self, ctx: &mut StepContext<'_>) {
        match self {
            Integrator::WHFast(i) => i.post_force(ctx),
            Integrator::Saba(i) => i.post_force(ctx),
            Integrator::Mercurius(i) => i.post_force(ctx),
            Integrator::Janus(i) => i.post_force(ctx),
            Integrator::Eos(i) => i.post_force(ctx),
            Integrator::Ias15(i) => i.post_force(ctx),
            Integrator::LeapFrog(i) => i.post_force(ctx),
            Integrator::Sei(i) => i.post_force(ctx),
            Integrator::Gbs(i) => i.post_force(ctx),
            Integrator::Trace(i) => i.post_force(ctx),
            Integrator::None => {}
        }
    }
}

pub struct SyncContext<'a> {
    pub particles: &'a mut [Particle],
}

pub struct StepContext<'a> {
    pub particles: &'a mut [Particle],
    pub t: &'a mut f64,
    pub dt: f64,
    pub dt_last_done: &'a mut f64,
    pub ignore_gravity_terms: &'a mut IgnoreGravityTerms,
}

pub trait Synchronizable {
    fn synchronize(&mut self, ctx: SyncContext<'_>);
}

pub trait ForceSplitIntegrator {
    fn pre_force(&mut self, ctx: &mut StepContext<'_>);
    fn post_force(&mut self, ctx: &mut StepContext<'_>);
}

pub trait Reset {
    fn reset(&mut self);
}
