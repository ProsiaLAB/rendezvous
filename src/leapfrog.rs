//! This is the standard leap frog integrator. It is second order
//! and symplectic. No configuration is available
//! (the timestep is set in the simulation structure).

use rayon::iter::ParallelIterator;

use crate::gravity::IgnoreGravityTerms;
use crate::integrator::ForceSplit;
use crate::integrator::StepContext;

pub struct LeapFrog;

impl ForceSplit for LeapFrog {
    fn pre_force(&mut self, ctx: &mut StepContext<'_>) {
        *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreAll;
        ctx.particles.par_iter_mut().for_each(|p| {
            p.x += 0.5 * ctx.dt * p.vx;
            p.y += 0.5 * ctx.dt * p.vy;
            p.z += 0.5 * ctx.dt * p.vz;
        });
        *ctx.t += ctx.dt / 2.0;
    }

    fn post_force(&mut self, ctx: &mut StepContext<'_>) {
        ctx.particles.par_iter_mut().for_each(|p| {
            p.vx += ctx.dt * p.ax;
            p.vy += ctx.dt * p.ay;
            p.vz += ctx.dt * p.az;

            p.x += 0.5 * ctx.dt * p.vx;
            p.y += 0.5 * ctx.dt * p.vy;
            p.z += 0.5 * ctx.dt * p.vz;
        });
        *ctx.t += ctx.dt / 2.0;
        *ctx.dt_last_done = ctx.dt;
    }
}
