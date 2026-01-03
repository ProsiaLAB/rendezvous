#![allow(clippy::excessive_precision)]

use std::f64::consts::PI;
use std::mem;

use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator,
};
use thiserror::Error;

use crate::gravity::{Gravity, IgnoreGravityTerms};
use crate::integrator::{ForceSplit, Synchronize};
use crate::integrator::{StepContext, SyncContext};
use crate::particle::{self, Particle, Particles, TestParticleKind, Transformations};

const CORRECTOR_A_1: f64 = 0.41833001326703777398908601289259374469640768464934;
const CORRECTOR_A_2: f64 = 0.83666002653407554797817202578518748939281536929867;
const CORRECTOR_A_3: f64 = 1.2549900398011133219672580386777812340892230539480;
const CORRECTOR_A_4: f64 = 1.6733200530681510959563440515703749787856307385973;
const CORRECTOR_A_5: f64 = 2.0916500663351888699454300644629687234820384232467;
const CORRECTOR_A_6: f64 = 2.5099800796022266439345160773555624681784461078960;
const CORRECTOR_A_7: f64 = 2.9283100928692644179236020902481562128748537925454;
const CORRECTOR_A_8: f64 = 3.3466401061363021919126881031407499575712614771947;
const CORRECTOR_B_31: f64 = -0.024900596027799867499350357910273437184309981229127;
const CORRECTOR_B_51: f64 = -0.0083001986759332891664501193034244790614366604097090;
const CORRECTOR_B_52: f64 = 0.041500993379666445832250596517122395307183302048545;
const CORRECTOR_B_71: f64 = 0.0024926811426922105779030593952776964450539008582219;
const CORRECTOR_B_72: f64 = -0.018270923246702131478062356884535264841652263842597;
const CORRECTOR_B_73: f64 = 0.053964399093127498721765893493510877532452806339655;
const CORRECTOR_B_111: f64 = 0.00020361579647854651301632818774633716473696537436847;
const CORRECTOR_B_112: f64 = -0.0023487215292295354188307328851055489876255097419754;
const CORRECTOR_B_113: f64 = 0.012309078592019946317544564763237909911330686448336;
const CORRECTOR_B_114: f64 = -0.038121613681288650508647613260247372125243616270670;
const CORRECTOR_B_115: f64 = 0.072593394748842738674253180742744961827622366521517;
const CORRECTOR_B_178: f64 = 0.093056103771425958591541059067553547100903397724386;
const CORRECTOR_B_177: f64 = -0.065192863576377893658290760803725762027864651086787;
const CORRECTOR_B_176: f64 = 0.032422198864713580293681523029577130832258806467604;
const CORRECTOR_B_175: f64 = -0.012071760822342291062449751726959664253913904872527;
const CORRECTOR_B_174: f64 = 0.0033132577069380655655490196833451994080066801611459;
const CORRECTOR_B_173: f64 = -0.00063599983075817658983166881625078545864140848560259;
const CORRECTOR_B_172: f64 = 0.000076436355227935738363241846979413475106795392377415;
const CORRECTOR_B_171: f64 = -0.0000043347415473373580190650223498124944896789841432241;
const CORRECTOR2_B: f64 = 0.03486083443891981449909050107438281205803;

const INVERSE_FACTORIALS: [f64; 35] = [
    1.,
    1.,
    1. / 2.,
    1. / 6.,
    1. / 24.,
    1. / 120.,
    1. / 720.,
    1. / 5040.,
    1. / 40320.,
    1. / 362880.,
    1. / 3628800.,
    1. / 39916800.,
    1. / 479001600.,
    1. / 6227020800.,
    1. / 87178291200.,
    1. / 1307674368000.,
    1. / 20922789888000.,
    1. / 355687428096000.,
    1. / 6402373705728000.,
    1. / 121645100408832000.,
    1. / 2432902008176640000.,
    1. / 51090942171709440000.,
    1. / 1124000727777607680000.,
    1. / 25852016738884976640000.,
    1. / 620448401733239439360000.,
    1. / 15511210043330985984000000.,
    1. / 403291461126605635584000000.,
    1. / 10888869450418352160768000000.,
    1. / 304888344611713860501504000000.,
    1. / 8841761993739701954543616000000.,
    1. / 265252859812191058636308480000000.,
    1. / 8222838654177922817725562880000000.,
    1. / 263130836933693530167218012160000000.,
    1. / 8683317618811886495518194401280000000.,
    1. / 295232799039604140847618609643520000000.,
];

const NMAX_QUARTIC: usize = 64;
const NMAX_NEWTON: usize = 32;

pub struct WHFast {
    recalculate_coordinates_this_time_step: bool,
    coordinates: Coordinates,
    kernel: Kernel,
    corrector: Corrector,
    use_secondary_corrector: bool,
    keep_unsynchoronized: bool,
    is_synchronized: bool,
    safe_mode: SafeMode,
    particles: Particles,
    time_step_warning: usize,
}

impl WHFast {
    fn init(&mut self, ctx: &mut SyncContext<'_>) -> Result<(), WHFastError> {
        if let Some(var_cfg) = ctx.var_cfg {
            for v in var_cfg.iter() {
                if v.order != 1 {
                    return Err(WHFastError::UnsupportedVariationalOrder);
                }
                if v.is_test_particle {
                    return Err(WHFastError::TestParticleVariationalNotSupported);
                }
            }
        }

        if ctx.var_cfg.is_some() && !matches!(self.coordinates, Coordinates::Jacobi) {
            return Err(WHFastError::VariationalJacobiOnly);
        }

        if !matches!(self.kernel, Kernel::Default)
            && !matches!(self.coordinates, Coordinates::Jacobi)
        {
            return Err(WHFastError::NonStandardKernelJacobiOnly);
        }

        if ctx.var_cfg.is_some() && !matches!(self.kernel, Kernel::Default) {
            return Err(WHFastError::VariationalStandardKernelOnly);
        }

        if !matches!(self.corrector, Corrector::None)
            && !matches!(
                self.coordinates,
                Coordinates::Jacobi | Coordinates::Barycentric
            )
        {
            return Err(WHFastError::SymplecticCorrectorJacobiOrBarycentricOnly);
        }

        if self.keep_unsynchoronized && matches!(self.safe_mode, SafeMode::Combine) {
            return Err(WHFastError::InvalidSafeModeCombination);
        }

        if matches!(self.kernel, Kernel::ModifiedKick | Kernel::Lazy) {
            *ctx.gravity = Gravity::Jacobi;
        } else {
            match self.coordinates {
                Coordinates::Jacobi => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreWHFastwithJacobi;
                }
                Coordinates::Barycentric => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreAll
                }
                _ => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreWHFastwithDHC;
                }
            }
        }

        self.recalculate_coordinates();

        Ok(())
    }

    pub fn recalculate_coordinates(&mut self) {
        self.recalculate_coordinates_this_time_step = true;
    }

    fn kepler_step(&mut self, ctx: &SyncContext<'_>, dt: f64) {
        let m0 = ctx.particles[0].m;
        let n_active = if ctx.particles.are_all_active()
            || matches!(ctx.test_particle_kind, TestParticleKind::Massive)
        {
            ctx.particles.n_real()
        } else {
            ctx.particles.active.len()
        };

        let mut eta = m0;

        match self.coordinates {
            Coordinates::Jacobi => {
                self.particles
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(i, p)| {});
            }

            // (1..ctx.particles.n_real()).into_par_iter().for_each(|i| {
            //     let mut eta = m0;
            //     let particles = self.particles.as_mut().unwrap();
            //     for j in 1..(n_active.min(i)) {
            //         let mass = particles[j].m;
            //         eta += mass;
            //     }
            //     let mut kpl = KeplerSystem::new(&mut particles[i], eta * ctx.g);
            //     kpl.solve(dt);
            // }),
            Coordinates::DemocraticHeliocentric => {}
            Coordinates::Whds => {}
            Coordinates::Barycentric => {}
        }
    }

    fn com_step(&mut self, dt: f64) {
        self.particles[0].x += self.particles[0].vx * dt;
        self.particles[0].y += self.particles[0].vy * dt;
        self.particles[0].z += self.particles[0].vz * dt;
    }

    fn apply_corrector(&mut self) {}

    fn apply_secondary_corrector(&mut self) {}
}

impl Synchronize for WHFast {
    fn synchronize(&mut self, ctx: &mut SyncContext<'_>) {
        if self.init(ctx).is_err() {
            return;
        }

        if !self.is_synchronized {
            let n_real = ctx.particles.n_real();
            let n_active = if ctx.particles.are_all_active()
                || matches!(ctx.test_particle_kind, TestParticleKind::Massive)
            {
                n_real
            } else {
                ctx.particles.active.len()
            };

            let sync_particles = if self.keep_unsynchoronized {
                Some(self.particles.clone())
            } else {
                None
            };

            match self.kernel {
                Kernel::Default | Kernel::ModifiedKick => {}
                Kernel::Lazy => {
                    self.kepler_step(ctx, ctx.dt / 2.0);
                    self.com_step(ctx.dt / 2.0);
                }
                Kernel::Composition => {
                    self.kepler_step(ctx, 3.0 * ctx.dt / 8.0);
                    self.com_step(3.0 * ctx.dt / 8.0);
                }
            }

            if self.use_secondary_corrector {
                self.apply_secondary_corrector();
            }

            if !matches!(self.corrector, Corrector::None) {
                self.apply_corrector();
            }

            match self.coordinates {
                Coordinates::Jacobi => {
                    let masses = ctx
                        .particles
                        .active
                        .iter()
                        .map(|p| p.m)
                        .collect::<Vec<f64>>();
                    ctx.particles.active.transform_jacobi_to_inertial_posvel(
                        &self.particles.active,
                        &masses,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::DemocraticHeliocentric => {
                    ctx.particles.active.transform_dhc_to_inertial_posvel(
                        &self.particles.active,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::Whds => {
                    ctx.particles.active.transform_whds_to_inertial_posvel(
                        &self.particles.active,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::Barycentric => {
                    ctx.particles
                        .active
                        .transform_barycentric_to_inertial_posvel(
                            &self.particles.active,
                            n_real,
                            n_active,
                        );
                }
            }

            if let Some(vcs) = ctx.var_cfg {
                todo!()
            }

            if let Some(sp) = sync_particles {
                self.particles = sp;
            }
        }
    }
}

impl ForceSplit for WHFast {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}

enum Coordinates {
    /// Jacobi coordinates (default)
    Jacobi,
    /// Democratic Heliocentric coordinates
    DemocraticHeliocentric,
    /// WHDS coordinates (Hernandez and Dehnen 2017)
    Whds,
    /// Barycentric coordinates
    Barycentric,
}

enum Kernel {
    Default,
    ModifiedKick,
    Composition,
    Lazy,
}

enum Corrector {
    None,
    Order3,
    Order5,
    Order7,
    Order11,
    Order17,
}

impl Corrector {
    fn order(self) -> Option<u32> {
        match self {
            Corrector::None => None,
            Corrector::Order3 => Some(3),
            Corrector::Order5 => Some(5),
            Corrector::Order7 => Some(7),
            Corrector::Order11 => Some(11),
            Corrector::Order17 => Some(17),
        }
    }
}

enum SafeMode {
    DriftKickDrift,
    Combine,
}

#[derive(Debug, Error)]
pub enum WHFastError {
    #[error("WHFast/MEGNO only supports first-order variational equations")]
    UnsupportedVariationalOrder,

    #[error("Test particle variational particles are not supported in WHFast")]
    TestParticleVariationalNotSupported,

    #[error("Variational particles require Jacobi coordinates in WHFast")]
    VariationalJacobiOnly,

    #[error("Non-standard kernels require Jacobi coordinates in WHFast")]
    NonStandardKernelJacobiOnly,

    #[error("Variational particles are only supported with the standard kernel in WHFast")]
    VariationalStandardKernelOnly,

    #[error("Symplectic correctors require Jacobi or Barycentric coordinates in WHFast")]
    SymplecticCorrectorJacobiOrBarycentricOnly,

    #[error("Cannot keep unsynchronized particles when using SafeMode::Combine in WHFast")]
    InvalidSafeModeCombination,
}

mod stiefel {
    use crate::whfast::stumpff;

    pub(crate) fn gs(gs: &mut [f64], beta: f64, x: f64) {}

    pub(crate) fn gs3(gs: &mut [f64], beta: f64, x: f64) {
        let x2 = x * x;
        stumpff::cs3(gs, beta * x2);
        gs[1] *= x;
        gs[2] *= x2;
        gs[3] *= x2 * x;
    }
}

mod stumpff {
    use crate::whfast::INVERSE_FACTORIALS;

    pub(crate) fn cs() {
        todo!()
    }

    pub(crate) fn cs3(cs: &mut [f64], mut z: f64) {
        let mut n = 0;

        while z.abs() > 0.1 {
            z /= 4.0;
            n += 1;
        }

        const NMAX: usize = 13;

        let mut c_odd = INVERSE_FACTORIALS[NMAX];
        let mut c_even = INVERSE_FACTORIALS[NMAX - 1];

        for np in (3..NMAX - 2).rev().step_by(2) {
            c_odd = INVERSE_FACTORIALS[np] - z * c_odd;
            c_even = INVERSE_FACTORIALS[np - 1] - z * c_even;
        }

        cs[3] = c_odd;
        cs[2] = c_even;
        cs[1] = INVERSE_FACTORIALS[1] - z * c_odd;
        cs[0] = INVERSE_FACTORIALS[0] - z * c_even;

        while n > 0 {
            cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
            cs[2] = cs[1] * cs[1] * 0.5;
            cs[1] *= cs[0];
            cs[0] = 2.0 * cs[0] * cs[0] - 1.0;
            n -= 1;
        }
    }
}

struct KeplerSystem<'a> {
    mass: f64,
    r0: f64,
    r0i: f64,
    v2: f64,
    beta: f64,
    eta0: f64,
    zeta0: f64,
    gs: [f64; 6],
    p: &'a mut Particle,
}

impl<'a> KeplerSystem<'a> {
    fn new(p: &'a mut Particle, mass: f64) -> Self {
        let r0 = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        let r0i = 1.0 / r0;
        let v2 = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
        let beta = 2.0 * mass * r0i - v2;
        let eta0 = p.x * p.vx + p.y * p.vy + p.z * p.vz;
        let zeta0 = mass - beta * r0;

        KeplerSystem {
            mass,
            r0,
            r0i,
            v2,
            beta,
            eta0,
            zeta0,
            gs: [0.0; 6],
            p,
        }
    }

    fn solve(&mut self, dt: f64) {
        struct BetaVals {
            x_per_period: f64,
            inv_period: f64,
        }

        let betvals = if self.beta > 0.0 {
            Some(BetaVals {
                x_per_period: 2.0 * PI * self.beta.sqrt(),
                inv_period: self.beta.sqrt() * self.beta / (2.0 * PI * self.mass),
            })
        } else {
            None
        };

        let mut nwarns = 0;

        let mut x = if let Some(bv) = &betvals {
            if dt.abs() * bv.inv_period > 1.0 && nwarns == 0 {
                nwarns += 1;
                eprintln!(
                    "WARNING: WHFast is having convergence issues because the time step is \
                     comparable to or larger than the orbital period (dt * inv_period = {}). \
                     Consider reducing the time step.",
                    dt.abs() * bv.inv_period
                );
            }

            let dtr0i = dt * self.r0i;
            dtr0i * (1.0 - dtr0i * self.eta0 * 0.5 * self.r0i)
        } else {
            0.0
        };

        let mut converged = false;
        let mut old_x = x;

        stiefel::gs3(&mut self.gs, self.beta, x);

        let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];

        let mut ri = 1.0 / (self.r0 + val);
        x = ri * (x * val - self.eta0 * self.gs[2] - self.zeta0 * self.gs[3] + dt);

        if let Some(bv) = &betvals {
            if (x - old_x).abs() > 0.01 * bv.x_per_period {
                // Quartic solver
                // Linear initial guess
                x = self.beta * dt / self.mass;
                let mut prev_x = [0.0; NMAX_QUARTIC + 1];
                for n_lag in 1..NMAX_QUARTIC {
                    stiefel::gs3(&mut self.gs, self.beta, x);
                    let f = self.r0 * x + self.eta0 * self.gs[2] + self.zeta0 * self.gs[3] - dt;
                    let fp = self.r0 + self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                    let fpp = self.eta0 * self.gs[0] + self.zeta0 * self.gs[1];
                    let denom = fp + (16.0 * fp * fp - 20.0 * f * fpp).abs().sqrt();
                    let x_new = (x * denom - 5.0 * f) / denom;

                    for &pxi in prev_x.iter().take(n_lag).skip(1) {
                        if (x_new - pxi).abs() < 1e-12 {
                            converged = true;
                            break;
                        }
                    }

                    if converged {
                        break;
                    }

                    prev_x[n_lag] = x_new;
                    x = x_new;
                }
                let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                ri = 1.0 / (self.r0 + val);
            } else {
                // Newton's method

                for _ in 1..NMAX_NEWTON {
                    let old_x2 = old_x;
                    old_x = x;
                    stiefel::gs3(&mut self.gs, self.beta, x);
                    let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                    let ri = 1.0 / (self.r0 + val);
                    x = ri * (x * val - self.eta0 * self.gs[2] - self.zeta0 * self.gs[3] + dt);
                    if x == old_x || x == old_x2 {
                        converged = true;
                        break;
                    }
                }
            }
        }

        if !converged {
            let (mut xmin, mut xmax) = if let Some(bv) = &betvals {
                let xmin = bv.x_per_period * (dt * bv.inv_period).floor();
                let xmax = xmin + bv.x_per_period;
                (xmin, xmax)
            } else {
                // Hyperbolic
                let h2 = self.r0 * self.r0 * self.v2 - self.eta0 * self.eta0;
                let q = h2
                    / self.mass
                    / (1.0 + (1.0 - h2 * self.beta / (self.mass * self.mass)).sqrt());
                let vq = (h2.sqrt() / q).copysign(dt);
                let mut xmin = dt / ((vq * dt).abs() + self.r0);
                let mut xmax = dt / q;
                if dt < 0.0 {
                    mem::swap(&mut xmin, &mut xmax);
                }
                (xmin, xmax)
            };
            x = 0.5 * (xmin + xmax);

            loop {
                stiefel::gs3(&mut self.gs, self.beta, x);
                let s = self.r0 * x + self.eta0 * self.gs[2] + self.zeta0 * self.gs[3] - dt;
                if s >= 0.0 {
                    xmax = x;
                } else {
                    xmin = x;
                }
                x = 0.5 * (xmin + xmax);

                if (xmax - xmin).abs() > ((xmax + xmin) * 1e-15).abs() {
                    break;
                }
            }
            let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
            ri = 1.0 / (self.r0 + val);
        }

        if ri.is_nan() {
            ri = 0.0;
            self.gs[1] = 0.0;
            self.gs[2] = 0.0;
            self.gs[3] = 0.0;
        }

        let f = -self.mass * self.gs[2] * self.r0i;
        let g = dt - self.mass * self.gs[3];
        let fd = -self.mass * self.gs[1] * self.r0i * ri;
        let gd = -self.mass * self.gs[2] * ri;
    }
}

struct KeplerResult {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    nwarns: usize,
    val: f64,
}

// fn kepler_solver_variations(
//     ctx: &SyncContext<'_>,
//     dp1: &mut Particle,
//     kpl: &mut KeplerSystem,
//     i: usize,
//     x: f64,
//     mass: f64,
//     dt: f64,
// ) {
//     // Variations
//     if let Some(vcs) = ctx.var_cfg {
//         for vc in vcs.iter() {
//             stiefel::gs(&mut kpl.gs, kpl.beta, x);
//             let dr0 = (dp1.x * x0 + dp1.y * y0 + dp1.z * z0) * r0i;
//             let dbeta =
//                 -2.0 * mass * dr0 * r0i * r0i - 2.0 * (dp1.vx * vx0 + dp1.vy * vy0 + dp1.vz * vz0);
//             let deta0 =
//                 dp1.x * vx0 + x0 * dp1.vx + dp1.y * vy0 + y0 * dp1.vy + dp1.z * vz0 + z0 * dp1.vz;
//             let dzeta0 = -beta * dr0 - r0 * dbeta;
//             let g3_beta = 0.5 * (3.0 * gs[5] - x * gs[4]);
//             let g2_beta = 0.5 * (2.0 * gs[4] - x * gs[3]);
//             let g1_beta = 0.5 * (gs[3] - x * gs[2]);
//             let tbeta = eta0 * g2_beta + zeta0 * g3_beta;
//             let dx = -ri * (x * dr0 + gs[2] * deta0 + gs[3] * dzeta0 + tbeta * dbeta);
//             let dg1 = gs[0] * dx + g1_beta * dbeta;
//             let dg2 = gs[1] * dx + g2_beta * dbeta;
//             let dg3 = gs[2] * dx + g3_beta * dbeta;
//             let dr = dr0 + gs[1] * deta0 + gs[2] * dzeta0 + eta0 * dg1 + zeta0 * dg2;
//             let df = mass * gs[2] * dr0 * r0i * r0i - mass * dg2 * r0i;
//             let dg = -mass * dg3;
//             let dfd = -mass * dg1 * r0i * ri + mass * gs[1] * (dr0 * r0i + dr * ri) * r0i * ri;
//             let dgd = -mass * dg2 * ri + mass * gs[2] * dr * ri * ri;

//             dp1.x += f * dp1.x + g * dp1.vx + df * x0 + dg * vx0;
//             dp1.y += f * dp1.y + g * dp1.vy + df * y0 + dg * vy0;
//             dp1.z += f * dp1.z + g * dp1.vz + df * z0 + dg * vz0;

//             dp1.vx += fd * dp1.x + gd * dp1.vx + dfd * x0 + dgd * vx0;
//             dp1.vy += fd * dp1.y + gd * dp1.vy + dfd * y0 + dgd * vy0;
//             dp1.vz += fd * dp1.z + gd * dp1.vz + dfd * z0 + dgd * vz0;
//         }
//     }
// }
