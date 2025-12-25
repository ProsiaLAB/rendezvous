#![allow(clippy::excessive_precision)]

use thiserror::Error;

use crate::gravity::{Gravity, IgnoreGravityTerms};
use crate::integrator::{ForceSplit, Synchronize};
use crate::integrator::{StepContext, SyncContext};
use crate::particle::{Particles, TestParticleKind, Transformations};

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

pub struct WHFast {
    recalculate_coordinates_this_time_step: bool,
    coordinates: Coordinates,
    kernel: Kernel,
    corrector: Corrector,
    use_secondary_corrector: bool,
    keep_unsynchoronized: bool,
    safe_mode: SafeMode,
    particles: Option<Particles>,
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

    pub fn kepler_step(&mut self) {}

    pub fn com_step(&mut self) {}

    pub fn apply_corrector(&mut self) {}

    pub fn apply_secondary_corrector(&mut self) {}
}

impl Synchronize for WHFast {
    fn synchronize(&mut self, ctx: &mut SyncContext<'_>) {
        if self.init(ctx).is_err() {
            return;
        }

        let n_real = ctx.particles.n_real();
        let n_active = if ctx.particles.are_all_active()
            || matches!(ctx.test_particle_kind, TestParticleKind::Massive)
        {
            n_real
        } else {
            ctx.particles.n_active()
        };

        let mut p_jh = match self.particles.take() {
            Some(p) => p,
            None => return,
        };

        match self.kernel {
            Kernel::Default | Kernel::ModifiedKick => {}
            Kernel::Lazy => {
                self.kepler_step();
                self.com_step();
            }
            Kernel::Composition => {
                self.kepler_step();
                self.com_step();
            }
        }

        if self.use_secondary_corrector {
            self.apply_secondary_corrector();
        }

        if !matches!(self.corrector, Corrector::None) {
            self.apply_corrector();
        }

        let masses = ctx
            .particles
            .as_slice()
            .iter()
            .map(|p| p.m)
            .collect::<Vec<f64>>();

        match self.coordinates {
            Coordinates::Jacobi => {
                ctx.particles
                    .as_mut_slice()
                    .transform_jacobi_to_inertial_posvel(
                        p_jh.as_slice(),
                        &masses,
                        n_real,
                        n_active,
                    );
            }
            Coordinates::DemocraticHeliocentric => {
                ctx.particles
                    .as_mut_slice()
                    .transform_dhc_to_inertial_posvel();
            }
            Coordinates::Whds => {
                ctx.particles
                    .as_mut_slice()
                    .transform_whds_to_inertial_posvel();
            }
            Coordinates::Barycentric => {
                ctx.particles
                    .as_mut_slice()
                    .transform_barycentric_to_inertial_posvel();
            }
        }

        if let Some(vcs) = ctx.var_cfg {
            // TODO
        }

        if self.keep_unsynchoronized {
            self.particles = Some(p_jh);
        } else {
            self.particles = None;
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
