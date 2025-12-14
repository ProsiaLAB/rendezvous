use crate::{ias15::Ias15, particle::Particle, rendezvous::Simulation, whfast::WHFast};

pub struct Mercurius {
    pub r_crit_hill: f64,
    pub recalculate_coordinates_this_time_step: bool,
    pub recalculate_r_crit_this_time_step: bool,
    pub safe_mode: bool,

    pub(crate) is_synchronized: bool,
    pub(crate) mode: MercuriusMode,
    pub(crate) n_encounter: usize,
    pub(crate) n_encounter_active: usize,
    pub(crate) tp_only_encounter: bool,
    pub(crate) dcrit: Vec<f64>,
    pub(crate) particles_backup: Vec<Particle>,
    pub(crate) particles_backup_additional_forces: Vec<Particle>,
    pub(crate) encounter_map: Vec<usize>,

    pub(crate) whfast: WHFast,
    pub(crate) ias15: Ias15,
}

impl Mercurius {
    /// Switching function between close-encounter (IAS15) and long-range (WHFast) integrators.
    pub fn switch(&self) {
        todo!()
    }

    pub fn set_dcrit(&mut self, p0: &Particle, pi: &Particle, g: f64, dt: f64, index: usize) {
        let m0 = p0.m;
        let dx = pi.x;
        let dy = pi.y;
        let dz = pi.z;
        let dvx = pi.vx - p0.vx;
        let dvy = pi.vy - p0.vy;
        let dvz = pi.vz - p0.vz;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let v2 = dvx * dvx + dvy * dvy + dvz * dvz;
        let gm = g * (m0 + pi.m);
        let a = gm * r / (2.0 * gm - r * v2);
        let vc = (gm / a.abs()).sqrt();
        let mut dcrit: f64 = 0.0;
        // Criteria 1: average velocity
        dcrit = dcrit.max(vc * 0.4 * dt);
        // Criteria 2: current velocity
        dcrit = dcrit.max(v2.sqrt() * 0.4 * dt);
        // Criteria 3: Hill radius
        let hill_radius = self.r_crit_hill * a * (pi.m / (3.0 * p0.m)).cbrt();
        dcrit = dcrit.max(hill_radius);
        // Criteria 4: physical length
        dcrit = dcrit.max(2.0 * pi.r);
        self.dcrit[index] = dcrit;
    }
}

pub(crate) enum MercuriusMode {
    LongRange,
    CloseEncounter,
}
