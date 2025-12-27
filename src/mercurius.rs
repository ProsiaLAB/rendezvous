use crate::integrator::{ForceSplit, Reset, Synchronize};
use crate::integrator::{StepContext, SyncContext};
use crate::particle::{Particle, Particles};

use crate::ias15::Ias15;
use crate::whfast::WHFast;

pub struct Mercurius {
    r_crit_hill: f64,
    recalculate_coordinates_this_time_step: bool,
    recalculate_r_crit_this_time_step: bool,
    safe_mode: bool,

    is_synchronized: bool,
    pub(crate) mode: MercuriusMode,
    pub(crate) n_encounter_active: usize,
    tp_only_encounter: bool,
    pub(crate) dcrit: Vec<f64>,
    particles_backup: Particles,
    particles_backup_additional_forces: Particles,
    pub(crate) encounter_map: Vec<usize>,
    switch_fn: SwitchFunction,

    whfast: WHFast,
    ias15: Ias15,
}

impl Mercurius {
    /// Switching function between close-encounter (IAS15) and long-range (WHFast) integrators.
    pub fn switch(&self, d: f64, dcrit_max: f64) -> f64 {
        match self.switch_fn {
            SwitchFunction::Mercury => Self::switch_mercury(d, dcrit_max),
            SwitchFunction::Hernandez19C4 => Self::switch_hernandez19_c4(d, dcrit_max),
            SwitchFunction::Hernandez19C5 => Self::switch_hernandez19_c5(d, dcrit_max),
            SwitchFunction::Infinity => Self::switch_infinity(d, dcrit_max),
        }
    }

    fn switch_mercury(d: f64, dcrit: f64) -> f64 {
        let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
        if y < 0.0 {
            0.0
        } else if y > 1.0 {
            1.0
        } else {
            10.0 * y * y * y - 15.0 * y * y * y * y + 6.0 * y * y * y * y * y
        }
    }

    fn switch_hernandez19_c4(d: f64, dcrit: f64) -> f64 {
        let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
        if y < 0.0 {
            0.0
        } else if y > 1.0 {
            1.0
        } else {
            (70.0 * y * y * y * y - 315.0 * y * y * y + 540.0 * y * y - 420.0 * y + 126.0)
                * y
                * y
                * y
                * y
                * y
        }
    }

    fn switch_hernandez19_c5(d: f64, dcrit: f64) -> f64 {
        let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
        if y < 0.0 {
            0.0
        } else if y > 1.0 {
            1.0
        } else {
            (-252.0 * y * y * y * y * y + 1386.0 * y * y * y * y - 3080.0 * y * y * y
                + 3465.0 * y * y
                - 1980.0 * y
                + 462.0)
                * y
                * y
                * y
                * y
                * y
                * y
        }
    }

    fn switch_infinity(d: f64, dcrit: f64) -> f64 {
        let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
        if y < 0.0 {
            0.0
        } else if y > 1.0 {
            1.0
        } else {
            f(y) / (f(y) + f(1.0 - y))
        }
    }

    pub fn add(&mut self, g: f64, dt: f64, particles: &Particles) {
        match self.mode {
            MercuriusMode::LongRange => {
                // WHFast mode
                self.recalculate_r_crit_this_time_step = true;
                self.recalculate_coordinates_this_time_step = true;
            }
            MercuriusMode::CloseEncounter => {
                // IAS15 mode
                self.ias15.reset();
                if self.dcrit.len() < particles.len() {
                    self.dcrit.resize(particles.len(), 0.0);
                }
                let new_index = particles.len() - 1;
                let p0 = &particles[0];
                let pi = &particles[new_index];
                self.set_dcrit(p0, pi, g, dt, new_index);
                if self.particles_backup.len() < particles.len() {
                    self.particles_backup.resize_as(&particles);
                    self.particles_backup_additional_forces
                        .resize_as(&particles);
                }
                self.encounter_map.push(new_index);
                if particles.are_all_active() {
                    self.n_encounter_active += 1;
                }
            }
        }
    }

    pub fn remove(&mut self, index: usize) -> bool {
        if !self.dcrit.is_empty() && index < self.dcrit.len() {
            self.dcrit.remove(index);
        }
        self.ias15.reset();
        if matches!(self.mode, MercuriusMode::CloseEncounter)
            && let Some(pos) = self.encounter_map.iter().position(|&x| x == index)
        {
            self.encounter_map.remove(pos);
            if pos < self.n_encounter_active {
                self.n_encounter_active -= 1;
            }
        }

        true
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

    pub fn recalculate_coordinates(&mut self) {
        self.recalculate_coordinates_this_time_step = true;
    }
}

impl Synchronize for Mercurius {
    fn synchronize(&mut self, _ctx: &mut SyncContext<'_>) {
        todo!()
    }
}

impl ForceSplit for Mercurius {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}

pub enum MercuriusMode {
    /// WHFast
    LongRange,
    /// IAS15
    CloseEncounter,
}

pub enum SwitchFunction {
    Mercury,
    Hernandez19C4,
    Hernandez19C5,
    Infinity,
}

fn f(x: f64) -> f64 {
    if x < 0.0 { 0.0 } else { (-1.0 / x).exp() }
}
