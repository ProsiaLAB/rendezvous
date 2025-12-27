use crate::gbs::Gbs;
use crate::ias15::Ias15;
use crate::integrator::StepContext;
use crate::integrator::{ForceSplit, Reset};
use crate::particle::Particles;
use crate::whfast::WHFast;

pub struct Trace {
    pub pericentric_mode: PericentricMode,
    pub r_crit_hill: f64,
    pub peri_crit_eta: f64,

    pub mode: TraceMode,
    pub n_encounter_active: usize,
    pub tp_only_encounter: bool,
    pub particles_backup: Particles,
    pub particles_backup_kepler: Particles,
    pub particles_backup_additional_forces: Particles,
    pub encounter_map: Vec<usize>,
    pub current_ks: Vec<usize>,

    pub gbs: Gbs,
    pub ias15: Ias15,
    pub whfast: WHFast,
}

impl Trace {
    pub fn switch(&mut self) {
        todo!()
    }

    pub fn switch_pericentric(&mut self) {
        todo!()
    }

    pub fn resize_current_ks(&mut self, old_n: usize, new_n: usize) {
        let mut new_ks = vec![0; new_n * new_n];

        // Copy old matrix into expanded matrix
        for i in 0..old_n {
            for j in 0..old_n {
                new_ks[i * new_n + j] = self.current_ks[i * old_n + j];
            }
        }

        self.current_ks = new_ks;
    }

    pub fn add(&mut self, particles: &Particles) {
        if matches!(self.mode, TraceMode::Kepler | TraceMode::Full) {
            // GBS part
            let new_n = particles.len();
            let old_n = new_n - 1;

            self.particles_backup.resize_as(&particles);
            self.particles_backup_kepler.resize_as(&particles);
            self.particles_backup_additional_forces
                .resize_as(&particles);
            self.resize_current_ks(old_n, new_n);
            for &p in self.encounter_map.iter().skip(1) {
                self.current_ks[p * new_n + old_n] = 1;
            }

            self.encounter_map.push(old_n);
            if particles.are_all_active() {
                self.n_encounter_active += 1;
            }
        }
    }

    pub fn remove(&mut self, n: usize, index: usize) -> bool {
        self.gbs.reset();
        if matches!(self.mode, TraceMode::Kepler | TraceMode::Full) {
            if let Some(pos) = self.encounter_map.iter().position(|&x| x == index) {
                self.encounter_map.remove(pos);
                if pos < self.n_encounter_active {
                    self.n_encounter_active -= 1;
                }
            }

            let new_n = n - 1;
            let old_n = n;

            for i in 0..new_n {
                let src_i = if i < index { i } else { i + 1 };
                for j in 0..new_n {
                    let src_j = if j < index { j } else { j + 1 };
                    self.current_ks[i * new_n + j] = self.current_ks[src_i * old_n + src_j];
                }
            }

            self.current_ks.truncate(new_n * new_n);
        }

        true
    }
}

impl ForceSplit for Trace {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}

pub enum PericentricMode {
    PartialGbs,
    FullGbs,
    FullIas15,
}

pub enum TraceMode {
    Interaction,
    Kepler,
    None,
    Full,
}
