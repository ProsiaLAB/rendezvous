use crate::integrator::ForceSplitIntegrator;
use crate::integrator::StepContext;
use crate::particle::Particle;

pub struct Trace {
    pub pericentric_mode: PericentricMode,
    pub r_crit_hill: f64,
    pub peri_crit_eta: f64,

    pub mode: TraceMode,
    pub n_encounter: usize,
    pub n_encounter_active: usize,
    pub tp_only_encounter: bool,
    pub particles_backup: Vec<Particle>,
    pub particles_backup_kepler: Vec<Particle>,
    pub particles_backup_additional_forces: Vec<Particle>,
    pub encounter_map: Vec<usize>,
    pub current_ks: Vec<usize>,
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

    pub fn update_encounters(&mut self, n_active: usize, particles: &[Particle]) {
        if matches!(self.mode, TraceMode::Kepler | TraceMode::Full) {
            // GBS part
            let new_n = particles.len();
            let old_n = new_n - 1;

            self.particles_backup.resize(new_n, Particle::default());
            self.particles_backup_kepler
                .resize(new_n, Particle::default());
            self.particles_backup_additional_forces
                .resize(new_n, Particle::default());
            self.resize_current_ks(old_n, new_n);
            for &p in self.encounter_map.iter().skip(1) {
                self.current_ks[p * new_n + old_n] = 1;
            }

            self.encounter_map.push(old_n);
            self.n_encounter += 1;
            if n_active == usize::MAX {
                self.n_encounter_active += 1;
            }
        }
    }
}

impl ForceSplitIntegrator for Trace {
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
