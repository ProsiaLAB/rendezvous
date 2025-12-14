use crate::particle::Particle;

pub struct Trace {
    pub pericentric_mode: PericentricMode,
    pub r_crit_hill: f64,
    pub peri_crit_eta: f64,

    pub(crate) mode: TraceMode,
    pub(crate) n_encounter: usize,
    pub(crate) n_encounter_active: usize,
    pub(crate) tp_only_encounter: bool,
    pub(crate) particles_backup: Vec<Particle>,
    pub(crate) particles_backup_kepler: Vec<Particle>,
    pub(crate) particles_backup_additional_forces: Vec<Particle>,
    pub(crate) encounter_map: Vec<usize>,
    pub(crate) current_ks: Vec<usize>,
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
