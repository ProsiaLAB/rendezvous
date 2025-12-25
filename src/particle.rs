use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::tree::NodeId;

#[derive(Default, Clone)]
pub struct Particle {
    pub x: f64,
    pub y: f64,
    pub z: f64,

    pub vx: f64,
    pub vy: f64,
    pub vz: f64,

    pub ax: f64,
    pub ay: f64,
    pub az: f64,

    pub m: f64,
    pub r: f64,

    pub last_collision_time: f64,
    pub c: Option<NodeId>,
    pub hash: u64,
    pub removed: bool,
}

pub struct Particles {
    data: Vec<Particle>,
    n_active: usize,
    n_test: usize,
    n_variational: usize,
}

impl Particles {
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            n_active: 0,
            n_test: 0,
            n_variational: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns an iterator over all real (non-variational) particles.
    pub fn real_iter(&self) -> impl Iterator<Item = &Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].iter()
    }

    /// Returns a mutable iterator over all real (non-variational) particles.
    pub fn real_iter_mut(&mut self) -> impl Iterator<Item = &mut Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].iter_mut()
    }

    /// Returns an iterator over all particles.
    pub fn iter(&self) -> impl Iterator<Item = &Particle> {
        self.data.iter()
    }

    /// Returns a mutable iterator over all particles.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Particle> {
        self.data.iter_mut()
    }

    /// Returns a parallel iterator over all real (non-variational) particles.
    pub fn real_par_iter(&self) -> impl ParallelIterator<Item = &Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].par_iter()
    }

    /// Returns a mutable parallel iterator over all real (non-variational) particles.
    pub fn real_par_iter_mut(&mut self) -> impl ParallelIterator<Item = &mut Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].par_iter_mut()
    }

    /// Returns a parallel iterator over all particles.
    pub fn par_iter(&self) -> impl ParallelIterator<Item = &Particle> {
        self.data.par_iter()
    }

    /// Returns a mutable parallel iterator over all particles.
    pub fn par_iter_mut(&mut self) -> impl ParallelIterator<Item = &mut Particle> {
        self.data.par_iter_mut()
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.n_active = 0;
        self.n_test = 0;
        self.n_variational = 0;
    }
}

pub struct ParticleData<T> {
    pub active: Vec<T>,
    pub test: Vec<T>,
    pub variational: Vec<T>,
}

/// Type of test particles in the simulation.
pub enum TestParticleKind {
    /// Active particles feel test particles
    Massive,
    /// Active particles do not feel test particles
    Passive,
}

pub enum RemoveKind {
    Active(usize),
    Test(usize),
}
