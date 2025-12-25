use std::ops::{Index, IndexMut};

use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator};

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

#[derive(Default)]
pub struct Particles {
    data: Vec<Particle>,
    n_active: usize,
    n_test: usize,
    n_variational: usize,
    all_active: bool,
}

impl Particles {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn n_real(&self) -> usize {
        self.n_active + self.n_test
    }

    pub fn n_active(&self) -> usize {
        self.n_active
    }

    pub fn active(&self) -> &[Particle] {
        &self.data[..self.n_active]
    }

    pub fn active_mut(&mut self) -> &mut [Particle] {
        &mut self.data[..self.n_active]
    }

    pub fn test(&self) -> &[Particle] {
        &self.data[self.n_active..self.n_active + self.n_test]
    }

    pub fn test_mut(&mut self) -> &mut [Particle] {
        &mut self.data[self.n_active..self.n_active + self.n_test]
    }

    pub fn variational(&self) -> &[Particle] {
        &self.data[self.n_active + self.n_test..]
    }

    pub fn variational_mut(&mut self) -> &mut [Particle] {
        &mut self.data[self.n_active + self.n_test..]
    }

    pub fn n_test(&self) -> usize {
        self.n_test
    }

    pub fn n_variational(&self) -> usize {
        self.n_variational
    }

    pub fn any_test(&self) -> bool {
        self.n_test != 0
    }

    pub fn any_variational(&self) -> bool {
        self.n_variational != 0
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn are_all_active(&self) -> bool {
        self.all_active
    }

    pub fn push(&mut self, p: Particle) {
        self.data.push(p);
    }

    pub fn as_slice(&self) -> &[Particle] {
        &self.data
    }

    pub fn as_mut_slice(&mut self) -> &mut [Particle] {
        &mut self.data
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
    pub fn real_par_iter(&self) -> impl IndexedParallelIterator<Item = &Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].par_iter()
    }

    /// Returns a mutable parallel iterator over all real (non-variational) particles.
    pub fn real_par_iter_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Particle> {
        let n_real = self.n_active + self.n_test;
        self.data[..n_real].par_iter_mut()
    }

    /// Returns a parallel iterator over all particles.
    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = &Particle> {
        self.data.par_iter()
    }

    /// Returns a mutable parallel iterator over all particles.
    pub fn par_iter_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Particle> {
        self.data.par_iter_mut()
    }

    pub fn resize(&mut self, new_size: usize) {
        self.data.resize(new_size, Particle::default());
    }

    pub fn remove(&mut self, index: usize) -> Particle {
        let removed = self.data.remove(index);

        if index < self.n_active {
            self.n_active -= 1;
        }

        removed
    }

    pub fn swap_remove(&mut self, index: usize) {
        self.data.swap_remove(index);
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        self.data.swap(i, j);
    }

    pub fn pop(&mut self) -> Option<Particle> {
        self.data.pop()
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.n_active = 0;
        self.n_test = 0;
        self.n_variational = 0;
    }
}

pub trait Transformations {
    fn transform_jacobi_to_inertial_posvel(
        &mut self,
        jacobi: &[Particle],
        masses: &[f64],
        n_real: usize,
        n_active: usize,
    );

    fn transform_dhc_to_inertial_posvel(&mut self);
    fn transform_whds_to_inertial_posvel(&mut self);
    fn transform_barycentric_to_inertial_posvel(&mut self);
}

impl Transformations for [Particle] {
    fn transform_jacobi_to_inertial_posvel(
        &mut self,
        jacobi: &[Particle],
        masses: &[f64],
        n_real: usize,
        n_active: usize,
    ) {
        let centre = &jacobi[0];
        let mut eta = centre.m;
        let mut s_x = centre.x * eta;
        let mut s_y = centre.y * eta;
        let mut s_z = centre.z * eta;
        let mut s_vx = centre.vx * eta;
        let mut s_vy = centre.vy * eta;
        let mut s_vz = centre.vz * eta;

        for i in (n_active..n_real).rev() {
            let p_j = &jacobi[i];
            let ei = 1.0 / eta;
            self[i].x = p_j.x + s_x * ei;
            self[i].y = p_j.y + s_y * ei;
            self[i].z = p_j.z + s_z * ei;
            self[i].vx = p_j.vx + s_vx * ei;
            self[i].vy = p_j.vy + s_vy * ei;
            self[i].vz = p_j.vz + s_vz * ei;
        }

        for i in (1..n_active).rev() {
            let p_j = &jacobi[i];
            let ei = 1.0 / eta;
            s_x = (s_x - masses[i] * p_j.x) * ei;
            s_y = (s_y - masses[i] * p_j.y) * ei;
            s_z = (s_z - masses[i] * p_j.z) * ei;
            s_vx = (s_vx - masses[i] * p_j.vx) * ei;
            s_vy = (s_vy - masses[i] * p_j.vy) * ei;
            s_vz = (s_vz - masses[i] * p_j.vz) * ei;
            self[i].x = p_j.x + s_x;
            self[i].y = p_j.y + s_y;
            self[i].z = p_j.z + s_z;
            self[i].vx = p_j.vx + s_vx;
            self[i].vy = p_j.vy + s_vy;
            self[i].vz = p_j.vz + s_vz;
            eta -= masses[i];
            s_x *= eta;
            s_y *= eta;
            s_z *= eta;
            s_vx *= eta;
            s_vy *= eta;
            s_vz *= eta;
        }

        let mi = 1.0 / eta;
        self[0].x = s_x * mi;
        self[0].y = s_y * mi;
        self[0].z = s_z * mi;
        self[0].vx = s_vx * mi;
        self[0].vy = s_vy * mi;
        self[0].vz = s_vz * mi;
    }

    fn transform_dhc_to_inertial_posvel(&mut self) {}
    fn transform_whds_to_inertial_posvel(&mut self) {}
    fn transform_barycentric_to_inertial_posvel(&mut self) {}
}

impl Index<usize> for Particles {
    type Output = Particle;

    fn index(&self, i: usize) -> &Self::Output {
        debug_assert!(i < self.n_active + self.n_test + self.n_variational);
        &self.data[i]
    }
}

impl IndexMut<usize> for Particles {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        debug_assert!(i < self.n_active + self.n_test + self.n_variational);
        &mut self.data[i]
    }
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
