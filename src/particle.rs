use std::ops::{Index, IndexMut};

use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

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

    fn transform_dhc_to_inertial_posvel(
        &mut self,
        dhc: &[Particle],
        n_real: usize,
        n_active: usize,
    );

    fn transform_dhc_to_inertial_pos(&mut self, dhc: &[Particle], n_real: usize, n_active: usize);

    fn transform_whds_to_inertial_pos(&mut self, whds: &[Particle], n_real: usize, n_active: usize);

    fn transform_whds_to_inertial_posvel(
        &mut self,
        whds: &[Particle],
        n_real: usize,
        n_active: usize,
    );

    fn transform_barycentric_to_inertial_posvel(
        &mut self,
        bary: &[Particle],
        n_real: usize,
        n_active: usize,
    );
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

    fn transform_dhc_to_inertial_pos(&mut self, dhc: &[Particle], n_real: usize, n_active: usize) {
        let mtot = dhc[0].m;

        let (x0, y0, z0) = self[1..n_active]
            .par_iter_mut()
            .zip(&dhc[1..n_active])
            .map(|(p, ph)| {
                let m = ph.m;
                p.m = m;
                (ph.x * m / mtot, ph.y * m / mtot, ph.z * m / mtot)
            })
            .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

        self[0].x = dhc[0].x - x0;
        self[0].y = dhc[0].y - y0;
        self[0].z = dhc[0].z - z0;

        let x0 = self[0].x;
        let y0 = self[0].y;
        let z0 = self[0].z;

        self[1..n_real]
            .par_iter_mut()
            .zip(&dhc[1..n_real])
            .for_each(|(p, ph)| {
                p.x = ph.x + x0;
                p.y = ph.y + y0;
                p.z = ph.z + z0;
            });
    }

    fn transform_dhc_to_inertial_posvel(
        &mut self,
        dhc: &[Particle],
        n_real: usize,
        n_active: usize,
    ) {
        self.transform_dhc_to_inertial_pos(dhc, n_real, n_active);

        let m0 = self[0].m;

        let vx0 = self[0].vx;
        let vy0 = self[0].vy;
        let vz0 = self[0].vz;

        self[1..n_real]
            .par_iter_mut()
            .zip(&dhc[1..n_real])
            .for_each(|(p, ph)| {
                p.vx = ph.vx + vx0;
                p.vy = ph.vy + vy0;
                p.vz = ph.vz + vz0;
            });

        let (vx0, vy0, vz0) = self[1..n_active]
            .par_iter_mut()
            .zip(&dhc[1..n_active])
            .map(|(p, ph)| {
                let m = p.m;
                (ph.vx * m / m0, ph.vy * m / m0, ph.vz * m / m0)
            })
            .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

        self[0].vx = dhc[0].vx - vx0;
        self[0].vy = dhc[0].vy - vy0;
        self[0].vz = dhc[0].vz - vz0;
    }

    fn transform_whds_to_inertial_pos(
        &mut self,
        whds: &[Particle],
        n_real: usize,
        n_active: usize,
    ) {
        self.transform_dhc_to_inertial_pos(whds, n_real, n_active);
    }

    fn transform_whds_to_inertial_posvel(
        &mut self,
        whds: &[Particle],
        n_real: usize,
        n_active: usize,
    ) {
        self.transform_whds_to_inertial_pos(whds, n_real, n_active);

        let m0 = self[0].m;

        let whds_vx0 = whds[0].vx;
        let whds_vy0 = whds[0].vy;
        let whds_vz0 = whds[0].vz;

        self[1..n_active]
            .par_iter_mut()
            .zip(&whds[1..n_active])
            .for_each(|(p, ph)| {
                let mi = p.m;
                let mf = (m0 + mi) / m0;
                p.vx = ph.vx / mf + whds_vx0;
                p.vy = ph.vy / mf + whds_vy0;
                p.vz = ph.vz / mf + whds_vz0;
            });

        self[n_active..n_real]
            .par_iter_mut()
            .zip(&whds[n_active..n_real])
            .for_each(|(p, ph)| {
                p.vx = ph.vx + whds_vx0;
                p.vy = ph.vy + whds_vy0;
                p.vz = ph.vz + whds_vz0;
            });

        let (vx0, vy0, vz0) = self[1..n_active]
            .par_iter_mut()
            .zip(&whds[1..n_active])
            .map(|(p, ph)| {
                let m = p.m;
                (
                    ph.vx * m / (m0 + m),
                    ph.vy * m / (m0 + m),
                    ph.vz * m / (m0 + m),
                )
            })
            .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

        self[0].vx = whds[0].vx - vx0;
        self[0].vy = whds[0].vy - vy0;
        self[0].vz = whds[0].vz - vz0;
    }

    fn transform_barycentric_to_inertial_posvel(
        &mut self,
        bary: &[Particle],
        n_real: usize,
        n_active: usize,
    ) {
        self[0].x = bary[0].m * bary[0].x;
        self[0].y = bary[0].m * bary[0].y;
        self[0].z = bary[0].m * bary[0].z;
        self[0].vx = bary[0].m * bary[0].vx;
        self[0].vy = bary[0].m * bary[0].vy;
        self[0].vz = bary[0].m * bary[0].vz;
        self[0].m = bary[0].m;

        let mut s_x = 0.0;
        let mut s_y = 0.0;
        let mut s_z = 0.0;
        let mut s_vx = 0.0;
        let mut s_vy = 0.0;
        let mut s_vz = 0.0;
        let mut s_m = 0.0;

        for i in 1..n_real {
            let p_bi = &bary[i];
            let p_b0 = &bary[0];
            let pi = &mut self[i];
            pi.x = p_bi.x + p_b0.x;
            pi.y = p_bi.y + p_b0.y;
            pi.z = p_bi.z + p_b0.z;
            pi.vx = p_bi.vx + p_b0.vx;
            pi.vy = p_bi.vy + p_b0.vy;
            pi.vz = p_bi.vz + p_b0.vz;
            if i < n_active {
                let m = p_bi.m;
                pi.m = m;
                s_x += pi.x * m;
                s_y += pi.y * m;
                s_z += pi.z * m;
                s_vx += pi.vx * m;
                s_vy += pi.vy * m;
                s_vz += pi.vz * m;
                s_m += m;
            }
        }

        self[0].x -= s_x;
        self[0].y -= s_y;
        self[0].z -= s_z;
        self[0].vx -= s_vx;
        self[0].vy -= s_vy;
        self[0].vz -= s_vz;
        self[0].m -= s_m;

        let m0i = 1.0 / self[0].m;
        self[0].x *= m0i;
        self[0].y *= m0i;
        self[0].z *= m0i;
        self[0].vx *= m0i;
        self[0].vy *= m0i;
        self[0].vz *= m0i;
    }
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
