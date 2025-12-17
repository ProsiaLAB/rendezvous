use itertools::iproduct;
use prosia_extensions::types::Vec3;
use rayon::iter::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::{
    boundary::{Boundary, BoundaryContext, GhostBox},
    integrator::Integrator,
    particle::{Particle, TestParticleType},
};

pub enum Gravity {
    None,
    Basic,
    Compensated,
    Tree,
    Mercurius,
    Jacobi,
    Trace,
}

/// Which gravity terms to ignore in the simulation.
pub enum IgnoreGravityTerms {
    /// Include all gravity terms.
    IgnoreAll,
    /// Ignore gravity terms not required for WHFast
    /// with Jacobi coordinates.
    IgnoreWHFastwithJacobi,
    /// Ignore gravity terms not required for WHFast
    /// with Democratic Heliocentric Coordinates.
    IgnoreWHFastwithDHC,
}

impl IgnoreGravityTerms {
    fn for_massive_test_particles(&self, i: usize, j: usize) -> bool {
        match self {
            IgnoreGravityTerms::IgnoreWHFastwithJacobi => j == 1 && i == 0,
            IgnoreGravityTerms::IgnoreWHFastwithDHC => j == 0 && i == 0,
            _ => false,
        }
    }

    fn for_active_particles(&self, i: usize, j: usize) -> bool {
        if i == j {
            return true;
        }

        match self {
            IgnoreGravityTerms::IgnoreWHFastwithJacobi => (j == 1 && i == 0) || (i == 1 && j == 0),
            IgnoreGravityTerms::IgnoreWHFastwithDHC => j == 0 && i == 0,
            _ => false,
        }
    }
}

pub struct GravityContext<'a> {
    pub particles: &'a mut [Particle],
    pub n: usize,
    pub n_active: usize,
    pub g: f64,
    pub t: f64,
    pub integrator: &'a Integrator,
    pub boundary: &'a Boundary,
    pub n_ghost_x: usize,
    pub n_ghost_y: usize,
    pub n_ghost_z: usize,
    pub box_size: &'a Vec3,
    pub ignore_gravity_terms: &'a IgnoreGravityTerms,
    pub softening: f64,
    pub test_particle_type: &'a TestParticleType,
}

impl GravityContext<'_> {
    pub fn get_begin_indices(&self) -> (usize, usize) {
        let start_i = if matches!(self.ignore_gravity_terms, IgnoreGravityTerms::IgnoreAll) {
            1
        } else {
            2
        };
        let start_j = if matches!(
            self.ignore_gravity_terms,
            IgnoreGravityTerms::IgnoreWHFastwithDHC
        ) {
            1
        } else {
            0
        };

        (start_i, start_j)
    }

    fn update_accels<F>(
        &mut self,
        i_range: std::ops::Range<usize>,
        j_range: std::ops::Range<usize>,
        soft2: f64,
        gb: &GhostBox,
        ignore: F,
    ) where
        F: Fn(usize, usize) -> bool + Sync,
    {
        let accels: Vec<_> = i_range
            .into_par_iter()
            .map(|i| {
                let mut ax = 0.0;
                let mut ay = 0.0;
                let mut az = 0.0;

                for j in j_range.clone() {
                    if ignore(i, j) {
                        continue;
                    }

                    let dx = gb.position.x + self.particles[i].x - self.particles[j].x;
                    let dy = gb.position.y + self.particles[i].y - self.particles[j].y;
                    let dz = gb.position.z + self.particles[i].z - self.particles[j].z;

                    let dr = (dx * dx + dy * dy + dz * dz + soft2).sqrt();
                    let prefactor = -self.g / (dr * dr * dr) * self.particles[i].m;

                    ax += prefactor * dx;
                    ay += prefactor * dy;
                    az += prefactor * dz;
                }

                (ax, ay, az)
            })
            .collect();

        for (p, (ax, ay, az)) in self.particles.iter_mut().zip(accels.into_iter()) {
            p.ax += ax;
            p.ay += ay;
            p.az += az;
        }
    }

    fn apply_none(&mut self) {
        self.particles.iter_mut().for_each(|p| {
            p.ax = 0.0;
            p.ay = 0.0;
            p.az = 0.0;
        });
    }

    fn apply_jacobi(&mut self) {
        if !matches!(self.integrator, Integrator::WHFast(_) | Integrator::Saba(_)) {
            eprintln!(
                "WARNING: Jacobi gravity is intended to be used with WHFast or Saba integrators."
            )
        }
        let mut rjx = 0.0;
        let mut rjy = 0.0;
        let mut rjz = 0.0;
        let mut mj = 0.0;
        for j in 0..self.particles.len() {
            self.particles[j].ax = 0.0;
            self.particles[j].ay = 0.0;
            self.particles[j].az = 0.0;
            for i in 0..j + 1 {
                if j > 1 {
                    // Jacobi term
                    let qjx = self.particles[j].x - rjx / mj;
                    let qjy = self.particles[j].y - rjy / mj;
                    let qjz = self.particles[j].z - rjz / mj;
                    let dr = (qjx * qjx + qjy * qjy + qjz * qjz).sqrt();
                    let dqjdri = if i < j { -self.particles[i].m } else { mj };
                    let prefactor = self.g * dqjdri / (dr * dr * dr);
                    self.particles[j].ax += prefactor * qjx;
                    self.particles[j].ay += prefactor * qjy;
                    self.particles[j].az += prefactor * qjz;
                }
                if i != j && (i != 0 || j != 1) {
                    // Direct term
                    let dx = self.particles[i].x - self.particles[j].x;
                    let dy = self.particles[i].y - self.particles[j].y;
                    let dz = self.particles[i].z - self.particles[j].z;
                    let dr = (dx * dx + dy * dy + dz * dz).sqrt();
                    let prefactor = self.g / (dr * dr * dr);
                    let prefactor_i = prefactor * self.particles[i].m;
                    let prefactor_j = prefactor * self.particles[j].m;

                    self.particles[j].ax -= prefactor_j * dx;
                    self.particles[j].ay -= prefactor_j * dy;
                    self.particles[j].az -= prefactor_j * dz;
                    self.particles[i].ax += prefactor_i * dx;
                    self.particles[i].ay += prefactor_i * dy;
                    self.particles[i].az += prefactor_i * dz;
                }
            }
            rjx += self.particles[j].m * self.particles[j].x;
            rjy += self.particles[j].m * self.particles[j].y;
            rjz += self.particles[j].m * self.particles[j].z;
            mj += self.particles[j].m;
        }
    }

    fn apply_basic(&mut self) {
        let n_active = if self.n_active == usize::MAX {
            self.n
        } else {
            self.n_active
        };

        let soft2 = self.softening * self.softening;

        self.particles.par_iter_mut().for_each(|p| {
            p.ax = 0.0;
            p.ay = 0.0;
            p.az = 0.0;
        });

        let ngx = self.n_ghost_x as isize;
        let ngy = self.n_ghost_y as isize;
        let ngz = self.n_ghost_z as isize;

        let boundary_ctx = BoundaryContext {
            t: self.t,
            box_size: self.box_size,
            integrator: self.integrator,
        };

        for (gbx, gby, gbz) in iproduct!(-ngx..=ngx, -ngy..=ngy, -ngz..=ngz) {
            let gb = self.boundary.get_ghost_box(&boundary_ctx, gbx, gby, gbz);

            self.update_accels(0..self.n, 0..n_active, soft2, &gb, |i, j| {
                self.ignore_gravity_terms.for_active_particles(i, j)
            });

            if matches!(self.test_particle_type, TestParticleType::Massive) {
                self.update_accels(0..n_active, n_active..self.n, soft2, &gb, |i, j| {
                    self.ignore_gravity_terms.for_massive_test_particles(i, j)
                });
            }
        }
    }

    fn apply_compensated(&mut self) {
        todo!()
    }

    fn apply_tree(&mut self) {
        todo!()
    }

    fn apply_mercurius(&mut self) {
        todo!()
    }

    fn apply_trace(&mut self) {
        todo!()
    }
}

impl Gravity {
    pub fn apply(&self, ctx: &mut GravityContext<'_>) {
        match self {
            Gravity::None => {
                ctx.apply_none();
            }
            Gravity::Jacobi => {
                ctx.apply_jacobi();
            }
            Gravity::Basic => {
                ctx.apply_basic();
            }
            Gravity::Compensated => {
                ctx.apply_compensated();
            }
            Gravity::Tree => {
                ctx.apply_tree();
            }
            Gravity::Mercurius => {
                ctx.apply_mercurius();
            }
            Gravity::Trace => {
                ctx.apply_trace();
            }
        }
    }
}
