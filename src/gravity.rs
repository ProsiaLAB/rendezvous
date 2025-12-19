use itertools::iproduct;
use prosia_extensions::types::Vec3;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator,
};

use crate::boundary::{Boundary, BoundaryContext, GhostBox};
use crate::integrator::Integrator;
use crate::mercurius::{Mercurius, MercuriusMode};
use crate::particle::{Particle, TestParticleType};
use crate::trace::{Trace, TraceMode};
use crate::tree::{NodeId, NodeKind, TreeType};

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
            IgnoreGravityTerms::IgnoreWHFastwithJacobi => (j == 1 && i == 0) || (i == 1 && j == 0),
            IgnoreGravityTerms::IgnoreWHFastwithDHC => i == 0 || j == 0,
            _ => false,
        }
    }

    fn for_active_particles(&self, i: usize, j: usize) -> bool {
        if i == j {
            return true;
        }

        match self {
            IgnoreGravityTerms::IgnoreWHFastwithJacobi => (j == 1 && i == 0) || (i == 1 && j == 0),
            IgnoreGravityTerms::IgnoreWHFastwithDHC => j == 0 || i == 0,
            _ => false,
        }
    }
}

pub struct GravityContext<'a> {
    pub particles: &'a mut [Particle],
    pub n_real: usize,
    pub n_active: usize,
    pub n_root: usize,
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
    pub gravity_cs: &'a mut [Vec3],
    pub tree: Option<&'a TreeType>,
    pub opening_angle: f64,
}

impl GravityContext<'_> {
    fn update_basic<F>(
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

    fn update_compensated_ij<F>(
        &mut self,
        i_range: std::ops::Range<usize>,
        j_range: std::ops::Range<usize>,
        soft2: f64,
        ignore: F,
    ) where
        F: Fn(usize, usize) -> bool + Sync,
    {
        let updates: Vec<_> = i_range
            .into_par_iter()
            .map(|i| {
                let mut cs = Vec3::new(0.0, 0.0, 0.0);
                let mut acc = Vec3::new(0.0, 0.0, 0.0);

                for j in j_range.clone() {
                    if ignore(i, j) {
                        continue;
                    }

                    let dx = self.particles[i].x - self.particles[j].x;
                    let dy = self.particles[i].y - self.particles[j].y;
                    let dz = self.particles[i].z - self.particles[j].z;

                    let r2 = dx * dx + dy * dy + dz * dz + soft2;
                    let r = r2.sqrt();
                    let prefactor = self.g / (r2 * r);
                    let prefactor_j = -prefactor * self.particles[j].m;

                    let ix = prefactor_j * dx;
                    let yx = ix - cs.x;
                    let tx = acc.x + yx;
                    cs.x = (tx - acc.x) - yx;
                    acc.x = tx;

                    let iy = prefactor_j * dy;
                    let yy = iy - cs.y;
                    let ty = acc.y + yy;
                    cs.y = (ty - acc.y) - yy;
                    acc.y = ty;

                    let iz = prefactor_j * dz;
                    let zy = iz - cs.z;
                    let tz = acc.z + zy;
                    cs.z = (tz - acc.z) - zy;
                    acc.z = tz;
                }
                (cs, acc)
            })
            .collect();

        for (i, (cs, acc)) in updates.into_iter().enumerate() {
            self.gravity_cs[i] = cs;
            self.particles[i].ax = acc.x;
            self.particles[i].ay = acc.y;
            self.particles[i].az = acc.z;
        }
    }

    fn update_compensated_ji<F>(
        &mut self,
        j_range: std::ops::Range<usize>,
        i_range: std::ops::Range<usize>,
        soft2: f64,
        ignore: F,
    ) where
        F: Fn(usize, usize) -> bool + Sync,
    {
        let updates: Vec<_> = j_range
            .into_par_iter()
            .map(|j| {
                let mut cs = Vec3::new(0.0, 0.0, 0.0);
                let mut acc = Vec3::new(0.0, 0.0, 0.0);

                for i in i_range.clone() {
                    if ignore(i, j) {
                        continue;
                    }

                    let dx = self.particles[i].x - self.particles[j].x;
                    let dy = self.particles[i].y - self.particles[j].y;
                    let dz = self.particles[i].z - self.particles[j].z;

                    let r2 = dx * dx + dy * dy + dz * dz + soft2;
                    let r = r2.sqrt();
                    let prefactor = self.g / (r2 * r);
                    let prefactor_i = prefactor * self.particles[i].m;

                    let ix = prefactor_i * dx;
                    let yx = ix - cs.x;
                    let tx = acc.x + yx;
                    cs.x = (tx - acc.x) - yx;
                    acc.x = tx;

                    let iy = prefactor_i * dy;
                    let yy = iy - cs.y;
                    let ty = acc.y + yy;
                    cs.y = (ty - acc.y) - yy;
                    acc.y = ty;

                    let iz = prefactor_i * dz;
                    let zy = iz - cs.z;
                    let tz = acc.z + zy;
                    cs.z = (tz - acc.z) - zy;
                    acc.z = tz;
                }
                (cs, acc)
            })
            .collect();

        for (i, (cs, acc)) in updates.into_iter().enumerate() {
            self.gravity_cs[i] = cs;
            self.particles[i].ax = acc.x;
            self.particles[i].ay = acc.y;
            self.particles[i].az = acc.z;
        }
    }

    fn calculate_acceleration_for_particle(
        &self,
        pi: usize,
        gb: &GhostBox,
        accum: &mut (f64, f64, f64),
    ) {
        for i in 0..self.n_root {
            self.calculate_acceleration_from_node(NodeId(i), pi, gb, accum);
        }
    }

    fn calculate_acceleration_from_node(
        &self,
        node_id: NodeId,
        pi: usize,
        gb: &GhostBox,
        accum: &mut (f64, f64, f64),
    ) {
        let node = self.tree.unwrap().get_node(node_id);

        let soft2 = self.softening * self.softening;
        let dx = gb.position.x - node.mx;
        let dy = gb.position.y - node.my;
        let dz = gb.position.z - node.mz;

        let r2 = dx * dx + dy * dy + dz * dz;

        match node.kind {
            NodeKind::NonLeaf => {
                if node.w * node.w > self.opening_angle * r2 {
                    node.children.iter().for_each(|child_opt| {
                        if let Some(child_id) = child_opt {
                            self.calculate_acceleration_from_node(*child_id, pi, gb, accum);
                        }
                    });
                } else {
                    let r = (r2 + soft2).sqrt();
                    let prefactor = -self.g / (r * r * r) * node.m;

                    #[cfg(feature = "quadrupole")]
                    {
                        let mut q_prefactor = self.g / (r * r * r * r * r);
                        accum.0 += q_prefactor
                            * (dx * node.moment.mxx + dy * node.moment.mxy + dz * node.moment.mxz);
                        accum.1 += q_prefactor
                            * (dx * node.moment.mxy + dy * node.moment.myy + dz * node.moment.myz);
                        accum.2 += q_prefactor
                            * (dx * node.moment.mxz + dy * node.moment.myz + dz * node.moment.mzz);

                        let mrr = node.moment.mxx * dx * dx
                            + 2.0 * node.moment.mxy * dx * dy
                            + 2.0 * node.moment.mxz * dx * dz
                            + node.moment.myy * dy * dy
                            + 2.0 * node.moment.myz * dy * dz
                            + node.moment.mzz * dz * dz;

                        q_prefactor *= -5.0 / (2.0 * r * r) * mrr;

                        accum.0 += (q_prefactor + prefactor) * dx;
                        accum.1 += (q_prefactor + prefactor) * dy;
                        accum.2 += (q_prefactor + prefactor) * dz;
                    }

                    #[cfg(not(feature = "quadrupole"))]
                    {
                        accum.0 += prefactor * dx;
                        accum.1 += prefactor * dy;
                        accum.2 += prefactor * dz;
                    }
                }
            }
            NodeKind::Leaf(pt) => {
                if !node.remote && pt == pi {
                    // Skip self-interaction
                } else {
                    let r = (r2 + soft2).sqrt();
                    let prefactor = -self.g / (r * r * r) * node.m;
                    accum.0 += prefactor * dx;
                    accum.1 += prefactor * dy;
                    accum.2 += prefactor * dz;
                }
            }
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

    fn apply_basic(&mut self, n_active: usize, soft2: f64) {
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

            self.update_basic(0..self.n_real, 0..n_active, soft2, &gb, |i, j| {
                self.ignore_gravity_terms.for_active_particles(i, j)
            });

            if matches!(self.test_particle_type, TestParticleType::Massive) {
                self.update_basic(0..n_active, n_active..self.n_real, soft2, &gb, |i, j| {
                    self.ignore_gravity_terms.for_massive_test_particles(i, j)
                });
            }
        }
    }

    fn apply_compensated(&mut self, n_active: usize, soft2: f64) {
        self.particles
            .par_iter_mut()
            .take(self.n_real)
            .zip(self.gravity_cs.par_iter_mut().take(self.n_real))
            .for_each(|(p, cs)| {
                p.ax = 0.0;
                p.ay = 0.0;
                p.az = 0.0;

                cs.x = 0.0;
                cs.y = 0.0;
                cs.z = 0.0;
            });

        self.update_compensated_ij(0..n_active, 0..n_active, soft2, |i, j| {
            self.ignore_gravity_terms.for_active_particles(i, j)
        });

        self.update_compensated_ij(n_active..self.n_real, 0..n_active, soft2, |i, j| {
            self.ignore_gravity_terms.for_massive_test_particles(i, j)
        });

        if matches!(self.test_particle_type, TestParticleType::Massive) {
            self.update_compensated_ji(0..n_active, n_active..self.n_real, soft2, |i, j| {
                self.ignore_gravity_terms.for_massive_test_particles(i, j)
            });
        }
    }

    fn apply_tree(&mut self) {
        self.particles.iter_mut().for_each(|p| {
            p.ax = 0.0;
            p.ay = 0.0;
            p.az = 0.0;
        });

        let ngx = self.n_ghost_x as isize;
        let ngy = self.n_ghost_y as isize;
        let ngz = self.n_ghost_z as isize;

        for (gbx, gby, gbz) in iproduct!(-ngx..=ngx, -ngy..=ngy, -ngz..=ngz) {
            let boundary_ctx = BoundaryContext {
                t: self.t,
                box_size: self.box_size,
                integrator: self.integrator,
            };

            let accum: Vec<(f64, f64, f64)> = (0..self.particles.len())
                .into_par_iter()
                .map(|i| {
                    let mut gb = self.boundary.get_ghost_box(&boundary_ctx, gbx, gby, gbz);
                    gb.position.x += self.particles[i].x;
                    gb.position.y += self.particles[i].y;
                    gb.position.z += self.particles[i].z;
                    let mut accum = (0.0, 0.0, 0.0);
                    self.calculate_acceleration_for_particle(i, &gb, &mut accum);
                    accum
                })
                .collect();

            for (i, (ax, ay, az)) in accum.into_iter().enumerate() {
                self.particles[i].ax += ax;
                self.particles[i].ay += ay;
                self.particles[i].az += az;
            }
        }
    }

    fn apply_mercurius(&mut self, n_active: usize, soft2: f64) {
        match self.integrator {
            Integrator::Mercurius(m) => match m.mode {
                MercuriusMode::LongRange => {
                    self.apply_mercurius_long_range(m, n_active, soft2);
                }
                MercuriusMode::CloseEncounter => {
                    self.apply_mercurius_close_encounter(m, soft2);
                }
            },
            _ => {
                eprintln!(
                    "WARNING: Mercurius gravity is intended to be used with the Mercurius integrator."
                )
            }
        }
    }

    fn apply_mercurius_long_range(&mut self, m: &Mercurius, n_active: usize, soft2: f64) {
        self.particles[0].ax = 0.0;
        self.particles[0].ay = 0.0;
        self.particles[0].az = 0.0;

        let accels: Vec<(f64, f64, f64)> = (1..self.n_real)
            .into_par_iter()
            .map(|i| {
                let mut ax = 0.0;
                let mut ay = 0.0;
                let mut az = 0.0;

                for j in 1..n_active {
                    if i == j {
                        continue;
                    }

                    let dx = self.particles[i].x - self.particles[j].x;
                    let dy = self.particles[i].y - self.particles[j].y;
                    let dz = self.particles[i].z - self.particles[j].z;

                    let dr = (dx * dx + dy * dy + dz * dz + soft2).sqrt();
                    let dcrit_max = m.dcrit[i].max(m.dcrit[j]);
                    let val = m.switch(dr, dcrit_max);
                    let prefactor = -self.g * self.particles[j].m * val / (dr * dr * dr);
                    ax += prefactor * dx;
                    ay += prefactor * dy;
                    az += prefactor * dz;
                }
                (ax, ay, az)
            })
            .collect();

        for (i, (ax, ay, az)) in accels.into_iter().enumerate() {
            self.particles[i + 1].ax = ax;
            self.particles[i + 1].ay = ay;
            self.particles[i + 1].az = az;
        }

        if matches!(self.test_particle_type, TestParticleType::Massive) {
            let accels_tp: Vec<(f64, f64, f64)> = (1..n_active)
                .into_par_iter()
                .map(|i| {
                    let mut ax = 0.0;
                    let mut ay = 0.0;
                    let mut az = 0.0;

                    for j in n_active..self.n_real {
                        let dx = self.particles[i].x - self.particles[j].x;
                        let dy = self.particles[i].y - self.particles[j].y;
                        let dz = self.particles[i].z - self.particles[j].z;

                        let dr = (dx * dx + dy * dy + dz * dz + soft2).sqrt();
                        let dcrit_max = m.dcrit[i].max(m.dcrit[j]);
                        let val = m.switch(dr, dcrit_max);
                        let prefactor = -self.g * self.particles[j].m * val / (dr * dr * dr);
                        ax += prefactor * dx;
                        ay += prefactor * dy;
                        az += prefactor * dz;
                    }
                    (ax, ay, az)
                })
                .collect();

            for (i, (ax, ay, az)) in accels_tp.into_iter().enumerate() {
                self.particles[i + 1].ax += ax;
                self.particles[i + 1].ay += ay;
                self.particles[i + 1].az += az;
            }
        }
    }

    fn apply_mercurius_close_encounter(&mut self, m: &Mercurius, soft2: f64) {
        self.particles[0].ax = 0.0;
        self.particles[0].ay = 0.0;
        self.particles[0].az = 0.0;

        // In heliocentric coordinates, the star feels no acceleration
        let accels: Vec<(f64, f64, f64)> = (1..m.n_encounter)
            .into_par_iter()
            .map(|i| {
                let mi = m.encounter_map[i];

                let mut ax = 0.0;
                let mut ay = 0.0;
                let mut az = 0.0;

                // Acceleration due to star
                let x = self.particles[mi].x;
                let y = self.particles[mi].y;
                let z = self.particles[mi].z;
                let r = (x * x + y * y + z * z + soft2).sqrt();
                let prefactor = -self.g / (r * r * r) * self.particles[0].m;
                ax += prefactor * x;
                ay += prefactor * y;
                az += prefactor * z;

                for j in 1..m.n_encounter_active {
                    if i == j {
                        continue;
                    }

                    let mj = m.encounter_map[j];

                    let dx = x - self.particles[mj].x;
                    let dy = y - self.particles[mj].y;
                    let dz = z - self.particles[mj].z;
                    let r = (dx * dx + dy * dy + dz * dz + soft2).sqrt();
                    let dcrit_max = m.dcrit[mi].max(m.dcrit[mj]);
                    let val = m.switch(r, dcrit_max);
                    let prefactor = -self.g * self.particles[mj].m * (1.0 - val) / (r * r * r);
                    ax += prefactor * dx;
                    ay += prefactor * dy;
                    az += prefactor * dz;
                }
                (ax, ay, az)
            })
            .collect();

        for (i, (ax, ay, az)) in accels.into_iter().enumerate() {
            let mi = m.encounter_map[i + 1];
            self.particles[mi].ax = ax;
            self.particles[mi].ay = ay;
            self.particles[mi].az = az;
        }

        if matches!(self.test_particle_type, TestParticleType::Massive) {
            let accels: Vec<(f64, f64, f64)> = (1..m.n_encounter_active)
                .into_par_iter()
                .map(|i| {
                    let mi = m.encounter_map[i];

                    let mut ax = 0.0;
                    let mut ay = 0.0;
                    let mut az = 0.0;

                    let x = self.particles[mi].x;
                    let y = self.particles[mi].y;
                    let z = self.particles[mi].z;

                    for j in m.n_encounter_active..m.n_encounter {
                        let mj = m.encounter_map[j];

                        let dx = x - self.particles[mj].x;
                        let dy = y - self.particles[mj].y;
                        let dz = z - self.particles[mj].z;
                        let r = (dx * dx + dy * dy + dz * dz + soft2).sqrt();
                        let dcrit_max = m.dcrit[mi].max(m.dcrit[mj]);
                        let val = m.switch(r, dcrit_max);
                        let prefactor = -self.g * self.particles[mj].m * (1.0 - val) / (r * r * r);
                        ax += prefactor * dx;
                        ay += prefactor * dy;
                        az += prefactor * dz;
                    }
                    (ax, ay, az)
                })
                .collect();

            for (i, (ax, ay, az)) in accels.into_iter().enumerate() {
                let mi = m.encounter_map[i + 1];
                self.particles[mi].ax += ax;
                self.particles[mi].ay += ay;
                self.particles[mi].az += az;
            }
        }
    }

    fn apply_trace(&mut self) {
        match self.integrator {
            Integrator::Trace(t) => match t.mode {
                TraceMode::Interaction => {
                    self.apply_trace_interaction(t);
                }
                TraceMode::Kepler => {
                    self.apply_trace_kepler(t);
                }
                _ => {}
            },
            _ => {
                eprintln!(
                    "WARNING: Trace gravity is intended to be used with the Trace integrator."
                )
            }
        }
    }

    fn apply_trace_interaction(&mut self, t: &Trace) {
        todo!()
    }

    fn apply_trace_kepler(&mut self, t: &Trace) {
        todo!()
    }
}

impl Gravity {
    pub fn apply(&self, ctx: &mut GravityContext<'_>) {
        let n_active = if ctx.n_active == usize::MAX {
            ctx.n_real
        } else {
            ctx.n_active
        };

        let soft2 = ctx.softening * ctx.softening;
        match self {
            Gravity::None => {
                ctx.apply_none();
            }
            Gravity::Jacobi => {
                ctx.apply_jacobi();
            }
            Gravity::Basic => {
                ctx.apply_basic(n_active, soft2);
            }
            Gravity::Compensated => {
                ctx.apply_compensated(n_active, soft2);
            }
            Gravity::Tree => {
                ctx.apply_tree();
            }
            Gravity::Mercurius => {
                ctx.apply_mercurius(n_active, soft2);
            }
            Gravity::Trace => {
                ctx.apply_trace();
            }
        }
    }
}
