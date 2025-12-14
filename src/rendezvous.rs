use crate::{
    collision::CollisionResolver, integrator::Integrator, mercurius::MercuriusMode,
    particle::Particle, trace::TraceMode,
};

#[allow(non_snake_case)]
pub struct Simulation {
    pub t: f64,
    pub G: f64,
    pub softening: f64,
    pub dt: f64,
    pub dt_last_done: f64,

    pub root_size: f64,
    pub root_x: usize,
    pub root_y: usize,
    pub root_z: usize,
    pub box_size: [f64; 3],
    pub n_root: usize,
    pub box_size_max: f64,

    pub gravity: Gravity,
    pub collision: Collision,
    pub collision_resolve: Option<CollisionResolver>,
    pub boundary: Boundary,

    pub n_active: usize,

    pub particles: Vec<Particle>,

    pub integrator: Integrator,
}

impl Simulation {
    pub fn init(integrator: Integrator) -> Self {
        Simulation {
            t: 0.0,
            G: 1.0,
            softening: 0.0,
            dt: 0.01,
            dt_last_done: 0.0,

            root_size: -1.0,
            root_x: 1,
            root_y: 1,
            root_z: 1,
            box_size: [-1.0, -1.0, -1.0],
            n_root: 1,
            box_size_max: -1.0,

            gravity: Gravity::Basic,
            collision: Collision::None,
            collision_resolve: None,
            boundary: Boundary::None,

            n_active: usize::MAX,

            particles: Vec::new(),

            integrator,
        }
    }

    pub fn configure_box(&mut self, root_size: f64, x: usize, y: usize, z: usize) {
        self.root_size = root_size;
        self.root_x = x;
        self.root_y = y;
        self.root_z = z;
        self.box_size = [
            root_size * x as f64,
            root_size * y as f64,
            root_size * z as f64,
        ];

        self.n_root = x * y * z;
        self.box_size_max = self.box_size.iter().cloned().fold(f64::NAN, f64::max);
        if self.root_x == 0 || self.root_y == 0 || self.root_z == 0 {
            panic!("Number of root cells in each dimension must be positive");
        }
    }

    pub fn add(&mut self, p: Particle) {
        if !self.is_particle_in_box(&p) {
            if self.box_size[0] == 0.0 && self.box_size[1] == 0.0 && self.box_size[2] == 0.0 {
                eprintln!(
                    "ERROR: Adding particle outside of box when box size is zero. \
                     Did you forget to configure the box?"
                );
            } else {
                // Particle has left the box
                eprintln!(
                    "ERROR: Adding particle outside of box. Particle position: ({}, {}, {}), \
                     Box size: ({}, {}, {})",
                    p.x, p.y, p.z, self.box_size[0], self.box_size[1], self.box_size[2]
                );
            }
            return;
        }

        if matches!(self.gravity, Gravity::Tree)
            || matches!(self.collision, Collision::Tree | Collision::LineTree)
        {
            if self.root_size == -1.0 {
                eprintln!("ERROR: Box not configured before adding particles.");
                return;
            }
            if p.x.abs() > self.box_size[0] / 2.0
                || p.y.abs() > self.box_size[1] / 2.0
                || p.z.abs() > self.box_size[2] / 2.0
            {
                eprintln!(
                    "ERROR: Particle position outside of box when adding particle. \
                     Did you forget to configure the box?"
                );
            }
        }

        self.particles.push(p);

        match &mut self.integrator {
            Integrator::Mercurius(rim) => match rim.mode {
                MercuriusMode::LongRange => {
                    // WHFast mode
                    rim.recalculate_r_crit_this_time_step = true;
                    rim.recalculate_coordinates_this_time_step = true;
                }
                MercuriusMode::CloseEncounter => {
                    // IAS15 mode
                    rim.ias15.reset();
                    if rim.dcrit.len() < self.particles.len() {
                        rim.dcrit.resize(self.particles.len(), 0.0);
                    }
                    let new_index = self.particles.len() - 1;
                    let p0 = &self.particles[0];
                    let pi = &self.particles[new_index];
                    let g = self.G;
                    let dt = self.dt;
                    rim.set_dcrit(p0, pi, g, dt, new_index);
                    if rim.particles_backup.len() < self.particles.len() {
                        rim.particles_backup
                            .resize(self.particles.len(), Particle::default());
                        rim.particles_backup_additional_forces
                            .resize(self.particles.len(), Particle::default());
                    }
                    rim.encounter_map.push(new_index);
                    rim.n_encounter += 1;
                    if self.n_active == usize::MAX {
                        rim.n_encounter_active += 1;
                    }
                }
            },
            Integrator::Trace(trace) => {
                if matches!(trace.mode, TraceMode::Kepler | TraceMode::Full) {
                    // GBS part
                    let new_n = self.particles.len();
                    let old_n = new_n - 1;

                    trace.particles_backup.resize(new_n, Particle::default());
                    trace
                        .particles_backup_kepler
                        .resize(new_n, Particle::default());
                    trace
                        .particles_backup_additional_forces
                        .resize(new_n, Particle::default());
                    trace.resize_current_ks(old_n, new_n);
                    for &p in trace.encounter_map.iter().skip(1) {
                        trace.current_ks[p * new_n + old_n] = 1;
                    }

                    trace.encounter_map.push(old_n);
                    trace.n_encounter += 1;
                    if self.n_active == usize::MAX {
                        trace.n_encounter_active += 1;
                    }
                }
            }
            _ => {}
        }
    }

    pub fn is_particle_in_box(&self, p: &Particle) -> bool {
        match self.boundary {
            Boundary::Open | Boundary::Shear | Boundary::Periodic => {
                if p.x > self.box_size[0] / 2.0 {
                    return false;
                }
                if p.x < -self.box_size[0] / 2.0 {
                    return false;
                }
                if p.y > self.box_size[1] / 2.0 {
                    return false;
                }
                if p.y < -self.box_size[1] / 2.0 {
                    return false;
                }
                if p.z > self.box_size[2] / 2.0 {
                    return false;
                }
                if p.z < -self.box_size[2] / 2.0 {
                    return false;
                }

                true
            }
            Boundary::None => true,
        }
    }

    pub fn additional_forces(&mut self) {
        unimplemented!()
    }

    pub fn pre_time_step(&mut self, dt: f64) {
        unimplemented!()
    }

    pub fn post_time_step(&mut self, dt: f64) {
        unimplemented!()
    }

    pub fn coefficient_of_restitution(&self, p1: &Particle, p2: &Particle) -> f64 {
        unimplemented!()
    }

    pub fn integrate(&mut self, t_end: f64) {
        unimplemented!()
    }
}

pub enum Gravity {
    None,
    Basic,
    Compensated,
    Tree,
    Mercurius,
    Jacobi,
    Trace,
}

pub enum Collision {
    None,
    Direct,
    Tree,
    Line,
    LineTree,
}

pub enum Boundary {
    None,
    Open,
    Periodic,
    Shear,
}

pub trait Heartbeat {
    fn heartbeat(&mut self);
}

impl Heartbeat for Simulation {
    fn heartbeat(&mut self) {
        // Default: do nothing
    }
}
