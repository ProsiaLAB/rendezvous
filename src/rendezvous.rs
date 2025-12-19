use std::ops::ControlFlow;

use prosia_extensions::types::Vec3;

use crate::{
    boundary::Boundary,
    collision::{Collision, CollisionResolver},
    gravity::{Gravity, GravityContext, IgnoreGravityTerms},
    integrator::{ForceSplitIntegrator, Integrator, StepContext, SyncContext, Synchronizable},
    mercurius::MercuriusMode,
    ode::OdeState,
    particle::{Particle, TestParticleType},
    trace::TraceMode,
};

#[cfg(feature = "quadrupole")]
use crate::tree::Quadrupole;

type StepDecision = ControlFlow<ExitStatus, ()>;

#[cfg(feature = "quadrupole")]
pub type TreeParam = Quadrupole;

#[cfg(not(feature = "quadrupole"))]
pub type TreeParam = ();

pub type Tree = crate::tree::Tree<TreeParam>;

pub struct Simulation {
    /// Current simulation time. Initialized to zero by default.
    pub t: f64,
    /// Gravitational constant. Initialized to 1.0 by default.
    pub g: f64,
    /// Gravitational softening. Initialized to 0.0 by default.
    pub softening: f64,
    /// Time step. Initialized to 0.001 by default.
    pub dt: f64,
    /// Last successful time step.
    pub dt_last_done: f64,
    /// Number of steps done.
    pub steps_done: usize,

    /// Number of *real* particles.
    /// Corresponds to `(r->N - r->N_var)` in the C version.
    pub n_real: usize,
    /// Number of variational particles, if any.
    pub n_var: Option<usize>,

    pub run_state: RunState,

    pub root_size: f64,
    pub root_x: usize,
    pub root_y: usize,
    pub root_z: usize,
    pub box_size: Vec3,
    pub n_root: usize,
    pub box_size_max: f64,

    pub n_ghost_x: usize,
    pub n_ghost_y: usize,
    pub n_ghost_z: usize,

    pub gravity: Gravity,
    pub collision: Collision,
    pub collision_resolve: Option<CollisionResolver>,
    pub boundary: Boundary,

    pub n_active: usize,
    pub exit_max_distance: Option<f64>,
    pub exit_min_distance: Option<f64>,
    pub exact_finish_time: bool,
    /// This variable determines if the gravity form the
    /// central object is included in the gravity calculation.
    ///
    /// In general, the integrators will set this variable
    /// automatically and nothing needs to be changed by the user.
    pub ignore_gravity_terms: IgnoreGravityTerms,
    pub test_particle_type: TestParticleType,

    pub particles: Vec<Particle>,
    pub odes: Vec<OdeState>,
    pub gravity_cs: Vec<Vec3>,

    pub integrator: Integrator,
    /// Executed at each timestep once.
    /// Use this to do extra output/work during a simulation.
    pub heartbeat: Option<fn(&mut Self)>,

    pub additional_forces: Option<fn(&mut Self)>,
    pub pre_time_step_modifications: Option<fn(&mut Self)>,
    pub post_time_step_modifications: Option<fn(&mut Self)>,

    pub tree: Option<Tree>,
    pub tree_needs_update: bool,
    pub opening_angle: f64,
}

impl Simulation {
    pub fn init(integrator: Integrator) -> Self {
        Simulation {
            t: 0.0,
            g: 1.0,
            softening: 0.0,
            dt: 0.01,
            dt_last_done: 0.0,
            steps_done: 0,

            root_size: -1.0,
            root_x: 1,
            root_y: 1,
            root_z: 1,
            box_size: Vec3::new(-1.0, -1.0, -1.0),
            n_root: 1,
            box_size_max: -1.0,

            n_ghost_x: 0,
            n_ghost_y: 0,
            n_ghost_z: 0,

            exit_max_distance: None,
            exit_min_distance: None,
            exact_finish_time: true,
            ignore_gravity_terms: IgnoreGravityTerms::IgnoreAll,
            run_state: RunState::Running,

            heartbeat: None,
            n_real: 0,
            n_var: None,

            gravity: Gravity::Basic,
            collision: Collision::None,
            collision_resolve: None,
            boundary: Boundary::None,

            n_active: usize::MAX,
            test_particle_type: TestParticleType::Massless,

            particles: Vec::new(),
            odes: Vec::new(),
            gravity_cs: Vec::new(),

            integrator,
            additional_forces: None,
            pre_time_step_modifications: None,
            post_time_step_modifications: None,

            tree: None,
            tree_needs_update: false,
            opening_angle: 0.25,
        }
    }

    pub fn configure_box(&mut self, root_size: f64, x: usize, y: usize, z: usize) {
        self.root_size = root_size;
        self.root_x = x;
        self.root_y = y;
        self.root_z = z;
        self.box_size = Vec3::new(
            root_size * x as f64,
            root_size * y as f64,
            root_size * z as f64,
        );

        self.n_root = x * y * z;
        self.box_size_max = self.box_size.x.max(self.box_size.y).max(self.box_size.z);
        if self.root_x == 0 || self.root_y == 0 || self.root_z == 0 {
            panic!("Number of root cells in each dimension must be positive");
        }
    }

    pub fn add(&mut self, p: Particle) {
        if !self.is_particle_in_box(&p) {
            if self.box_size.x == 0.0 && self.box_size.y == 0.0 && self.box_size.z == 0.0 {
                eprintln!(
                    "ERROR: Adding particle outside of box when box size is zero. \
                     Did you forget to configure the box?"
                );
            } else {
                // Particle has left the box
                eprintln!(
                    "ERROR: Adding particle outside of box. Particle position: ({}, {}, {}), \
                     Box size: ({}, {}, {})",
                    p.x, p.y, p.z, self.box_size.x, self.box_size.y, self.box_size.z
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
            if p.x.abs() > self.box_size.x / 2.0
                || p.y.abs() > self.box_size.y / 2.0
                || p.z.abs() > self.box_size.z / 2.0
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
                    let g = self.g;
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
                if p.x > self.box_size.x / 2.0 {
                    return false;
                }
                if p.x < -self.box_size.x / 2.0 {
                    return false;
                }
                if p.y > self.box_size.y / 2.0 {
                    return false;
                }
                if p.y < -self.box_size.y / 2.0 {
                    return false;
                }
                if p.z > self.box_size.z / 2.0 {
                    return false;
                }
                if p.z < -self.box_size.z / 2.0 {
                    return false;
                }

                true
            }
            Boundary::None => true,
        }
    }

    pub fn pre_time_step(&mut self, _dt: f64) {
        unimplemented!()
    }

    pub fn post_time_step(&mut self, _dt: f64) {
        unimplemented!()
    }

    pub fn coefficient_of_restitution(&self, _p1: &Particle, _p2: &Particle) -> f64 {
        unimplemented!()
    }

    pub fn run_heartbeat(&mut self) -> StepDecision {
        if let Some(hb) = self.heartbeat {
            hb(self);
        };

        if let Some(max_distance) = self.exit_max_distance {
            let max2 = max_distance * max_distance;
            for p in self.particles.iter().take(self.n_real) {
                let r2 = p.x * p.x + p.y * p.y + p.z * p.z;
                if r2 > max2 {
                    return ControlFlow::Break(ExitStatus::Escape);
                }
            }
        }

        if let Some(min_distance) = self.exit_min_distance {
            let min2 = min_distance * min_distance;
            for i in 0..self.n_real {
                let p = &self.particles[i];
                for j in 0..i {
                    let q = &self.particles[j];
                    let x = p.x - q.x;
                    let y = p.y - q.y;
                    let z = p.z - q.z;
                    let r2 = x * x + y * y + z * z;
                    if r2 < min2 {
                        return ControlFlow::Break(ExitStatus::Encounter);
                    }
                }
            }
        }

        ControlFlow::Continue(())
    }

    pub fn sigint_received(&self) -> bool {
        todo!()
    }

    pub fn check_exit(&mut self, tmax: f64, last_full_dt: &mut f64) -> StepDecision {
        if let RunState::SingleStep { remaining } = self.run_state {
            if remaining > 0 {
                self.run_state = RunState::SingleStep {
                    remaining: remaining - 1,
                };
            } else {
                self.run_state = RunState::Paused;
            }
        }

        // TODO: Pause handling
        // while matches!(self.run_state, RunState::Paused | RunState::Screenshot) {
        //     thread::sleep(Duration::from_millis(1));

        //     // TODO: The SIGINT stuff
        //     if self.sigint_received() {
        //         return ControlFlow::Break(ExitStatus::SigInt);
        //     }

        //     return ControlFlow::Continue(());
        // }

        let dt_sign = 1.0f64.copysign(self.dt);

        // TODO: Display pending messages

        if tmax != f64::INFINITY {
            if self.exact_finish_time {
                if (self.t + self.dt) * dt_sign >= tmax * dt_sign {
                    if self.t == tmax {
                        return ControlFlow::Break(ExitStatus::Success);
                    } else if matches!(self.run_state, RunState::LastStep) {
                        let mut tscale = 1e-12 * tmax.abs();
                        if tscale < f64::MIN_POSITIVE {
                            tscale = 1e-12;
                        }
                        if (self.t - tmax).abs() < tscale {
                            return ControlFlow::Break(ExitStatus::Success);
                        } else {
                            self.synchronize();
                            self.dt = tmax - self.t;
                        }
                    } else {
                        self.run_state = RunState::LastStep;
                        self.synchronize();
                        if self.dt_last_done != 0.0 {
                            *last_full_dt = self.dt_last_done;
                        }
                        self.dt = tmax - self.t;
                    }
                } else if matches!(self.run_state, RunState::LastStep) {
                    self.run_state = RunState::Running;
                }
            } else if (self.t * dt_sign) >= (tmax * dt_sign) {
                return ControlFlow::Break(ExitStatus::Success);
            }
        }

        if self.particles.is_empty() {
            if self.odes.is_empty() {
                eprintln!("WARNING: No particles or ODE systems in simulation.");
                return ControlFlow::Break(ExitStatus::NoParticles);
            } else if matches!(self.integrator, Integrator::Gbs(_)) {
                eprintln!(
                    "WARNING: No particles in simulation with GBS integrator. \
                     Use GBS integrator to integrate user-defined ODEs \
                     without any particles present."
                );
                return ControlFlow::Break(ExitStatus::NoParticles);
            }
        }

        ControlFlow::Continue(())
    }

    pub fn check_boundary(&mut self) {
        todo!()
    }

    pub fn update_tree(&mut self) {
        todo!()
    }

    pub fn update_tree_gravity_data(&mut self) {
        todo!()
    }

    pub fn update_accelerations(&mut self) {
        if !matches!(self.integrator, Integrator::Mercurius(_))
            && matches!(self.gravity, Gravity::Mercurius)
        {
            eprintln!(
                "WARNING: You are using Mercurius gravity with a non-Mercurius integrator. \
                 This will likely lead to unexpected behavior. \
                 Switching to basic gravity."
            );
            self.gravity = Gravity::Basic;
        }

        if matches!(self.gravity, Gravity::Compensated)
            && self.gravity_cs.len() < self.particles.len()
        {
            self.gravity_cs
                .resize(self.particles.len(), Vec3::new(0.0, 0.0, 0.0));
        }

        let mut ctx = GravityContext {
            particles: &mut self.particles,
            n_real: self.n_real,
            n_active: self.n_active,
            n_root: self.n_root,
            g: self.g,
            t: self.t,
            integrator: &self.integrator,
            boundary: &self.boundary,
            n_ghost_x: self.n_ghost_x,
            n_ghost_y: self.n_ghost_y,
            n_ghost_z: self.n_ghost_z,
            box_size: &self.box_size,
            ignore_gravity_terms: &self.ignore_gravity_terms,
            softening: self.softening,
            test_particle_type: &self.test_particle_type,
            gravity_cs: &mut self.gravity_cs,
            tree: self.tree.as_ref(),
            opening_angle: self.opening_angle,
        };

        self.gravity.apply(&mut ctx);
    }

    pub fn update_accelerations_for_variational_particles(&mut self) {
        todo!()
    }

    pub fn rescale_variational_particles(&mut self) {
        todo!()
    }

    pub fn search_collisions(&mut self) {
        todo!()
    }

    pub fn get_sync_context(&mut self) -> SyncContext<'_> {
        SyncContext {
            particles: &mut self.particles,
        }
    }

    pub fn synchronize(&mut self) {
        let ctx = SyncContext {
            particles: &mut self.particles,
            // TODO: add other fields as needed
        };

        match &mut self.integrator {
            Integrator::WHFast(w) => w.synchronize(ctx),
            Integrator::Saba(s) => s.synchronize(ctx),
            Integrator::Mercurius(m) => m.synchronize(ctx),
            Integrator::Janus(j) => j.synchronize(ctx),
            Integrator::Eos(e) => e.synchronize(ctx),

            // No-ops
            Integrator::Ias15(_)
            | Integrator::LeapFrog(_)
            | Integrator::Sei(_)
            | Integrator::Gbs(_)
            | Integrator::Trace(_)
            | Integrator::None => {}
        }
    }

    pub fn integrate(&mut self, t_end: f64) -> ExitStatus {
        if t_end != self.t {
            let dt_sign = (t_end - self.t).signum();
            self.dt = self.dt.copysign(dt_sign);
        }

        let mut last_full_dt = self.dt;
        self.dt_last_done = 0.0;

        // TODO: add the warning stuff here

        if !matches!(self.run_state, RunState::Paused | RunState::Screenshot) {
            self.run_state = RunState::Running;
        }

        match self.run_heartbeat() {
            ControlFlow::Continue(()) => {}
            ControlFlow::Break(exit) => return exit,
        }

        let exit_status = loop {
            // Check time-based exit / exact-finish logic
            match self.check_exit(t_end, &mut last_full_dt) {
                ControlFlow::Continue(()) => {}
                ControlFlow::Break(exit) => break exit,
            }

            // TODO: Archive heartbeat (optional)
            // if self.simulationarchive_filename.is_some() {
            //     self.simulationarchive_heartbeat();
            // }

            self.step();

            match self.run_heartbeat() {
                ControlFlow::Continue(()) => {}
                ControlFlow::Break(exit) => break exit,
            }

            if self.sigint_received() {
                break ExitStatus::SigInt;
            }

            // TODO: Throttling (optional) usleep stuff
        };

        self.synchronize();

        if self.exact_finish_time {
            self.dt = last_full_dt;
        }

        // TODO: final archive stuff
        // if self.simulationarchive_filename.is_some() {
        //     self.simulationarchive_heartbeat();
        // }

        exit_status
    }

    pub fn step(&mut self) {
        if let Some(pre_mod) = self.pre_time_step_modifications {
            self.synchronize();
            pre_mod(self);
            match &mut self.integrator {
                Integrator::WHFast(w) => {
                    w.recalculate_coordinates_this_time_step = true;
                }
                Integrator::Mercurius(m) => {
                    m.recalculate_coordinates_this_time_step = true;
                }
                _ => {}
            }
        }

        self.pre_force_step();

        if self.tree_needs_update
            || matches!(self.gravity, Gravity::Tree)
            || matches!(self.collision, Collision::Tree | Collision::LineTree)
        {
            self.check_boundary();
            self.update_tree();
        }

        if self.tree.is_some() && matches!(self.gravity, Gravity::Tree) {
            self.update_tree_gravity_data();
        }

        self.update_accelerations();

        if self.n_var.is_some() {
            self.update_accelerations_for_variational_particles();
        }

        if let Some(addf) = self.additional_forces {
            addf(self);
        }

        self.post_force_step();

        if let Some(post_mod) = self.post_time_step_modifications {
            self.synchronize();
            post_mod(self);
            match &mut self.integrator {
                Integrator::WHFast(w) => {
                    w.recalculate_coordinates_this_time_step = true;
                }
                Integrator::Mercurius(m) => {
                    m.recalculate_coordinates_this_time_step = true;
                }
                _ => {}
            }
        }

        if self.n_var.is_some() {
            self.rescale_variational_particles();
        }

        self.check_boundary();

        if self.tree_needs_update {
            self.update_tree();
        }

        self.search_collisions();

        self.steps_done += 1;
    }

    pub fn pre_force_step(&mut self) {
        let mut ctx = StepContext {
            particles: &mut self.particles,
            t: &mut self.t,
            dt: self.dt,
            dt_last_done: &mut self.dt_last_done,
            ignore_gravity_terms: &mut self.ignore_gravity_terms,
        };

        self.integrator.pre_force(&mut ctx);
    }

    pub fn post_force_step(&mut self) {
        let mut ctx = StepContext {
            particles: &mut self.particles,
            t: &mut self.t,
            dt: self.dt,
            dt_last_done: &mut self.dt_last_done,
            ignore_gravity_terms: &mut self.ignore_gravity_terms,
        };

        self.integrator.post_force(&mut ctx);
    }
}

pub enum RunState {
    Running,
    Paused,
    Screenshot,
    ScreenshotReady,
    LastStep,
    SingleStep { remaining: u8 },
}

#[derive(Copy, Clone)]
pub enum ExitStatus {
    Success,
    GenericError,
    NoParticles,
    Encounter,
    Escape,
    User,
    SigInt,
    Collision,
}
