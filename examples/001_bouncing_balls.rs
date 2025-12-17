use rendezvous::{
    collision::{Collision, CollisionResolver},
    gravity::Gravity,
    integrator::Integrator,
    leapfrog::LeapFrog,
    particle::Particle,
    rendezvous::Simulation,
};

fn main() {
    let mut sim = Simulation::init(Integrator::LeapFrog(LeapFrog));

    sim.dt = 1e-2;
    sim.gravity = Gravity::Basic;
    sim.collision = Collision::Direct;
    sim.collision_resolve = CollisionResolver::HardSphere.into();

    sim.configure_box(3.0, 1, 1, 1);

    // Initial conditions
    let p = Particle {
        x: 1.0,
        y: 1.0,
        z: 1.0,
        m: 1.0,
        r: 0.1,
        ..Default::default()
    };
    sim.add(p);

    let p = Particle {
        x: -1.0,
        y: -1.0,
        z: -1.0,
        m: 1.0,
        r: 0.1,
        ..Default::default()
    };
    sim.add(p);

    sim.integrate(f64::INFINITY);
}
