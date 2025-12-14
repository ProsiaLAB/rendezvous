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
    pub hash: u64,
}
