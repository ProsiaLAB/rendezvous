use prosia_extensions::types::Vec3;

use crate::integrator::Integrator;

pub struct BoundaryContext<'a> {
    pub t: f64,
    pub box_size: &'a Vec3,
    pub integrator: &'a Integrator,
}

pub enum Boundary {
    None,
    Open,
    Periodic,
    Shear,
}

impl Boundary {
    pub fn get_ghost_box(
        &self,
        ctx: &BoundaryContext<'_>,
        i: isize,
        j: isize,
        k: isize,
    ) -> GhostBox {
        match self {
            Boundary::Open => {
                let x = ctx.box_size.x * (i as f64);
                let y = ctx.box_size.y * (j as f64);
                let z = ctx.box_size.z * (k as f64);
                GhostBox {
                    position: Vec3::new(x, y, z),
                    velocity: Vec3::zero(),
                }
            }
            Boundary::Shear => {
                let omega = match ctx.integrator {
                    Integrator::Sei(i) => i.omega,
                    _ => {
                        panic!("Shear boundary condition requires Sei integrator");
                    }
                };
                let vx = 0.0;
                let vy = -1.5 * (i as f64) * omega * ctx.box_size.x;
                let vz = 0.0;

                let shift = if i == 0 {
                    -(vy * ctx.t) % ctx.box_size.y
                } else if i > 0 {
                    let rem = -(vy * ctx.t - ctx.box_size.y / 2.0) % ctx.box_size.y;
                    rem - ctx.box_size.y / 2.0
                } else {
                    let rem = -(vy * ctx.t + ctx.box_size.y / 2.0) % ctx.box_size.y;
                    rem + ctx.box_size.y / 2.0
                };

                let x = ctx.box_size.x * (i as f64);
                let y = ctx.box_size.y * (j as f64) - shift;
                let z = ctx.box_size.z * (k as f64);
                GhostBox {
                    position: Vec3::new(x, y, z),
                    velocity: Vec3::new(vx, vy, vz),
                }
            }
            Boundary::Periodic => {
                let x = ctx.box_size.x * (i as f64);
                let y = ctx.box_size.y * (j as f64);
                let z = ctx.box_size.z * (k as f64);
                GhostBox {
                    position: Vec3::new(x, y, z),
                    velocity: Vec3::zero(),
                }
            }
            Boundary::None => GhostBox {
                position: Vec3::zero(),
                velocity: Vec3::zero(),
            },
        }
    }
}

pub struct GhostBox {
    pub position: Vec3,
    pub velocity: Vec3,
}
