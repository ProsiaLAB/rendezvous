use std::{
    iter::Sum,
    ops::{Add, AddAssign, Mul, Sub, SubAssign},
};

pub enum IntegratorError {
    KeplerStepFailed,
}

const G: f64 = 2.9591220828559093e-4; // AU^3 / (day^2 * solar mass) 

#[derive(Copy, Clone, Debug)]
pub struct Vector {
    data: [f64; 3],
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vector { data: [x, y, z] }
    }

    pub fn zero() -> Self {
        Vector {
            data: [0.0, 0.0, 0.0],
        }
    }

    pub fn norm(&self) -> f64 {
        (self.data[0] * self.data[0] + self.data[1] * self.data[1] + self.data[2] * self.data[2])
            .sqrt()
    }

    pub fn dot(&self, other: &Vector) -> f64 {
        self.data[0] * other.data[0] + self.data[1] * other.data[1] + self.data[2] * other.data[2]
    }
}

// Vector * scalar
impl Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, rhs: f64) -> Vector {
        Vector {
            data: [self.data[0] * rhs, self.data[1] * rhs, self.data[2] * rhs],
        }
    }
}

// scalar * Vector (optional but often useful)
impl Mul<Vector> for f64 {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Vector {
        rhs * self
    }
}

// Vector + Vector
impl Add for Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Vector {
        Vector {
            data: [
                self.data[0] + rhs.data[0],
                self.data[1] + rhs.data[1],
                self.data[2] + rhs.data[2],
            ],
        }
    }
}

impl AddAssign for Vector {
    fn add_assign(&mut self, rhs: Vector) {
        self.data[0] += rhs.data[0];
        self.data[1] += rhs.data[1];
        self.data[2] += rhs.data[2];
    }
}

impl Sub for Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Vector {
        Vector {
            data: [
                self.data[0] - rhs.data[0],
                self.data[1] - rhs.data[1],
                self.data[2] - rhs.data[2],
            ],
        }
    }
}

impl SubAssign for Vector {
    fn sub_assign(&mut self, rhs: Vector) {
        self.data[0] -= rhs.data[0];
        self.data[1] -= rhs.data[1];
        self.data[2] -= rhs.data[2];
    }
}

// Enable iterator .sum()
impl Sum for Vector {
    fn sum<I: Iterator<Item = Vector>>(iter: I) -> Vector {
        iter.fold(Vector::zero(), |a, b| a + b)
    }
}

pub struct Bodies {
    pub q: Vec<Vector>,
    pub v: Vec<Vector>,
}

impl Bodies {
    pub fn momentum(&self, masses: &[f64]) -> Vector {
        masses
            .iter()
            .zip(self.v.iter())
            .map(|(&m, &vel)| vel * m)
            .sum()
    }

    pub fn map_sun(&mut self, h: f64, masses: &[f64]) {
        let mut sp = Vector::zero();

        for (&m, &v) in masses.iter().skip(1).zip(self.v.iter().skip(1)) {
            sp += m * v;
        }

        let fac = h / masses[0] * sp;

        for i in 1..masses.len() {
            self.q[i] += fac;
        }
    }

    pub fn adjust_sun(&mut self, masses: &[f64], p0: Vector) {
        // adjustSun
    }
}

impl Bodies {
    pub fn inter_pair(&mut self, h: f64, masses: &[f64], j: usize, k: usize) {
        // Pair interaction
        let q_vec = self.q[j] - self.q[k];
        let q_mod = q_vec.norm();
        let fac = G / q_mod.powi(3) * q_vec;
        self.v[j] -= h * masses[k] * fac;
        self.v[k] += h * masses[j] * fac;
    }

    pub fn kepler_pair(
        &mut self,
        h: f64,
        masses: &[f64],
        i: usize,
        j: usize,
    ) -> Result<(), IntegratorError> {
        let gm = G * masses[0];
        let a_i = &mut self.q[i];
        let b_i = &mut self.v[i];

        let mut kpl_i = KeplerSystem::new(a_i, b_i, gm, h);
        kpl_i.step(0)?;

        let a_j = &mut self.q[j];
        let b_j = &mut self.v[j];
        let mut kpl_j = KeplerSystem::new(a_j, b_j, gm, h);
        kpl_j.step(0)?;

        Ok(())
    }
}

pub struct SystemParams {
    pub masses: Vec<f64>,
    pub r_levels: Vec<f64>,
    pub h_levels: Vec<f64>,
    pub steps: Vec<usize>,
}

/// Symba: Symplectic Massive Body Algorithm
pub struct Symba<'a> {
    pub params: &'a SystemParams,
    pub level_constraint: usize,
}

impl<'a> Symba<'a> {
    pub fn step(&self, bodies: &mut Bodies) -> Result<(), IntegratorError> {
        let n = self.params.masses.len();
        let r_levels = self.params.r_levels.len();

        let p0 = bodies.momentum(&self.params.masses);

        bodies.map_sun(self.params.h_levels[0] / 2.0, &self.params.masses);

        for level in 0..r_levels {
            if self.level_constraint != level + 1 {
                continue;
            }

            let stepsn = self.params.steps[level];
            let hlevn = self.params.h_levels[level];

            for j in 1..n {
                for k in (j + 1)..n {
                    for _ in 0..stepsn {
                        bodies.inter_pair(hlevn / 2.0, &self.params.masses, j, k);
                        bodies.kepler_pair(hlevn, &self.params.masses, j, k)?;
                        bodies.inter_pair(hlevn / 2.0, &self.params.masses, j, k);
                    }
                }
            }
        }

        bodies.map_sun(self.params.h_levels[0] / 2.0, &self.params.masses);

        bodies.adjust_sun(&self.params.masses, p0);

        Ok(())
    }
}

struct KeplerSystem<'a> {
    pub v: &'a mut Vector,
    pub vd: &'a mut Vector,
    pub kc: f64,
    pub dt: f64,
    pub r0: f64,
    pub eta: f64,
    pub beta: f64,
    pub zeta: f64,
    pub b: f64,
}

impl<'a> KeplerSystem<'a> {
    fn new(v: &'a mut Vector, vd: &'a mut Vector, kc: f64, dt: f64) -> Self {
        let v2 = vd.norm().powi(2);

        let r0 = v.norm();
        let eta = v.dot(vd);
        let beta = 2.0 * kc / r0 - v2;
        let zeta = kc - beta * r0;
        let b = beta.abs().sqrt();
        Self {
            v,
            vd,
            kc,
            dt,
            r0,
            eta,
            beta,
            zeta,
            b,
        }
    }

    fn initial_guess(&self) -> f64 {
        if self.zeta != 0.0 {
            let a = 3.0 * self.eta / self.zeta;
            let b = 6.0 * self.r0 / self.zeta;
            let c = -6.0 * self.dt / self.zeta;
            cubic(a, b, c)
        } else if self.eta != 0.0 {
            let r_eta = self.r0 / self.eta;
            let disc = r_eta * r_eta + 8.0 * self.dt / self.eta;
            if disc >= 0.0 {
                -r_eta + disc.sqrt()
            } else {
                self.dt / self.r0
            }
        } else {
            self.dt / self.r0
        }
    }

    fn step(&mut self, depth: u32) -> Result<(), IntegratorError> {
        if depth > 30 {
            panic!("Kepler step did not converge");
        }

        self.step_internal().or_else(|_| {
            self.dt /= 4.0;
            self.step(depth + 1)
        })
    }

    fn step_internal(&mut self) -> Result<(), IntegratorError> {
        let (g1, g2, r, ca, bsa) = if self.beta < 0.0 {
            let x0 = self.initial_guess();
            let mut x = x0;

            let result = self.solve_universal_hyperbolic_newton(&mut x).or_else(|_| {
                x = x0;
                self.solve_universal_hyperbolic_laguerre(&mut x)
            });

            match result {
                Ok((s2, c2)) => {
                    let a = self.kc / (-self.beta);
                    let g1 = 2.0 * s2 * c2 / self.b;
                    let c = 2.0 * s2 * s2;
                    let g2 = c / (-self.beta);
                    let ca = c * a;
                    let r = self.r0 + self.eta * g1 + self.zeta * g2;
                    let bsa = (a / r) * (self.b / self.r0) * 2.0 * s2 * c2;
                    (g1, g2, r, ca, bsa)
                }
                Err(_) => return Err(IntegratorError::KeplerStepFailed),
            }
        } else if self.beta > 0.0 {
            let mut x0 = self.dt / self.r0;
            let ff = self.zeta * x0 * x0 * x0 + 3.0 * self.eta * x0 * x0;
            let fp = 3.0 * self.zeta * x0 * x0 + 6.0 * self.eta * x0 + 6.0 * self.r0;
            x0 -= ff / fp;

            let mut x = x0;
            let result = self.solve_universal_hyperbolic_newton(&mut x).or_else(|_| {
                x = x0;
                self.solve_universal_hyperbolic_laguerre(&mut x)
            });

            match result {
                Ok((s2, c2)) => {
                    let a = self.kc / self.beta;
                    let g1 = 2.0 * s2 * c2 / self.b;
                    let c = 2.0 * s2 * s2;
                    let g2 = c / self.beta;
                    let ca = c * a;
                    let r = self.r0 + self.eta * g1 + self.zeta * g2;
                    let bsa = (a / r) * (self.b / self.r0) * 2.0 * s2 * c2;
                    (g1, g2, r, ca, bsa)
                }
                Err(_) => return Err(IntegratorError::KeplerStepFailed),
            }
        } else {
            let mut x = self.dt / self.r0;
            let result = self.solve_universal_parabolic(&mut x);
            match result {
                Ok((_, _)) => {
                    let g1 = x;
                    let g2 = x * x / 2.0;
                    let ca = self.kc * g2;
                    let r = self.r0 + self.eta * g1 + self.zeta * g2;
                    let bsa = (self.kc / x) / (r * self.r0);
                    (g1, g2, r, ca, bsa)
                }
                Err(_) => return Err(IntegratorError::KeplerStepFailed),
            }
        };

        let fhat = -(ca / self.r0);
        let g = self.eta * g2 + self.r0 * g1;
        let fdot = -bsa;
        let gdothat = -(ca / r);

        let s1 = self.v.data[0] + fhat * self.v.data[0] + g * self.vd.data[0];
        let s2 = self.v.data[1] + fhat * self.v.data[1] + g * self.vd.data[1];
        let s3 = self.v.data[2] + fhat * self.v.data[2] + g * self.vd.data[2];

        let s4 = self.vd.data[0] + fdot * self.v.data[0] + gdothat * self.vd.data[0];
        let s5 = self.vd.data[1] + fdot * self.v.data[1] + gdothat * self.vd.data[1];
        let s6 = self.vd.data[2] + fdot * self.v.data[2] + gdothat * self.vd.data[2];

        self.v.data[0] = s1;
        self.v.data[1] = s2;
        self.v.data[2] = s3;

        self.vd.data[0] = s4;
        self.vd.data[1] = s5;
        self.vd.data[2] = s6;

        Ok(())
    }

    fn solve_universal_parabolic(&self, x: &mut f64) -> Result<(f64, f64), ()> {
        let a = 3.0 * self.eta / self.zeta;
        let b = 6.0 * self.r0 / self.zeta;
        let c = -6.0 * self.dt / self.zeta;
        *x = cubic(a, b, c);
        Ok((0.0, 1.0))
    }

    fn solve_universal_hyperbolic_newton(&self, x: &mut f64) -> Result<(f64, f64), ()> {
        let mut xnew = *x;
        let mut count = 0;
        let mut flag = false;

        while !flag {
            let x0 = xnew;
            let arg = self.b * x0 / 2.0;
            if arg.abs() > 200.0 {
                return Err(());
            }
            let s2 = arg.sinh();
            let c2 = arg.cosh();
            let g1 = 2.0 * s2 * c2 / self.b;
            let g2 = 2.0 * s2 * s2 / (-self.beta);
            let g3 = -(x0 - g1) / (-self.beta);
            let g = self.eta * g1 + self.zeta * g2;
            xnew = (x0 * g - self.eta * g2 - self.zeta * g3 + self.dt) / (self.r0 + g);
            count += 1;
            if count > 10 {
                return Err(());
            }
            if (x0 - xnew) <= 1e-9 * xnew.abs() {
                flag = true;
            }
        }

        *x = xnew;
        let arg = self.b * xnew / 2.0;
        let s2 = arg.sinh();
        let c2 = arg.cosh();
        Ok((s2, c2))
    }

    fn solve_universal_hyperbolic_laguerre(&self, x: &mut f64) -> Result<(f64, f64), ()> {
        let mut xnew = *x;
        let mut count = 0;
        let mut flag = false;

        while !flag {
            let x0 = xnew;
            let arg = self.b * x0 / 2.0;
            if arg.abs() > 50.0 {
                return Err(());
            }
            let s2 = arg.sinh();
            let c2 = arg.cosh();
            let g1 = 2.0 * s2 * c2 / self.b;
            let g2 = 2.0 * s2 * s2 / (-self.beta);
            let g3 = -(x0 - g1) / (-self.beta);
            let f = self.r0 * x0 + self.eta * g2 + self.zeta * g3 - self.dt;
            let fp = self.r0 + self.eta * g1 + self.zeta * g2;
            let g0 = 1.0 + (-self.beta * g2);
            let fpp = self.eta * g0 + self.zeta * g1;
            let denom = fp + (16.0 * fp * fp - 20.0 * f * fpp).abs().sqrt();
            if denom == 0.0 {
                return Err(());
            }
            let dx = -5.0 * f / denom;
            xnew = x0 + dx;
            count += 1;
            if count > 20 {
                return Err(());
            }
            if dx.abs() <= 1e-9 * xnew.abs() {
                flag = true;
            }
        }
        *x = xnew;
        let arg = self.b * xnew / 2.0;
        if arg.abs() > 200.0 {
            return Err(());
        }
        let s2 = arg.sinh();
        let c2 = arg.cosh();
        let g1 = 2.0 * s2 * c2 / self.b;
        let g2 = 2.0 * s2 * s2 / (-self.beta);
        let g3 = -(*x - g1) / (-self.beta);
        let g = self.eta * g1 + self.zeta * g2;
        xnew = (*x * g - self.eta * g2 - self.zeta * g3 + self.dt) / (self.r0 + g);

        *x = xnew;
        let arg = self.b * xnew / 2.0;
        let s2 = arg.sinh();
        let c2 = arg.cosh();
        Ok((s2, c2))
    }

    #[allow(dead_code)]
    fn solve_universal_hyperbolic_bisection(&self, x: &mut f64) -> Result<(f64, f64), ()> {
        todo!()
    }
}

fn cubic(a: f64, b: f64, c: f64) -> f64 {
    let q = (a * a - 3.0 * b) / 9.0;
    let r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

    let r2 = r * r;
    let q3 = q * q * q;

    if r2 < q3 {
        0.0
    } else {
        let sgn = r.signum();
        let a = -sgn * (r.abs() + (r2 - q3).sqrt()).powf(1.0 / 3.0);
        let b = if a == 0.0 { 0.0 } else { q / a };
        (a + b) - a / 3.0
    }
}
