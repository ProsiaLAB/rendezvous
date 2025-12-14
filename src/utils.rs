pub mod vector {
    use std::{
        iter::Sum,
        ops::{Add, AddAssign, Mul, Sub, SubAssign},
    };

    #[derive(Copy, Clone, Debug)]
    pub struct Vector {
        data: [f64; 3],
    }

    impl Vector {
        pub fn new(x: f64, y: f64, z: f64) -> Self {
            Vector { data: [x, y, z] }
        }

        #[inline]
        pub fn x(&self) -> f64 {
            self.data[0]
        }

        #[inline]
        pub fn y(&self) -> f64 {
            self.data[1]
        }

        #[inline]
        pub fn z(&self) -> f64 {
            self.data[2]
        }

        pub fn set(&mut self, x: f64, y: f64, z: f64) {
            self.data[0] = x;
            self.data[1] = y;
            self.data[2] = z;
        }

        pub fn zero() -> Self {
            Vector {
                data: [0.0, 0.0, 0.0],
            }
        }

        pub fn norm(&self) -> f64 {
            (self.data[0] * self.data[0]
                + self.data[1] * self.data[1]
                + self.data[2] * self.data[2])
                .sqrt()
        }

        pub fn dot(&self, other: &Vector) -> f64 {
            self.data[0] * other.data[0]
                + self.data[1] * other.data[1]
                + self.data[2] * other.data[2]
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
}
