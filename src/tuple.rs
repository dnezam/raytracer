use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::utils;

#[derive(Debug, Copy, Clone)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl PartialEq for Tuple {
    fn eq(&self, other: &Self) -> bool {
        utils::eq(self.x, other.x)
            && utils::eq(self.y, other.y)
            && utils::eq(self.z, other.z)
            && utils::eq(self.w, other.w)
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }
}

impl Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

// Tuple * f64
impl Mul<f64> for Tuple {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

// f64 * Tuple
impl Mul<Tuple> for f64 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        Self::Output {
            x: rhs.x * self,
            y: rhs.y * self,
            z: rhs.z * self,
            w: rhs.w * self,
        }
    }
}

// Tuple / f64
impl Div<f64> for Tuple {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Tuple { x, y, z, w }
    }

    pub fn point(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 1.0 }
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 0.0 }
    }

    pub fn is_point(self) -> bool {
        utils::eq(self.w, 1.0)
    }

    pub fn is_vector(self) -> bool {
        utils::eq(self.w, 0.0)
    }

    pub fn magnitude(self) -> f64 {
        let xx = self.x.powi(2);
        let yy = self.y.powi(2);
        let zz = self.z.powi(2);
        let ww = self.w.powi(2);
        (xx + yy + zz + ww).sqrt()
    }

    pub fn normalize(self) -> Self {
        self / self.magnitude()
    }

    pub fn dot(self, other: Self) -> f64 {
        // Only implement dot product for vectors: Here we deviate from the book,
        // because the reasons of generalizing to all dimensions don't convince me.
        // Although it is possible that I'll consider the w coordinate in the future.
        assert!(self.is_vector());
        assert!(other.is_vector());

        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(self, other: Self) -> Self {
        // Only implement cross-product for vectors, not tuples: four-dimensional
        // cross product is more complicated and not needed
        assert!(self.is_vector());
        assert!(other.is_vector());

        let x = self.y * other.z - self.z * other.y;
        let y = self.z * other.x - self.x * other.z;
        let z = self.x * other.y - self.y * other.x;
        Tuple::vector(x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_point() {
        let a = Tuple::new(4.3, -4.2, 3.1, 1.0);
        assert!(utils::eq(a.x, 4.3));
        assert!(utils::eq(a.y, -4.2));
        assert!(utils::eq(a.z, 3.1));
        assert!(utils::eq(a.w, 1.0));
        assert!(a.is_point());
        assert!(!a.is_vector());
    }

    #[test]
    fn tuple_vector() {
        let a = Tuple::new(4.3, -4.2, 3.1, 0.0);
        assert!(utils::eq(a.x, 4.3));
        assert!(utils::eq(a.y, -4.2));
        assert!(utils::eq(a.z, 3.1));
        assert!(utils::eq(a.w, 0.0));
        assert!(!a.is_point());
        assert!(a.is_vector());
    }

    #[test]
    fn new_point() {
        let p = Tuple::point(4.0, -4.0, 3.0);
        assert_eq!(p, Tuple::new(4.0, -4.0, 3.0, 1.0));
    }

    #[test]
    fn new_vector() {
        let v = Tuple::vector(4.0, -4.0, 3.0);
        assert_eq!(v, Tuple::new(4.0, -4.0, 3.0, 0.0));
    }

    #[test]
    fn add_tuples() {
        let a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);
        assert_eq!(a1 + a2, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }

    #[test]
    fn subtract_points() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::point(5.0, 6.0, 7.0);
        assert_eq!(p1 - p2, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn subtract_vector_from_point() {
        let p = Tuple::point(3.0, 2.0, 1.0);
        let v = Tuple::vector(5.0, 6.0, 7.0);
        assert_eq!(p - v, Tuple::point(-2.0, -4.0, -6.0));
    }

    #[test]
    fn substract_vectors() {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);
        assert_eq!(v1 - v2, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn subtract_vector_from_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);
        assert_eq!(zero - v, Tuple::vector(-1.0, 2.0, -3.0));
    }

    #[test]
    fn negate_tuple() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(-a, Tuple::new(-1.0, 2.0, -3.0, 4.0));
    }

    #[test]
    fn multiply_tuple_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 3.5, Tuple::new(3.5, -7.0, 10.5, -14.0));
    }

    #[test]
    fn multiply_tuple_by_fraction() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 0.5, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn multiply_commutative() {
        let a = Tuple::new(1.0, 0.2, -4.66, 253.0);
        assert_eq!(-1.3 * a, a * -1.3);
    }

    #[test]
    fn divide_tuple_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a / 2.0, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn magnitude_1() {
        let v = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_2() {
        let v = Tuple::vector(0.0, 1.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_3() {
        let v = Tuple::vector(0.0, 0.0, 1.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_4() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        assert_eq!(v.magnitude(), f64::sqrt(14.0));
    }

    #[test]
    fn magnitude_5() {
        let v = Tuple::vector(-1.0, -2.0, -3.0);
        assert_eq!(v.magnitude(), f64::sqrt(14.0));
    }

    #[test]
    fn normalize_1() {
        let v = Tuple::vector(4.0, 0.0, 0.0);
        assert_eq!(v.normalize(), Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normalize_2() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        // vector(1/√14, 2/√14, 3/√14)​
        assert_eq!(v.normalize(), Tuple::vector(0.26726, 0.53452, 0.80178));
    }

    #[test]
    fn normalize_magnitude() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let norm = v.normalize();
        assert_eq!(norm.magnitude(), 1.0);
    }

    #[test]
    fn dot_vectors() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b), 20.0);
    }

    #[test]
    fn cross_vectors() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b), Tuple::vector(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(a), Tuple::vector(1.0, -2.0, 1.0));
    }
}
