use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::{errors::TupleError, utils};

type Result<T> = std::result::Result<T, TupleError>;

/// Used to encapsulate position, direction and distance.
///
/// We will use a left-handed coordinate system. This is reflected by the z-axis
/// pointing *away* from us.
#[derive(Debug, Copy, Clone)]
pub struct Tuple {
    /// The x-axis points to the right.
    pub x: f64,
    /// The y-axis points up.
    pub y: f64,
    /// The z-axis points away from us.
    pub z: f64,
    /// Determines whether the tuple is a point or vector.
    ///
    /// If w is 1, we have a point. If w is 0, we have a vector.
    pub w: f64,
}

impl PartialEq for Tuple {
    /// Implements approximate equality for tuples.
    fn eq(&self, other: &Self) -> bool {
        utils::eq(self.x, other.x)
            && utils::eq(self.y, other.y)
            && utils::eq(self.z, other.z)
            && utils::eq(self.w, other.w)
    }
}

impl Add for Tuple {
    type Output = Self;

    /// Implements the addition of two tuples.
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

    /// Implements the subtraction of two tuples.
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

    /// Implements the negation of a tuple.
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

/// Scalar multiplication of the form Tuple * f64.
impl Mul<f64> for Tuple {
    type Output = Self;

    /// Implements scalar multiplication of the form: Tuple * f64.
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

/// Scalar multiplication of the form f64 * Tuple.
impl Mul<Tuple> for f64 {
    type Output = Tuple;

    /// Implements scalar multiplication of the form: Tuple * f64.
    fn mul(self, rhs: Tuple) -> Self::Output {
        Self::Output {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
            w: self * rhs.w,
        }
    }
}

/// Scalar division of the form Tuple / f64.
impl Div<f64> for Tuple {
    type Output = Self;

    /// Implements scalar division of the form Tuple / f64.
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
    /// Create a new tuple.
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Tuple { x, y, z, w }
    }

    /// Create a new point.
    pub fn point(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 1.0 }
    }

    /// Create a new vector.
    pub fn vector(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 0.0 }
    }

    /// Returns true if and only if this tuple is a point.
    pub fn is_point(self) -> bool {
        utils::eq(self.w, 1.0)
    }

    /// Returns true if and only if this tuple is a vector.
    pub fn is_vector(self) -> bool {
        utils::eq(self.w, 0.0)
    }

    /// Returns the magnitude (length).
    ///
    /// # Returns
    /// Returns `Ok(f64)` containing the magnitude if the tuple is a vector.
    /// Returns `Err(TupleError::NotAVector)` if the tuple is not a vector.
    pub fn magnitude(self) -> Result<f64> {
        if !self.is_vector() {
            return Err(TupleError::NotAVector);
        }

        // We now know that self is a vector, hence w must be 0 and we don't need
        // to consider it in the calculation.
        let xx = self.x.powi(2);
        let yy = self.y.powi(2);
        let zz = self.z.powi(2);
        Ok((xx + yy + zz).sqrt())
    }

    /// Returns the normalized tuple.
    ///
    /// # Returns
    /// Returns `Ok(Tuple)` containing the normalized vector if the tuple is a vector.
    /// Returns `Err(TupleError::NotAVector)` if the tuple is not a vector.
    pub fn normalize(self) -> Result<Self> {
        Ok(self / self.magnitude()?)
    }

    /// Returns the dot product.
    ///
    /// Intuitively, the smaller the result, the larger the angle between the vectors.
    ///
    /// # Returns
    /// Returns `Ok(Tuple)` containing the dot product if both tuples are vectors.
    /// Returns `Err(TupleError::NotAVector)` if at least one tuple is not a vector.
    pub fn dot(self, other: Self) -> Result<f64> {
        if !self.is_vector() || !other.is_vector() {
            return Err(TupleError::NotAVector);
        }

        Ok(self.x * other.x + self.y * other.y + self.z * other.z)
    }

    /// Returns the cross product.
    ///
    /// The cross produt is a vector that is perpendicular to both of the original vectors.
    ///
    /// # Returns
    /// Returns `Ok(Tuple)` containing the cross product if both tuples are vectors.
    /// Returns `Err(TupleError::NotAVector)` if at least one tuple is not a vector.
    pub fn cross(self, other: Self) -> Result<Self> {
        // Only implement cross-product for vectors: four-dimensional
        // cross product is more complicated and not needed
        if !self.is_vector() || !other.is_vector() {
            return Err(TupleError::NotAVector);
        }

        let x = self.y * other.z - self.z * other.y;
        let y = self.z * other.x - self.x * other.z;
        let z = self.x * other.y - self.y * other.x;
        Ok(Tuple::vector(x, y, z))
    }

    /// Returns the vector resulting from reflecting the current vector around `normal`
    ///
    /// # Returns
    /// Returns `Ok(Tuple)` containing the reflected vector if `self` is a vector and
    /// `normal` a normalized vector.
    /// Returns `Err(TupleError::NotAVector)` if at least one tuple is not a vector.
    /// Returns `Err(TupleError::NotNormalized)` if `normal` isn't normalized.
    pub fn reflect(self, normal: Self) -> Result<Self> {
        if !self.is_vector() || !normal.is_vector() {
            Err(TupleError::NotAVector)
        // unwrap(): Both self and normal must be vectors here.
        } else if !utils::eq(normal.magnitude().unwrap(), 1.0) {
            Err(TupleError::NotNormalized)
        } else {
            Ok(self - normal * 2.0 * self.dot(normal).unwrap())
        }
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
    fn magnitude_point() {
        let p = Tuple::point(1.0, 0.0, 0.0);
        assert_eq!(p.magnitude().unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn magnitude_1() {
        let v = Tuple::vector(1.0, 0.0, 0.0);
        assert!(utils::eq(v.magnitude().unwrap(), 1.0));
    }

    #[test]
    fn magnitude_2() {
        let v = Tuple::vector(0.0, 1.0, 0.0);
        assert!(utils::eq(v.magnitude().unwrap(), 1.0));
    }

    #[test]
    fn magnitude_3() {
        let v = Tuple::vector(0.0, 0.0, 1.0);
        assert!(utils::eq(v.magnitude().unwrap(), 1.0));
    }

    #[test]
    fn magnitude_4() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        assert!(utils::eq(v.magnitude().unwrap(), f64::sqrt(14.0)));
    }

    #[test]
    fn magnitude_5() {
        let v = Tuple::vector(-1.0, -2.0, -3.0);
        assert!(utils::eq(v.magnitude().unwrap(), f64::sqrt(14.0)));
    }

    #[test]
    fn normalize_point() {
        let p = Tuple::point(4.0, 0.0, 0.0);
        assert_eq!(p.normalize().unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn normalize_1() {
        let v = Tuple::vector(4.0, 0.0, 0.0);
        assert_eq!(v.normalize().unwrap(), Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normalize_2() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        // vector(1/√14, 2/√14, 3/√14)​
        assert_eq!(
            v.normalize().unwrap(),
            Tuple::vector(0.26726, 0.53452, 0.80178)
        );
    }

    #[test]
    fn normalize_magnitude() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let norm = v.normalize().unwrap();
        assert!(utils::eq(norm.magnitude().unwrap(), 1.0));
    }

    #[test]
    fn dot_vector_point() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn dot_point_vector() {
        let a = Tuple::point(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn dot_point_point() {
        let a = Tuple::point(1.0, 2.0, 3.0);
        let b = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn dot_vectors() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert!(utils::eq(a.dot(b).unwrap(), 20.0));
    }

    #[test]
    fn cross_vector_point() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn cross_point_vector() {
        let a = Tuple::point(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn cross_point_point() {
        let a = Tuple::point(1.0, 2.0, 3.0);
        let b = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b).unwrap_err(), TupleError::NotAVector);
    }

    #[test]
    fn cross_vectors() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b).unwrap(), Tuple::vector(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(a).unwrap(), Tuple::vector(1.0, -2.0, 1.0));
    }

    #[test]
    fn reflect_45() {
        let v = Tuple::vector(1.0, -1.0, 0.0);
        let n = Tuple::vector(0.0, 1.0, 0.0);
        let r = v.reflect(n).unwrap();

        assert_eq!(r, Tuple::vector(1.0, 1.0, 0.0));
    }

    #[test]
    fn reflect_slanted() {
        let v = Tuple::vector(0.0, -1.0, 0.0);
        let frc_sqrt_2_2 = (2.0_f64).sqrt() / 2.0;
        let n = Tuple::vector(frc_sqrt_2_2, frc_sqrt_2_2, 0.0);
        let r = v.reflect(n).unwrap();

        assert_eq!(r, Tuple::vector(1.0, 0.0, 0.0));
    }
}
