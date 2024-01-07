use crate::{errors::RayError, tuple::Tuple};

type Result<T> = std::result::Result<T, RayError>;

/// Used to represent a ray.
#[derive(Debug, Copy, Clone)]
pub struct Ray {
    // We make these private and implement getters, since we have invariants to uphold.
    origin: Tuple,    // is a point
    direction: Tuple, // is a vector
}

impl Ray {
    /// Returns a new ray.
    ///
    /// # Returns
    /// Returns Ok(Ray) if origin is a point and direction is a vector.
    /// Returns Err(RayError::NotAPoint) if origin is not a point.
    /// Returns Err(RayError:NotAVector) if direction is not a vector.
    pub fn new(origin: Tuple, direction: Tuple) -> Result<Self> {
        if !origin.is_point() {
            return Err(RayError::NotAPoint);
        }
        if !direction.is_vector() {
            return Err(RayError::NotAVector);
        }

        Ok(Self { origin, direction })
    }

    /// Returns the origin of the ray.
    pub fn origin(self) -> Tuple {
        self.origin
    }

    /// Returns the direction of the ray.
    pub fn direction(self) -> Tuple {
        self.direction
    }

    /// Returns the position of the ray after `t` time units.
    pub fn position(self, t: f64) -> Tuple {
        self.origin + self.direction * t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let origin = Tuple::point(1.0, 2.0, 3.0);
        let direction = Tuple::vector(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction).unwrap();

        assert_eq!(r.origin, origin);
        assert_eq!(r.direction, direction);
    }

    #[test]
    fn origin() {
        let origin = Tuple::point(1.0, 2.0, 3.0);
        let direction = Tuple::vector(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction).unwrap();

        assert_eq!(r.origin(), r.origin);
    }

    #[test]
    fn direction() {
        let origin = Tuple::point(1.0, 2.0, 3.0);
        let direction = Tuple::vector(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction).unwrap();

        assert_eq!(r.direction(), r.direction);
    }

    #[test]
    fn position() {
        let r = Ray::new(Tuple::point(2.0, 3.0, 4.0), Tuple::vector(1.0, 0.0, 0.0)).unwrap();

        assert_eq!(r.position(0.0), Tuple::point(2.0, 3.0, 4.0));
        assert_eq!(r.position(1.0), Tuple::point(3.0, 3.0, 4.0));
        assert_eq!(r.position(-1.0), Tuple::point(1.0, 3.0, 4.0));
        assert_eq!(r.position(2.5), Tuple::point(4.5, 3.0, 4.0));
    }
}
