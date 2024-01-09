use crate::intersection::{Intersection, Intersections};
use crate::utils::{self, Object};
use crate::{errors::MatrixError, ray::Ray, tuple::Tuple, Matrix};
use std::result::Result;

/// Represents a sphere located at the (world) origin.
#[derive(Debug, Copy, Clone)]
pub struct Sphere {
    id: usize,
    transform: Matrix<4>,
}

impl PartialEq for Sphere {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Sphere {
    /// Returns a new Sphere.
    pub fn new() -> Self {
        Self {
            id: utils::get_id(),
            transform: Matrix::<4>::IDENTITY,
        }
    }

    /// Returns the current transformation.
    pub fn transform(self) -> Matrix<4> {
        self.transform
    }

    /// Sets the transformation matrix.
    ///
    /// # Returns
    /// Returns Ok(()) if `transform` is an invertible matrix.
    /// Returns Err(MatrixError::NotInvertible) if `transform` is not an invertible matrix.
    pub fn set_transform(&mut self, transform: Matrix<4>) -> Result<(), MatrixError> {
        if !transform.invertible() {
            Err(MatrixError::NotInvertible)
        } else {
            self.transform = transform;
            Ok(())
        }
    }

    /// Returns the t values where the ray with the sphere.
    pub fn intersect(self, ray: Ray) -> Intersections {
        let ray = ray.transform(self.transform().inverse().unwrap());

        // Vector from sphere's center to the ray origin: Sphere is located at world origin
        let sphere_to_ray = ray.origin() - Tuple::point(0.0, 0.0, 0.0);

        // ray.direction() must be a vector by the invariants in Ray
        let a = ray.direction().dot(ray.direction()).unwrap();
        // sphere_to_ray must be a vector since ray.origin() is a point (invariant in Ray)
        // and subtracting a point from a point results in a vector.
        let b = 2.0 * ray.direction().dot(sphere_to_ray).unwrap();
        let c = sphere_to_ray.dot(sphere_to_ray).unwrap() - 1.0;

        let discriminant = b.powi(2) - 4.0 * a * c;
        if discriminant < 0.0 {
            return Intersections::empty();
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

        let i1 = Intersection::new(t1, Object::Sphere(self));
        let i2 = Intersection::new(t2, Object::Sphere(self));

        // (Hopefully) t1 and t2 are both numbers, hence this unwrap should always succeed.
        Intersections::new(&[i1, i2]).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_at_two_points() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let xs = sphere.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), 4.0);
        assert_eq!(xs[1].t(), 6.0);
    }

    #[test]
    fn intersect_tangent() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let xs = sphere.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), 5.0);
        assert_eq!(xs[1].t(), 5.0);
    }

    #[test]
    fn intersect_miss() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        assert!(sphere.intersect(r).is_empty());
    }

    #[test]
    fn intersect_ray_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let xs = sphere.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), -1.0);
        assert_eq!(xs[1].t(), 1.0);
    }

    #[test]
    fn intersect_sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let xs = sphere.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), -6.0);
        assert_eq!(xs[1].t(), -4.0);
    }

    #[test]
    fn intersect_sets_intersection_object() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let xs = sphere.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].object(), Object::Sphere(sphere));
        assert_eq!(xs[1].object(), Object::Sphere(sphere));
    }

    #[test]
    fn default_transformation() {
        let sphere = Sphere::new();
        assert_eq!(sphere.transform(), Matrix::<4>::IDENTITY);
    }

    #[test]
    fn change_transformation() {
        let mut sphere = Sphere::new();
        let t = Matrix::<4>::translation(2.0, 3.0, 4.0);
        sphere.set_transform(t).unwrap();

        assert_eq!(sphere.transform(), t);
    }

    #[test]
    fn change_transformation_invalid() {
        let mut sphere = Sphere::new();
        let t = Matrix::<4>::default();
        assert_eq!(
            sphere.set_transform(t).unwrap_err(),
            MatrixError::NotInvertible
        );
    }

    #[test]
    fn intersect_scaled_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let mut s = Sphere::new();
        s.transform = Matrix::<4>::scaling(2.0, 2.0, 2.0);
        let xs = s.intersect(r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), 3.0);
        assert_eq!(xs[1].t(), 7.0);
    }

    #[test]
    fn intersect_translated_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let mut s = Sphere::new();
        s.transform = Matrix::<4>::translation(5.0, 0.0, 0.0);
        let xs = s.intersect(r);

        assert_eq!(xs.len(), 0);
    }
}
