use crate::intersection::{Intersection, Intersections};
use crate::utils::{self, Object};
use crate::{ray::Ray, tuple::Tuple};

/// Represents a sphere located at the (world) origin.
#[derive(Debug, Copy, Clone)]
pub struct Sphere {
    id: usize,
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
        }
    }

    /// Returns the t values where the ray with the sphere.
    pub fn intersect(self, ray: Ray) -> Intersections {
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
}
