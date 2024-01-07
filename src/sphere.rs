use crate::{ray::Ray, tuple::Tuple, utils};

/// Represents a sphere located at the (world) origin.
pub struct Sphere {
    id: usize,
}

impl Sphere {
    /// Returns a new Sphere.
    pub fn new() -> Self {
        Self {
            id: utils::get_id(),
        }
    }

    /// Returns the t values where the ray with the sphere.
    ///
    /// # Returns
    /// Returns Some((f64, f64)) containing the t values if the ray intersects the
    /// sphere in 1 or 2 points. The t values are sorted in increasing order. If the
    /// sphere is intersected in one point only, the same point is returned twice.
    /// Returns None if the ray does not intersect the sphere.
    pub fn intersect(self, ray: Ray) -> Option<(f64, f64)> {
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
            return None;
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

        Some((t1, t2))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_at_two_points() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let (xs0, xs1) = sphere.intersect(r).unwrap();

        assert_eq!(xs0, 4.0);
        assert_eq!(xs1, 6.0);
    }

    #[test]
    fn intersect_tangent() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let (xs0, xs1) = sphere.intersect(r).unwrap();

        assert_eq!(xs0, 5.0);
        assert_eq!(xs1, 5.0);
    }

    #[test]
    fn intersect_miss() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        assert!(sphere.intersect(r).is_none());
    }

    #[test]
    fn intersect_ray_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let (xs0, xs1) = sphere.intersect(r).unwrap();

        assert_eq!(xs0, -1.0);
        assert_eq!(xs1, 1.0);
    }

    #[test]
    fn intersect_sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0)).unwrap();
        let sphere = Sphere::new();
        let (xs0, xs1) = sphere.intersect(r).unwrap();

        assert_eq!(xs0, -6.0);
        assert_eq!(xs1, -4.0);
    }
}
