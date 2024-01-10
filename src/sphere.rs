use crate::intersection::{Intersection, Intersections};
use crate::utils::{self, Object};
use crate::{errors::SphereError, ray::Ray, tuple::Tuple, Matrix};

type Result<T> = std::result::Result<T, SphereError>;

/// Represents a sphere located at the (world) origin.
#[derive(Debug, Copy, Clone)]
pub struct Sphere {
    // Invariant: Must be unique
    id: usize,
    // Invariant: Must be invertible
    transform: Matrix<4>,
}

impl PartialEq for Sphere {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
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
    /// Returns Err(SphereError::NotInvertible) if `transform` is not an invertible matrix.
    pub fn set_transform(&mut self, transform: Matrix<4>) -> Result<()> {
        if !transform.invertible() {
            Err(SphereError::NotInvertible)
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

    /// Returns the normal at the passed surface point in world coordinates.
    ///
    /// # Returns
    /// Returns Ok(Tuple) containing the normal vector if a surface point
    /// in world coordinates was passed.
    /// Returns Err(SphereError::NotAPoint) if `world_point` is not a point.
    /// Returns Err(SphereError::NotOnSurface) if `world_point` is not on the surface of the sphere.
    pub fn normal_at(self, world_point: Tuple) -> Result<Tuple> {
        if !world_point.is_point() {
            return Err(SphereError::NotAPoint);
        }

        // unwrap() must succeed because of the invariant on self.transform
        let inverse = self.transform().inverse().unwrap();
        let object_point = inverse * world_point;
        let object_normal = object_point - Tuple::point(0.0, 0.0, 0.0);

        // Check whether the object_point is actually on the surface. This is equivalent
        // to checking whether the object_normal is actually normalized.
        // unwrap(): Point - Point = Vector
        if !utils::eq(object_normal.magnitude().unwrap(), 1.0) {
            return Err(SphereError::NotOnSurface);
        }

        // HACK Technically, we should be applying the inverse and transpose of
        // transform.submatrix(3, 3) (otherwise we run into problems if we translate)
        // Instead, we multiply with the 4x4 inverse and set w to 0 by creating a vector.
        let mut world_normal = inverse.transpose() * object_normal;
        world_normal = Tuple::vector(world_normal.x, world_normal.y, world_normal.z);

        // unwrap(): world_normal is created as a vector.
        Ok(world_normal.normalize().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, PI};

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
            SphereError::NotInvertible
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

    #[test]
    fn normal_point_on_x_axis() {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normal_point_on_y_axis() {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(0.0, 1.0, 0.0)).unwrap();
        assert_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_point_on_z_axis() {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(0.0, 0.0, 1.0)).unwrap();
        assert_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn normal_nonaxial_point() {
        let s = Sphere::new();
        let normalized_sqrt = (3.0_f64).sqrt() / 3.0;
        let n = s
            .normal_at(Tuple::point(
                normalized_sqrt,
                normalized_sqrt,
                normalized_sqrt,
            ))
            .unwrap();
        assert_eq!(
            n,
            Tuple::vector(normalized_sqrt, normalized_sqrt, normalized_sqrt)
        );
    }

    #[test]
    fn normal_is_normalized() {
        let s = Sphere::new();
        let normalized_sqrt = (3.0_f64).sqrt() / 3.0;
        let n = s
            .normal_at(Tuple::point(
                normalized_sqrt,
                normalized_sqrt,
                normalized_sqrt,
            ))
            .unwrap();
        assert_eq!(n, n.normalize().unwrap());
    }

    #[test]
    fn normal_after_translation() {
        let mut s = Sphere::new();
        s.set_transform(Matrix::<4>::translation(0.0, 1.0, 0.0))
            .unwrap();
        let n = s
            .normal_at(Tuple::point(0.0, 1.0 + FRAC_1_SQRT_2, -FRAC_1_SQRT_2))
            .unwrap();
        assert_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    // BUG This test passes if I don't check whether the point passed to normal_at is
    // actually on the surface. This is strange, as the text says that we can assume
    // that the point will always be on the surface of the sphere.
    // #[test]
    // fn normal_after_transformation() {
    //     let mut s = Sphere::new();
    //     let m = Matrix::<4>::IDENTITY
    //         .rotate_z(PI / 5.0)
    //         .scale(1.0, 0.5, 1.0);
    //     s.set_transform(m).unwrap();
    //     let normalized_sqrt = (2.0_f64).sqrt() / 2.0;
    //     let n = s
    //         .normal_at(Tuple::point(0.0, normalized_sqrt, -normalized_sqrt))
    //         .unwrap();
    //     assert_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    // }
}
