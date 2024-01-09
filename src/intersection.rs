use std::ops::Index;

use crate::utils::{self, Object};

/// Represent where an object is intersected.
#[derive(Debug, Copy, Clone)]
pub struct Intersection {
    t: f64,
    object: Object,
}

impl PartialEq for Intersection {
    fn eq(&self, other: &Self) -> bool {
        utils::eq(self.t, other.t) && self.object == other.object
    }
}

impl Intersection {
    /// Create a new intersection.
    pub fn new(t: f64, object: Object) -> Self {
        Self { t, object }
    }

    /// Returns at what distance the object was intersected.
    pub fn t(self) -> f64 {
        self.t
    }

    /// Returns the object that was intersected.
    pub fn object(self) -> Object {
        self.object
    }
}

/// Represents a list of intersections.
#[derive(PartialEq, Debug)]
pub struct Intersections {
    elements: Vec<Intersection>,
}

impl Index<usize> for Intersections {
    type Output = Intersection;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl Intersections {
    /// Create a new list of intersections.
    ///
    /// # Returns
    /// Returns Some(Intersections), if none of the distances of the intersections is NaN.
    /// Returns None, if some distance of the intersection is NaN.
    pub fn new(elements: &[Intersection]) -> Option<Self> {
        let mut elements = elements.to_vec();

        if elements
            .iter()
            .any(|intersection| intersection.t().is_nan())
        {
            None
        } else {
            // elements will not contain NaN here
            elements.sort_by(|i1, i2| i1.t().partial_cmp(&i2.t()).unwrap());
            Some(Self { elements })
        }
    }

    // Create an empty list of intersections.
    pub fn empty() -> Self {
        // An empty list is always sortable.
        Self::new(&[]).unwrap()
    }

    /// Returns the number of intersections.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// Returns whether the list of intersections is empty.
    pub fn is_empty(&self) -> bool {
        self.elements.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sphere::Sphere;

    #[test]
    fn new_intersection() {
        let s = Sphere::new();
        let i = Intersection::new(3.5, Object::Sphere(s));
        assert_eq!(i.t(), 3.5);
        assert_eq!(i.object(), Object::Sphere(s));
    }

    #[test]
    fn new_intersections() {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, Object::Sphere(s));
        let i2 = Intersection::new(2.0, Object::Sphere(s));
        let xs = Intersections::new(&[i1, i2]).unwrap();
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t(), 1.0);
        assert_eq!(xs[1].t(), 2.0);
    }

    #[test]
    fn empty() {
        assert_eq!(Intersections::empty(), Intersections::new(&[]).unwrap());
    }

    #[test]
    fn is_empty() {
        assert!(Intersections::empty().is_empty());
    }

    #[test]
    fn is_not_empty() {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, Object::Sphere(s));
        let xs = Intersections::new(&[i1]).unwrap();
        assert!(!xs.is_empty());
    }
}
