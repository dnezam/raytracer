//! Raytracer based on the book "The Ray Tracer Challenge" by Jamis Buck.
#![deny(missing_docs)]

// TODO Check if constants like sqrt(2)/2 are in our code and can be replaced
// using the standard math library.

mod canvas;
mod color;
mod errors;
mod intersection;
mod light;
mod material;
mod matrix;
mod ray;
mod sphere;
mod tuple;
mod utils;

pub use canvas::Canvas;
pub use color::Color;
pub use matrix::Matrix;
pub use ray::Ray;
pub use sphere::Sphere;
pub use tuple::Tuple;
