//! Raytracer based on the book "The Ray Tracer Challenge" by Jamis Buck.
#![deny(missing_docs)]

mod canvas;
mod color;
mod errors;
mod matrix;
mod ray;
mod sphere;
mod tuple;
mod utils;

pub use canvas::Canvas;
pub use color::Color;
pub use matrix::Matrix;
pub use tuple::Tuple;
