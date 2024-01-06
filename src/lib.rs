//! Raytracer based on the book "The Ray Tracer Challenge" by Jamis Buck.
#![deny(missing_docs)]

// TODO Replace as with (try_)into

mod canvas;
mod color;
mod errors;
mod matrix;
mod tuple;
mod utils;

pub use canvas::Canvas;
pub use color::Color;
pub use errors::CanvasError;
pub use errors::MatrixError;
pub use errors::TupleError;
pub use matrix::Matrix;
pub use tuple::Tuple;
