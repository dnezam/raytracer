/// An error type for matrix operations.
#[derive(Debug, PartialEq)]
pub enum MatrixError {
    /// Indicates that an invertible matrix was expected, but not passed.
    NotInvertible,
}

/// An error type for canvas operations.
#[derive(Debug, PartialEq)]
pub enum CanvasError {
    /// Indicates that at least one of the coordinates is out-of-bounds.
    OutOfBounds,
}

/// An error type for tuple operations.
#[derive(Debug, PartialEq)]
pub enum TupleError {
    /// Indicates that a vector was expected, but not passed.
    NotAVector,
}
