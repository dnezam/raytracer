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

/// An error type for ray operations.
#[derive(Debug, PartialEq)]
pub enum RayError {
    /// Indicates that a point was expected, but not passed.
    NotAPoint,
    /// Indicates that a vector was expected, but not passed.
    NotAVector,
}

/// An error type for sphere operations.
#[derive(Debug, PartialEq)]
pub enum SphereError {
    /// Indicates that an invertible matrix was expected, but not passed.
    NotInvertible,
    /// Indicates that a point was expected, but not passed.
    NotAPoint,
    /// Indicates that a point on a surface was expected, but not passed.
    NotOnSurface,
}
