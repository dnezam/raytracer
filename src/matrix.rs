use crate::{errors::MatrixError, tuple::Tuple};
use std::ops::{Index, IndexMut, Mul};

use crate::utils;

type Result<T> = std::result::Result<T, MatrixError>;

/// An N x N matrix.
#[derive(Debug, Copy, Clone)]
pub struct Matrix<const N: usize> {
    elements: [[f64; N]; N],
}

impl<const N: usize> Index<[usize; 2]> for Matrix<N> {
    type Output = f64;

    /// Implements indexing of the form matrix[[row, column]] in immutable contexts.
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.elements[index[0]][index[1]]
    }
}

impl<const N: usize> IndexMut<[usize; 2]> for Matrix<N> {
    /// Implements indexing of the form matrix[[row, column]] in mutable contexts.
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.elements[index[0]][index[1]]
    }
}

impl<const N: usize> PartialEq for Matrix<N> {
    /// Implements approximate equality for matrices.
    fn eq(&self, other: &Self) -> bool {
        let flat_self = self.elements.iter().flatten();
        let flat_other = other.elements.iter().flatten();
        flat_self
            .zip(flat_other)
            .all(|(e1, e2)| utils::eq(*e1, *e2))
    }
}

impl<const N: usize> Default for Matrix<N> {
    /// Create a new matrix initialized with zeros.
    fn default() -> Self {
        Self {
            elements: [[0.0; N]; N],
        }
    }
}

impl<const N: usize> Matrix<N> {
    /// Create a new matrix.
    pub fn new(elements: [[f64; N]; N]) -> Self {
        Matrix { elements }
    }

    /// Returns a transposed matrix.
    pub fn transpose(self) -> Self {
        let mut transposed = Self::default();

        for i in 0..N {
            for j in 0..N {
                transposed[[j, i]] = self[[i, j]];
            }
        }

        transposed
    }
}

impl Matrix<2> {
    /// Returns the determinant.
    fn determinant(self) -> f64 {
        self[[0, 0]] * self[[1, 1]] - self[[0, 1]] * self[[1, 0]]
    }
}

impl Matrix<3> {
    /// Returns a copy of the given matrix with the given row and column removed.
    fn submatrix(self, row: usize, column: usize) -> Matrix<2> {
        // Cannot be implemented more generically, as we cannot return Matrix<N-1>
        // (no const operations allowed)

        // submatrix is a private function; hence, if this assert fires we made a mistake somewhere.
        assert!(row < 3 && column < 3);

        let mut submatrix = Matrix::<2>::default();

        for (sub_i, i) in (0..3).filter(|&x| x != row).enumerate() {
            for (sub_j, j) in (0..3).filter(|&x| x != column).enumerate() {
                submatrix[[sub_i, sub_j]] = self[[i, j]]
            }
        }

        submatrix
    }

    /// Returns the (row, column) minor.
    fn minor(self, row: usize, column: usize) -> f64 {
        self.submatrix(row, column).determinant()
    }

    /// Returns the (row, column) cofactor.
    fn cofactor(self, row: usize, column: usize) -> f64 {
        let minor = self.minor(row, column);
        if (row + column) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }

    /// Returns the determinant.
    fn determinant(self) -> f64 {
        (0..3).map(|j| self[[0, j]] * self.cofactor(0, j)).sum()
    }
}

// We will only need to multiply 4x4 matrices.
impl Mul<Matrix<4>> for Matrix<4> {
    type Output = Matrix<4>;

    /// Implements matrix-matrix multiplication.
    fn mul(self, rhs: Matrix<4>) -> Self::Output {
        let mut result = Self::Output::default();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    result[[i, j]] += self[[i, k]] * rhs[[k, j]];
                }
            }
        }
        result
    }
}

// 4x4 Matrix * Tuple
impl Mul<Tuple> for Matrix<4> {
    type Output = Tuple;

    /// Implements matrix-tuple multiplication.
    fn mul(self, rhs: Tuple) -> Self::Output {
        let tuple = [rhs.x, rhs.y, rhs.z, rhs.w];
        let tuple_dot_row = |row: usize| {
            tuple
                .iter()
                .zip(self.elements[row].iter())
                .map(|(x, y)| x * y)
                .sum::<f64>()
        };
        Self::Output {
            x: tuple_dot_row(0),
            y: tuple_dot_row(1),
            z: tuple_dot_row(2),
            w: tuple_dot_row(3),
        }
    }
}

impl Matrix<4> {
    /// We only need the identity matrix for 4x4 matrices.
    pub const IDENTITY: Self = Self {
        elements: [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    };

    /// Returns a copy of the given matrix with the given row and column removed.
    fn submatrix(self, row: usize, column: usize) -> Matrix<3> {
        // Cannot be implemented more generically, as we cannot return Matrix<N-1>
        // (no const operations allowed)

        // submatrix is a private function; hence, if this assert fires we made a mistake somewhere.
        assert!(row < 4 && column < 4);

        let mut submatrix = Matrix::<3>::default();

        for (sub_i, i) in (0..4).filter(|&x| x != row).enumerate() {
            for (sub_j, j) in (0..4).filter(|&x| x != column).enumerate() {
                submatrix[[sub_i, sub_j]] = self[[i, j]]
            }
        }

        submatrix
    }

    /// Returns the (row, column) minor.
    fn minor(self, row: usize, column: usize) -> f64 {
        self.submatrix(row, column).determinant()
    }

    /// Returns the (row, column) cofactor.
    fn cofactor(self, row: usize, column: usize) -> f64 {
        let minor = self.minor(row, column);
        if (row + column) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }

    /// Returns the determinant.
    fn determinant(self) -> f64 {
        (0..4).map(|j| self[[0, j]] * self.cofactor(0, j)).sum()
    }

    /// Returns whether the matrix is invertible.
    pub fn invertible(self) -> bool {
        self.determinant() != 0.0
    }

    /// Returns the inverse.
    ///
    /// # Returns
    /// Returns `Ok(Matrix)` containing the inverse of the given matrix if it is invertible.
    /// Returns `Err(MatrixError::NotInvertible)` if the matrix is not invertible.
    pub fn inverse(self) -> Result<Self> {
        if !self.invertible() {
            return Err(MatrixError::NotInvertible);
        }

        let mut inverse = Self::default();

        for i in 0..4 {
            for j in 0..4 {
                let cofactor = self.cofactor(i, j);
                inverse[[j, i]] = cofactor / self.determinant();
            }
        }

        Ok(inverse)
    }

    /// Returns a translation matrix.
    pub fn translation(x: f64, y: f64, z: f64) -> Self {
        Self {
            elements: [
                [1.0, 0.0, 0.0, x],
                [0.0, 1.0, 0.0, y],
                [0.0, 0.0, 1.0, z],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating translation.
    pub fn translate(self, x: f64, y: f64, z: f64) -> Self {
        Self::translation(x, y, z) * self
    }

    /// Returns a scaling matrix.
    pub fn scaling(x: f64, y: f64, z: f64) -> Self {
        Self {
            elements: [
                [x, 0.0, 0.0, 0.0],
                [0.0, y, 0.0, 0.0],
                [0.0, 0.0, z, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating scaling.
    pub fn scale(self, x: f64, y: f64, z: f64) -> Self {
        Self::scaling(x, y, z) * self
    }

    /// Returns a rotation matrix around the x-axis.
    pub fn rotation_x(radians: f64) -> Self {
        let (sin, cos) = radians.sin_cos();
        Self {
            elements: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, cos, -sin, 0.0],
                [0.0, sin, cos, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating rotation around the x-axis.
    pub fn rotate_x(self, radians: f64) -> Self {
        Self::rotation_x(radians) * self
    }

    /// Returns a rotation matrix around the y-axis.
    pub fn rotation_y(radians: f64) -> Self {
        let (sin, cos) = radians.sin_cos();
        Self {
            elements: [
                [cos, 0.0, sin, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [-sin, 0.0, cos, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating rotation around the y-axis.
    pub fn rotate_y(self, radians: f64) -> Self {
        Self::rotation_y(radians) * self
    }

    /// Returns a rotation matrix around the z-axis.
    pub fn rotation_z(radians: f64) -> Self {
        let (sin, cos) = radians.sin_cos();
        Self {
            elements: [
                [cos, -sin, 0.0, 0.0],
                [sin, cos, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating rotation around the z-axis.
    pub fn rotate_z(self, radians: f64) -> Self {
        Self::rotation_z(radians) * self
    }

    /// Returns a shearing matrix.
    ///
    /// # Parameters
    /// - `xy`: Shearing of the x-axis in proportion to the y-axis.
    /// - `xz`: Shearing of the x-axis in proportion to the z-axis.
    /// - `yx`: Shearing of the y-axis in proportion to the x-axis.
    /// - `yz`: Shearing of the y-axis in proportion to the z-axis.
    /// - `zx`: Shearing of the z-axis in proportion to the x-axis.
    /// - `zy`: Shearing of the z-axis in proportion to the y-axis.
    pub fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        Self {
            elements: [
                [1.0, xy, xz, 0.0],
                [yx, 1.0, yz, 0.0],
                [zx, zy, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Returns a transformation by concatenating shearing.
    ///
    /// # Parameters
    /// - `xy`: Shearing of the x-axis in proportion to the y-axis.
    /// - `xz`: Shearing of the x-axis in proportion to the z-axis.
    /// - `yx`: Shearing of the y-axis in proportion to the x-axis.
    /// - `yz`: Shearing of the y-axis in proportion to the z-axis.
    /// - `zx`: Shearing of the z-axis in proportion to the x-axis.
    /// - `zy`: Shearing of the z-axis in proportion to the y-axis.
    pub fn shear(self, xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        Self::shearing(xy, xz, yx, yz, zx, zy) * self
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    #[test]
    fn new_4x4() {
        let matrix = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);

        assert_eq!(matrix[[0, 0]], 1.0);
        assert_eq!(matrix[[0, 3]], 4.0);
        assert_eq!(matrix[[1, 0]], 5.5);
        assert_eq!(matrix[[1, 2]], 7.5);
        assert_eq!(matrix[[2, 2]], 11.0);
        assert_eq!(matrix[[3, 0]], 13.5);
        assert_eq!(matrix[[3, 2]], 15.5);
    }

    #[test]
    fn new_2x2() {
        let matrix = Matrix::<2>::new([[-3.0, 5.0], [1.0, -2.0]]);
        assert_eq!(matrix[[0, 0]], -3.0);
        assert_eq!(matrix[[0, 1]], 5.0);
        assert_eq!(matrix[[1, 0]], 1.0);
        assert_eq!(matrix[[1, 1]], -2.0);
    }

    #[test]
    fn new_3x3() {
        let matrix = Matrix::<3>::new([[-3.0, 5.0, 0.0], [1.0, -2.0, -7.0], [0.0, 1.0, 1.0]]);
        assert_eq!(matrix[[0, 0]], -3.0);
        assert_eq!(matrix[[1, 1]], -2.0);
        assert_eq!(matrix[[2, 2]], 1.0);
    }

    #[test]
    fn equal() {
        let matrix1 = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);

        let matrix2 = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        assert_eq!(matrix1, matrix2);
    }

    #[test]
    fn not_equal() {
        let matrix1 = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);

        let matrix2 = Matrix::<4>::new([
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [8.0, 7.0, 6.0, 5.0],
            [4.0, 3.0, 2.0, 1.0],
        ]);
        assert_ne!(matrix1, matrix2);
    }

    #[test]
    fn multiply_4x4_matrix_matrix() {
        let matrix1 = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let matrix2 = Matrix::<4>::new([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);

        let expected_matrix = Matrix::<4>::new([
            [20.0, 22.0, 50.0, 48.0],
            [44.0, 54.0, 114.0, 108.0],
            [40.0, 58.0, 110.0, 102.0],
            [16.0, 26.0, 46.0, 42.0],
        ]);

        assert_eq!(matrix1 * matrix2, expected_matrix);
    }

    #[test]
    fn multiply_4x4_matrix_tuple() {
        let matrix = Matrix::<4>::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let tuple = Tuple::new(1.0, 2.0, 3.0, 1.0);

        let expected_tuple = Tuple::new(18.0, 24.0, 33.0, 1.0);

        assert_eq!(matrix * tuple, expected_tuple);
    }

    #[test]
    fn multiply_4x4_matrix_identity() {
        let matrix = Matrix::<4>::new([
            [0.0, 1.0, 2.0, 3.0],
            [1.0, 2.0, 4.0, 8.0],
            [2.0, 4.0, 8.0, 16.0],
            [4.0, 8.0, 16.0, 32.0],
        ]);

        assert_eq!(matrix * Matrix::<4>::IDENTITY, matrix)
    }

    #[test]
    fn transpose() {
        let matrix = Matrix::<4>::new([
            [0.0, 9.0, 3.0, 0.0],
            [9.0, 8.0, 0.0, 8.0],
            [1.0, 8.0, 5.0, 3.0],
            [0.0, 0.0, 5.0, 8.0],
        ]);

        let expected_matrix = Matrix::<4>::new([
            [0.0, 9.0, 1.0, 0.0],
            [9.0, 8.0, 8.0, 0.0],
            [3.0, 0.0, 5.0, 5.0],
            [0.0, 8.0, 3.0, 8.0],
        ]);

        assert_eq!(matrix.transpose(), expected_matrix);
    }

    #[test]
    fn transpose_identity() {
        assert_eq!(Matrix::<4>::IDENTITY.transpose(), Matrix::<4>::IDENTITY);
    }

    #[test]
    fn determinant_2x2() {
        let matrix = Matrix::<2>::new([[1.0, 5.0], [-3.0, 2.0]]);
        assert_eq!(matrix.determinant(), 17.0);
    }

    #[test]
    fn submatrix_3x3() {
        let matrix = Matrix::<3>::new([[1.0, 5.0, 0.0], [-3.0, 2.0, 7.0], [0.0, 6.0, -3.0]]);
        let expected_submatrix = Matrix::<2>::new([[-3.0, 2.0], [0.0, 6.0]]);

        assert_eq!(matrix.submatrix(0, 2), expected_submatrix);
    }

    #[test]
    fn submatrix_4x4() {
        let matrix = Matrix::<4>::new([
            [-6.0, 1.0, 1.0, 6.0],
            [-8.0, 5.0, 8.0, 6.0],
            [-1.0, 0.0, 8.0, 2.0],
            [-7.0, 1.0, -1.0, 1.0],
        ]);

        let expected_submatrix =
            Matrix::<3>::new([[-6.0, 1.0, 6.0], [-8.0, 8.0, 6.0], [-7.0, -1.0, 1.0]]);

        assert_eq!(matrix.submatrix(2, 1), expected_submatrix);
    }

    #[test]
    fn minor_3x3() {
        let matrix = Matrix::<3>::new([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);
        let submatrix = matrix.submatrix(1, 0);
        assert_eq!(matrix.minor(1, 0), submatrix.determinant());
    }

    #[test]
    fn cofactor_3x3() {
        let matrix = Matrix::<3>::new([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);
        assert_eq!(matrix.minor(0, 0), -12.0);
        assert_eq!(matrix.cofactor(0, 0), -12.0);
        assert_eq!(matrix.minor(1, 0), 25.0);
        assert_eq!(matrix.cofactor(1, 0), -25.0);
    }

    #[test]
    fn determinant_3x3() {
        let matrix = Matrix::<3>::new([[1.0, 2.0, 6.0], [-5.0, 8.0, -4.0], [2.0, 6.0, 4.0]]);

        assert_eq!(matrix.cofactor(0, 0), 56.0);
        assert_eq!(matrix.cofactor(0, 1), 12.0);
        assert_eq!(matrix.cofactor(0, 2), -46.0);
        assert_eq!(matrix.determinant(), -196.0);
    }

    #[test]
    fn determinant_4x4() {
        let matrix = Matrix::<4>::new([
            [-2.0, -8.0, 3.0, 5.0],
            [-3.0, 1.0, 7.0, 3.0],
            [1.0, 2.0, -9.0, 6.0],
            [-6.0, 7.0, 7.0, -9.0],
        ]);

        assert_eq!(matrix.cofactor(0, 0), 690.0);
        assert_eq!(matrix.cofactor(0, 1), 447.0);
        assert_eq!(matrix.cofactor(0, 2), 210.0);
        assert_eq!(matrix.cofactor(0, 3), 51.0);
        assert_eq!(matrix.determinant(), -4071.0);
    }

    #[test]
    fn not_invertible_4x4() {
        // Last row is a duplicate
        let matrix = Matrix::<4>::new([
            [6.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 7.0, 6.0],
            [4.0, -9.0, 3.0, -7.0],
            [4.0, -9.0, 3.0, -7.0],
        ]);

        assert_eq!(matrix.inverse().unwrap_err(), MatrixError::NotInvertible)
    }

    #[test]
    fn invertible_4x4() {
        let matrix = Matrix::<4>::new([
            [6.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 7.0, 6.0],
            [4.0, -9.0, 3.0, -7.0],
            [9.0, 1.0, 7.0, -6.0],
        ]);

        assert_eq!(matrix.determinant(), -2120.0);
        assert!(matrix.invertible());
    }

    #[test]
    fn noninvertible_4x4() {
        let matrix = Matrix::<4>::new([
            [-4.0, 2.0, -2.0, -3.0],
            [9.0, 6.0, 2.0, 6.0],
            [0.0, -5.0, 1.0, -5.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        assert_eq!(matrix.determinant(), 0.0);
        assert!(!matrix.invertible());
    }

    #[test]
    fn inverse_4x4_1() {
        let matrix = Matrix::<4>::new([
            [-5.0, 2.0, 6.0, -8.0],
            [1.0, -5.0, 1.0, 8.0],
            [7.0, 7.0, -6.0, -7.0],
            [1.0, -3.0, 7.0, 4.0],
        ]);
        let inverse = matrix.inverse().unwrap();

        let expected_inverse = Matrix::<4>::new([
            [0.21805, 0.45113, 0.24060, -0.04511],
            [-0.80827, -1.45677, -0.44361, 0.52068],
            [-0.07895, -0.22368, -0.05263, 0.19737],
            [-0.52256, -0.81391, -0.30075, 0.30639],
        ]);

        assert_eq!(matrix.determinant(), 532.0);
        assert_eq!(matrix.cofactor(2, 3), -160.0);
        assert_eq!(inverse[[3, 2]], -160.0 / 532.0);
        assert_eq!(matrix.cofactor(3, 2), 105.0);
        assert_eq!(inverse[[2, 3]], 105.0 / 532.0);
        assert_eq!(inverse, expected_inverse);
    }

    #[test]
    fn inverse_4x4_2() {
        let matrix = Matrix::<4>::new([
            [8.0, -5.0, 9.0, 2.0],
            [7.0, 5.0, 6.0, 1.0],
            [-6.0, 0.0, 9.0, 6.0],
            [-3.0, 0.0, -9.0, -4.0],
        ]);
        let inverse = matrix.inverse().unwrap();

        let expected_inverse = Matrix::<4>::new([
            [-0.15385, -0.15385, -0.28205, -0.53846],
            [-0.07692, 0.12308, 0.02564, 0.03077],
            [0.35897, 0.35897, 0.43590, 0.92308],
            [-0.69231, -0.69231, -0.76923, -1.92308],
        ]);
        assert_eq!(inverse, expected_inverse);
    }

    #[test]
    fn inverse_4x4_3() {
        let matrix = Matrix::<4>::new([
            [9.0, 3.0, 0.0, 9.0],
            [-5.0, -2.0, -6.0, -3.0],
            [-4.0, 9.0, 6.0, 4.0],
            [-7.0, 6.0, 6.0, 2.0],
        ]);
        let inverse = matrix.inverse().unwrap();

        let expected_inverse = Matrix::<4>::new([
            [-0.04074, -0.07778, 0.14444, -0.22222],
            [-0.07778, 0.03333, 0.36667, -0.33333],
            [-0.02901, -0.14630, -0.10926, 0.12963],
            [0.17778, 0.06667, -0.26667, 0.33333],
        ]);
        assert_eq!(inverse, expected_inverse);
    }

    #[test]
    fn multiply_by_inverse_4x4() {
        let matrix1 = Matrix::<4>::new([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);
        let matrix2 = Matrix::<4>::new([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);

        let product = matrix1 * matrix2;
        assert_eq!(product * matrix2.inverse().unwrap(), matrix1);
    }

    #[test]
    fn multiply_by_translation_4x4() {
        let transform = Matrix::<4>::translation(5.0, -3.0, 2.0);
        let p = Tuple::point(-3.0, 4.0, 5.0);
        assert_eq!(transform * p, Tuple::point(2.0, 1.0, 7.0));
    }

    #[test]
    fn multiply_by_inverse_translation_4x4() {
        let transform = Matrix::<4>::translation(5.0, -3.0, 2.0);
        let inv = transform.inverse().unwrap();
        let p = Tuple::point(-3.0, 4.0, 5.0);
        assert_eq!(inv * p, Tuple::point(-8.0, 7.0, 3.0));
    }

    #[test]
    fn vector_translation_invariant() {
        let transform = Matrix::<4>::translation(5.0, -3.0, 2.0);
        let v = Tuple::vector(-3.0, 4.0, 5.0);
        assert_eq!(transform * v, v);
    }

    #[test]
    fn translation_translate_equal() {
        assert_eq!(
            Matrix::<4>::translation(5.0, -3.0, 2.0),
            Matrix::<4>::IDENTITY.translate(5.0, -3.0, 2.0)
        );
    }

    #[test]
    fn scaling_point() {
        let transform = Matrix::<4>::scaling(2.0, 3.0, 4.0);
        let p = Tuple::point(-4.0, 6.0, 8.0);
        assert_eq!(transform * p, Tuple::point(-8.0, 18.0, 32.0));
    }

    #[test]
    fn scaling_vector() {
        let transform = Matrix::<4>::scaling(2.0, 3.0, 4.0);
        let p = Tuple::vector(-4.0, 6.0, 8.0);
        assert_eq!(transform * p, Tuple::vector(-8.0, 18.0, 32.0));
    }

    #[test]
    fn inverse_scaling() {
        let transform = Matrix::<4>::scaling(2.0, 3.0, 4.0);
        let inv = transform.inverse().unwrap();
        let v = Tuple::vector(-4.0, 6.0, 8.0);
        assert_eq!(inv * v, Tuple::vector(-2.0, 2.0, 2.0));
    }

    #[test]
    fn reflection_is_negative_scaling() {
        let transform = Matrix::<4>::scaling(-1.0, 1.0, 1.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(-2.0, 3.0, 4.0));
    }

    #[test]
    fn scaling_scale_equal() {
        assert_eq!(
            Matrix::<4>::scaling(2.0, 3.0, 4.0),
            Matrix::<4>::IDENTITY.scale(2.0, 3.0, 4.0)
        );
    }

    #[test]
    fn rotation_x() {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix::<4>::rotation_x(PI / 4.0);
        let full_quarter = Matrix::<4>::rotation_x(PI / 2.0);

        assert_eq!(
            half_quarter * p,
            Tuple::point(0.0, (2.0_f64).sqrt() / 2.0, (2.0_f64).sqrt() / 2.0)
        );
        assert_eq!(full_quarter * p, Tuple::point(0.0, 0.0, 1.0));
    }

    #[test]
    fn inverse_rotation_x() {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix::<4>::rotation_x(PI / 4.0);
        let inv = half_quarter.inverse().unwrap();
        assert_eq!(
            inv * p,
            Tuple::point(0.0, (2.0_f64).sqrt() / 2.0, -(2.0_f64).sqrt() / 2.0)
        )
    }

    #[test]
    fn rotation_x_rotate_x_equal() {
        assert_eq!(
            Matrix::<4>::rotation_x(PI),
            Matrix::<4>::IDENTITY.rotate_x(PI)
        );
    }

    #[test]
    fn rotation_y() {
        let p = Tuple::point(0.0, 0.0, 1.0);
        let half_quarter = Matrix::<4>::rotation_y(PI / 4.0);
        let full_quarter = Matrix::<4>::rotation_y(PI / 2.0);

        assert_eq!(
            half_quarter * p,
            Tuple::point((2.0_f64).sqrt() / 2.0, 0.0, (2.0_f64).sqrt() / 2.0)
        );
        assert_eq!(full_quarter * p, Tuple::point(1.0, 0.0, 0.0));
    }

    #[test]
    fn rotation_y_rotate_y_equal() {
        assert_eq!(
            Matrix::<4>::rotation_y(PI),
            Matrix::<4>::IDENTITY.rotate_y(PI)
        );
    }

    #[test]
    fn rotation_z() {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix::<4>::rotation_z(PI / 4.0);
        let full_quarter = Matrix::<4>::rotation_z(PI / 2.0);

        assert_eq!(
            half_quarter * p,
            Tuple::point(-((2.0_f64).sqrt() / 2.0), (2.0_f64).sqrt() / 2.0, 0.0)
        );
        assert_eq!(full_quarter * p, Tuple::point(-1.0, 0.0, 0.0));
    }

    #[test]
    fn rotation_z_rotate_z_equal() {
        assert_eq!(
            Matrix::<4>::rotation_z(PI),
            Matrix::<4>::IDENTITY.rotate_z(PI)
        );
    }

    #[test]
    fn shearing_x_proportional_y() {
        let transform = Matrix::<4>::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(5.0, 3.0, 4.0));
    }

    #[test]
    fn shearing_x_proportional_z() {
        let transform = Matrix::<4>::shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(6.0, 3.0, 4.0));
    }

    #[test]
    fn shearing_y_proportional_x() {
        let transform = Matrix::<4>::shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(2.0, 5.0, 4.0));
    }

    #[test]
    fn shearing_y_proportional_z() {
        let transform = Matrix::<4>::shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(2.0, 7.0, 4.0));
    }

    #[test]
    fn shearing_z_proportional_x() {
        let transform = Matrix::<4>::shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(2.0, 3.0, 6.0));
    }

    #[test]
    fn shearing_z_proportional_y() {
        let transform = Matrix::<4>::shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Tuple::point(2.0, 3.0, 7.0));
    }

    #[test]
    fn shearing_shear_equal() {
        assert_eq!(
            Matrix::<4>::shearing(0.0, 2.0, 1.0, 0.5, 1.0, 0.0),
            Matrix::<4>::IDENTITY.shear(0.0, 2.0, 1.0, 0.5, 1.0, 0.0)
        );
    }

    #[test]
    fn sequentially_applied_transformations() {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let rotation_x = Matrix::<4>::rotation_x(PI / 2.0);
        let scaling = Matrix::<4>::scaling(5.0, 5.0, 5.0);
        let translation = Matrix::<4>::translation(10.0, 5.0, 7.0);

        let p2 = rotation_x * p;
        assert_eq!(p2, Tuple::point(1.0, -1.0, 0.0));

        let p3 = scaling * p2;
        assert_eq!(p3, Tuple::point(5.0, -5.0, 0.0));

        let p4 = translation * p3;
        assert_eq!(p4, Tuple::point(15.0, 0.0, 7.0))
    }

    #[test]
    fn chained_transformations() {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let rotation_x = Matrix::<4>::rotation_x(PI / 2.0);
        let scaling = Matrix::<4>::scaling(5.0, 5.0, 5.0);
        let translation = Matrix::<4>::translation(10.0, 5.0, 7.0);
        let transform = translation * scaling * rotation_x;

        assert_eq!(transform * p, Tuple::point(15.0, 0.0, 7.0));
    }

    // Make sure we don't concatenate in a wrong manner (A.B() should B * A, not A * B)
    // This test is probably weak though, since I don't think it does too well of a job
    // of exercising these behaviors.
    #[test]
    fn chained_transforms_fluent() {
        let (translation_x, translation_y, translation_z) = (1.2, 0.4, -3.3);
        let (scaling_x, scaling_y, scaling_z) = (0.2, -4.2, 1.2);
        let radians_x = 13.0;
        let radians_y = -12.0;
        let radians_z = -0.5;
        let (shearing_xy, shearing_xz, shearing_yx, shearing_yz, shearing_zx, shearing_zy) =
            (0.0, 1.2, -0.9, 3.5, 2.1, 15.0);

        let translation = Matrix::<4>::translation(translation_x, translation_y, translation_z);
        let scaling = Matrix::<4>::scaling(scaling_x, scaling_y, scaling_z);
        let rotation_x = Matrix::<4>::rotation_x(radians_x);
        let rotation_y = Matrix::<4>::rotation_y(radians_y);
        let rotation_z = Matrix::<4>::rotation_z(radians_z);
        let shearing = Matrix::<4>::shearing(
            shearing_xy,
            shearing_xz,
            shearing_yx,
            shearing_yz,
            shearing_zx,
            shearing_zy,
        );

        assert_eq!(
            translation * scaling * rotation_x * rotation_y * rotation_z * shearing,
            Matrix::<4>::IDENTITY
                .shear(
                    shearing_xy,
                    shearing_xz,
                    shearing_yx,
                    shearing_yz,
                    shearing_zx,
                    shearing_zy
                )
                .rotate_z(radians_z)
                .rotate_y(radians_y)
                .rotate_x(radians_x)
                .scale(scaling_x, scaling_y, scaling_z)
                .translate(translation_x, translation_y, translation_z)
        )
    }
}
