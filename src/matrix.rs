use crate::Tuple;
use std::ops::{Index, IndexMut, Mul};

use crate::utils;

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

/// Functions applicable to all matrices
impl<const N: usize> Matrix<N> {
    /// Create a new matrix.
    pub fn new(elements: [[f64; N]; N]) -> Self {
        Matrix { elements }
    }

    /// Create a new matrix initialized with zeros.
    pub fn default() -> Self {
        Self::new([[0.0; N]; N])
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
        let mut submatrix = Matrix::<2>::default();

        // Assert that row and column are in-bounds
        assert!(row < 3 && column < 3);

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
    // Cannot be implemented more generically, as we cannot return Matrix<N-1>
    // (no const operations allowed)
    fn submatrix(self, row: usize, column: usize) -> Matrix<3> {
        let mut submatrix = Matrix::<3>::default();

        // Assert that row and column are in-bounds
        assert!(row < 4 && column < 4);

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

    /// Returns whether a function is invertible.
    // TODO Instead of panicking in inverse, remove invertible and make inverse return a result type
    pub fn invertible(self) -> bool {
        self.determinant() != 0.0
    }

    /// Returns the inverse.
    pub fn inverse(self) -> Self {
        assert!(self.invertible());

        let mut inverse = Self::default();

        for i in 0..4 {
            for j in 0..4 {
                let cofactor = self.cofactor(i, j);
                inverse[[j, i]] = cofactor / self.determinant();
            }
        }

        inverse
    }
}

#[cfg(test)]
mod tests {
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
        let inverse = matrix.inverse();

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
        let inverse = matrix.inverse();

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
        let inverse = matrix.inverse();

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
        assert_eq!(product * matrix2.inverse(), matrix1);
    }
}
