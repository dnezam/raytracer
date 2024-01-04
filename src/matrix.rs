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
        let mut result = Self::Output::new([[0.0; 4]; 4]);
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

// We will only need to multiply tuples with 4x4 matrices.
impl Mul<Tuple> for Matrix<4> {
    type Output = Tuple;

    /// Implements matrix-tuple multiplication.
    fn mul(self, rhs: Tuple) -> Self::Output {
        // If I ever implement these operations again, I will make everything matrices.
        Self::Output {
            x: self[[0, 0]] * rhs.x
                + self[[0, 1]] * rhs.y
                + self[[0, 2]] * rhs.z
                + self[[0, 3]] * rhs.w,
            y: self[[1, 0]] * rhs.x
                + self[[1, 1]] * rhs.y
                + self[[1, 2]] * rhs.z
                + self[[1, 3]] * rhs.w,
            z: self[[2, 0]] * rhs.x
                + self[[2, 1]] * rhs.y
                + self[[2, 2]] * rhs.z
                + self[[2, 3]] * rhs.w,
            w: self[[3, 0]] * rhs.x
                + self[[3, 1]] * rhs.y
                + self[[3, 2]] * rhs.z
                + self[[3, 3]] * rhs.w,
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
        let mut transposed = Self::new([[0.0; N]; N]);

        for i in 0..N {
            for j in 0..N {
                transposed[[j, i]] = self[[i, j]];
            }
        }

        transposed
    }
}

impl Matrix<4> {
    // We only need the identity matrix for 4x4 matrices.
    const IDENTITY: Self = Self {
        elements: [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    };
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
}
