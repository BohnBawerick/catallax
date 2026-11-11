// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: a square matrix with only the LU factorisation and
// solve that SMLib's Newton solver uses; statuses become `Fault` and the
// near-singular warning becomes a returned flag.
// Adapted from rust-brep crates/kernel/src/matrix.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, a translation of SMLib
// SmMatrix.cpp (LUDecompose, SolveWithLUMatrix, SolveLinearSystem).
//! Dense square matrices and SMLib's partially pivoted LU solve.

use crate::Fault;
use crate::scalar::EFF_ZERO;
use std::ops::{Index, IndexMut};

/// A dense row-major matrix (`SmMatrix`), zero on construction.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    rows: usize,
    columns: usize,
    data: Vec<f64>,
}

/// How an LU solve went (`SolveLinearSystem`'s success and warning).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LuQuality {
    /// Every pivot was at least `SM_EFF_ZERO` in magnitude.
    Regular,
    /// A pivot below `SM_EFF_ZERO` was replaced by `±SM_EFF_ZERO`
    /// (`SM_WARN`); the solution is still returned.
    NearSingular,
}

impl Matrix {
    /// A `rows` by `columns` zero matrix.
    pub fn new(rows: usize, columns: usize) -> Self {
        Self {
            rows,
            columns,
            data: vec![0.0; rows * columns],
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn columns(&self) -> usize {
        self.columns
    }

    fn swap_rows(&mut self, a: usize, b: usize) {
        for column in 0..self.columns {
            self.data
                .swap(a * self.columns + column, b * self.columns + column);
        }
    }

    /// `LUDecompose`: Crout's method with implicit row scaling. `indices`
    /// receives the row permutation.
    pub fn lu_decompose(&mut self, indices: &mut [usize]) -> Result<LuQuality, Fault> {
        if self.rows != self.columns || indices.len() < self.rows {
            return Err(Fault::InvalidInput);
        }
        let n = self.rows;
        let mut scale = vec![0.0; n];
        for (row, scale) in scale.iter_mut().enumerate() {
            let mut big = 0.0;
            for col in 0..n {
                let value = self[row][col].abs();
                if value > big {
                    big = value;
                }
            }
            if big < EFF_ZERO {
                return Err(Fault::Degenerate);
            }
            *scale = 1.0 / big;
        }
        let mut quality = LuQuality::Regular;
        for col in 0..n {
            for row in 0..col {
                let mut sum = self[row][col];
                for k in 0..row {
                    sum -= self[row][k] * self[k][col];
                }
                self[row][col] = sum;
            }
            let mut big = 0.0;
            let mut maximum = 0;
            for row in col..n {
                let mut sum = self[row][col];
                for k in 0..col {
                    sum -= self[row][k] * self[k][col];
                }
                self[row][col] = sum;
                let value = scale[row] * sum.abs();
                if value >= big {
                    big = value;
                    maximum = row;
                }
            }
            if col != maximum {
                self.swap_rows(maximum, col);
                scale[maximum] = scale[col];
            }
            indices[col] = maximum;
            if self[col][col].abs() < EFF_ZERO {
                self[col][col] = if self[col][col] > 0.0 {
                    EFF_ZERO
                } else {
                    -EFF_ZERO
                };
                quality = LuQuality::NearSingular;
            }
            if col + 1 != n {
                let reciprocal = 1.0 / self[col][col];
                for row in col + 1..n {
                    self[row][col] *= reciprocal;
                }
            }
        }
        Ok(quality)
    }

    /// `SolveWithLUMatrix`: forward and back substitution on an LU factor.
    pub fn solve_with_lu(
        &self,
        indices: &[usize],
        rhs: &[f64],
        solution: &mut [f64],
    ) -> Result<(), Fault> {
        let n = self.rows;
        if n != self.columns
            || rhs.len() != n
            || solution.len() != n
            || indices.len() < n
            || indices[..n].iter().any(|&i| i >= n)
        {
            return Err(Fault::InvalidInput);
        }
        solution.copy_from_slice(rhs);
        let mut first = None;
        for row in 0..n {
            let permutation = indices[row];
            let mut sum = solution[permutation];
            solution[permutation] = solution[row];
            if let Some(start) = first {
                for col in start..row {
                    sum -= self[row][col] * solution[col];
                }
            } else if sum != 0.0 {
                first = Some(row);
            }
            solution[row] = sum;
        }
        for row in (0..n).rev() {
            let mut sum = solution[row];
            for col in row + 1..n {
                sum -= self[row][col] * solution[col];
            }
            let denominator = self[row][row];
            if denominator.abs() < EFF_ZERO {
                return Err(Fault::Degenerate);
            }
            solution[row] = sum / denominator;
        }
        Ok(())
    }

    /// `SolveLinearSystem`: solves `self * x = rhs` on a copy of `self`.
    pub fn solve_linear_system(
        &self,
        rhs: &[f64],
        solution: &mut Vec<f64>,
    ) -> Result<LuQuality, Fault> {
        if self.rows != self.columns || rhs.len() != self.rows {
            return Err(Fault::InvalidInput);
        }
        solution.resize(rhs.len(), 0.0);
        let mut factor = self.clone();
        let mut indices = vec![0; self.rows];
        let quality = factor
            .lu_decompose(&mut indices)
            .map_err(|_| Fault::Degenerate)?;
        factor.solve_with_lu(&indices, rhs, solution)?;
        Ok(quality)
    }
}

impl Index<usize> for Matrix {
    type Output = [f64];
    fn index(&self, row: usize) -> &[f64] {
        &self.data[row * self.columns..(row + 1) * self.columns]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, row: usize) -> &mut [f64] {
        &mut self.data[row * self.columns..(row + 1) * self.columns]
    }
}
