//! Rational elimination candidates for algebraic coefficient polynomials.
//!
//! Every shared expression leaf receives one auxiliary variable. Elimination
//! supplies candidates only; callers certify them in the original expression.

use crate::{MPoly, Real, RealPoly, UPoly, decider::Work, real::Expr};
use std::{collections::HashMap, sync::Arc};

enum Definition<'a> {
    Expression(&'a Real),
    System(&'a crate::IsolatedRoot, usize),
}

struct Lowering<'a> {
    definitions: Vec<Definition<'a>>,
    indices: HashMap<*const crate::real::Node, usize>,
    systems: Vec<(&'a crate::IsolatedRoot, usize)>,
}

impl<'a> Lowering<'a> {
    fn visit(&mut self, value: &'a Real, work: &mut Work) -> Result<(), &'static str> {
        work.charge(1)?;
        if value.as_rational().is_some() || self.indices.contains_key(&Arc::as_ptr(&value.node)) {
            return Ok(());
        }
        match &value.node.expr {
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                self.visit(a, work)?;
                self.visit(b, work)?;
            }
            Expr::Neg(a) | Expr::Sqrt(a) | Expr::Square(a) | Expr::Pow(a, _) => {
                self.visit(a, work)?
            }
            Expr::Coordinate(root, _) => {
                if !self.systems.iter().any(|(existing, _)| {
                    Arc::ptr_eq(existing.system(), root.system())
                        && existing.root_box() == root.root_box()
                }) {
                    let start = self.definitions.len() + 1;
                    self.systems.push((root, start));
                    for axis in 0..root.root_box().dimension() {
                        work.charge(1)?;
                        self.definitions.push(Definition::System(root, axis));
                    }
                }
            }
            Expr::Rational(_) | Expr::Algebraic(_) => {}
        }
        let index = self.definitions.len() + 1;
        self.definitions.push(Definition::Expression(value));
        self.indices.insert(Arc::as_ptr(&value.node), index);
        Ok(())
    }

    fn expression(&self, value: &Real, nvars: usize) -> MPoly {
        if let Some(rational) = value.as_rational() {
            MPoly::constant(nvars, rational)
        } else {
            MPoly::var(nvars, self.indices[&Arc::as_ptr(&value.node)])
        }
    }

    fn system_start(&self, root: &crate::IsolatedRoot) -> usize {
        self.systems
            .iter()
            .find(|(existing, _)| {
                Arc::ptr_eq(existing.system(), root.system())
                    && existing.root_box() == root.root_box()
            })
            .unwrap()
            .1
    }

    fn equations(&self, work: &mut Work) -> Result<Vec<MPoly>, &'static str> {
        let nvars = self.definitions.len() + 1;
        let mut equations = Vec::with_capacity(self.definitions.len());
        for (index, definition) in self.definitions.iter().enumerate() {
            work.charge(1)?;
            let variable = MPoly::var(nvars, index + 1);
            let equation = match definition {
                Definition::System(root, axis) => {
                    let start = self.system_start(root);
                    MPoly::from_terms(
                        nvars,
                        root.system().equations()[*axis]
                            .terms()
                            .into_iter()
                            .map(|(m, c)| {
                                let mut exponents = vec![0; nvars];
                                exponents[start..start + root.root_box().dimension()]
                                    .copy_from_slice(m.exponents());
                                (exponents, c)
                            })
                            .collect(),
                    )
                }
                Definition::Expression(value) => {
                    let operand = |v| self.expression(v, nvars);
                    match &value.node.expr {
                        Expr::Rational(c) => &variable - &MPoly::constant(nvars, c.clone()),
                        Expr::Algebraic(a) => {
                            let polynomial = a.defining_polynomial();
                            let terms = polynomial
                                .coefficients()
                                .iter()
                                .enumerate()
                                .map(|(power, c)| {
                                    let mut exponents = vec![0; nvars];
                                    exponents[index + 1] =
                                        power.try_into().expect("polynomial degree fits u32");
                                    (exponents, c.clone())
                                })
                                .collect();
                            MPoly::from_terms(nvars, terms)
                        }
                        Expr::Coordinate(root, axis) => {
                            &variable - &MPoly::var(nvars, self.system_start(root) + axis)
                        }
                        Expr::Add(a, b) => &variable - &(&operand(a) + &operand(b)),
                        Expr::Sub(a, b) => &variable - &(&operand(a) - &operand(b)),
                        Expr::Mul(a, b) => {
                            let (a, b) = (operand(a), operand(b));
                            work.mpoly(&a, &b)?;
                            &variable - &(&a * &b)
                        }
                        Expr::Div(a, b) => {
                            let b = operand(b);
                            work.mpoly(&variable, &b)?;
                            &(&variable * &b) - &operand(a)
                        }
                        Expr::Neg(a) => &variable + &operand(a),
                        Expr::Sqrt(a) => &(&variable * &variable) - &operand(a),
                        Expr::Square(a) => {
                            let a = operand(a);
                            work.mpoly(&a, &a)?;
                            &variable - &(&a * &a)
                        }
                        Expr::Pow(a, power) => {
                            let a = operand(a);
                            work.charge(u64::from(*power).saturating_mul(a.terms().len() as u64))?;
                            &variable - &a.power(*power)
                        }
                    }
                }
            };
            equations.push(equation);
        }
        Ok(equations)
    }
}

pub(crate) fn eliminant(polynomial: &RealPoly, work: &mut Work) -> Result<UPoly, &'static str> {
    let terms = polynomial.terms();
    let mut lowering = Lowering {
        definitions: Vec::new(),
        indices: HashMap::new(),
        systems: Vec::new(),
    };
    for (_, coefficient) in &terms {
        lowering.visit(coefficient, work)?;
    }
    let nvars = lowering.definitions.len() + 1;
    let mut equations = lowering.equations(work)?;
    let x = MPoly::var(nvars, 0);
    let mut target = MPoly::zero(nvars);
    for (monomial, coefficient) in &terms {
        let coefficient = lowering.expression(coefficient, nvars);
        let power = x.power(monomial.exponents()[0]);
        work.mpoly(&coefficient, &power)?;
        target = &target + &(&coefficient * &power);
    }
    equations.push(target);
    for variable in (1..nvars).rev() {
        work.charge(1)?;
        let Some(pivot) = equations
            .iter()
            .position(|p| p.terms().iter().any(|(m, _)| m.exponents()[variable] != 0))
        else {
            continue;
        };
        let defining = equations.remove(pivot);
        for equation in &mut equations {
            if equation
                .terms()
                .iter()
                .any(|(m, _)| m.exponents()[variable] != 0)
            {
                work.mpoly(equation, &defining)?;
                *equation =
                    equation.resultant_checked(&defining, variable, &|| work.check_deadline())?;
            }
        }
        equations.retain(|p| !p.is_zero());
    }
    let mut result: Option<UPoly> = None;
    for equation in equations {
        let candidate = equation.to_upoly(0).ok_or("numbers.root-not-certified")?;
        if candidate.is_zero() {
            continue;
        }
        result = Some(match result {
            None => candidate,
            Some(previous) => {
                work.poly(&previous, &candidate)?;
                previous.gcd_checked(&candidate, &|| work.check_deadline())?
            }
        });
    }
    result.ok_or("numbers.root-not-certified")
}

pub(crate) fn resultant(
    left: &RealPoly,
    right: &RealPoly,
    var: usize,
    work: &mut Work,
) -> Result<RealPoly, &'static str> {
    let nvars = left.nvars();
    if left.stored_len() == 0 || right.stored_len() == 0 {
        return Ok(RealPoly::zero(nvars));
    }
    let coefficients = |polynomial: &RealPoly| {
        let terms = polynomial.terms();
        let degree = terms
            .iter()
            .map(|(m, _)| m.exponents()[var] as usize)
            .max()
            .unwrap_or(0);
        let mut coefficients = vec![Vec::new(); degree + 1];
        for (m, c) in terms {
            let mut exponents = m.exponents().to_vec();
            let power = exponents[var] as usize;
            exponents[var] = 0;
            coefficients[power].push((exponents, c));
        }
        coefficients
            .into_iter()
            .map(|terms| RealPoly::from_terms(nvars, terms))
            .collect::<Vec<_>>()
    };
    // Each coefficient above is nonzero by the caller's prior certificates.
    // Reserve before allocating the Sylvester matrix or expanding a determinant.
    let left_degree = left
        .terms()
        .iter()
        .map(|(m, _)| m.exponents()[var] as usize)
        .max()
        .unwrap();
    let right_degree = right
        .terms()
        .iter()
        .map(|(m, _)| m.exponents()[var] as usize)
        .max()
        .unwrap();
    let size = left_degree
        .checked_add(right_degree)
        .ok_or("numbers.budget-exhausted")?;
    let mut bound = 1u64;
    for i in 1..=size {
        work.charge(1)?;
        bound = bound.saturating_mul(i as u64);
        if bound > work.policy.max_work.saturating_sub(work.used) {
            return Err("numbers.budget-exhausted");
        }
    }
    work.charge(bound.saturating_add((size as u64).saturating_pow(2)))?;
    let a = coefficients(left);
    let b = coefficients(right);
    let mut matrix = vec![vec![RealPoly::zero(nvars); size]; size];
    for i in 0..right_degree {
        for j in 0..=left_degree {
            matrix[i][i + j] = a[left_degree - j].clone();
        }
    }
    for i in 0..left_degree {
        for j in 0..=right_degree {
            matrix[right_degree + i][i + j] = b[right_degree - j].clone();
        }
    }
    fn determinant(
        matrix: &[Vec<RealPoly>],
        row: usize,
        columns: &mut [usize],
        nvars: usize,
        work: &mut Work,
    ) -> Result<RealPoly, &'static str> {
        work.charge(1)?;
        if columns.is_empty() {
            return Ok(RealPoly::constant(nvars, Real::one()));
        }
        let mut result = RealPoly::zero(nvars);
        for i in 0..columns.len() {
            columns.swap(0, i);
            let column = columns[0];
            let minor = determinant(matrix, row + 1, &mut columns[1..], nvars, work)?;
            columns.swap(0, i);
            work.charge(
                (matrix[row][column].stored_len() as u64).saturating_mul(minor.stored_len() as u64),
            )?;
            let term = &matrix[row][column] * &minor;
            // Swapping the selected column to the front changes determinant
            // sign once when i is nonzero; the remaining column order is kept.
            result = if i == 0 {
                &result + &term
            } else {
                &result - &term
            };
        }
        Ok(result)
    }
    determinant(&matrix, 0, &mut (0..size).collect::<Vec<_>>(), nvars, work)
}
