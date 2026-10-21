//! Validated Taylor integration. All coefficients are balls, including the remainder.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::ops::{Add, Div, Mul, Neg, Sub};

use catallax_intersect::TracePiece;
use catallax_numbers::{Ball, BudgetPolicy, Decider, Rational, Refusal};

mod analytic;
pub(crate) use analytic::Analytic;

const ORDER: usize = 20;
const SIZE: usize = ORDER + 1;

#[derive(Clone, Debug)]
pub(crate) enum Stop {
    Budget,
    Split,
    Refusal(Refusal),
}

impl From<Refusal> for Stop {
    fn from(refusal: Refusal) -> Self {
        Self::Refusal(refusal)
    }
}

type Checked<T> = Result<T, Stop>;
type Jet = [Ball; SIZE];

pub(crate) struct Work<'a> {
    pub(crate) precision: u32,
    pub(crate) dec: &'a Decider,
    remaining: u64,
    zero: Ball,
    pub(crate) quadrature_dominates: bool,
}

impl<'a> Work<'a> {
    pub(crate) fn new(
        policy: &BudgetPolicy,
        precision: u32,
        remaining: u64,
        dec: &'a Decider,
    ) -> Checked<Self> {
        dec.check_deadline()?;
        assert!(precision > 0 && precision <= policy.max_precision_bits);
        let limbs = u64::from(precision.div_ceil(64));
        let remaining = remaining.checked_sub(limbs * limbs).ok_or(Stop::Budget)?;
        Ok(Self {
            precision,
            dec,
            remaining,
            zero: Ball::from_rational(&Rational::zero(), precision),
            quadrature_dominates: false,
        })
    }

    pub(crate) fn remaining(&self) -> u64 {
        self.remaining
    }

    fn charge(&mut self) -> Checked<()> {
        self.dec.check_deadline()?;
        let limbs = u64::from(self.precision.div_ceil(64));
        let cost = limbs * limbs;
        self.remaining = self.remaining.checked_sub(cost).ok_or(Stop::Budget)?;
        Ok(())
    }

    pub(crate) fn rational(&mut self, r: &Rational) -> Checked<Ball> {
        if r.is_zero() {
            return Ok(self.zero.clone());
        }
        self.charge()?;
        Ok(Ball::from_rational(r, self.precision))
    }

    pub(crate) fn pi(&mut self) -> Checked<Ball> {
        self.charge()?;
        Ok(Ball::pi_checked(self.precision, self.dec)?)
    }

    pub(crate) fn add(&mut self, a: &Ball, b: &Ball) -> Checked<Ball> {
        if a.is_exact_zero() {
            return Ok(b.clone());
        }
        if b.is_exact_zero() {
            return Ok(a.clone());
        }
        self.charge()?;
        Ok(a + b)
    }

    pub(crate) fn sub(&mut self, a: &Ball, b: &Ball) -> Checked<Ball> {
        if b.is_exact_zero() {
            return Ok(a.clone());
        }
        self.charge()?;
        Ok(a - b)
    }

    pub(crate) fn mul(&mut self, a: &Ball, b: &Ball) -> Checked<Ball> {
        if a.is_exact_zero() || b.is_exact_zero() {
            return Ok(self.zero.clone());
        }
        self.charge()?;
        Ok(a * b)
    }

    pub(crate) fn div(&mut self, a: &Ball, b: &Ball) -> Checked<Ball> {
        if b.contains_zero() {
            return Err(Stop::Split);
        }
        if a.is_exact_zero() {
            return Ok(self.zero.clone());
        }
        self.charge()?;
        self.charge()?;
        a.checked_div(b).ok_or(Stop::Split)
    }

    pub(crate) fn sqrt(&mut self, a: &Ball) -> Checked<Ball> {
        self.charge()?;
        a.sqrt().ok_or(Stop::Split)
    }

    pub(crate) fn atan(&mut self, a: &Ball) -> Checked<Ball> {
        self.charge()?;
        Ok(a.atan_checked(self.dec)?)
    }

    pub(crate) fn hull(&mut self, a: &Ball, b: &Ball) -> Checked<Ball> {
        self.charge()?;
        Ok(a.hull(b))
    }

    pub(crate) fn interval(&mut self, lo: &Rational, hi: &Rational) -> Checked<Ball> {
        let a = self.rational(lo)?;
        let b = self.rational(hi)?;
        self.hull(&a, &b)
    }

    fn zeros<const N: usize>(&self) -> [Ball; N] {
        std::array::from_fn(|_| self.zero.clone())
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
enum Node {
    Constant(Rational),
    Input(usize),
    Add(usize, usize),
    Mul(usize, usize),
    Div(usize, usize),
    Neg(usize),
    Sqrt(usize),
}

#[derive(Default)]
struct Nodes {
    entries: Vec<Node>,
    defined: Vec<bool>,
    nonzero: Vec<bool>,
    lookup: BTreeMap<Node, usize>,
}

#[derive(Default)]
pub(crate) struct Builder(RefCell<Nodes>);

#[derive(Clone, Copy)]
pub(crate) struct Expr<'a> {
    index: usize,
    builder: &'a Builder,
}

impl Builder {
    fn node(&self, node: Node) -> Expr<'_> {
        let mut nodes = self.0.borrow_mut();
        let constant = |i: usize| match &nodes.entries[i] {
            Node::Constant(r) => Some(r),
            _ => None,
        };
        // Fold exact rational identities before interval evaluation. Identities
        // that erase a child are used only when that child is everywhere defined.
        let node = match node {
            Node::Add(a, b) => {
                if let (Some(x), Some(y)) = (constant(a), constant(b)) {
                    Node::Constant(x + y)
                } else if constant(a).is_some_and(Rational::is_zero) {
                    return Expr {
                        index: b,
                        builder: self,
                    };
                } else if constant(b).is_some_and(Rational::is_zero) {
                    return Expr {
                        index: a,
                        builder: self,
                    };
                } else if (nodes.entries[a] == Node::Neg(b) && nodes.defined[b])
                    || (nodes.entries[b] == Node::Neg(a) && nodes.defined[a])
                {
                    Node::Constant(Rational::zero())
                } else {
                    Node::Add(a.min(b), a.max(b))
                }
            }
            Node::Mul(a, b) => {
                if let (Some(x), Some(y)) = (constant(a), constant(b)) {
                    Node::Constant(x * y)
                } else if (constant(a).is_some_and(Rational::is_zero) && nodes.defined[b])
                    || (constant(b).is_some_and(Rational::is_zero) && nodes.defined[a])
                {
                    Node::Constant(Rational::zero())
                } else if constant(a).is_some_and(|x| *x == Rational::one()) {
                    return Expr {
                        index: b,
                        builder: self,
                    };
                } else if constant(b).is_some_and(|x| *x == Rational::one()) {
                    return Expr {
                        index: a,
                        builder: self,
                    };
                } else {
                    Node::Mul(a.min(b), a.max(b))
                }
            }
            Node::Div(a, b) => {
                if let (Some(x), Some(y)) = (constant(a), constant(b))
                    && let Some(value) = x.checked_div(y)
                {
                    Node::Constant(value)
                } else if constant(a).is_some_and(Rational::is_zero) && nodes.nonzero[b] {
                    Node::Constant(Rational::zero())
                } else if constant(b).is_some_and(|x| *x == Rational::one()) {
                    return Expr {
                        index: a,
                        builder: self,
                    };
                } else {
                    Node::Div(a, b)
                }
            }
            Node::Neg(a) => {
                if let Some(value) = constant(a) {
                    Node::Constant(-value)
                } else if let Node::Neg(child) = nodes.entries[a] {
                    return Expr {
                        index: child,
                        builder: self,
                    };
                } else {
                    Node::Neg(a)
                }
            }
            Node::Sqrt(a) => {
                if constant(a).is_some_and(|x| x.is_zero() || *x == Rational::one()) {
                    return Expr {
                        index: a,
                        builder: self,
                    };
                } else {
                    Node::Sqrt(a)
                }
            }
            other => other,
        };
        if let Some(index) = nodes.lookup.get(&node) {
            return Expr {
                index: *index,
                builder: self,
            };
        }
        let (defined, nonzero) = match &node {
            Node::Constant(value) => (true, !value.is_zero()),
            Node::Input(_) => (true, false),
            Node::Add(a, b) => (nodes.defined[*a] && nodes.defined[*b], false),
            Node::Mul(a, b) => (
                nodes.defined[*a] && nodes.defined[*b],
                nodes.nonzero[*a] && nodes.nonzero[*b],
            ),
            Node::Div(a, b) => (
                nodes.defined[*a] && nodes.nonzero[*b],
                nodes.nonzero[*a] && nodes.nonzero[*b],
            ),
            Node::Neg(a) => (nodes.defined[*a], nodes.nonzero[*a]),
            Node::Sqrt(a) => (
                constant(*a).is_some_and(|x| *x >= Rational::zero()),
                constant(*a).is_some_and(|x| *x > Rational::zero()),
            ),
        };
        let index = nodes.entries.len();
        nodes.lookup.insert(node.clone(), index);
        nodes.entries.push(node);
        nodes.defined.push(defined);
        nodes.nonzero.push(nonzero);
        Expr {
            index,
            builder: self,
        }
    }

    pub(crate) fn r(&self, r: Rational) -> Expr<'_> {
        self.node(Node::Constant(r))
    }

    pub(crate) fn k(&self, k: i64) -> Expr<'_> {
        self.r(Rational::from_i64(k))
    }

    pub(crate) fn input(&self, i: usize) -> Expr<'_> {
        self.node(Node::Input(i))
    }

    pub(crate) fn finish(&self, outputs: &[Expr<'_>]) -> Program {
        let nodes = self.0.borrow();
        let source = &nodes.entries;
        let mut needed = vec![false; source.len()];
        for output in outputs {
            needed[output.index] = true;
        }
        for i in (0..source.len()).rev() {
            if needed[i] {
                match source[i] {
                    Node::Add(a, b) | Node::Mul(a, b) | Node::Div(a, b) => {
                        needed[a] = true;
                        needed[b] = true;
                    }
                    Node::Neg(a) | Node::Sqrt(a) => needed[a] = true,
                    _ => {}
                }
            }
        }
        let mut indices = vec![0; source.len()];
        let mut nodes = Vec::new();
        for (i, node) in source.iter().enumerate() {
            if needed[i] {
                indices[i] = nodes.len();
                nodes.push(match node {
                    Node::Add(a, b) => Node::Add(indices[*a], indices[*b]),
                    Node::Mul(a, b) => Node::Mul(indices[*a], indices[*b]),
                    Node::Div(a, b) => Node::Div(indices[*a], indices[*b]),
                    Node::Neg(a) => Node::Neg(indices[*a]),
                    Node::Sqrt(a) => Node::Sqrt(indices[*a]),
                    other => other.clone(),
                });
            }
        }
        Program {
            nodes,
            outputs: outputs.iter().map(|e| indices[e.index]).collect(),
        }
    }
}

impl Expr<'_> {
    pub(crate) fn sqrt(self) -> Self {
        self.builder.node(Node::Sqrt(self.index))
    }

    pub(crate) fn square(self) -> Self {
        self * self
    }
}

impl<'a> Add for Expr<'a> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.builder.node(Node::Add(self.index, rhs.index))
    }
}
impl<'a> Sub for Expr<'a> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}
impl<'a> Mul for Expr<'a> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.builder.node(Node::Mul(self.index, rhs.index))
    }
}
impl<'a> Div for Expr<'a> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.builder.node(Node::Div(self.index, rhs.index))
    }
}
impl Neg for Expr<'_> {
    type Output = Self;
    fn neg(self) -> Self {
        self.builder.node(Node::Neg(self.index))
    }
}

pub(crate) struct Program {
    nodes: Vec<Node>,
    outputs: Vec<usize>,
}

impl Program {
    fn evaluate<const N: usize>(
        &self,
        input: &[[Ball; N]],
        work: &mut Work,
    ) -> Checked<Vec<[Ball; N]>> {
        let mut values: Vec<[Ball; N]> = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            work.dec.check_deadline()?;
            if let Node::Input(i) = node {
                values.push(input[*i].clone());
                continue;
            }
            let mut value = work.zeros::<N>();
            match node {
                Node::Constant(r) => value[0] = work.rational(r)?,
                Node::Input(_) => unreachable!(),
                Node::Add(a, b) => {
                    for (k, coefficient) in value.iter_mut().enumerate() {
                        *coefficient = work.add(&values[*a][k], &values[*b][k])?;
                    }
                }
                Node::Neg(a) => {
                    for (k, coefficient) in value.iter_mut().enumerate() {
                        if !values[*a][k].is_exact_zero() {
                            work.charge()?;
                            *coefficient = -&values[*a][k];
                        }
                    }
                }
                Node::Mul(a, b) if a == b && N > 1 => {
                    let order = values[*a]
                        .iter()
                        .rposition(|b| !b.is_exact_zero())
                        .unwrap_or(0);
                    let mut products: [Option<Ball>; N] = std::array::from_fn(|_| None);
                    for (k, coefficient) in
                        value.iter_mut().enumerate().take((2 * order + 1).min(N))
                    {
                        for j in k.saturating_sub(order)..=k.min(order) {
                            let index = j.min(k - j);
                            if j == index {
                                products[index] = if !values[*a][j].is_exact_zero()
                                    && !values[*a][k - j].is_exact_zero()
                                {
                                    Some(work.mul(&values[*a][j], &values[*a][k - j])?)
                                } else {
                                    None
                                };
                            }
                            // Reuse products without reordering rounded additions.
                            if let Some(product) = &products[index] {
                                *coefficient = work.add(coefficient, product)?;
                            }
                        }
                    }
                }
                Node::Mul(a, b) => {
                    let a_order = values[*a]
                        .iter()
                        .rposition(|b| !b.is_exact_zero())
                        .unwrap_or(0);
                    let b_order = values[*b]
                        .iter()
                        .rposition(|b| !b.is_exact_zero())
                        .unwrap_or(0);
                    for (k, coefficient) in value
                        .iter_mut()
                        .enumerate()
                        .take((a_order + b_order + 1).min(N))
                    {
                        for j in k.saturating_sub(b_order)..=k.min(a_order) {
                            if !values[*a][j].is_exact_zero() && !values[*b][k - j].is_exact_zero()
                            {
                                let product = work.mul(&values[*a][j], &values[*b][k - j])?;
                                *coefficient = work.add(coefficient, &product)?;
                            }
                        }
                    }
                }
                Node::Div(a, b) => {
                    if values[*b][0].contains_zero() {
                        return Err(Stop::Split);
                    }
                    if N == 1 {
                        value[0] = work.div(&values[*a][0], &values[*b][0])?;
                    } else if !values[*a].iter().all(Ball::is_exact_zero) {
                        let one = work.rational(&Rational::one())?;
                        let inverse = work.div(&one, &values[*b][0])?;
                        let order = values[*b]
                            .iter()
                            .rposition(|b| !b.is_exact_zero())
                            .unwrap_or(0);
                        for k in 0..N {
                            let mut numerator = values[*a][k].clone();
                            for j in 1..=k.min(order) {
                                if !values[*b][j].is_exact_zero() && !value[k - j].is_exact_zero() {
                                    let product = work.mul(&values[*b][j], &value[k - j])?;
                                    numerator = work.sub(&numerator, &product)?;
                                }
                            }
                            value[k] = work.mul(&numerator, &inverse)?;
                        }
                    }
                }
                Node::Sqrt(a) => {
                    value[0] = work.sqrt(&values[*a][0])?;
                    if value[0].contains_zero() {
                        return Err(Stop::Split);
                    }
                    if N > 1 {
                        let two = work.rational(&Rational::from_i64(2))?;
                        let denominator = work.mul(&two, &value[0])?;
                        let one = work.rational(&Rational::one())?;
                        let inverse = work.div(&one, &denominator)?;
                        let mut products: [Option<Ball>; N] = std::array::from_fn(|_| None);
                        for k in 1..N {
                            let mut numerator = values[*a][k].clone();
                            for j in 1..k {
                                let index = j.min(k - j);
                                if j == index {
                                    products[index] = if !value[j].is_exact_zero()
                                        && !value[k - j].is_exact_zero()
                                    {
                                        Some(work.mul(&value[j], &value[k - j])?)
                                    } else {
                                        None
                                    };
                                }
                                if let Some(product) = &products[index] {
                                    numerator = work.sub(&numerator, product)?;
                                }
                            }
                            value[k] = work.mul(&numerator, &inverse)?;
                        }
                    }
                }
            }
            values.push(value);
        }
        Ok(self.outputs.iter().map(|i| values[*i].clone()).collect())
    }

    pub(crate) fn point(
        &self,
        position: &[Ball; 3],
        extras: &[Ball],
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        let mut input: Vec<[Ball; 1]> = position.iter().map(|b| [b.clone()]).collect();
        input.extend((0..3).map(|_| work.zeros()));
        input.extend(extras.iter().map(|b| [b.clone()]));
        Ok(self
            .evaluate(&input, work)?
            .into_iter()
            .map(|[value]| value)
            .collect())
    }

    pub(crate) fn panel(
        &self,
        piece: &TracePiece,
        a: &Rational,
        b: &Rational,
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        let midpoint = (a + b).mul_pow2(-1);
        let center = work.rational(&midpoint)?;
        let range = work.interval(&a.clone().min(b.clone()), &a.clone().max(b.clone()))?;
        // Symmetric integration uses center orders 0, 2, ..., ORDER-2.
        // Only the whole-panel remainder needs the full order.
        let at_center = input::<{ ORDER - 1 }>(piece, &center, work)?;
        let on_range = input::<SIZE>(piece, &range, work)?;
        let center = self.evaluate(&at_center, work)?;
        let range = self.evaluate(&on_range, work)?;
        center
            .iter()
            .zip(&range)
            .map(|(c, r)| integrate(c, r, a, b, &midpoint, work))
            .collect()
    }

    pub(crate) fn rectangle(
        &self,
        piece: &TracePiece,
        a: &Rational,
        b: &Rational,
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        let range = work.interval(&a.clone().min(b.clone()), &a.clone().max(b.clone()))?;
        let coefficients = piece.taylor(&range, 1, work.dec)?;
        let mut input = Vec::with_capacity(6);
        for row in coefficients.iter().take(2) {
            for coefficient in row {
                input.push([coefficient.clone()]);
            }
        }
        let width = work.rational(&(b - a))?;
        self.evaluate(&input, work)?
            .iter()
            .map(|j| work.mul(&j[0], &width))
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn lifted_panel(
        &self,
        angular: &Program,
        piece: &TracePiece,
        a: &Rational,
        b: &Rational,
        start: &Ball,
        rectangle: bool,
        work: &mut Work,
    ) -> Checked<(Vec<Ball>, Ball)> {
        let midpoint = (a + b).mul_pow2(-1);
        let center_argument = work.rational(&midpoint)?;
        let range_argument = work.interval(&a.clone().min(b.clone()), &a.clone().max(b.clone()))?;
        let mut range_input = if rectangle {
            let coefficients = piece.taylor(&range_argument, 1, work.dec)?;
            coefficients
                .iter()
                .flat_map(|row| row.iter())
                .map(|b| {
                    let mut jet = work.zeros();
                    jet[0] = b.clone();
                    jet
                })
                .collect::<Vec<_>>()
        } else {
            input::<SIZE>(piece, &range_argument, work)?
        };
        let range_angular = angular.evaluate(&range_input, work)?.remove(0);
        let displacement = b - a;
        let offsets = work.interval(
            &displacement.clone().min(Rational::zero()),
            &displacement.clone().max(Rational::zero()),
        )?;
        let variation = work.mul(&range_angular[0], &offsets)?;
        let value_range = work.add(start, &variation)?;
        let mut angle_range = work.zeros();
        angle_range[0] = value_range;
        for k in 1..SIZE {
            let divisor = work.rational(&Rational::from_i64(k as i64))?;
            angle_range[k] = work.div(&range_angular[k - 1], &divisor)?;
        }
        range_input.push(angle_range);
        let range = self.evaluate(&range_input, work)?;
        if rectangle {
            let width = work.rational(&displacement)?;
            let increment = work.mul(&range_angular[0], &width)?;
            let end = work.add(start, &increment)?;
            let values = range
                .iter()
                .map(|j| work.mul(&j[0], &width))
                .collect::<Checked<_>>()?;
            return Ok((values, end));
        }
        let mut center_input = input::<SIZE>(piece, &center_argument, work)?;
        let center_angular = angular.evaluate(&center_input, work)?.remove(0);
        let increment = integrate(&center_angular, &range_angular, a, b, &midpoint, work)?;
        let half = integrate(
            &center_angular,
            &range_angular,
            a,
            &midpoint,
            &midpoint,
            work,
        )?;
        let mut angle_center = work.zeros();
        angle_center[0] = work.add(start, &half)?;
        for k in 1..SIZE {
            let divisor = work.rational(&Rational::from_i64(k as i64))?;
            angle_center[k] = work.div(&center_angular[k - 1], &divisor)?;
        }
        center_input.push(angle_center);
        let center = self.evaluate(&center_input, work)?;
        let values = center
            .iter()
            .zip(&range)
            .map(|(c, r)| integrate(c, r, a, b, &midpoint, work))
            .collect::<Checked<_>>()?;
        Ok((values, work.add(start, &increment)?))
    }
}

fn input<const N: usize>(piece: &TracePiece, t: &Ball, work: &mut Work) -> Checked<Vec<[Ball; N]>> {
    // Piece 4's internal operations are excluded from the measurement work count.
    let coefficients = piece.taylor(t, N, work.dec)?;
    let mut inputs = Vec::with_capacity(6);
    for (i, _) in coefficients[0].iter().enumerate() {
        inputs.push(std::array::from_fn(|k| coefficients[k][i].clone()));
    }
    for (i, _) in coefficients[0].iter().enumerate() {
        let mut derivative = work.zeros();
        for (k, coefficient) in derivative.iter_mut().enumerate() {
            let factor = work.rational(&Rational::from_i64((k + 1) as i64))?;
            *coefficient = work.mul(&coefficients[k + 1][i], &factor)?;
        }
        inputs.push(derivative);
    }
    Ok(inputs)
}

fn integrate(
    center: &[Ball],
    range: &Jet,
    a: &Rational,
    b: &Rational,
    midpoint: &Rational,
    work: &mut Work,
) -> Checked<Ball> {
    let left = a - midpoint;
    let right = b - midpoint;
    let mut result = work.zero.clone();
    for (k, coefficient) in center.iter().enumerate().take(ORDER) {
        let power = i32::try_from(k + 1).unwrap();
        let factor = (right.pow(power).unwrap() - left.pow(power).unwrap())
            / Rational::from_i64(i64::from(power));
        let factor = work.rational(&factor)?;
        let term = work.mul(coefficient, &factor)?;
        result = work.add(&result, &term)?;
    }
    // ORDER is even, so the integral of the nonnegative remainder weight
    // multiplies the whole-range Taylor coefficient without changing its sign.
    let power = (ORDER + 1) as i32;
    let factor = (right.pow(power).unwrap() - left.pow(power).unwrap())
        / Rational::from_i64(i64::from(power));
    let factor = work.rational(&factor)?;
    let remainder = work.mul(&range[ORDER], &factor)?;
    let truncation = remainder.mid().abs() + remainder.rad();
    if truncation > result.hi() - result.lo() {
        work.quadrature_dominates = true;
    }
    work.add(&result, &remainder)
}

pub(crate) fn angle(cosine: &Ball, sine: &Ball, work: &mut Work) -> Checked<Ball> {
    if !cosine.contains_zero() {
        let ratio = work.div(sine, cosine)?;
        let angle = work.atan(&ratio)?;
        if cosine.hi() < Rational::zero() {
            let pi = work.pi()?;
            work.add(&angle, &pi)
        } else {
            Ok(angle)
        }
    } else if !sine.contains_zero() {
        let ratio = work.div(cosine, sine)?;
        let angle = work.atan(&ratio)?;
        let pi = work.pi()?;
        let two = work.rational(&Rational::from_i64(2))?;
        let half = work.div(&pi, &two)?;
        let half = if sine.hi() < Rational::zero() {
            let minus = work.rational(&Rational::from_i64(-1))?;
            work.mul(&half, &minus)?
        } else {
            half
        };
        work.sub(&half, &angle)
    } else {
        Err(Stop::Split)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_square_and_positive_root_preserve_zero_coefficients() {
        let policy = BudgetPolicy::default();
        let dec = Decider::new(policy);
        let mut work = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
        let builder = Builder::default();
        let x = builder.input(0);
        let program = builder.finish(&[x.square(), x.sqrt()]);
        // Input is (1 + t - t^3)^2. The outputs are its square and
        // its positive square root near t=0, through degree six.
        let input =
            [1, 2, 1, -2, -2, 0, 1].map(|n| Ball::from_rational(&Rational::from_i64(n), 128));
        let actual = program.evaluate(&[input], &mut work).unwrap();
        for (row, expected) in [[1, 4, 6, 0, -11, -12, 2], [1, 1, 0, -1, 0, 0, 0]]
            .iter()
            .enumerate()
        {
            for (k, coefficient) in expected.iter().enumerate() {
                assert!(actual[row][k].contains(&Rational::from_i64(*coefficient)));
            }
        }
    }

    #[test]
    fn symmetric_panel_retains_the_signed_order_twenty_remainder() {
        let policy = BudgetPolicy::default();
        let dec = Decider::new(policy);
        let mut work = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
        let center: [Ball; ORDER - 1] = std::array::from_fn(|k| {
            let coefficient = (k as i64 + 1) * if k % 2 == 0 { 1 } else { -1 };
            Ball::from_rational(&Rational::from_i64(coefficient), 128)
        });
        // The ten even center terms each integrate to 2 on [-1,1].
        // The order-20 term contributes +2 or -2, not a zero odd moment.
        for sign in [-1, 1] {
            let mut range = work.zeros();
            range[ORDER] = Ball::from_rational(&Rational::from_i64(sign * 21), 128);
            let value = integrate(
                &center,
                &range,
                &Rational::from_i64(-1),
                &Rational::one(),
                &Rational::zero(),
                &mut work,
            )
            .unwrap();
            assert!(value.contains(&Rational::from_i64(20 + sign * 2)));
        }
    }
}
