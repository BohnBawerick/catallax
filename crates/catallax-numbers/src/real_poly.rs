//! Polynomials over exact algebraic expressions, with a rational arithmetic path.

use crate::{
    Ball, BoxRoots, Decider, Integer, MPoly, Monomial, Outcome, PolySystem, Rational, Real,
    Refusal, RootBox, Sign, StrictSign,
};
use std::{collections::BTreeMap, sync::Arc};

/// A fixed-variable polynomial whose coefficient definitions retain shared roots.
///
/// Construction order is graded lexicographic monomial order, then immutable
/// expression order within each coefficient sum. Neither construction identity
/// nor the number of stored terms decides mathematical equality or degree.
#[derive(Clone, Debug)]
pub struct RealPoly {
    coefficients: Coefficients,
}

#[derive(Clone, Debug)]
enum Coefficients {
    Rational(MPoly),
    Algebraic {
        nvars: usize,
        terms: BTreeMap<Monomial, Real>,
    },
}

impl RealPoly {
    /// Zero in the specified number of variables.
    pub fn zero(nvars: usize) -> Self {
        MPoly::zero(nvars).into()
    }

    /// A constant, without validating the coefficient's expression domain.
    pub fn constant(nvars: usize, coefficient: Real) -> Self {
        if let Some(value) = coefficient.rational_ref() {
            return MPoly::constant(nvars, value.clone()).into();
        }
        Self::from_terms(nvars, vec![(vec![0; nvars], coefficient)])
    }

    /// A variable. Panics when the index is out of range.
    pub fn var(nvars: usize, index: usize) -> Self {
        MPoly::var(nvars, index).into()
    }

    /// Combines like terms in a fixed construction order. Only syntactically
    /// rational zeros are removed; other zeros require a decision.
    pub fn from_terms(nvars: usize, mut terms: Vec<(Vec<u32>, Real)>) -> Self {
        assert!(terms.iter().all(|(e, _)| e.len() == nvars));
        if terms.iter().all(|(_, c)| c.rational_ref().is_some()) {
            let rational = terms
                .into_iter()
                .map(|(exponents, coefficient)| {
                    (
                        exponents,
                        coefficient
                            .rational_ref()
                            .expect("rational coefficient")
                            .clone(),
                    )
                })
                .collect();
            return MPoly::from_terms(nvars, rational).into();
        }
        terms.sort_by(|(a, x), (b, y)| {
            a.iter()
                .map(|e| u64::from(*e))
                .sum::<u64>()
                .cmp(&b.iter().map(|e| u64::from(*e)).sum::<u64>())
                .then_with(|| a.cmp(b))
                .then_with(|| x.construction_cmp(y))
        });
        let mut combined = BTreeMap::new();
        for (exponents, coefficient) in terms {
            let monomial = Monomial::new(exponents);
            combined
                .entry(monomial)
                .and_modify(|sum: &mut Real| *sum = &*sum + &coefficient)
                .or_insert(coefficient);
        }
        combined.retain(|_, c| !c.rational_ref().is_some_and(Rational::is_zero));
        Self {
            coefficients: Coefficients::Algebraic {
                nvars,
                terms: combined,
            },
        }
    }

    /// Number of variables.
    pub fn nvars(&self) -> usize {
        match &self.coefficients {
            Coefficients::Rational(p) => p.nvars(),
            Coefficients::Algebraic { nvars, .. } => *nvars,
        }
    }

    /// Coefficient definitions in graded lexicographic monomial order.
    pub fn terms(&self) -> Vec<(Monomial, Real)> {
        match &self.coefficients {
            Coefficients::Rational(p) => {
                p.terms().into_iter().map(|(m, c)| (m, c.into())).collect()
            }
            Coefficients::Algebraic { terms, .. } => {
                terms.iter().map(|(m, c)| (m.clone(), c.clone())).collect()
            }
        }
    }

    /// Cheap exact conversion. Does not refine, eliminate, or round coefficients.
    pub fn as_rational(&self) -> Option<MPoly> {
        match &self.coefficients {
            Coefficients::Rational(p) => Some(p.clone()),
            Coefficients::Algebraic { nvars, terms } => {
                let terms: Option<Vec<_>> = terms
                    .iter()
                    .map(|(m, c)| Some((m.exponents().to_vec(), c.as_rational()?)))
                    .collect();
                Some(MPoly::from_terms(*nvars, terms?))
            }
        }
    }

    /// Borrowed rational fast path, without cloning a rational polynomial.
    pub fn rational(&self) -> Option<&MPoly> {
        match &self.coefficients {
            Coefficients::Rational(p) => Some(p),
            Coefficients::Algebraic { .. } => None,
        }
    }

    /// Consumes a cached rational polynomial without cloning its coefficients.
    pub fn into_rational(self) -> Option<MPoly> {
        match self.coefficients {
            Coefficients::Rational(polynomial) => Some(polynomial),
            Coefficients::Algebraic { nvars, terms } => {
                let terms: Option<Vec<_>> = terms
                    .into_iter()
                    .map(|(m, c)| Some((m.exponents().to_vec(), c.as_rational()?)))
                    .collect();
                Some(MPoly::from_terms(nvars, terms?))
            }
        }
    }

    /// Exact evaluation in monomial order. Builds expressions without decisions.
    pub fn eval_real(&self, x: &[Real]) -> Real {
        assert_eq!(x.len(), self.nvars());
        match &self.coefficients {
            Coefficients::Rational(p) => p.eval_real(x),
            Coefficients::Algebraic { terms, .. } => {
                terms.iter().fold(Real::zero(), |sum, (m, c)| {
                    let value = m
                        .exponents()
                        .iter()
                        .zip(x)
                        .filter(|(e, _)| **e != 0)
                        .fold(c.clone(), |term, (e, v)| &term * &v.pow_u32(*e));
                    &sum + &value
                })
            }
        }
    }

    /// Certified enclosure including both coefficient and input uncertainty.
    pub fn eval_ball(&self, x: &[Ball], dec: &Decider) -> Result<Ball, Refusal> {
        assert_eq!(x.len(), self.nvars());
        dec.check_deadline()?;
        match &self.coefficients {
            Coefficients::Rational(p) => Ok(p.eval_ball(x)),
            Coefficients::Algebraic { terms, .. } => {
                let precision = x.iter().map(Ball::precision).max().unwrap_or(1);
                let mut sum = Ball::from_rational(&Rational::zero(), precision);
                for (m, c) in terms {
                    let mut value = dec.enclose_at(c, precision)?;
                    for (e, v) in m.exponents().iter().zip(x).filter(|(e, _)| **e != 0) {
                        value = &value * &v.pow_u32(*e);
                    }
                    sum = &sum + &value;
                }
                Ok(sum)
            }
        }
    }

    /// Partial derivative. Does not decide whether resulting coefficients vanish.
    pub fn partial(&self, var: usize) -> Self {
        assert!(var < self.nvars());
        if let Some(p) = self.rational() {
            return p.partial(var).into();
        }
        Self::from_terms(
            self.nvars(),
            self.terms()
                .into_iter()
                .filter(|(m, _)| m.exponents()[var] != 0)
                .map(|(m, c)| {
                    let mut exponents = m.exponents().to_vec();
                    let factor = Real::from(Rational::from_integer(Integer::from_u64(u64::from(
                        exponents[var],
                    ))));
                    exponents[var] -= 1;
                    (exponents, &c * &factor)
                })
                .collect(),
        )
    }

    /// Multiplies every coefficient by the exact scalar definition.
    pub fn scale(&self, scalar: &Real) -> Self {
        if let (Some(p), Some(c)) = (self.rational(), scalar.rational_ref()) {
            return p.scale(c).into();
        }
        Self::from_terms(
            self.nvars(),
            self.terms()
                .into_iter()
                .map(|(m, c)| (m.exponents().to_vec(), &c * scalar))
                .collect(),
        )
    }

    /// Substitutes a polynomial for one variable without eliminating coefficients.
    pub fn substitute(&self, var: usize, value: &Self) -> Self {
        assert!(var < self.nvars());
        assert_eq!(self.nvars(), value.nvars());
        if let (Some(p), Some(v)) = (self.rational(), value.rational()) {
            return p.substitute(var, v).into();
        }
        self.terms()
            .into_iter()
            .fold(Self::zero(self.nvars()), |sum, (m, c)| {
                let mut exponents = m.exponents().to_vec();
                let power = exponents[var];
                exponents[var] = 0;
                &sum + &(&Self::from_terms(self.nvars(), vec![(exponents, c)])
                    * &value.pow_u32(power))
            })
    }

    /// Nonnegative integer power; exponent zero is the constant polynomial one.
    pub fn pow_u32(&self, mut exponent: u32) -> Self {
        let mut factor = self.clone();
        let mut result = Self::constant(self.nvars(), Real::one());
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = &result * &factor;
            }
            exponent >>= 1;
            if exponent != 0 {
                factor = &factor * &factor;
            }
        }
        result
    }

    /// Mathematical degree, after certifying the coefficients in descending order.
    pub fn total_degree(&self, dec: &mut Decider) -> Result<Option<u32>, Refusal> {
        if let Some(p) = self.rational() {
            return Ok(p.total_degree());
        }
        for (m, c) in self.terms().into_iter().rev() {
            if nonzero(&c, dec)? {
                return Ok(Some(m.degree()));
            }
        }
        Ok(None)
    }

    /// Distinct real roots of a univariate polynomial in increasing order.
    ///
    /// Rational coefficients use the existing rational solver. Linear and
    /// quadratic algebraic equations retain their exact closed definitions.
    /// Higher degrees take the rational norm when every coefficient lies in one
    /// quadratic tower, and otherwise lower shared coefficient definitions to
    /// rational equations. Either way they isolate candidates and certify
    /// every candidate in the original equation. Root domains and
    /// conjugate choices are never inferred from an eliminated equation alone.
    pub fn real_roots(&self, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        assert_eq!(self.nvars(), 1);
        if let Some(p) = self.rational() {
            return root_values(dec.real_roots(&p.to_upoly(0).unwrap()));
        }
        for (_, coefficient) in self.terms() {
            dec.enclose(&coefficient, &Rational::one())?;
        }
        let Some(degree) = self.total_degree(dec)? else {
            return Err(Refusal::new(
                "numbers.zero-polynomial",
                "every real number is a root",
            ));
        };
        if degree <= 2 {
            let mut coefficients = [Real::zero(), Real::zero(), Real::zero()];
            for (m, coefficient) in self.terms() {
                if m.exponents()[0] <= degree {
                    coefficients[m.exponents()[0] as usize] = coefficient;
                }
            }
            let [c, b, a] = coefficients;
            if degree == 0 {
                return Ok(vec![]);
            }
            if degree == 1 {
                return Ok(vec![-c / b]);
            }
            let discriminant = &b.square() - &(&Real::from_i64(4) * &(&a * &c));
            let sign = |value: &Real, dec: &mut Decider| match dec.sign(value) {
                Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
                Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
                Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
                Outcome::NotExact(n) => Err(Refusal::from(n)),
                Outcome::Refused(r) => Err(r),
            };
            let denominator = &Real::from_i64(2) * &a;
            return match sign(&discriminant, dec)? {
                Sign::Negative => Ok(vec![]),
                Sign::Zero => Ok(vec![-b / denominator]),
                Sign::Positive => {
                    let root = discriminant.sqrt();
                    let lower = &(-&b) - &root;
                    let upper = &(-&b) + &root;
                    let mut roots = vec![&lower / &denominator, &upper / &denominator];
                    if sign(&a, dec)? == Sign::Negative {
                        roots.swap(0, 1);
                    }
                    Ok(roots)
                }
            };
        }
        // Coefficients in one quadratic tower have a rational norm, taken one
        // generator at a time as `(A + sqrt(r) B)(A - sqrt(r) B)`. Every real
        // root is a root of the norm, and the sign tree decides exactly which
        // of those the polynomial keeps. Only other coefficients need the
        // general eliminant, whose auxiliary variables grow with every
        // expression node.
        if let Some(split) = dec.real_poly_tower_split(self)? {
            let mut roots = Vec::new();
            for candidate in root_values(dec.real_roots(&split.norm))? {
                if tree_sign(&split.sign, &candidate, dec)? == Sign::Zero {
                    roots.push(candidate);
                }
            }
            return Ok(roots);
        }
        let eliminant = dec.real_poly_eliminant(self)?;
        let candidates = root_values(dec.real_roots(&eliminant))?;
        let mut roots = Vec::new();
        for candidate in candidates {
            if !nonzero(&self.eval_real(std::slice::from_ref(&candidate)), dec)? {
                roots.push(candidate);
            }
        }
        Ok(roots)
    }

    /// The one simple common root of `n` polynomials in `n` variables inside a
    /// closed rational box, exactly, where every coefficient lies in one quadratic
    /// tower. The tower's generators join the variables, so the coordinates are
    /// leaves of a small rational system and no eliminant of the tower is built.
    /// The caller establishes the root; a box that certifies none refuses, and
    /// `None` only means a coefficient leaves a quadratic tower.
    ///
    /// Returned coordinates retain the original coefficient definitions and
    /// certified equations. `Decider::sign` can reuse those equations as exact
    /// identities, independently of enclosure width or a separation bound.
    pub fn simple_root_in(
        equations: &[Self],
        bounds: &[[Rational; 2]],
        precision: u32,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Real>>, Refusal> {
        assert_eq!(equations.len(), bounds.len());
        assert!(equations.iter().all(|e| e.nvars() == bounds.len()));
        let Some((lowered, generators)) = dec.tower_root_system(equations, precision)? else {
            return Ok(None);
        };
        // Room above the arithmetic's rounding, or the box cannot map into itself.
        let slack = Rational::one().mul_pow2(-i64::from(precision / 2));
        let (mut lo, mut hi): (Vec<_>, Vec<_>) = bounds
            .iter()
            .map(|[lo, hi]| (lo.clone(), hi.clone()))
            .unzip();
        for generator in &generators {
            lo.push(generator.lo() - slack.clone());
            hi.push(generator.hi() + slack.clone());
        }
        let not_certified = || {
            Refusal::new(
                "numbers.root-not-certified",
                "the box does not hold one simple root in the coefficient tower",
            )
        };
        let system = PolySystem::new(lowered).ok_or_else(not_certified)?;
        let candidate = RootBox::new(lo, hi).ok_or_else(not_certified)?;
        match dec.certify_root_at(&Arc::new(system), &candidate, precision) {
            Outcome::Proven(BoxRoots::One(mut root)) => {
                // A lowered generator denotes the original positive square root,
                // never another conjugate selected only by the remaining equations.
                for lo in &root.root_box().lo()[bounds.len()..] {
                    dec.check_deadline_if_set()?;
                    if lo.sign() != Sign::Positive {
                        return Err(not_certified());
                    }
                }
                root.original_equations = Some(Arc::from(equations));
                let root = Arc::new(root);
                Ok(Some(
                    (0..bounds.len())
                        .map(|axis| Real::from_shared_root_coordinate(root.clone(), axis))
                        .collect(),
                ))
            }
            Outcome::Proven(BoxRoots::None) | Outcome::ProvenCoincident(_) => Err(not_certified()),
            Outcome::NotExact(n) => Err(Refusal::from(n)),
            Outcome::Refused(refusal) => Err(refusal),
        }
    }

    /// The simple root of a square system with retained algebraic coefficients.
    /// The caller establishes existence for the actual coefficients in `bounds`.
    /// Shared coefficient expressions and pinned system coordinates are lifted
    /// together, then a rational Krawczyk certificate selects one joint root.
    /// This does not prove completeness outside that caller-certified graph.
    /// Principal square-root signs and original leaf/system isolators constrain
    /// every algebraic choice. The joint certificate binds lifted coefficients
    /// to their original DAGs before retaining the original equations as exact
    /// identity authority.
    pub fn simple_algebraic_root_in(
        equations: &[Self],
        bounds: &[[Rational; 2]],
        precision: u32,
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        Self::simple_algebraic_root_in_mode(equations, bounds, precision, None, dec)
    }

    /// Requested original-system axes, in order, as univariate algebraic values.
    ///
    /// The caller retains the existence obligation of
    /// [`Self::simple_algebraic_root_in`]. Each requested coordinate must have
    /// exactly one projected value in its certified closed box; otherwise this
    /// view refuses. The joint root and coefficient authority are unchanged.
    ///
    /// # Panics
    ///
    /// Panics if a requested axis is outside the original system.
    pub fn simple_algebraic_root_in_reconstructed(
        equations: &[Self],
        bounds: &[[Rational; 2]],
        axes: &[usize],
        precision: u32,
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        assert!(axes.iter().all(|&axis| axis < bounds.len()));
        Self::simple_algebraic_root_in_mode(equations, bounds, precision, Some(axes), dec)
    }

    fn simple_algebraic_root_in_mode(
        equations: &[Self],
        bounds: &[[Rational; 2]],
        precision: u32,
        reconstruct: Option<&[usize]>,
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let refused = || {
            Refusal::new(
                "numbers.root-not-certified",
                "the joint box does not certify the caller's simple algebraic-system root",
            )
        };
        if equations.is_empty()
            || equations.len() != bounds.len()
            || equations.iter().any(|e| e.nvars() != bounds.len())
        {
            return Err(refused());
        }
        let (lowered, constants) = dec.algebraic_root_system(equations, reconstruct.is_some())?;
        let slack = Rational::one().mul_pow2(-i64::from(precision / 2));
        let (mut lo, mut hi): (Vec<_>, Vec<_>) = bounds
            .iter()
            .map(|[lo, hi]| (lo.clone(), hi.clone()))
            .unzip();
        for constant in constants {
            let enclosure = dec.enclose_budgeted_at(&constant, precision)?;
            let enclosed_lo = enclosure.lo();
            let enclosed_hi = enclosure.hi();
            let mut lower = &enclosed_lo - &slack;
            let mut upper = &enclosed_hi + &slack;
            // Arithmetic definitions are functional once their inputs are
            // fixed. Retain exact selectors for the non-functional leaves,
            // rather than identifying them from overlapping enclosures.
            match &constant.node.expr {
                crate::real::Expr::Sqrt(_) => {
                    if lower.sign() == Sign::Negative {
                        lower = Rational::zero();
                    }
                }
                crate::real::Expr::Algebraic(value) => {
                    let (selected_lo, selected_hi) = value.isolating_interval();
                    // Algebraic isolators are open, whereas Krawczyk uses
                    // closed containment, including roots on a box face.
                    if lower <= selected_lo {
                        if enclosed_lo <= selected_lo {
                            return Err(refused());
                        }
                        lower = (&enclosed_lo + &selected_lo).mul_pow2(-1);
                    }
                    if upper >= selected_hi {
                        if enclosed_hi >= selected_hi {
                            return Err(refused());
                        }
                        upper = (&enclosed_hi + &selected_hi).mul_pow2(-1);
                    }
                }
                crate::real::Expr::Coordinate(root, axis) => {
                    if lower < root.root_box().lo()[*axis] {
                        lower = root.root_box().lo()[*axis].clone();
                    }
                    if upper > root.root_box().hi()[*axis] {
                        upper = root.root_box().hi()[*axis].clone();
                    }
                }
                _ => {}
            }
            lo.push(lower);
            hi.push(upper);
        }
        let system = Arc::new(PolySystem::new(lowered).ok_or_else(refused)?);
        let candidate = RootBox::new(lo, hi).ok_or_else(refused)?;
        match dec.certify_root_at(&system, &candidate, precision) {
            Outcome::Proven(BoxRoots::One(mut root)) => {
                // The lifted equations retain every coefficient dependency.
                // Their exact selectors bind each non-functional choice to
                // its original root; the remaining arithmetic is functional.
                // Only this certified correspondence admits original f/g.
                root.original_equations = Some(equations.to_vec().into());
                if let Some(axes) = reconstruct {
                    let mut coordinates = Vec::with_capacity(axes.len());
                    for &axis in axes {
                        let value = dec.reconstruct_coordinate(&root, axis)?.ok_or_else(|| {
                            Refusal::new(
                                "numbers.coordinate-reconstruction-unproven",
                                format!(
                                    "coordinate {axis} has no unique proved projection in the certified root box"
                                ),
                            )
                        })?;
                        coordinates.push(value);
                    }
                    return Ok(coordinates);
                }
                let root = Arc::new(root);
                Ok((0..bounds.len())
                    .map(|axis| Real::from_shared_joint_root_coordinate(root.clone(), axis))
                    .collect())
            }
            Outcome::Proven(BoxRoots::None) | Outcome::ProvenCoincident(_) => Err(refused()),
            Outcome::NotExact(n) => Err(Refusal::from(n)),
            Outcome::Refused(refusal) => Err(refusal),
        }
    }

    /// Monic greatest common divisor of two univariate polynomials over the
    /// shared exact algebraic coefficient field. Both zero inputs give zero.
    /// Coefficient zeros and every leading divisor are certified, not inferred
    /// from expression keys or interval overlap. This does not eliminate the
    /// coefficients into a larger rational polynomial.
    pub fn gcd(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        assert_eq!(self.nvars(), 1);
        assert_eq!(other.nvars(), 1);
        dec.check_deadline()?;
        let mut a = univariate_coefficients(self, dec)?;
        let mut b = univariate_coefficients(other, dec)?;
        while let Some((&degree, _)) = b.last_key_value() {
            dec.check_deadline()?;
            if degree == 0 {
                return Ok(Self::constant(1, Real::one()));
            }
            let remainder = univariate_remainder(a, &b, dec)?;
            a = b;
            b = remainder;
        }
        let Some((&degree, leading)) = a.last_key_value() else {
            return Ok(Self::zero(1));
        };
        let leading = leading.clone();
        let mut terms = Vec::with_capacity(a.len());
        for (power, coefficient) in a {
            dec.check_deadline()?;
            terms.push((
                vec![power],
                if power == degree {
                    Real::one()
                } else {
                    &coefficient / &leading
                },
            ));
        }
        Ok(Self::from_terms(1, terms))
    }

    /// Certifies whether every coefficient is zero.
    pub fn is_zero(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.total_degree(dec)?.is_none())
    }

    /// Certifies polynomial equality, rather than comparing expression identities.
    pub fn equals(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        assert_eq!(self.nvars(), other.nvars());
        if let (Some(a), Some(b)) = (self.rational(), other.rational()) {
            return Ok(a == b);
        }
        (self - other).is_zero(dec)
    }

    /// Sylvester resultant after certifying all coefficient zeros and degrees.
    /// Rational inputs retain the existing fraction-free determinant algorithm.
    pub fn resultant(&self, other: &Self, var: usize, dec: &mut Decider) -> Result<Self, Refusal> {
        assert_eq!(self.nvars(), other.nvars());
        assert!(var < self.nvars());
        if let (Some(a), Some(b)) = (self.rational(), other.rational()) {
            return Ok(a.resultant(b, var).into());
        }
        let mut clean = |polynomial: &Self| -> Result<Self, Refusal> {
            let mut terms = Vec::new();
            for (m, coefficient) in polynomial.terms() {
                if nonzero(&coefficient, dec)? {
                    terms.push((m.exponents().to_vec(), coefficient));
                }
            }
            Ok(Self::from_terms(polynomial.nvars(), terms))
        };
        let left = clean(self)?;
        let right = clean(other)?;
        dec.real_poly_resultant(&left, &right, var)
    }

    pub(crate) fn stored_len(&self) -> usize {
        match &self.coefficients {
            Coefficients::Rational(p) => p.term_count(),
            Coefficients::Algebraic { terms, .. } => terms.len(),
        }
    }

    /// Safe construction-key agreement. Different keys do not prove different
    /// polynomials; use `equals` for a mathematical equality decision.
    pub fn same_definition(&self, other: &Self) -> bool {
        self.match_definition(other, |a, b| a.same_definition(b))
    }
    /// Matches storage and monomials, delegating coefficient graph agreement.
    pub(crate) fn match_definition<'a>(
        &'a self,
        other: &'a Self,
        mut coefficient: impl FnMut(&'a Real, &'a Real) -> bool,
    ) -> bool {
        match (&self.coefficients, &other.coefficients) {
            (Coefficients::Rational(a), Coefficients::Rational(b)) => a == b,
            (
                Coefficients::Algebraic {
                    nvars: a,
                    terms: left,
                },
                Coefficients::Algebraic {
                    nvars: b,
                    terms: right,
                },
            ) => {
                a == b
                    && left.len() == right.len()
                    && left
                        .iter()
                        .zip(right)
                        .all(|((m, a), (n, b))| m == n && coefficient(a, b))
            }
            _ => false,
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        assert_eq!(self.nvars(), rhs.nvars());
        if let (Some(a), Some(b)) = (self.rational(), rhs.rational()) {
            return (a + b).into();
        }
        Self::from_terms(
            self.nvars(),
            self.terms()
                .into_iter()
                .chain(rhs.terms())
                .map(|(m, c)| (m.exponents().to_vec(), c))
                .collect(),
        )
    }

    fn sub(&self, rhs: &Self) -> Self {
        if let (Some(a), Some(b)) = (self.rational(), rhs.rational()) {
            return (a - b).into();
        }
        self + &(-rhs)
    }

    fn mul(&self, rhs: &Self) -> Self {
        assert_eq!(self.nvars(), rhs.nvars());
        if let (Some(a), Some(b)) = (self.rational(), rhs.rational()) {
            return (a * b).into();
        }
        let left = self.terms();
        let right = rhs.terms();
        let mut terms = Vec::with_capacity(
            left.len()
                .checked_mul(right.len())
                .expect("polynomial term count overflow"),
        );
        for (a, x) in &left {
            for (b, y) in &right {
                let exponents = a
                    .exponents()
                    .iter()
                    .zip(b.exponents())
                    .map(|(a, b)| a.checked_add(*b).expect("monomial exponent overflow"))
                    .collect();
                terms.push((exponents, x * y));
            }
        }
        Self::from_terms(self.nvars(), terms)
    }

    fn neg(&self) -> Self {
        if let Some(p) = self.rational() {
            return (-p).into();
        }
        self.scale(&Real::from_i64(-1))
    }
}

type Univariate = BTreeMap<u32, Real>;

fn univariate_coefficients(
    polynomial: &RealPoly,
    dec: &mut Decider,
) -> Result<Univariate, Refusal> {
    let mut result = BTreeMap::new();
    for (monomial, coefficient) in polynomial.terms() {
        dec.check_deadline()?;
        // Validate every coefficient before a constant gcd could short-circuit.
        if nonzero(&coefficient, dec)? {
            result.insert(monomial.exponents()[0], coefficient);
        }
    }
    Ok(result)
}

fn univariate_remainder(
    mut a: Univariate,
    b: &Univariate,
    dec: &mut Decider,
) -> Result<Univariate, Refusal> {
    let (&degree, leading) = b.last_key_value().expect("nonzero certified divisor");
    while let Some((&power, coefficient)) = a.last_key_value() {
        dec.check_deadline()?;
        if power < degree {
            break;
        }
        let shift = power - degree;
        let scale = coefficient / leading;
        // The leading term cancels by field arithmetic and the certified
        // nonzero divisor. Do not build and re-prove that known-zero expression.
        a.pop_last();
        for (&power, coefficient) in b.range(..degree) {
            dec.check_deadline()?;
            let term = &scale * coefficient;
            a.entry(power + shift)
                .and_modify(|old| *old = &*old - &term)
                .or_insert_with(|| -term);
        }
        while let Some((_, coefficient)) = a.last_key_value() {
            if nonzero(coefficient, dec)? {
                break;
            }
            a.pop_last();
        }
    }
    Ok(a)
}

// The sign at `at` of a polynomial over a quadratic tower.
fn tree_sign(tree: &crate::tower::SignTree, at: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match tree {
        crate::tower::SignTree::Rational(poly) if poly.is_zero() => Ok(Sign::Zero),
        crate::tower::SignTree::Rational(poly) => strict_sign(&poly.eval_real(at), dec),
        crate::tower::SignTree::Pair(parts) => {
            let [a, b, norm] = &**parts;
            let a = tree_sign(a, at, dec)?;
            let b = tree_sign(b, at, dec)?;
            if a == Sign::Zero || a == b {
                return Ok(b);
            }
            Ok(match tree_sign(norm, at, dec)? {
                Sign::Positive => a,
                Sign::Negative => b,
                Sign::Zero => Sign::Zero,
            })
        }
    }
}

fn strict_sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

fn nonzero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(_) => Ok(true),
        Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

fn root_values(outcome: Outcome<Vec<crate::RealRoot>>) -> Result<Vec<Real>, Refusal> {
    match outcome {
        Outcome::Proven(roots) => Ok(roots
            .into_iter()
            .map(|root| Real::from_algebraic(root.value))
            .collect()),
        Outcome::Refused(r) => Err(r),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        // Root isolation does not answer with a coincidence; refuse rather
        // than panic if it ever does.
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.root-isolation",
            "root isolation answered a coincidence instead of a root list",
        )),
    }
}

impl From<MPoly> for RealPoly {
    fn from(value: MPoly) -> Self {
        Self {
            coefficients: Coefficients::Rational(value),
        }
    }
}

binops!(RealPoly; Add::add, Sub::sub, Mul::mul);
negop!(RealPoly);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Coincidence, ZeroProof};

    // Root isolation answers with a list or a refusal; a coincidence must
    // refuse, never panic.
    #[test]
    fn a_coincident_root_list_refuses_instead_of_panicking() {
        let refusal = root_values(Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ExactArithmetic,
        }))
        .expect_err("a coincidence is not a root list");
        assert_eq!(refusal.code, "numbers.root-isolation");
    }
}
