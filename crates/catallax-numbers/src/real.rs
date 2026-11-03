//! Immutable definitions with eager hardware bounds and demand-evaluated construction metadata.

use crate::{Algebraic, Ball, Interval, IsolatedRoot, Rational};
use std::sync::{Arc, LazyLock, OnceLock};

/// Exact real expression. Cloning shares its immutable definition and caches.
#[derive(Clone, Debug)]
pub struct Real {
    pub(crate) node: Arc<Node>,
}
#[derive(Clone, Debug)]
pub(crate) enum Expr {
    Rational(Rational),
    Algebraic(Algebraic),
    Coordinate(Arc<IsolatedRoot>, usize),
    Add(Real, Real),
    Sub(Real, Real),
    Mul(Real, Real),
    Div(Real, Real),
    Neg(Real),
    Sqrt(Real),
    Square(Real),
    Pow(Real, u32),
}
pub(crate) struct Node {
    pub expr: Expr,
    pub approx: Interval,
    coarse: OnceLock<Option<Ball>>,
    rational: OnceLock<Option<Rational>>,
    pub fast_defined: bool,
    // Construction metadata for routing an optional proof method: bit 0 means
    // an algebraic leaf outside the native quadratic/quartic routes; bit 1
    // means a pinned system coordinate. This
    // records definition kinds only and never proves a value or a domain.
    pub identity_atoms: u8,
    rational_possible: bool,
}

impl Node {
    fn new(
        expr: Expr,
        approx: Interval,
        fast_defined: bool,
        frozen: Option<(Option<Rational>, Option<Ball>)>,
    ) -> Self {
        let (identity_atoms, rational_possible) = match &expr {
            Expr::Algebraic(a) => (
                u8::from(matches!(&a.repr,
                    crate::algebraic::AlgebraicRepr::Root(root)
                        if root.polynomial.degree().is_some_and(|degree| degree == 3 || degree > 4))),
                frozen
                    .as_ref()
                    .is_some_and(|(rational, _)| rational.is_some()),
            ),
            Expr::Coordinate(..) => (2, false),
            Expr::Sqrt(a) => (a.node.identity_atoms, false),
            Expr::Neg(a) | Expr::Square(a) | Expr::Pow(a, _) => {
                (a.node.identity_atoms, a.node.rational_possible)
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => (
                a.node.identity_atoms | b.node.identity_atoms,
                a.node.rational_possible && b.node.rational_possible,
            ),
            Expr::Rational(_) => (0, true),
        };
        let (rational, coarse) = match frozen {
            Some((rational, coarse)) => (OnceLock::from(rational), OnceLock::from(coarse)),
            None => (OnceLock::new(), OnceLock::new()),
        };
        Self {
            expr,
            approx,
            coarse,
            rational,
            fast_defined,
            identity_atoms,
            rational_possible,
        }
    }

    fn rational_ready(&self) -> bool {
        !self.rational_possible
            || matches!(&self.expr, Expr::Rational(_))
            || self.rational.get().is_some()
    }

    fn rational_value(&self) -> Option<&Rational> {
        if !self.rational_possible {
            return None;
        }
        if let Expr::Rational(rational) = &self.expr {
            return Some(rational);
        }
        self.rational
            .get()
            .expect("rational metadata dependencies must be initialized")
            .as_ref()
    }

    fn coarse_value(&self) -> Option<&Ball> {
        self.coarse
            .get()
            .expect("coarse metadata dependencies must be initialized")
            .as_ref()
    }

    fn rational_ref(&self) -> Option<&Rational> {
        CacheTask::Rational(self).ensure();
        self.rational_value()
    }

    pub(crate) fn coarse_ref(&self) -> Option<&Ball> {
        CacheTask::Coarse(self).ensure();
        self.coarse_value()
    }

    fn calculate_rational(&self) -> Option<Rational> {
        match &self.expr {
            Expr::Add(a, b) => a
                .node
                .rational_value()
                .and_then(|a| b.node.rational_value().map(|b| a + b)),
            Expr::Sub(a, b) => a
                .node
                .rational_value()
                .and_then(|a| b.node.rational_value().map(|b| a - b)),
            Expr::Mul(a, b) => a
                .node
                .rational_value()
                .and_then(|a| b.node.rational_value().map(|b| a * b)),
            Expr::Div(a, b) => a
                .node
                .rational_value()
                .and_then(|a| b.node.rational_value().and_then(|b| a.checked_div(b))),
            Expr::Neg(a) => a.node.rational_value().map(|a| -a),
            Expr::Square(a) => a.node.rational_value().map(|a| a * a),
            Expr::Pow(a, exp) => a.node.rational_value().map(|a| {
                let mut value = Rational::one();
                let mut n = *exp;
                if n != 0 {
                    let mut a = a.clone();
                    while n != 0 {
                        if n & 1 != 0 {
                            value = &value * &a;
                        }
                        n >>= 1;
                        if n != 0 {
                            a = &a * &a;
                        }
                    }
                }
                value
            }),
            Expr::Rational(_) | Expr::Algebraic(_) | Expr::Coordinate(..) | Expr::Sqrt(_) => {
                unreachable!("leaf and nonrational metadata is already available")
            }
        }
    }

    fn calculate_coarse(&self) -> Option<Ball> {
        if let Some(rational) = self.rational_value() {
            return Some(Ball::from_rational(rational, 1));
        }
        let coarse = match &self.expr {
            Expr::Add(a, b) => a
                .node
                .coarse_value()
                .and_then(|a| b.node.coarse_value().map(|b| a + b)),
            Expr::Sub(a, b) => a
                .node
                .coarse_value()
                .and_then(|a| b.node.coarse_value().map(|b| a - b)),
            Expr::Mul(a, b) => a
                .node
                .coarse_value()
                .and_then(|a| b.node.coarse_value().map(|b| a * b)),
            Expr::Div(a, b) => a
                .node
                .coarse_value()
                .and_then(|a| b.node.coarse_value().and_then(|b| a.checked_div(b))),
            Expr::Neg(a) => a.node.coarse_value().map(|a| -a),
            Expr::Sqrt(a) => a.node.coarse_value().and_then(Ball::sqrt),
            Expr::Square(a) => a.node.coarse_value().map(Ball::square),
            Expr::Pow(a, exp) => a.node.coarse_value().map(|a| a.pow_u32(*exp)),
            Expr::Rational(_) | Expr::Algebraic(_) | Expr::Coordinate(..) => {
                unreachable!("leaf construction bounds are frozen")
            }
        };
        coarse.or_else(|| {
            if self.fast_defined && self.approx.is_finite() {
                Some(Ball::bounds(
                    Rational::from_f64_exact(self.approx.lo()).unwrap(),
                    Rational::from_f64_exact(self.approx.hi()).unwrap(),
                    1,
                ))
            } else {
                None
            }
        })
    }
}

impl std::fmt::Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Node")
            .field("expr", &self.expr)
            .field("approx", &self.approx)
            .field("coarse", &self.coarse_ref())
            .field("rational", &self.rational_ref())
            .field("fast_defined", &self.fast_defined)
            .field("identity_atoms", &self.identity_atoms)
            .finish()
    }
}

#[derive(Clone, Copy)]
enum CacheTask<'a> {
    Rational(&'a Node),
    Coarse(&'a Node),
}

impl<'a> CacheTask<'a> {
    fn node(self) -> &'a Node {
        match self {
            Self::Rational(node) | Self::Coarse(node) => node,
        }
    }

    fn child(self, value: &'a Real) -> Self {
        match self {
            Self::Rational(_) => Self::Rational(&value.node),
            Self::Coarse(_) => Self::Coarse(&value.node),
        }
    }

    fn ready(self) -> bool {
        match self {
            Self::Rational(node) => node.rational_ready(),
            Self::Coarse(node) => node.coarse.get().is_some(),
        }
    }

    fn present(self) -> bool {
        match self {
            Self::Rational(node) => node.rational_value().is_some(),
            Self::Coarse(node) => node.coarse_value().is_some(),
        }
    }

    fn dependency(self) -> Option<Self> {
        let node = self.node();
        if let Self::Coarse(_) = self {
            if !node.rational_ready() {
                return Some(Self::Rational(node));
            }
            if node.rational_value().is_some() {
                return None;
            }
        }
        let (first, second) = match &node.expr {
            Expr::Neg(a) | Expr::Sqrt(a) | Expr::Square(a) | Expr::Pow(a, _) => (a, None),
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => (a, Some(b)),
            Expr::Rational(_) | Expr::Algebraic(_) | Expr::Coordinate(..) => return None,
        };
        let first = self.child(first);
        if !first.ready() {
            return Some(first);
        }
        if !first.present() {
            return None;
        }
        let second = self.child(second?);
        (!second.ready()).then_some(second)
    }

    fn materialize(self) {
        match self {
            Self::Rational(node) => {
                node.rational.get_or_init(|| node.calculate_rational());
            }
            Self::Coarse(node) => {
                node.coarse.get_or_init(|| node.calculate_coarse());
            }
        }
    }

    fn ensure(self) {
        if self.ready() {
            return;
        }
        let Some(dependency) = self.dependency() else {
            self.materialize();
            return;
        };
        // Resolve dependencies before locking: no recursion or nested cache locks.
        let mut pending = vec![self, dependency];
        while let Some(&task) = pending.last() {
            if task.ready() {
                pending.pop();
            } else if let Some(dependency) = task.dependency() {
                pending.push(dependency);
            } else {
                task.materialize();
                pending.pop();
            }
        }
    }
}

impl Real {
    /// Zero.
    pub fn zero() -> Self {
        Self::from_i64(0)
    }
    /// One.
    pub fn one() -> Self {
        Self::from_i64(1)
    }
    /// Exact rational leaf.
    pub fn from_rational(r: Rational) -> Self {
        if r.is_integer()
            && let Some(shared) = r.numerator().to_i64().and_then(Self::shared_integer)
        {
            return shared;
        }
        Self::rational_leaf(r)
    }
    fn shared_integer(value: i64) -> Option<Self> {
        // Intern only literal values. Arithmetic definitions and their domain
        // checks remain intact, even when their rational result is -1, 0 or 1.
        static ZERO: LazyLock<Real> = LazyLock::new(|| Real::rational_leaf(Rational::zero()));
        static ONE: LazyLock<Real> = LazyLock::new(|| Real::rational_leaf(Rational::one()));
        static NEGATIVE_ONE: LazyLock<Real> =
            LazyLock::new(|| Real::rational_leaf(Rational::from_i64(-1)));
        match value {
            0 => Some(ZERO.clone()),
            1 => Some(ONE.clone()),
            -1 => Some(NEGATIVE_ONE.clone()),
            _ => None,
        }
    }
    fn rational_leaf(r: Rational) -> Self {
        let approx = Interval::from_rational(&r);
        Self::make(Expr::Rational(r), approx, true)
    }
    /// Exact integer leaf.
    pub fn from_i64(v: i64) -> Self {
        Self::shared_integer(v).unwrap_or_else(|| Self::rational_leaf(Rational::from_i64(v)))
    }
    /// Exact algebraic leaf.
    pub fn from_algebraic(a: Algebraic) -> Self {
        let (lo, hi) = a.isolating_interval();
        let approx = Interval::from_rational(&lo).hull(&Interval::from_rational(&hi));
        let coarse = || Ball::bounds(lo, hi, 1);
        let rational = a.as_rational().cloned();
        Self::frozen_leaf(Expr::Algebraic(a), approx, rational, coarse)
    }
    /// Coordinate of a certified system root; panics for an invalid axis.
    pub fn from_root_coordinate(root: &IsolatedRoot, axis: usize) -> Self {
        assert!(axis < root.root_box().dimension());
        let lo = root.root_box().lo()[axis].clone();
        let hi = root.root_box().hi()[axis].clone();
        let approx = Interval::from_rational(&lo).hull(&Interval::from_rational(&hi));
        Self::frozen_leaf(
            Expr::Coordinate(Arc::new(root.clone()), axis),
            approx,
            None,
            || Ball::bounds(lo, hi, 1),
        )
    }
    /// Lazy square root. Domain validation belongs to the decider.
    pub fn sqrt(&self) -> Self {
        let approx = self.node.approx.sqrt().unwrap_or(Interval::entire());
        Self::make(
            Expr::Sqrt(self.clone()),
            approx,
            self.node.fast_defined && self.node.approx.lo() >= 0.0,
        )
    }
    /// Lazy square, retaining the entire operand definition.
    pub fn square(&self) -> Self {
        Self::make(
            Expr::Square(self.clone()),
            self.node.approx.square(),
            self.node.fast_defined,
        )
    }
    /// Lazy non-negative integer power, retaining definedness even for exponent zero.
    pub fn pow_u32(&self, exp: u32) -> Self {
        let mut a = self.node.approx;
        let mut r = Interval::point(1.0).unwrap();
        let mut n = exp;
        while n != 0 {
            if n & 1 != 0 {
                r = r * a;
            }
            n >>= 1;
            if n != 0 {
                a = a.square();
            }
        }
        Self::make(Expr::Pow(self.clone(), exp), r, self.node.fast_defined)
    }
    /// Exact value of a structurally rational expression; does not refine.
    pub fn as_rational(&self) -> Option<Rational> {
        self.rational_ref().cloned()
    }
    /// Borrows a syntactically rational value, computing exact metadata without refinement.
    pub fn rational_ref(&self) -> Option<&Rational> {
        self.node.rational_ref()
    }

    /// Whether immutable construction definitions agree. Different definitions
    /// do not prove different values. This does not preserve expression-node sharing.
    pub fn same_definition(&self, other: &Self) -> bool {
        if let Some(a) = self.rational_ref()
            && let Some(b) = other.rational_ref()
        {
            return a == b;
        }
        self.construction_cmp(other).is_eq()
    }
    /// Whether this definition is the negation of `other`, without constructing
    /// an expression or refining either value. A false result is inconclusive.
    pub fn same_negated_definition(&self, other: &Self) -> bool {
        if let Some(a) = self.rational_ref()
            && let Some(b) = other.rational_ref()
        {
            return a.d == b.d
                && a.n.0.magnitude() == b.n.0.magnitude()
                && (a.is_zero() || a.sign() != b.sign());
        }
        match (&self.node.expr, &other.node.expr) {
            (Expr::Neg(a), _) => a.same_definition(other),
            (_, Expr::Neg(b)) => self.same_definition(b),
            _ => false,
        }
    }
    /// Whether two values share the same immutable definition node.
    /// Use only for caches whose reuse must preserve expression-DAG sharing;
    /// distinct nodes do not prove distinct mathematical values.
    pub fn same_node(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.node, &other.node)
    }
    /// Cached, conservative double enclosure.
    pub fn approx(&self) -> Interval {
        self.node.approx
    }
    /// Double approximation for display only.
    pub fn to_f64(&self) -> f64 {
        self.rational_ref()
            .map_or_else(|| self.approx().mid(), Rational::to_f64_nearest)
    }

    /// Orders immutable construction records, not their mathematical values.
    pub(crate) fn construction_cmp(&self, other: &Self) -> std::cmp::Ordering {
        use crate::algebraic::AlgebraicRepr;
        use std::cmp::Ordering;
        if Arc::ptr_eq(&self.node, &other.node) {
            return Ordering::Equal;
        }
        let tag = |e: &Expr| match e {
            Expr::Rational(_) => 0,
            Expr::Algebraic(_) => 1,
            Expr::Coordinate(_, _) => 2,
            Expr::Add(_, _) => 3,
            Expr::Sub(_, _) => 4,
            Expr::Mul(_, _) => 5,
            Expr::Div(_, _) => 6,
            Expr::Neg(_) => 7,
            Expr::Sqrt(_) => 8,
            Expr::Square(_) => 9,
            Expr::Pow(_, _) => 10,
        };
        tag(&self.node.expr)
            .cmp(&tag(&other.node.expr))
            .then_with(|| match (&self.node.expr, &other.node.expr) {
                (Expr::Rational(a), Expr::Rational(b)) => a.cmp(b),
                (Expr::Algebraic(a), Expr::Algebraic(b)) => match (&a.repr, &b.repr) {
                    (AlgebraicRepr::Rational(a), AlgebraicRepr::Rational(b)) => a.cmp(b),
                    (AlgebraicRepr::Rational(_), AlgebraicRepr::Root(_)) => Ordering::Less,
                    (AlgebraicRepr::Root(_), AlgebraicRepr::Rational(_)) => Ordering::Greater,
                    (AlgebraicRepr::Root(a), AlgebraicRepr::Root(b)) => a
                        .polynomial
                        .coefficients()
                        .cmp(b.polynomial.coefficients())
                        .then_with(|| a.lo.cmp(&b.lo))
                        .then_with(|| a.hi.cmp(&b.hi)),
                },
                (Expr::Coordinate(a, i), Expr::Coordinate(b, j)) => a
                    .system()
                    .dimension()
                    .cmp(&b.system().dimension())
                    .then_with(|| {
                        a.system()
                            .equations()
                            .iter()
                            .zip(b.system().equations())
                            .map(|(a, b)| a.construction_cmp(b))
                            .find(|o| *o != Ordering::Equal)
                            .unwrap_or(Ordering::Equal)
                    })
                    .then_with(|| a.root_box().lo().cmp(b.root_box().lo()))
                    .then_with(|| a.root_box().hi().cmp(b.root_box().hi()))
                    .then_with(|| i.cmp(j)),
                (Expr::Add(a, b), Expr::Add(c, d))
                | (Expr::Sub(a, b), Expr::Sub(c, d))
                | (Expr::Mul(a, b), Expr::Mul(c, d))
                | (Expr::Div(a, b), Expr::Div(c, d)) => {
                    a.construction_cmp(c).then_with(|| b.construction_cmp(d))
                }
                (Expr::Neg(a), Expr::Neg(b))
                | (Expr::Sqrt(a), Expr::Sqrt(b))
                | (Expr::Square(a), Expr::Square(b)) => a.construction_cmp(b),
                (Expr::Pow(a, n), Expr::Pow(b, m)) => a.construction_cmp(b).then_with(|| n.cmp(m)),
                _ => unreachable!("matching construction tags"),
            })
    }
    fn make(expr: Expr, approx: Interval, fast_defined: bool) -> Self {
        Self {
            node: Arc::new(Node::new(expr, approx, fast_defined, None)),
        }
    }

    fn frozen_leaf(
        expr: Expr,
        approx: Interval,
        rational: Option<Rational>,
        coarse: impl FnOnce() -> Ball,
    ) -> Self {
        // Preserve construction bounds independently of later refinement caches.
        let coarse = rational
            .as_ref()
            .map(|r| Ball::from_rational(r, 1))
            .unwrap_or_else(coarse);
        Self {
            node: Arc::new(Node::new(
                expr,
                approx,
                true,
                Some((rational, Some(coarse))),
            )),
        }
    }
    fn add(&self, rhs: &Self) -> Self {
        Self::make(
            Expr::Add(self.clone(), rhs.clone()),
            self.approx() + rhs.approx(),
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn sub(&self, rhs: &Self) -> Self {
        Self::make(
            Expr::Sub(self.clone(), rhs.clone()),
            self.approx() - rhs.approx(),
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn mul(&self, rhs: &Self) -> Self {
        Self::make(
            Expr::Mul(self.clone(), rhs.clone()),
            self.approx() * rhs.approx(),
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn div(&self, rhs: &Self) -> Self {
        Self::make(
            Expr::Div(self.clone(), rhs.clone()),
            self.approx() / rhs.approx(),
            self.node.fast_defined && rhs.node.fast_defined && rhs.approx().sign().is_some(),
        )
    }
    fn neg(&self) -> Self {
        Self::make(
            Expr::Neg(self.clone()),
            -self.approx(),
            self.node.fast_defined,
        )
    }
}
binops!(Real; Add::add, Sub::sub, Mul::mul, Div::div);
negop!(Real);
impl From<Rational> for Real {
    fn from(r: Rational) -> Self {
        Self::from_rational(r)
    }
}
impl From<i64> for Real {
    fn from(v: i64) -> Self {
        Self::from_i64(v)
    }
}
impl From<Algebraic> for Real {
    fn from(a: Algebraic) -> Self {
        Self::from_algebraic(a)
    }
}
