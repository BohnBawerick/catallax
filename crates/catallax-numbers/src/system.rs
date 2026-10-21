//! Polynomial systems and rational Krawczyk certificates.

use crate::decider::Work;
use crate::{Ball, MPoly, Rational, Real, UPoly};
use std::sync::Arc;

/// Square system of at least one rational polynomial equation.
#[derive(Clone, Debug)]
pub struct PolySystem {
    equations: Vec<MPoly>,
}
impl PolySystem {
    /// Requires `n >= 1` equations, all in `n` variables.
    pub fn new(equations: Vec<MPoly>) -> Option<Self> {
        let n = equations.len();
        if n != 0 && equations.iter().all(|p| p.nvars() == n) {
            Some(Self { equations })
        } else {
            None
        }
    }
    /// Equations in input order.
    pub fn equations(&self) -> &[MPoly] {
        &self.equations
    }
    /// Number of variables and equations.
    pub fn dimension(&self) -> usize {
        self.equations.len()
    }
}
/// Closed rational box.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RootBox {
    lo: Vec<Rational>,
    hi: Vec<Rational>,
}
impl RootBox {
    /// Requires equal nonzero lengths and componentwise ordered corners.
    pub fn new(lo: Vec<Rational>, hi: Vec<Rational>) -> Option<Self> {
        if !lo.is_empty() && lo.len() == hi.len() && lo.iter().zip(&hi).all(|(a, b)| a <= b) {
            Some(Self { lo, hi })
        } else {
            None
        }
    }
    /// Lower corner.
    pub fn lo(&self) -> &[Rational] {
        &self.lo
    }
    /// Upper corner.
    pub fn hi(&self) -> &[Rational] {
        &self.hi
    }
    /// Number of axes.
    pub fn dimension(&self) -> usize {
        self.lo.len()
    }
}
/// System and box certified to contain exactly one nonsingular solution.
#[derive(Clone, Debug)]
pub struct IsolatedRoot {
    pub(crate) system: Arc<PolySystem>,
    pub(crate) root_box: RootBox,
}
impl IsolatedRoot {
    /// Original system.
    pub fn system(&self) -> &Arc<PolySystem> {
        &self.system
    }
    /// Certified closed box.
    pub fn root_box(&self) -> &RootBox {
        &self.root_box
    }
    /// Exact coordinate leaf; panics for an invalid axis.
    pub fn coordinate(&self, axis: usize) -> Real {
        Real::from_root_coordinate(self, axis)
    }
    /// Proven coordinate enclosures at the requested working precision.
    pub fn enclose(&self, precision: u32) -> Vec<Ball> {
        self.enclose_impl(precision, None)
            .expect("unbudgeted root enclosure")
    }
    pub(crate) fn enclose_budgeted(
        &self,
        precision: u32,
        work: &mut Work,
    ) -> Result<Vec<Ball>, &'static str> {
        self.enclose_impl(precision, Some(work))
    }
    fn enclose_impl(
        &self,
        precision: u32,
        mut work: Option<&mut Work>,
    ) -> Result<Vec<Ball>, &'static str> {
        let mut bounds = self.root_box.clone();
        let p = precision.max(1);
        let g = work.as_ref().map_or(p.saturating_add(32), |w| {
            p.saturating_add(32).min(w.policy.max_precision_bits)
        });
        let target = Rational::one().mul_pow2(-i64::from(p));
        loop {
            if bounds
                .lo
                .iter()
                .zip(&bounds.hi)
                .all(|(a, b)| b - a <= target)
            {
                break;
            }
            if let Some(work) = work.as_deref_mut() {
                work.system(&self.system, g)?;
            }
            let Some((next, _, _)) = krawczyk(&self.system, &bounds, g) else {
                break;
            };
            let Some(next) = intersection(&bounds, &next) else {
                break;
            };
            if next == bounds {
                break;
            }
            bounds = next;
        }
        bounds
            .lo
            .iter()
            .zip(&bounds.hi)
            .map(|(a, b)| {
                if let Some(work) = work.as_deref_mut() {
                    work.ball(p, true)?;
                }
                Ok(Ball::bounds(a.clone(), b.clone(), p))
            })
            .collect()
    }
}
/// Certified contents of a candidate box.
#[derive(Clone, Debug)]
pub enum BoxRoots {
    /// No solution.
    None,
    /// Exactly one nonsingular solution.
    One(IsolatedRoot),
}

pub(crate) fn certify(
    system: &Arc<PolySystem>,
    candidate: &RootBox,
    p: u32,
    work: &mut Work,
) -> Result<Option<BoxRoots>, &'static str> {
    assert_eq!(system.dimension(), candidate.dimension());
    let mut bounds = candidate.clone();
    loop {
        work.system(system, p)?;
        let balls = box_balls(&bounds, p);
        if system
            .equations
            .iter()
            .any(|f| f.eval_ball(&balls).sign().is_some())
        {
            return Ok(Some(BoxRoots::None));
        }
        let Some((image, inclusion, contraction)) = krawczyk(system, &bounds, p) else {
            return Ok(None);
        };
        if inclusion && contraction {
            return Ok(Some(BoxRoots::One(IsolatedRoot {
                system: system.clone(),
                root_box: bounds,
            })));
        }
        let Some(next) = intersection(&bounds, &image) else {
            return Ok(Some(BoxRoots::None));
        };
        if next == bounds {
            return Ok(None);
        }
        bounds = next;
    }
}
fn box_balls(bounds: &RootBox, p: u32) -> Vec<Ball> {
    bounds
        .lo
        .iter()
        .zip(&bounds.hi)
        .map(|(a, b)| Ball::bounds(a.clone(), b.clone(), p))
        .collect()
}
fn intersection(a: &RootBox, b: &RootBox) -> Option<RootBox> {
    RootBox::new(
        a.lo.iter()
            .zip(&b.lo)
            .map(|(a, b)| a.clone().max(b.clone()))
            .collect(),
        a.hi.iter()
            .zip(&b.hi)
            .map(|(a, b)| a.clone().min(b.clone()))
            .collect(),
    )
}
fn inverse(matrix: Vec<Vec<Rational>>) -> Option<Vec<Vec<Rational>>> {
    let n = matrix.len();
    let mut a = matrix;
    let mut b = vec![vec![Rational::zero(); n]; n];
    for (i, row) in b.iter_mut().enumerate() {
        row[i] = Rational::one();
    }
    for k in 0..n {
        let pivot = (k..n).find(|i| !a[*i][k].is_zero())?;
        a.swap(k, pivot);
        b.swap(k, pivot);
        let d = a[k][k].clone();
        for j in 0..n {
            a[k][j] = &a[k][j] / &d;
            b[k][j] = &b[k][j] / &d;
        }
        for i in 0..n {
            if i == k {
                continue;
            }
            let c = a[i][k].clone();
            for j in 0..n {
                a[i][j] = &a[i][j] - &(&c * &a[k][j]);
                b[i][j] = &b[i][j] - &(&c * &b[k][j]);
            }
        }
    }
    Some(b)
}
// K(X) = m-C F(m) + (I-C J(X))(X-m), with C the exact inverse of J(m).
// Every root in X lies in K(X), so disjoint boxes prove absence of roots.
// T(x) = x-C F(x) lies in K(X) for every x in X. Closed containment
// K(X) ⊆ X makes T a self-map of the complete box, even if an axis has
// collapsed or a root is on a face. Together with the separately required
// strict bound ||I-C J(X)|| < 1, Banach gives a unique fixed point.
// Since C is invertible, it is a root; the same norm bound makes J invertible.
// Strict containment would lose exact pinned coordinates after intersection.
// Closed containment avoids inflation outside the candidate. Intersections
// retain every original root, including boundary roots, so exclusion is unchanged.
fn krawczyk(system: &PolySystem, bounds: &RootBox, p: u32) -> Option<(RootBox, bool, bool)> {
    let n = system.dimension();
    let x = box_balls(bounds, p);
    let mid: Vec<_> = bounds
        .lo
        .iter()
        .zip(&bounds.hi)
        .map(|(a, b)| (a + b).mul_pow2(-1))
        .collect();
    let jac: Vec<Vec<_>> = system
        .equations
        .iter()
        .map(|f| (0..n).map(|j| f.partial(j)).collect())
        .collect();
    let center_jac = jac
        .iter()
        .map(|row| row.iter().map(|f| f.eval(&mid)).collect())
        .collect();
    let c = inverse(center_jac)?;
    let zero = Ball::from_rational(&Rational::zero(), p);
    let one = Ball::from_rational(&Rational::one(), p);
    let f: Vec<_> = system
        .equations
        .iter()
        .map(|f| Ball::from_rational(&f.eval(&mid), p))
        .collect();
    let jac: Vec<Vec<_>> = jac
        .iter()
        .map(|row| row.iter().map(|f| f.eval_ball(&x)).collect())
        .collect();
    let mut lo = Vec::new();
    let mut hi = Vec::new();
    let mut contraction = true;
    for i in 0..n {
        let mut value = Ball::from_rational(&mid[i], p);
        for k in 0..n {
            value = &value - &(&Ball::from_rational(&c[i][k], p) * &f[k]);
        }
        let mut norm = Rational::zero();
        for j in 0..n {
            let mut r = if i == j { one.clone() } else { zero.clone() };
            for k in 0..n {
                r = &r - &(&Ball::from_rational(&c[i][k], p) * &jac[k][j]);
            }
            norm = &norm + &r.magnitude();
            value = &value + &(&r * &(&x[j] - &Ball::from_rational(&mid[j], p)));
        }
        contraction &= norm < Rational::one();
        lo.push(value.lo());
        hi.push(value.hi());
    }
    let inclusion = lo.iter().zip(&bounds.lo).all(|(a, b)| a >= b)
        && hi.iter().zip(&bounds.hi).all(|(a, b)| a <= b);
    Some((RootBox::new(lo, hi).unwrap(), inclusion, contraction))
}

pub(crate) fn coordinate_polynomial(
    system: &PolySystem,
    axis: usize,
    work: &mut Work,
) -> Result<Option<UPoly>, &'static str> {
    let n = system.dimension();
    let mut polys = system.equations.clone();
    // Linear combinations often expose affine constraints before elimination.
    for i in 0..system.equations.len() {
        for j in i + 1..system.equations.len() {
            let p = &system.equations[i] - &system.equations[j];
            insert_poly(&mut polys, p);
        }
    }
    for var in (0..n).filter(|v| *v != axis) {
        if let Some(p) = polys
            .iter()
            .filter_map(|p| p.to_upoly(axis))
            .find(|p| !p.is_zero() && p.degree().is_some_and(|d| d > 0))
        {
            return Ok(Some(p));
        }
        let linear = polys.iter().find_map(|f| {
            let c = f.coefficients_in(var);
            if c.len() != 2 {
                return None;
            }
            let terms = c[1].terms();
            if terms.len() == 1 && terms[0].0.degree() == 0 {
                Some((-&c[0]).scale(&terms[0].1.recip().unwrap()))
            } else {
                None
            }
        });
        if let Some(value) = linear {
            let mut next = Vec::new();
            for f in polys {
                work.mpoly(&f, &value)?;
                insert_poly(&mut next, f.substitute(var, &value));
            }
            polys = next;
        } else {
            let (independent, dependent): (Vec<_>, Vec<_>) = polys
                .into_iter()
                .partition(|f| f.coefficients_in(var).len() == 1);
            let mut next = independent;
            for i in 0..dependent.len() {
                for j in i + 1..dependent.len() {
                    work.mpoly(&dependent[i], &dependent[j])?;
                    insert_poly(
                        &mut next,
                        dependent[i]
                            .resultant_checked(&dependent[j], var, &|| work.check_deadline())?,
                    );
                }
            }
            polys = next;
        }
    }
    Ok(polys
        .iter()
        .filter_map(|p| p.to_upoly(axis))
        .find(|p| p.degree().is_some_and(|d| d > 0)))
}
fn insert_poly(polys: &mut Vec<MPoly>, p: MPoly) {
    if p.is_zero() {
        return;
    }
    let p = p.scale(&p.terms().last().unwrap().1.recip().unwrap());
    if !polys.contains(&p) {
        polys.push(p);
    }
}
