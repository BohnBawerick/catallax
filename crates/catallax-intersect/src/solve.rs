use crate::math::*;
use catallax_numbers::{Decider, MPoly, Outcome, Rational, Real, RealPoly, Refusal, UPoly};
use catallax_shapes::Point3;
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet, VecDeque},
    sync::Arc,
};

// Lexicographic elimination order y > z > t. The separating coordinate is
// t = x + a*y + b*z. All reductions are over the exact rational field.
#[derive(Clone, Debug, PartialEq, Eq)]
struct P(BTreeMap<[u32; 3], Rational>);
impl P {
    fn zero() -> Self {
        Self(BTreeMap::new())
    }
    fn constant(c: Rational) -> Self {
        let mut p = Self::zero();
        p.add([0; 3], c);
        p
    }
    fn var(k: usize) -> Self {
        let mut e = [0; 3];
        e[k] = 1;
        let mut p = Self::zero();
        p.add(e, r(1));
        p
    }
    fn add(&mut self, e: [u32; 3], c: Rational) {
        let v = self.0.get(&e).map_or_else(|| c.clone(), |v| v + &c);
        if v.is_zero() {
            self.0.remove(&e);
        } else {
            self.0.insert(e, v);
        }
    }
    fn leading(&self) -> Option<([u32; 3], Rational)> {
        self.0.last_key_value().map(|(e, c)| (*e, c.clone()))
    }
    fn monic(mut self) -> Self {
        if let Some((_, c)) = self.leading() {
            for v in self.0.values_mut() {
                *v = &*v / &c;
            }
        }
        self
    }
    fn multiple(&self, e: [u32; 3], c: &Rational) -> Self {
        let mut p = Self::zero();
        for (a, v) in &self.0 {
            p.add(std::array::from_fn(|i| a[i] + e[i]), v * c);
        }
        p
    }
    fn minus(&self, b: &Self) -> Self {
        let mut p = self.clone();
        for (e, c) in &b.0 {
            p.add(*e, -c);
        }
        p
    }
    fn plus(&self, b: &Self) -> Self {
        self.minus(&b.multiple([0; 3], &r(-1)))
    }
    fn times(&self, b: &Self) -> Self {
        let mut p = Self::zero();
        for (e, c) in &b.0 {
            p = p.plus(&self.multiple(*e, c));
        }
        p
    }
    fn pow(&self, n: u32) -> Self {
        let mut p = Self::constant(r(1));
        for _ in 0..n {
            p = p.times(self);
        }
        p
    }
    fn from_mpoly(p: &MPoly, a: &Rational, b: &Rational) -> Self {
        let x = Self::var(2)
            .minus(&Self::var(0).multiple([0; 3], a))
            .minus(&Self::var(1).multiple([0; 3], b));
        let vars = [x, Self::var(0), Self::var(1)];
        let mut out = Self::zero();
        for (m, c) in p.terms() {
            let mut term = Self::constant(c);
            for (i, e) in m.exponents().iter().enumerate() {
                term = term.times(&vars[i].pow(*e));
            }
            out = out.plus(&term);
        }
        out
    }
    fn univariate(&self) -> Option<UPoly> {
        if self.0.keys().any(|e| e[0] != 0 || e[1] != 0) {
            return None;
        }
        let mut c = vec![r(0); self.0.keys().map(|e| e[2] as usize + 1).max().unwrap_or(0)];
        for (e, v) in &self.0 {
            c[e[2] as usize] = v.clone();
        }
        Some(UPoly::new(c))
    }
    fn from_upoly(p: &UPoly) -> Self {
        let mut out = Self::zero();
        for (i, c) in p.coefficients().iter().enumerate() {
            out.add([0, 0, i as u32], c.clone());
        }
        out
    }
    fn coordinate(&self, k: usize) -> Option<UPoly> {
        let mut e = [0; 3];
        e[k] = 1;
        if self.leading()?.0 != e {
            return None;
        }
        let rest = self.minus(&Self::var(k));
        rest.univariate().map(|p| -p)
    }
}
fn divides(a: [u32; 3], b: [u32; 3]) -> bool {
    (0..3).all(|i| a[i] <= b[i])
}
fn reduce(mut p: P, basis: &[P]) -> P {
    let mut remainder = P::zero();
    while let Some((e, c)) = p.leading() {
        if let Some((g, (f, d))) = basis
            .iter()
            .filter_map(|g| g.leading().map(|l| (g, l)))
            .find(|(_, (f, _))| divides(*f, e))
        {
            p = p.minus(&g.multiple(std::array::from_fn(|i| e[i] - f[i]), &(c / d)));
        } else {
            p.0.remove(&e);
            remainder.add(e, c);
        }
    }
    remainder
}
struct CachedBasis {
    generators: Vec<P>,
    basis: Arc<[P]>,
}

thread_local! {
    // Basis reduction is pure rational arithmetic. Ordered polynomial values
    // keep independent stores and different shears distinct.
    static BASES: RefCell<VecDeque<CachedBasis>> = const { RefCell::new(VecDeque::new()) };
}

fn groebner(generators: Vec<P>) -> Arc<[P]> {
    if let Some(basis) = BASES.with(|bases| {
        bases
            .borrow()
            .iter()
            .find(|entry| entry.generators == generators)
            .map(|entry| entry.basis.clone())
    }) {
        return basis;
    }
    let basis: Arc<[P]> = uncached_groebner(generators.clone()).into();
    BASES.with(|bases| {
        let mut bases = bases.borrow_mut();
        if bases.len() == 32 {
            bases.pop_front();
        }
        bases.push_back(CachedBasis {
            generators,
            basis: basis.clone(),
        });
    });
    basis
}

fn critical_pair(a: &P, b: &P, i: usize, j: usize) -> (u64, [u32; 3], usize, usize) {
    let a = a.0.last_key_value().unwrap().0;
    let b = b.0.last_key_value().unwrap().0;
    let l = std::array::from_fn(|k| a[k].max(b[k]));
    (l.iter().map(|&e| u64::from(e)).sum(), l, i, j)
}

fn uncached_groebner(generators: Vec<P>) -> Vec<P> {
    let mut basis: Vec<P> = vec![];
    let mut pairs = BTreeSet::new();
    let mut resolved = BTreeSet::new();
    for p in generators {
        let p = reduce(p, &basis).monic();
        if !p.0.is_empty() {
            let j = basis.len();
            for i in 0..j {
                pairs.insert(critical_pair(&basis[i], &p, i, j));
            }
            basis.push(p);
        }
    }
    // Lowest LCM degree first limits coefficient swell; the remaining key
    // makes ties independent of hashing or allocation order.
    while let Some((_, l, i, j)) = pairs.pop_first() {
        let a = basis[i].0.last_key_value().unwrap().0;
        let b = basis[j].0.last_key_value().unwrap().0;
        // A processed pair has a standard representation in the current
        // basis, including when its nonzero remainder became a new generator.
        // The chain criterion may therefore reuse two such representations.
        let redundant = (0..3).all(|k| a[k] == 0 || b[k] == 0)
            || basis.iter().enumerate().any(|(k, g)| {
                k != i
                    && k != j
                    && divides(*g.0.last_key_value().unwrap().0, l)
                    && resolved.contains(&(i.min(k), i.max(k)))
                    && resolved.contains(&(j.min(k), j.max(k)))
            });
        resolved.insert((i, j));
        if redundant {
            continue;
        }
        let s = basis[i]
            .multiple(std::array::from_fn(|k| l[k] - a[k]), &r(1))
            .minus(&basis[j].multiple(std::array::from_fn(|k| l[k] - b[k]), &r(1)));
        let s = reduce(s, &basis).monic();
        if !s.0.is_empty() {
            let j = basis.len();
            for i in 0..j {
                pairs.insert(critical_pair(&basis[i], &s, i, j));
            }
            basis.push(s);
        }
    }
    let mut i = 0;
    while i < basis.len() {
        let p = basis.remove(i);
        let p = reduce(p, &basis).monic();
        if !p.0.is_empty() {
            basis.insert(i, p);
            i += 1;
        }
    }
    basis.sort_by_key(|p| p.leading().unwrap().0);
    basis
}
pub(super) enum Solutions {
    Contained,
    Points(Vec<Point3>),
}
pub(super) fn intersections(
    equations: &[RealPoly],
    surface: &RealPoly,
    dec: &mut Decider,
) -> Result<Solutions, Refusal> {
    if let Some(surface) = surface.rational()
        && equations.iter().all(|p| p.rational().is_some())
    {
        return rational_intersections(equations, surface, dec);
    }
    crate::algebraic_solve::intersections(equations, surface, dec)
}

fn rational_intersections(
    equations: &[RealPoly],
    surface: &MPoly,
    dec: &mut Decider,
) -> Result<Solutions, Refusal> {
    let mut shear = 0i64;
    loop {
        let a = r(shear);
        let b = r(shear).pow(2).unwrap();
        let curve: Vec<_> = equations
            .iter()
            .map(|p| P::from_mpoly(p.rational().expect("rational solver dispatch"), &a, &b))
            .collect();
        let basis = groebner(curve.clone());
        let surface = P::from_mpoly(surface, &a, &b);
        if reduce(surface.clone(), &basis).0.is_empty() {
            return Ok(Solutions::Contained);
        }
        let mut equations = curve;
        equations.push(surface);
        let basis = groebner(equations.clone());
        if basis
            .iter()
            .any(|p| p.0.len() == 1 && p.0.contains_key(&[0; 3]))
        {
            return Ok(Solutions::Points(vec![]));
        }
        let Some(poly) = basis
            .iter()
            .filter_map(P::univariate)
            .find(|p| p.degree().is_some_and(|d| d > 0))
        else {
            panic!("non-contained supported curve has finite meetings");
        };
        let square_free = poly.div_rem(&poly.gcd(&poly.derivative())).unwrap().0;
        equations.push(P::from_upoly(&square_free));
        // Square-free annihilators for every coordinate make the finite
        // quotient reduced, including roots with more than one tangent.
        for axis in 0..2 {
            let permute = |p: &P| {
                let mut out = P::zero();
                for (e, c) in &p.0 {
                    let mut e = *e;
                    e.swap(axis, 2);
                    out.add(e, c.clone());
                }
                out
            };
            let projected = groebner(equations.iter().map(permute).collect());
            let p = projected
                .iter()
                .filter_map(P::univariate)
                .find(|p| p.degree().is_some_and(|d| d > 0))
                .expect("finite coordinate projection");
            let p = p.div_rem(&p.gcd(&p.derivative())).unwrap().0;
            equations.push(permute(&P::from_upoly(&p)));
        }
        let radical = groebner(equations);
        let y = radical.iter().find_map(|p| p.coordinate(0));
        let z = radical.iter().find_map(|p| p.coordinate(1));
        if let (Some(y), Some(z)) = (y, z) {
            let roots = match dec.real_roots(&square_free) {
                Outcome::Proven(r) => r,
                Outcome::Refused(e) => return Err(e),
                _ => unreachable!("complete root list"),
            };
            let points = roots
                .into_iter()
                .map(|root| {
                    let t = Real::from_algebraic(root.value);
                    let y = y.eval_real(&t);
                    let z = z.eval_real(&t);
                    let x = t - &real(&a) * &y - &real(&b) * &z;
                    Point3::new(x, y, z)
                })
                .collect();
            return Ok(Solutions::Points(points));
        }
        // Only finitely many linear projections identify two distinct complex
        // solutions. This rational moment-curve family eventually separates them.
        shear += 1;
    }
}
