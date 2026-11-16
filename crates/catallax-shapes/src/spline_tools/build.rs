// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: exact rational interpolation and Boolean sums, immutable
// inputs, exact network intersections instead of averaging or moving curves.
// Adapted from SMLib NL_SrfAdv.cpp, N_CreateGordonSrf and N_CreateCoonsSrf,
// and NL_SrfFit.cpp, tensor interpolation, at
// c9c5979f04cea8c38213004d4f04df0dcfb40d05.

use super::edit::bezier_knots;
use super::*;

#[derive(Clone)]
struct Poly {
    degree: [usize; 2],
    values: Vec<Rational>,
}

fn size(dimensions: &[usize]) -> Result<usize, Refusal> {
    dimensions.iter().try_fold(1usize, |n, d| {
        n.checked_mul(*d)
            .ok_or_else(|| invalid("surface construction dimension product overflows"))
    })
}

fn capacity<T>(length: usize, dec: &mut Decider) -> Result<Vec<T>, Refusal> {
    let units =
        u64::try_from(length).map_err(|_| invalid("surface construction work count overflows"))?;
    dec.charge_work(units)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|e| Refusal::new("shapes.spline-construction-memory", e.to_string()))?;
    Ok(values)
}

fn filled<T: Clone>(length: usize, value: T, dec: &mut Decider) -> Result<Vec<T>, Refusal> {
    let mut values = capacity(length, dec)?;
    values.resize(length, value);
    Ok(values)
}

fn choose(n: usize, k: usize, dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut value = Rational::one();
    for i in 0..k.min(n - k) {
        dec.charge_work(3)?;
        value = value * Rational::from_ratio((n - i) as i64, (i + 1) as i64).unwrap();
    }
    Ok(value)
}

impl Poly {
    fn constant(x: Rational) -> Self {
        Self {
            degree: [0, 0],
            values: vec![x],
        }
    }

    fn elevated(&self, target: [usize; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        let mut out = self.clone();
        for axis in 0..2 {
            while out.degree[axis] < target[axis] {
                let mut degree = out.degree;
                degree[axis] += 1;
                let mut values = capacity(size(&[degree[0] + 1, degree[1] + 1])?, dec)?;
                for u in 0..=degree[0] {
                    for v in 0..=degree[1] {
                        dec.charge_work(6)?;
                        let i = if axis == 0 { u } else { v };
                        let alpha = Rational::from_ratio(i as i64, degree[axis] as i64).unwrap();
                        let mut value = Rational::zero();
                        if i < degree[axis] {
                            value = &out.values[u * (out.degree[1] + 1) + v]
                                * &(&Rational::one() - &alpha);
                        }
                        if i > 0 {
                            let index = if axis == 0 {
                                (u - 1) * (out.degree[1] + 1) + v
                            } else {
                                u * (out.degree[1] + 1) + v - 1
                            };
                            value = value + &out.values[index] * &alpha;
                        }
                        values.push(value);
                    }
                }
                out = Self { degree, values };
            }
        }
        Ok(out)
    }

    fn plus(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let degree = std::array::from_fn(|i| self.degree[i].max(other.degree[i]));
        let a = self.elevated(degree, dec)?;
        let b = other.elevated(degree, dec)?;
        let mut values = capacity(a.values.len(), dec)?;
        for (a, b) in a.values.iter().zip(&b.values) {
            dec.charge_work(1)?;
            values.push(a + b);
        }
        Ok(Self { degree, values })
    }

    fn times(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let degree = [
            self.degree[0]
                .checked_add(other.degree[0])
                .ok_or_else(|| invalid("surface product degree overflows"))?,
            self.degree[1]
                .checked_add(other.degree[1])
                .ok_or_else(|| invalid("surface product degree overflows"))?,
        ];
        let _products = size(&[self.values.len(), other.values.len()])?;
        let mut values = filled(
            size(&[degree[0] + 1, degree[1] + 1])?,
            Rational::zero(),
            dec,
        )?;
        for i in 0..=self.degree[0] {
            for j in 0..=self.degree[1] {
                for k in 0..=other.degree[0] {
                    for l in 0..=other.degree[1] {
                        dec.charge_work(9)?;
                        let scale = choose(self.degree[0], i, dec)?
                            * choose(other.degree[0], k, dec)?
                            / choose(degree[0], i + k, dec)?
                            * choose(self.degree[1], j, dec)?
                            * choose(other.degree[1], l, dec)?
                            / choose(degree[1], j + l, dec)?;
                        let index = (i + k) * (degree[1] + 1) + j + l;
                        values[index] = &values[index]
                            + &(&self.values[i * (self.degree[1] + 1) + j]
                                * &other.values[k * (other.degree[1] + 1) + l]
                                * scale);
                    }
                }
            }
        }
        Ok(Self { degree, values })
    }

    fn scaled(&self, factor: &Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut values = capacity(self.values.len(), dec)?;
        for value in &self.values {
            dec.charge_work(1)?;
            values.push(value * factor);
        }
        Ok(Self {
            degree: self.degree,
            values,
        })
    }
}

struct Patch {
    xyz: [Poly; 3],
    weight: Poly,
}

impl Patch {
    fn zero() -> Self {
        Self {
            xyz: std::array::from_fn(|_| Poly::constant(Rational::zero())),
            weight: Poly::constant(Rational::one()),
        }
    }

    fn add(&self, other: &Self, sign: i64, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut xyz = Vec::with_capacity(3);
        let sign = Rational::from_i64(sign);
        for axis in 0..3 {
            xyz.push(
                self.xyz[axis].times(&other.weight, dec)?.plus(
                    &other.xyz[axis]
                        .times(&self.weight, dec)?
                        .scaled(&sign, dec)?,
                    dec,
                )?,
            );
        }
        Ok(Self {
            xyz: xyz.try_into().unwrap_or_else(|_| unreachable!()),
            weight: self.weight.times(&other.weight, dec)?,
        })
    }

    fn curve(
        curve: &RationalSplineCurve,
        axis: usize,
        domain: [&Rational; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(
            u64::try_from(size(&[
                curve.controls.len(),
                curve.degree + 1,
                curve.degree + 1,
                4,
            ])?)
            .map_err(|_| invalid("surface extraction work count overflows"))?,
        )?;
        let constant = curve.weights.iter().all(|w| w == &curve.weights[0]);
        let curve = curve.extracted(domain, dec)?;
        let degree = if axis == 0 {
            [curve.degree, 0]
        } else {
            [0, curve.degree]
        };
        let xyz = std::array::from_fn(|i| Poly {
            degree,
            values: curve
                .controls
                .iter()
                .zip(&curve.weights)
                .map(|(p, w)| if constant { p[i].clone() } else { &p[i] * w })
                .collect(),
        });
        let weight = if constant {
            Poly::constant(Rational::one())
        } else {
            Poly {
                degree,
                values: curve.weights,
            }
        };
        Ok(Self { xyz, weight })
    }

    fn multiplied(&self, basis: &Poly, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut xyz = Vec::with_capacity(3);
        for x in &self.xyz {
            xyz.push(x.times(basis, dec)?);
        }
        Ok(Self {
            xyz: xyz.try_into().unwrap_or_else(|_| unreachable!()),
            weight: self.weight.clone(),
        })
    }

    fn net(&self, degrees: [usize; 2], dec: &mut Decider) -> Result<Vec<Homogeneous>, Refusal> {
        let xyz = self
            .xyz
            .iter()
            .map(|p| p.elevated(degrees, dec))
            .collect::<Result<Vec<_>, _>>()?;
        let w = self.weight.elevated(degrees, dec)?;
        let mut net = capacity(w.values.len(), dec)?;
        for i in 0..w.values.len() {
            net.push([
                xyz[0].values[i].clone(),
                xyz[1].values[i].clone(),
                xyz[2].values[i].clone(),
                w.values[i].clone(),
            ]);
        }
        Ok(net)
    }

    fn degree(&self) -> [usize; 2] {
        std::array::from_fn(|a| {
            self.xyz
                .iter()
                .map(|p| p.degree[a])
                .max()
                .unwrap()
                .max(self.weight.degree[a])
                .max(1)
        })
    }
}

// Cox-de Boor basis values in the active span, including the right endpoint.
fn basis(
    degree: usize,
    knots: &[Rational],
    count: usize,
    t: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Rational>, Refusal> {
    let span = if t == &knots[count] {
        count - 1
    } else {
        knots.partition_point(|k| k <= t) - 1
    };
    let mut values = filled(degree + 1, Rational::zero(), dec)?;
    values[0] = Rational::one();
    for j in 1..=degree {
        let mut saved = Rational::zero();
        for r in 0..j {
            dec.charge_work(7)?;
            let right = &knots[span + r + 1] - t;
            let left = t - &knots[span + r + 1 - j];
            let temp = &values[r] / &(&right + &left);
            values[r] = &saved + &(&right * &temp);
            saved = &left * &temp;
        }
        values[j] = saved;
    }
    let mut all = filled(count, Rational::zero(), dec)?;
    all[span - degree..=span].clone_from_slice(&values);
    Ok(all)
}

// Inverse interpolation matrix. Row i is the controls of cardinal spline i.
fn cardinals(
    nodes: &[Rational],
    degree: usize,
    domain: [&Rational; 2],
    dec: &mut Decider,
) -> Result<Vec<RationalSplineCurve>, Refusal> {
    let n = nodes.len();
    if n == 0
        || degree >= n
        || nodes.windows(2).any(|p| p[0] >= p[1])
        || nodes[0] < *domain[0]
        || nodes[n - 1] > *domain[1]
        || domain[0] >= domain[1]
    {
        return Err(invalid(
            "interpolation requires increasing parameters inside a nonempty domain and degree below point count",
        ));
    }
    if n == 1 {
        return Ok(vec![RationalSplineCurve::new(
            1,
            vec![
                domain[0].clone(),
                domain[0].clone(),
                domain[1].clone(),
                domain[1].clone(),
            ],
            vec![[Rational::one(), Rational::zero(), Rational::zero()]; 2],
            vec![Rational::one(); 2],
            dec,
        )?]);
    }
    if degree == 0 {
        return Err(invalid("interpolating degree must be positive"));
    }
    let knot_count = n
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("interpolation knot count overflows"))?;
    let mut knots = capacity(knot_count, dec)?;
    knots.resize(degree + 1, domain[0].clone());
    for j in 1..n - degree {
        dec.charge_work(
            u64::try_from(degree + 1)
                .map_err(|_| invalid("interpolation degree work overflows"))?,
        )?;
        knots.push(
            nodes[j..j + degree]
                .iter()
                .fold(Rational::zero(), |a, b| &a + b)
                / Rational::from_i64(degree as i64),
        );
    }
    knots.extend(std::iter::repeat_n(domain[1].clone(), degree + 1));
    cardinals_on_knots(nodes, degree, knots, dec)
}

fn cardinals_on_knots(
    nodes: &[Rational],
    degree: usize,
    knots: Vec<Rational>,
    dec: &mut Decider,
) -> Result<Vec<RationalSplineCurve>, Refusal> {
    let n = nodes.len();
    let width = size(&[2, n])?;
    let entries = size(&[n, width])?
        .checked_add(knots.len())
        .ok_or_else(|| invalid("interpolation work count overflows"))?;
    dec.charge_work(
        u64::try_from(entries).map_err(|_| invalid("interpolation work count overflows"))?,
    )?;
    axis(degree, &knots, n, dec)?;
    if knots.first() != Some(&knots[degree]) || knots.last() != Some(&knots[n]) {
        return Err(invalid("prescribed interpolation requires clamped knot ends"));
    }
    if nodes.windows(2).any(|pair| pair[0] >= pair[1])
        || nodes[0] < knots[degree]
        || nodes[n - 1] > knots[n]
    {
        return Err(invalid(
            "interpolation nodes must increase inside the active knot domain",
        ));
    }
    let mut matrix = capacity(n, dec)?;
    for (i, t) in nodes.iter().enumerate() {
        let mut row = basis(degree, &knots, n, t, dec)?;
        row.try_reserve_exact(n)
            .map_err(|e| Refusal::new("shapes.spline-construction-memory", e.to_string()))?;
        row.extend((0..n).map(|j| Rational::from_i64(i64::from(i == j))));
        matrix.push(row);
    }
    for pivot in 0..n {
        dec.charge_work(
            u64::try_from(n - pivot).map_err(|_| invalid("interpolation pivot work overflows"))?,
        )?;
        let source = (pivot..n)
            .find(|i| !matrix[*i][pivot].is_zero())
            .ok_or_else(|| invalid("singular interpolation parameters"))?;
        matrix.swap(pivot, source);
        let divisor = matrix[pivot][pivot].clone();
        for j in pivot..width {
            dec.charge_work(1)?;
            matrix[pivot][j] = &matrix[pivot][j] / &divisor;
        }
        for i in 0..n {
            if i == pivot {
                continue;
            }
            let factor = matrix[i][pivot].clone();
            if factor.is_zero() {
                continue;
            }
            for j in pivot..width {
                dec.charge_work(2)?;
                matrix[i][j] = &matrix[i][j] - &(&factor * &matrix[pivot][j]);
            }
        }
    }
    let mut result = capacity(n, dec)?;
    for i in 0..n {
        let mut controls = capacity(n, dec)?;
        for row in &matrix {
            controls.push([row[n + i].clone(), Rational::zero(), Rational::zero()]);
        }
        let mut copied_knots = capacity(knots.len(), dec)?;
        copied_knots.extend_from_slice(&knots);
        result.push(RationalSplineCurve::new(
            degree,
            copied_knots,
            controls,
            filled(n, Rational::one(), dec)?,
            dec,
        )?);
    }
    Ok(result)
}

fn scalar_piece(
    curve: &RationalSplineCurve,
    axis: usize,
    domain: [&Rational; 2],
    dec: &mut Decider,
) -> Result<Poly, Refusal> {
    Ok(Patch::curve(curve, axis, domain, dec)?.xyz[0].clone())
}

fn assemble(
    breaks: [&[Rational]; 2],
    patches: Vec<Patch>,
    dec: &mut Decider,
) -> Result<RationalSplineSurface, Refusal> {
    let degrees = std::array::from_fn(|a| patches.iter().map(|p| p.degree()[a]).max().unwrap());
    let cells = [breaks[0].len() - 1, breaks[1].len() - 1];
    let count = [
        size(&[cells[0], degrees[0]])?
            .checked_add(1)
            .ok_or_else(|| invalid("surface control count overflows"))?,
        size(&[cells[1], degrees[1]])?
            .checked_add(1)
            .ok_or_else(|| invalid("surface control count overflows"))?,
    ];
    let mut net = filled(
        size(&count)?,
        std::array::from_fn(|_| Rational::zero()),
        dec,
    )?;
    for (index, patch) in patches.into_iter().enumerate() {
        let cell = [index / cells[1], index % cells[1]];
        let local = patch.net(degrees, dec)?;
        for u in 0..=degrees[0] {
            for v in 0..=degrees[1] {
                let target = (cell[0] * degrees[0] + u) * count[1] + cell[1] * degrees[1] + v;
                let value = &local[u * (degrees[1] + 1) + v];
                dec.charge_work(10)?;
                if !net[target][3].is_zero() && &net[target] != value {
                    // Equivalent homogeneous scales are allowed, but never a geometric gap.
                    if (0..3).any(|a| &net[target][a] * &value[3] != &value[a] * &net[target][3]) {
                        return Err(invalid(
                            "surface construction produced incompatible patch seams",
                        ));
                    }
                    // The construction uses the same denominators on both sides of a seam.
                    return Err(invalid(
                        "surface construction produced inconsistent homogeneous seam scales",
                    ));
                }
                net[target] = value.clone();
            }
        }
    }
    let (controls, weights) = euclidean(net, dec)?;
    RationalSplineSurface::new(
        degrees,
        [
            bezier_knots(breaks[0], degrees[0]),
            bezier_knots(breaks[1], degrees[1]),
        ],
        count,
        controls,
        weights,
        dec,
    )
}

fn normalized(
    curves: &[RationalSplineCurve],
    dec: &mut Decider,
) -> Result<Vec<RationalSplineCurve>, Refusal> {
    curves
        .iter()
        .map(|c| c.reparameterized([&Rational::zero(), &Rational::one()], dec))
        .collect()
}

fn breaks(curves: &[RationalSplineCurve], bases: &[RationalSplineCurve]) -> Vec<Rational> {
    let mut result: Vec<_> = curves
        .iter()
        .chain(bases)
        .flat_map(|c| c.knots.iter().cloned())
        .collect();
    result.sort();
    result.dedup();
    result
}

impl RationalSplineSurface {
    /// Interpolates a rectangular u-major grid at the supplied parameter values.
    /// Tensor B-spline interpolation uses averaged knots and exact rational solves.
    /// This constructs a definition, not a regularity or solid certificate.
    pub fn through_grid(
        parameters: [&[Rational]; 2],
        degrees: [usize; 2],
        points: &[RationalPoint],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::interpolate_grid(parameters, degrees, points, None, dec)
    }

    /// Interpolates a u-major grid in a prescribed clamped expanded knot basis.
    /// Parameters must increase inside each active domain; singular systems refuse.
    pub fn through_grid_on_knots(
        parameters: [&[Rational]; 2],
        degrees: [usize; 2],
        knots: [Vec<Rational>; 2],
        points: &[RationalPoint],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::interpolate_grid(parameters, degrees, points, Some(knots), dec)
    }

    fn interpolate_grid(
        parameters: [&[Rational]; 2],
        degrees: [usize; 2],
        points: &[RationalPoint],
        knots: Option<[Vec<Rational>; 2]>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let count = parameters.map(<[Rational]>::len);
        if count[0] < 2 || count[1] < 2 || count[0].checked_mul(count[1]) != Some(points.len()) {
            return Err(invalid(
                "grid must contain at least two rows and columns and match its parameter counts",
            ));
        }
        let _tensor_work = size(&[points.len(), points.len(), 7])?;
        let bases = if let Some([u, v]) = knots {
            [
                cardinals_on_knots(parameters[0], degrees[0], u, dec)?,
                cardinals_on_knots(parameters[1], degrees[1], v, dec)?,
            ]
        } else {
            [
                cardinals(
                    parameters[0],
                    degrees[0],
                    [&parameters[0][0], &parameters[0][count[0] - 1]],
                    dec,
                )?,
                cardinals(
                    parameters[1],
                    degrees[1],
                    [&parameters[1][0], &parameters[1][count[1] - 1]],
                    dec,
                )?,
            ]
        };
        let mut controls = filled(points.len(), std::array::from_fn(|_| Rational::zero()), dec)?;
        for u in 0..count[0] {
            for v in 0..count[1] {
                for i in 0..count[0] {
                    for j in 0..count[1] {
                        dec.charge_work(7)?;
                        let scale = &bases[0][i].controls[u][0] * &bases[1][j].controls[v][0];
                        for a in 0..3 {
                            controls[u * count[1] + v][a] = &controls[u * count[1] + v][a]
                                + &(&points[i * count[1] + j][a] * &scale);
                        }
                    }
                }
            }
        }
        Self::new(
            degrees,
            [bases[0][0].knots.clone(), bases[1][0].knots.clone()],
            count,
            controls,
            filled(points.len(), Rational::one(), dec)?,
            dec,
        )
    }

    /// Builds a bilinearly blended Coons patch. `u_curves` are the v-min/v-max
    /// boundaries; `v_curves` are u-min/u-max. Each pair has the same direction.
    /// Input domains are mapped to [0,1]. All four corners must agree exactly.
    /// Rational boundaries retain their parameterization, with no polynomial fit.
    pub fn coons(
        u_curves: &[RationalSplineCurve; 2],
        v_curves: &[RationalSplineCurve; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let nodes = [Rational::zero(), Rational::one()];
        Self::gordon(u_curves, &nodes, v_curves, &nodes, dec)
    }

    /// Gordon Boolean sum of two interpolating skins minus their intersection grid.
    /// Profiles run in u at `v_parameters`, guides run in v at `u_parameters`.
    /// Curve domains map affinely to [0,1]. Intersections at the given parameters
    /// must be exactly equal. Interpolation degree is cubic, or lower for fewer curves.
    /// Unlike SMLib's polynomial approximation, rational networks remain exact.
    pub fn gordon(
        profiles: &[RationalSplineCurve],
        v_parameters: &[Rational],
        guides: &[RationalSplineCurve],
        u_parameters: &[Rational],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if profiles.is_empty()
            || guides.is_empty()
            || profiles.len() != v_parameters.len()
            || guides.len() != u_parameters.len()
        {
            return Err(invalid(
                "Gordon curve and parameter counts must match and be nonempty",
            ));
        }
        let profiles = normalized(profiles, dec)?;
        let guides = normalized(guides, dec)?;
        let zero = Rational::zero();
        let one = Rational::one();
        let bu = cardinals(u_parameters, (guides.len() - 1).min(3), [&zero, &one], dec)?;
        let bv = cardinals(
            v_parameters,
            (profiles.len() - 1).min(3),
            [&zero, &one],
            dec,
        )?;
        let mut grid = capacity(size(&[profiles.len(), guides.len()])?, dec)?;
        for (j, guide) in guides.iter().enumerate() {
            for (i, profile) in profiles.iter().enumerate() {
                let point = profile.point_at(&u_parameters[j], dec)?;
                if point != guide.point_at(&v_parameters[i], dec)? {
                    return Err(invalid(
                        "Gordon network curves do not intersect at the supplied parameters",
                    ));
                }
                grid.push(point);
            }
        }
        let u_breaks = breaks(&profiles, &bu);
        let v_breaks = breaks(&guides, &bv);
        let mut patches = capacity(size(&[u_breaks.len() - 1, v_breaks.len() - 1])?, dec)?;
        for u in u_breaks.windows(2) {
            for v in v_breaks.windows(2) {
                let ub = bu
                    .iter()
                    .map(|c| scalar_piece(c, 0, [&u[0], &u[1]], dec))
                    .collect::<Result<Vec<_>, _>>()?;
                let vb = bv
                    .iter()
                    .map(|c| scalar_piece(c, 1, [&v[0], &v[1]], dec))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut patch = Patch::zero();
                for (i, profile) in profiles.iter().enumerate() {
                    patch = patch.add(
                        &Patch::curve(profile, 0, [&u[0], &u[1]], dec)?.multiplied(&vb[i], dec)?,
                        1,
                        dec,
                    )?;
                }
                for (j, guide) in guides.iter().enumerate() {
                    patch = patch.add(
                        &Patch::curve(guide, 1, [&v[0], &v[1]], dec)?.multiplied(&ub[j], dec)?,
                        1,
                        dec,
                    )?;
                }
                for (j, ub) in ub.iter().enumerate() {
                    for (i, vb) in vb.iter().enumerate() {
                        let product = ub.times(vb, dec)?;
                        let point = &grid[j * profiles.len() + i];
                        let term = Patch {
                            xyz: [
                                product.scaled(&point[0], dec)?,
                                product.scaled(&point[1], dec)?,
                                product.scaled(&point[2], dec)?,
                            ],
                            weight: Poly::constant(one.clone()),
                        };
                        patch = patch.add(&term, -1, dec)?;
                    }
                }
                patches.push(patch);
            }
        }
        assemble([&u_breaks, &v_breaks], patches, dec)
    }

    /// Interpolates cross sections and one or more guide curves using the Gordon
    /// construction. Explicit intersection parameters avoid hidden curve snapping.
    pub fn guided_skin(
        profiles: &[RationalSplineCurve],
        v_parameters: &[Rational],
        guides: &[RationalSplineCurve],
        u_parameters: &[Rational],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if profiles.len() < 2 {
            return Err(invalid("guided skin requires at least two cross sections"));
        }
        Self::gordon(profiles, v_parameters, guides, u_parameters, dec)
    }

    /// Fills a directed closed N-sided boundary with N Coons quadrilaterals.
    /// Each boundary splits at its parameter midpoint. Straight midpoint-to-centre
    /// spokes are shared exactly. The caller chooses the centre; only C0 joins are
    /// claimed, not tangent continuity, injectivity or a sewn surface body.
    pub fn n_sided_fill(
        boundaries: &[RationalSplineCurve],
        centre: &RationalPoint,
        dec: &mut Decider,
    ) -> Result<Vec<Self>, Refusal> {
        if boundaries.len() < 3 {
            return Err(invalid("N-sided fill requires at least three boundaries"));
        }
        let curves = normalized(boundaries, dec)?;
        let zero = Rational::zero();
        let one = Rational::one();
        let half = Rational::from_ratio(1, 2).unwrap();
        let mut halves = Vec::with_capacity(curves.len());
        let mut spokes = Vec::with_capacity(curves.len());
        for (i, c) in curves.iter().enumerate() {
            if c.point_at(&one, dec)? != curves[(i + 1) % curves.len()].point_at(&zero, dec)? {
                return Err(invalid(
                    "N-sided boundary is not an exactly closed directed chain",
                ));
            }
            halves.push(c.split(&half, dec)?);
            spokes.push(RationalSplineCurve::new(
                1,
                vec![zero.clone(), zero.clone(), one.clone(), one.clone()],
                vec![c.point_at(&half, dec)?, centre.clone()],
                vec![one.clone(); 2],
                dec,
            )?);
        }
        (0..curves.len())
            .map(|i| {
                let previous = (i + curves.len() - 1) % curves.len();
                Self::coons(
                    &[halves[i][0].clone(), spokes[previous].clone()],
                    &[halves[previous][1].reversed(dec)?, spokes[i].clone()],
                    dec,
                )
            })
            .collect()
    }
}
#[cfg(test)]
mod work_tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn bernstein_product_debits_inner_rational_work() {
        let polynomial = Poly {
            degree: [3, 3],
            values: vec![Rational::one(); 16],
        };
        let (result, work) = with_operation_budget(1000, || {
            polynomial.times(&polynomial, &mut Decider::default())
        });
        assert_eq!(
            result
                .err()
                .expect("product must exhaust its allowance")
                .code,
            "numbers.operation-budget"
        );
        assert!(work.exhausted);
    }
}
