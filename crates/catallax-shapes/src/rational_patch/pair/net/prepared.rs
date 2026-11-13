use super::*;

/// Unrounded polynomial nets on a bounding rectangle. Only a restricted
/// chart's positive-weight check admits geometry; the bounding rectangle need
/// not itself be an authored or admitted chart.
pub(in crate::rational_patch::pair) struct SpanNets {
    domain: Domain,
    h: [Grid; 4],
    tu: [Grid; 3],
    tv: [Grid; 3],
    normal: [Grid; 3],
}

impl SpanNets {
    pub(in crate::rational_patch::pair) fn new(
        domain: &Domain,
        homogeneous: &[MPoly; 4],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut local = Vec::with_capacity(4);
        for poly in homogeneous {
            local.push(Dense::local(poly, domain, dec)?);
        }
        // The four homogeneous coordinates share one net, so one degree.
        let deg = local
            .iter()
            .fold([0, 0], |d, p| [d[0].max(p.deg[0]), d[1].max(p.deg[1])]);
        let local: Vec<Dense> = local.iter().map(|p| p.padded(deg)).collect();
        let local: [Dense; 4] = local.try_into().ok().expect("four coordinates");
        let mut tangents: Vec<[Dense; 3]> = Vec::with_capacity(2);
        for axis in 0..2 {
            let dw = local[3].derivative(axis, dec)?;
            let mut parts = Vec::with_capacity(3);
            for coordinate in local.iter().take(3) {
                let dx = coordinate.derivative(axis, dec)?;
                let left = dx.mul(&local[3], dec)?;
                let right = coordinate.mul(&dw, dec)?;
                parts.push(left.sub(&right, dec)?);
            }
            tangents.push(parts.try_into().ok().expect("three components"));
        }
        let [tu, tv]: [[Dense; 3]; 2] = tangents.try_into().ok().expect("two tangents");
        let mut normal = Vec::with_capacity(3);
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let left = tu[j].mul(&tv[k], dec)?;
            let right = tu[k].mul(&tv[j], dec)?;
            normal.push(left.sub(&right, dec)?);
        }
        let mut grids = |polys: &[Dense]| -> Result<Vec<Grid>, Refusal> {
            polys.iter().map(|p| p.bernstein(dec)).collect()
        };
        let h: [Grid; 4] = grids(&local)?.try_into().expect("four");
        let tu: [Grid; 3] = grids(&tu)?.try_into().expect("three");
        let tv: [Grid; 3] = grids(&tv)?.try_into().expect("three");
        let normal: [Grid; 3] = grids(&normal)?.try_into().expect("three");
        Ok(Self {
            domain: domain.clone(),
            h,
            tu,
            tv,
            normal,
        })
    }

    pub(in crate::rational_patch::pair) fn into_cell(
        self,
        index: usize,
        dec: &mut Decider,
    ) -> Result<Cell, Refusal> {
        let Self {
            domain,
            h,
            tu,
            tv,
            normal,
        } = self;
        if !h[3].all_positive() {
            return Err(Refusal::new(
                "intersect.spline-pair-weights",
                "a knot span has a homogeneous weight control that is not positive",
            ));
        }
        let unit_weight = h[3].v.iter().all(|w| *w == Rational::one());
        let unit = [Rational::zero(), Rational::one()];
        let cell = Cell::finish(
            index,
            domain,
            [unit.clone(), unit],
            0,
            unit_weight,
            h,
            None,
            None,
            dec,
        )?;
        let derived = Derived::new(tu, tv, normal, None, None, dec)?;
        let _ = cell.derived.set(derived);
        Ok(cell)
    }

    pub(in crate::rational_patch::pair) fn restricted_cell(
        &self,
        index: usize,
        domain: &Domain,
        dec: &mut Decider,
    ) -> Result<Cell, Refusal> {
        self.restricted(domain, dec)?.into_cell(index, dec)
    }

    pub(in crate::rational_patch::pair) fn restricted_u(
        &self,
        interval: &[Rational; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        self.restricted(&[interval.clone(), self.domain[1].clone()], dec)
    }

    fn restricted(&self, domain: &Domain, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut cuts: [[Rational; 2]; 2] =
            std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
        let mut ratios = [Rational::one(), Rational::one()];
        for axis in 0..2 {
            if less(&domain[axis][0], &self.domain[axis][0], dec)?
                || less(&self.domain[axis][1], &domain[axis][1], dec)?
                || !less(&domain[axis][0], &domain[axis][1], dec)?
            {
                return Err(Refusal::new(
                    "intersect.spline-pair-prepared-domain",
                    "a prepared cell must be a nonempty restriction of its polynomial rectangle",
                ));
            }
            let width = sub(&self.domain[axis][1], &self.domain[axis][0], dec)?;
            let low = sub(&domain[axis][0], &self.domain[axis][0], dec)?;
            let high = sub(&domain[axis][1], &self.domain[axis][0], dec)?;
            let child_width = sub(&domain[axis][1], &domain[axis][0], dec)?;
            dec.reserve_exact_arithmetic(&low, &high, true)?;
            dec.reserve_exact_arithmetic(&high, &width, true)?;
            dec.reserve_exact_arithmetic(&child_width, &width, true)?;
            cuts[axis] = [&low / &high, &high / &width];
            ratios[axis] = &child_width / &width;
        }
        let mut child = Self {
            domain: domain.clone(),
            h: self.h.clone(),
            tu: self.tu.clone(),
            tv: self.tv.clone(),
            normal: self.normal.clone(),
        };
        for (axis, [low, high]) in cuts.iter().enumerate() {
            for grid in child
                .h
                .iter_mut()
                .chain(&mut child.tu)
                .chain(&mut child.tv)
                .chain(&mut child.normal)
            {
                if grid.deg[axis] == 0 {
                    continue;
                }
                if *high != Rational::one() {
                    let [next, _] = grid.split(axis, high, dec)?;
                    *grid = next;
                }
                if !low.is_zero() {
                    let [_, next] = grid.split(axis, low, dec)?;
                    *grid = next;
                }
            }
        }
        // Cached derivatives use the bounding rectangle's unit parameters.
        // Each new root cell instead uses its own unit square.
        let normal_scale = mul(&ratios[0], &ratios[1], dec)?;
        for (grids, scale) in [
            (&mut child.tu, &ratios[0]),
            (&mut child.tv, &ratios[1]),
            (&mut child.normal, &normal_scale),
        ] {
            if *scale == Rational::one() {
                continue;
            }
            for grid in grids {
                for value in &mut grid.v {
                    *value = mul(value, scale, dec)?;
                }
            }
        }
        Ok(child)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Rational {
        Rational::from_ratio(n, d).unwrap()
    }

    #[test]
    fn reused_cells_enclose_positions_and_correctly_rescaled_tangents_and_normals() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let u = MPoly::var(2, 0);
            let v = MPoly::var(2, 1);
            let w = &MPoly::constant(2, Rational::one()) + &u;
            let h = [u.clone(), &w * &v, &u * &v, w];
            let base = [[q(2, 1), q(4, 1)], [q(-3, 1), q(1, 1)]];
            let prepared = SpanNets::new(&base, &h, &mut dec)?;
            for (rect, staged) in [
                (base.clone(), false),
                ([[q(5, 2), q(3, 1)], [q(-2, 1), q(0, 1)]], false),
                (base, true),
                ([[q(5, 2), q(3, 1)], [q(-2, 1), q(0, 1)]], true),
            ] {
                let cell = if staged {
                    prepared
                        .restricted_u(&rect[0], &mut dec)?
                        .restricted_cell(0, &rect, &mut dec)?
                } else {
                    prepared.restricted_cell(0, &rect, &mut dec)?
                };
                let derived = cell.derived.get().unwrap();
                let widths: [Rational; 2] =
                    std::array::from_fn(|axis| &rect[axis][1] - &rect[axis][0]);
                let parameters = rect
                    .each_ref()
                    .map(|[lo, hi]| [lo.clone(), (lo + hi).mul_pow2(-1), hi.clone()]);
                for u in &parameters[0] {
                    for v in &parameters[1] {
                        let at = [u.clone(), v.clone()];
                        let weight = h[3].eval(&at);
                        let point: [Rational; 3] =
                            std::array::from_fn(|i| &h[i].eval(&at) / &weight);
                        let tangent: [[Rational; 3]; 2] = std::array::from_fn(|axis| {
                            std::array::from_fn(|i| {
                                let numerator = &weight * &h[i].partial(axis).eval(&at)
                                    - h[i].eval(&at) * h[3].partial(axis).eval(&at);
                                &numerator * &widths[axis]
                            })
                        });
                        let normal: [Rational; 3] = std::array::from_fn(|i| {
                            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                            &tangent[0][j] * &tangent[1][k] - &tangent[0][k] * &tangent[1][j]
                        });
                        for (ranges, values) in [
                            (&cell.bbox, &point),
                            (&derived.tu_range, &tangent[0]),
                            (&derived.tv_range, &tangent[1]),
                            (&derived.normal_range, &normal),
                        ] {
                            for ([lo, hi], value) in ranges.iter().zip(values) {
                                assert!(
                                    lo <= value && value <= hi,
                                    "{value:?} leaves [{lo:?}, {hi:?}]"
                                );
                            }
                        }
                    }
                }
            }
            Ok(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }
}
