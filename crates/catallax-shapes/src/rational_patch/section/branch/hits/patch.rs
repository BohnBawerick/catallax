//! Exact clipping through hulls, retained cuts, or a jointly embedded boundary.

use super::*;
use crate::SplineBoundary;
use catallax_numbers::Interval;

fn hull(patch: &RationalPatch, dec: &mut Decider) -> Result<[Interval; 3], Refusal> {
    let mut points = patch.carrier().controls().iter();
    let first = points.next().expect("an admitted patch has controls");
    let mut bounds = [first.x.approx(), first.y.approx(), first.z.approx()];
    for point in points {
        dec.charge_work(1)?;
        let next = [point.x.approx(), point.y.approx(), point.z.approx()];
        for axis in 0..3 {
            bounds[axis] = bounds[axis].hull(&next[axis]);
        }
    }
    Ok(bounds)
}

fn cut_misses(
    surface: &SectionSurface,
    patch: &RationalPatch,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for (_, controls) in patch.bezier_spans(dec)? {
        if !surface.retains_span(controls, dec)? {
            continue;
        }
        let (values, _) = surface.net(controls, patch.carrier().degrees(), dec)?;
        if uniform(&values, dec)?.is_none() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn bezier(patch: &RationalPatch) -> bool {
    let surface = patch.carrier();
    let counts = surface.control_count();
    if counts != surface.degrees().map(|degree| degree + 1) {
        return false;
    }
    (0..2).all(|axis| {
        let knots = surface.knots()[axis];
        knots.len() == 2 * counts[axis]
            && knots[..counts[axis]]
                .iter()
                .all(|knot| knot == &patch.domain()[axis][0])
            && knots[counts[axis]..]
                .iter()
                .all(|knot| knot == &patch.domain()[axis][1])
    })
}

fn same_boundary(
    first: &RationalPatch,
    first_side: SplineBoundary,
    second: &RationalPatch,
    second_side: SplineBoundary,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let first = first.boundary(first_side);
    let second = second.boundary(second_side);
    if first.degree() != second.degree() || first.control_count() != second.control_count() {
        return Ok(false);
    }
    for index in 0..first.control_count() {
        dec.charge_work(1)?;
        let (Some(a), Some(b)) = (
            first.weight(index).rational_ref(),
            second.weight(index).rational_ref(),
        ) else {
            return Ok(false);
        };
        if a != b {
            return Ok(false);
        }
        let (a, b) = (first.control(index), second.control(index));
        for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
            let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) else {
                return Ok(false);
            };
            if a != b {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn positive_cover(
    bound: &ProjectedJacobian,
    domain: Domain,
    widths: &[Rational; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    // An optional exact certificate, not a new numerical tolerance or root rule.
    let mut pending = vec![(domain, 0u32)];
    let mut visited = 0usize;
    while let Some((cell, depth)) = pending.pop() {
        dec.check_deadline()?;
        visited += 1;
        if visited > 128 {
            return Ok(false);
        }
        if bound.margin(&cell, dec)? > Rational::zero() {
            continue;
        }
        if depth == 8 {
            return Ok(false);
        }
        let [low, high] = split(&cell, widths);
        pending.push((high, depth + 1));
        pending.push((low, depth + 1));
    }
    Ok(true)
}

fn jointly_embedded(
    branch: &PatchSectionBranch,
    other: &RationalPatch,
    axis: usize,
    high: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let source = &branch.charts[0];
    let mut support = source.domain().clone();
    if high {
        support[axis][0] = branch
            .cells
            .iter()
            .map(|cell| &cell.domain[axis][0])
            .min()
            .unwrap()
            .clone();
    } else {
        support[axis][1] = branch
            .cells
            .iter()
            .map(|cell| &cell.domain[axis][1])
            .max()
            .unwrap()
            .clone();
    }
    // Normalize each chart separately to a unit rectangle and offset its fixed
    // coordinate to the common seam. Equal Bezier boundary data makes this
    // piecewise map continuous. A common projection with a positive symmetric
    // Jacobian on both rectangles is strongly monotone across the seam too:
    // the line integral between two points splits once at that seam.
    let projection = other.embedding_projection();
    for (patch, domain) in [(source, support), (other, other.domain().clone())] {
        let [Some(x), Some(y), Some(z), Some(w)] = patch.data.spans[0]
            .coordinates
            .each_ref()
            .map(RealPoly::as_rational)
        else {
            return Ok(false);
        };
        let mut frame = projection.clone();
        let bound = ProjectedJacobian::new(
            &[x, y, z, w],
            &mut frame,
            patch.domain().clone(),
            &Rational::one(),
            dec,
        )?;
        if !positive_cover(&bound, domain, &patch.data.widths, dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn boundary_points(
    branch: &PatchSectionBranch,
    axis: usize,
    high: bool,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let source = &branch.charts[0];
    let span = &source.data.spans[0];
    let fixed = Real::from_rational(source.domain()[axis][usize::from(high)].clone());
    let restricted = branch
        .surface
        .equation(span)
        .substitute(axis, &RealPoly::constant(2, fixed.clone()));
    let polynomial = axis_polynomial(&restricted, 1 - axis, dec)?;
    if polynomial.is_zero(dec)? {
        return Ok(None);
    }
    let roots = polynomial.real_roots_with_coarse_brackets(dec)?;
    let mut points = Vec::new();
    let mut previous = None;
    for cell in branch.cells.iter() {
        let free = 1 - cell.solved;
        let forward =
            exact_sign(&(&cell.ends[1][free] - &cell.ends[0][free]), dec)? == Sign::Positive;
        for index in 0..roots.len() {
            let index = if forward {
                index
            } else {
                roots.len() - 1 - index
            };
            let mut at = [fixed.clone(), fixed.clone()];
            at[1 - axis] = roots[index].clone();
            if !contains(&cell.domain, &at, dec)? {
                continue;
            }
            if let Some(previous) = &previous {
                if same_uv(previous, &at, dec)? {
                    continue;
                }
            }
            points.push(span.point(&at));
            previous = Some(at);
        }
    }
    Ok(Some(points))
}

impl PatchSectionBranch {
    /// A complete exact finite meeting with a native patch when cheap certificates apply.
    ///
    /// Positive spline weights make disjoint Cartesian control hull boxes an
    /// exclusion proof. A strictly signed Bernstein net for the retained cut
    /// also excludes a target patch. Opposite Bezier sides with equal rational
    /// controls and weights identify the same normalized boundary. A common
    /// strictly monotone projection proves the supports meet only there; its
    /// scalar cut roots retain the section solver's exact representation.
    /// Finite side contacts need not be transverse. No numerical meeting
    /// becomes an incidence fact. `None` means no exact certificate applied.
    pub fn patch_hits(
        &self,
        other: &RationalPatch,
        dec: &mut Decider,
    ) -> Result<Option<SectionHits>, Refusal> {
        if let [source] = self.charts.as_ref() {
            let (a, b) = (hull(source, dec)?, hull(other, dec)?);
            if (0..3).any(|axis| a[axis].hi() < b[axis].lo() || b[axis].hi() < a[axis].lo()) {
                return Ok(Some(SectionHits::Points(Vec::new())));
            }
        }
        if cut_misses(&self.surface, other, dec)? {
            return Ok(Some(SectionHits::Points(Vec::new())));
        }
        let [source] = self.charts.as_ref() else {
            return Ok(None);
        };
        if !bezier(source) || !bezier(other) {
            return Ok(None);
        }
        for (first, second, axis, high) in [
            (SplineBoundary::UMin, SplineBoundary::UMax, 0, false),
            (SplineBoundary::UMax, SplineBoundary::UMin, 0, true),
            (SplineBoundary::VMin, SplineBoundary::VMax, 1, false),
            (SplineBoundary::VMax, SplineBoundary::VMin, 1, true),
        ] {
            if same_boundary(source, first, other, second, dec)?
                && jointly_embedded(self, other, axis, high, dec)?
            {
                return Ok(boundary_points(self, axis, high, dec)?.map(SectionHits::Points));
            }
        }
        Ok(None)
    }
}
