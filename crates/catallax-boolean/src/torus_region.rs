//! Exact cells for circular torus boundaries, including trimmed parallels and
//! meridians. Complete circles use angular labels; finite arcs use the connected
//! cells of their periodic arrangement, without introducing a parameter seam.

use crate::geometry::{refusal, sign};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{CurveId, Point3, PointId, Torus, Vector3};
use catallax_topology::{Loop, LoopId, Model};
use std::{cmp::Ordering, collections::BTreeSet};

mod trimmed;

struct Phase {
    x: Real,
    y: Real,
    half: bool,
}
impl Phase {
    fn new(x: Real, y: Real, dec: &mut Decider) -> Result<Self, Refusal> {
        let half = match sign(&y, dec, "boolean.classify-ambiguous")? {
            Sign::Positive => false,
            Sign::Negative => true,
            Sign::Zero => sign(&x, dec, "boolean.classify-ambiguous")? == Sign::Negative,
        };
        Ok(Self { x, y, half })
    }
    fn compare(&self, other: &Self, dec: &mut Decider) -> Result<Ordering, Refusal> {
        if self.half != other.half {
            return Ok(self.half.cmp(&other.half));
        }
        let cross = &(&self.x * &other.y) - &(&self.y * &other.x);
        Ok(match sign(&cross, dec, "boolean.classify-ambiguous")? {
            Sign::Positive => Ordering::Less,
            Sign::Negative => Ordering::Greater,
            Sign::Zero => Ordering::Equal,
        })
    }
}
fn add_cut(cuts: &mut Vec<Phase>, phase: Phase, dec: &mut Decider) -> Result<usize, Refusal> {
    for (i, old) in cuts.iter().enumerate() {
        if phase.compare(old, dec)? == Ordering::Equal {
            return Ok(i);
        }
    }
    let index = cuts.len();
    cuts.push(phase);
    Ok(index)
}
fn cell(cuts: &[Phase], query: &Phase, dec: &mut Decider) -> Result<i64, Refusal> {
    let mut before: Option<usize> = None;
    let mut last: Option<usize> = None;
    for (i, cut) in cuts.iter().enumerate() {
        let order = cut.compare(query, dec)?;
        if order == Ordering::Equal {
            return Err(refusal(
                "boolean.classify-ambiguous",
                "torus cell probe is on a cut",
            ));
        }
        let latest = match last {
            None => true,
            Some(previous) => cuts[previous].compare(cut, dec)? == Ordering::Less,
        };
        if latest {
            last = Some(i);
        }
        let nearer = match before {
            None => true,
            Some(previous) => cuts[previous].compare(cut, dec)? == Ordering::Less,
        };
        if order == Ordering::Less && nearer {
            before = Some(i);
        }
    }
    Ok(before.or(last).map_or(-1, |i| i as i64))
}
fn radius(squared: &Real) -> Real {
    if let Some(q) = squared.as_rational()
        && let Some(root) = crate::points::sqrt_rational(&q)
    {
        return Real::from(root);
    }
    squared.sqrt()
}

struct CircleCut {
    dimension: usize,
    phase: Phase,
    normal: Vector3,
    offset: Vector3,
}

fn circle_cut(
    model: &Model,
    torus: &Torus,
    curve: CurveId,
    dec: &mut Decider,
) -> Result<Option<CircleCut>, Refusal> {
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = curve_form(model.geom(), curve, dec)?
    else {
        return Ok(None);
    };
    let axis = &torus.axis;
    let offset = center.minus(&torus.center);
    if sign(
        &normal.cross(axis).norm_squared(),
        dec,
        "boolean.classify-ambiguous",
    )? == Sign::Zero
    {
        if sign(
            &offset.cross(axis).norm_squared(),
            dec,
            "boolean.classify-ambiguous",
        )? != Sign::Zero
        {
            return Ok(None);
        }
        Ok(Some(CircleCut {
            dimension: 1,
            phase: Phase::new(&radius(&radius_sq) - &torus.major, offset.dot(axis), dec)?,
            normal,
            offset,
        }))
    } else if sign(&normal.dot(axis), dec, "boolean.classify-ambiguous")? == Sign::Zero {
        for zero in [
            offset.dot(axis),
            offset.dot(&normal),
            offset.norm_squared() - torus.major.square(),
            &radius_sq - &torus.minor.square(),
        ] {
            if sign(&zero, dec, "boolean.classify-ambiguous")? != Sign::Zero {
                return Ok(None);
            }
        }
        Ok(Some(CircleCut {
            dimension: 0,
            phase: Phase::new(
                offset.dot(&torus.ref_dir),
                offset.dot(&axis.cross(&torus.ref_dir)),
                dec,
            )?,
            normal,
            offset,
        }))
    } else {
        Ok(None)
    }
}

fn point_phase(
    torus: &Torus,
    point: &Point3,
    dimension: usize,
    dec: &mut Decider,
) -> Result<Phase, Refusal> {
    let offset = point.minus(&torus.center);
    if dimension == 0 {
        Phase::new(
            offset.dot(&torus.ref_dir),
            offset.dot(&torus.axis.cross(&torus.ref_dir)),
            dec,
        )
    } else {
        let height = offset.dot(&torus.axis);
        let radial2 = offset.norm_squared() - height.square() / torus.axis.norm_squared();
        Phase::new(&radius(&radial2) - &torus.major, height, dec)
    }
}

pub(crate) fn signature(
    model: &Model,
    torus: &Torus,
    groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Vec<i64>>, Refusal> {
    if let Some(value) = whole_signature(model, torus, groups, point, dec)? {
        return Ok(Some(value));
    }
    trimmed::signature(model, torus, groups, point, dec)
}

fn whole_signature(
    model: &Model,
    torus: &Torus,
    groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Vec<i64>>, Refusal> {
    let mut azimuth = Vec::new();
    let mut tube = Vec::new();
    let mut seen = BTreeSet::new();
    for &loop_id in groups.iter().flatten() {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            return Ok(None);
        };
        let [coedge] = coedges.as_slice() else {
            return Ok(None);
        };
        let edge = model.edge(model.coedge(*coedge).edge);
        if edge.start.is_some() || edge.end.is_some() {
            return Ok(None);
        }
        if !seen.insert(edge.curve) {
            continue;
        }
        let Some(cut) = circle_cut(model, torus, edge.curve, dec)? else {
            return Ok(None);
        };
        add_cut(
            if cut.dimension == 0 {
                &mut azimuth
            } else {
                &mut tube
            },
            cut.phase,
            dec,
        )?;
    }
    let position = model.geom().position(point);
    let u = if azimuth.is_empty() {
        -1
    } else {
        let query = point_phase(torus, &position, 0, dec)?;
        cell(&azimuth, &query, dec)?
    };
    let v = if tube.is_empty() {
        -1
    } else {
        let query = point_phase(torus, &position, 1, dec)?;
        cell(&tube, &query, dec)?
    };
    Ok(Some(vec![u, v]))
}
