//! Circular noncontractible cuts partition either angular circle of a torus.
//! Label a cell by its preceding cut, not by a winding jump across an arbitrary
//! seam. The labels therefore agree on both sides of an uncut parameter seam.

use crate::geometry::{refusal, sign};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{PointId, Torus};
use catallax_topology::{Loop, LoopId, Model};
use std::{cmp::Ordering, collections::BTreeSet};

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
fn add_cut(cuts: &mut Vec<Phase>, phase: Phase, dec: &mut Decider) -> Result<(), Refusal> {
    for old in cuts.iter() {
        if phase.compare(old, dec)? == Ordering::Equal {
            return Ok(());
        }
    }
    cuts.push(phase);
    Ok(())
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

pub(crate) fn signature(
    model: &Model,
    torus: &Torus,
    groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Vec<i64>>, Refusal> {
    let axis = &torus.axis;
    let first = &torus.ref_dir;
    let second = torus.axis.cross(&torus.ref_dir);
    let center = &torus.center;
    let major = &torus.major;
    let major2 = torus.major.square();
    let minor2 = torus.minor.square();
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
        let CurveForm::Circle {
            center: circle_center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), edge.curve, dec)?
        else {
            return Ok(None);
        };
        let offset = circle_center.minus(center);
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
            add_cut(
                &mut tube,
                Phase::new(&radius(&radius_sq) - major, offset.dot(axis), dec)?,
                dec,
            )?;
        } else if sign(&normal.dot(axis), dec, "boolean.classify-ambiguous")? == Sign::Zero {
            for zero in [
                offset.dot(axis),
                offset.dot(&normal),
                offset.norm_squared() - major2.clone(),
                &radius_sq - &minor2,
            ] {
                if sign(&zero, dec, "boolean.classify-ambiguous")? != Sign::Zero {
                    return Ok(None);
                }
            }
            add_cut(
                &mut azimuth,
                Phase::new(offset.dot(first), offset.dot(&second), dec)?,
                dec,
            )?;
        } else {
            return Ok(None);
        }
    }
    let offset = model.geom().position(point).minus(center);
    let u = if azimuth.is_empty() {
        -1
    } else {
        let query = Phase::new(offset.dot(first), offset.dot(&second), dec)?;
        cell(&azimuth, &query, dec)?
    };
    let v = if tube.is_empty() {
        -1
    } else {
        let height = offset.dot(axis);
        let radial2 = offset.norm_squared() - height.square() / torus.axis.norm_squared();
        let query = Phase::new(&radius(&radial2) - major, height, dec)?;
        cell(&tube, &query, dec)?
    };
    Ok(Some(vec![u, v]))
}
