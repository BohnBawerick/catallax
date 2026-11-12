//! Finite queries over every retained chart, including identified rows.

use super::*;
use crate::rational_patch::exact_sign;
use crate::{Isometry, Point3, Pose, Vector3};
use catallax_numbers::{Ball, Real, Sign};

impl RationalAtlas {
    /// World axes giving an injective projection when every chart has the
    /// same nonsingular affine coordinate pair in original source parameters.
    /// The axes preserve the source's parameter orientation for winding.
    /// Other atlases return `None`; a chart-local projection alone is not enough.
    pub fn affine_projection(&self, dec: &mut Decider) -> Result<Option<[Vector3; 2]>, Refusal> {
        let patches = self
            .charts()
            .iter()
            .map(|chart| chart.patch())
            .collect::<Vec<_>>();
        // Every chart restricts the same single-valued source in its original
        // parameters. The common nonsingular pair recovers those parameters,
        // so equal projected positions cannot select different source points.
        RationalPatch::common_affine_projection(&patches, dec)
    }

    /// Exact parameter evaluation on the union of this logical face's charts.
    pub fn point_at(&self, at: &[Real; 2], dec: &mut Decider) -> Result<Point3, Refusal> {
        for chart in self.charts() {
            let mut inside = true;
            for (value, [lo, hi]) in at.iter().zip(chart.patch().domain()) {
                if exact_sign(&(value - &Real::from_rational(lo.clone())), dec)? == Sign::Negative
                    || exact_sign(&(&Real::from_rational(hi.clone()) - value), dec)?
                        == Sign::Negative
                {
                    inside = false;
                    break;
                }
            }
            if inside {
                return chart.patch().point_at(at, dec);
            }
        }
        Err(Refusal::new(
            "shapes.rational-atlas-domain",
            "parameter is outside this logical face",
        ))
    }

    /// Natural normal only at a proved finite point. Overlapping branches with
    /// incompatible normals refuse instead of selecting an arbitrary chart.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let mut normal: Option<Vector3> = None;
        for chart in self.charts() {
            if let Some(at) = chart.patch().parameters_at(point, dec)? {
                let [u, v] = chart.patch().derivatives_at(&at, dec)?;
                let next = u.cross(&v);
                if let Some(previous) = &normal {
                    if exact_sign(&previous.cross(&next).norm_squared(), dec)? != Sign::Zero
                        || exact_sign(&previous.dot(&next), dec)? != Sign::Positive
                    {
                        return Err(unproven("the point has incompatible atlas normals"));
                    }
                } else {
                    normal = Some(next);
                }
            }
        }
        normal.ok_or_else(|| {
            Refusal::new(
                "shapes.rational-atlas-domain",
                "normal point is not on the finite atlas",
            )
        })
    }

    /// A finite atlas has no globally defined implicit side. Incidence is
    /// answered exactly; an off-carrier side requires a trimmed body solver.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        for chart in self.charts() {
            if chart.patch().parameters_at(point, dec)?.is_some() {
                return Ok(Sign::Zero);
            }
        }
        Err(Refusal::new(
            "shapes.parametric-query",
            "a joined atlas needs a trimmed body side query",
        ))
    }

    /// Complete transverse line intersections, deduplicated exactly at chart
    /// seams and returned in increasing line parameter order.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        let mut hits: Vec<Point3> = Vec::new();
        for chart in self.charts() {
            for point in chart.patch().line_hits(origin, direction, dec)? {
                let mut position = hits.len();
                let mut duplicate = false;
                for (index, previous) in hits.iter().enumerate() {
                    if exact_sign(&point.minus(previous).norm_squared(), dec)? == Sign::Zero {
                        duplicate = true;
                        break;
                    }
                    if exact_sign(&point.minus(previous).dot(direction), dec)? == Sign::Negative {
                        position = index;
                        break;
                    }
                }
                if !duplicate {
                    hits.insert(position, point);
                }
            }
        }
        Ok(hits)
    }

    /// Complete conservative parameter cover. Separated chart branches refuse
    /// rather than turning a periodic seam into an arbitrary single inverse.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        let mut cover: Option<[Ball; 2]> = None;
        for chart in self.charts() {
            if let Some(next) = chart.patch().parameters_of_box(point, dec)? {
                if let Some(previous) = &mut cover {
                    if (0..2).any(|axis| {
                        previous[axis].hi() < next[axis].lo()
                            || next[axis].hi() < previous[axis].lo()
                    }) {
                        return Err(unproven("the box has separated atlas parameter branches"));
                    }
                    for axis in 0..2 {
                        previous[axis] = previous[axis].hull(&next[axis]);
                    }
                } else {
                    cover = Some(next);
                }
            }
        }
        Ok(cover)
    }

    /// Exact isometric transport retains the chart/gluing proofs.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }

    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        Self {
            data: Arc::new(Data {
                source: self.source().mapped(map),
                charts: self
                    .charts()
                    .iter()
                    .map(|chart| AtlasChart {
                        patch: chart.patch.mapped(map),
                        neighbors: chart.neighbors,
                    })
                    .collect(),
                periodic: [self.periodic_u(), self.periodic_v()],
            }),
        }
    }

    /// Positive identity for policy-independent definition cache reuse.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
}
