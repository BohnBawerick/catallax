//! A rolling-ball fillet face between two rational spline faces (M11).
//!
//! A ball of radius `r` rests on both supports. Its centre follows the spine,
//! where the supports' offsets `P + side r N/|N|` meet, and the spine is a
//! graph over one parameter `tau` of the first support. At each `tau` the
//! face's section is the circular arc round the centre from the first contact
//! to the second, the shorter way, on the side of `a + b` for the unit
//! directions `a` and `b` from the centre to the contacts. That face is not
//! rational, so the represented face is a rational B-spline fitted to it:
//!
//! - its two boundary rows are the supports composed with cubic Hermite
//!   curves of their contact parameters, so each rail lies exactly on its
//!   support and is the trim curve there;
//! - each section is a rational quadratic whose middle row and weight are
//!   cubic Hermite fits.
//!
//! Every proof is recomputed from exact definitions. On each span of `tau`
//! a parametric Krawczyk certificate proves a unique smooth spine point for
//! every `tau` in the span, and interval Taylor jets of that implicit solution
//! bound its derivatives to fourth order. The two-sided position bound
//! follows from the Hermite remainders, the node data and the weight errors,
//! as for the rolling envelope. Both surfaces share their parameters, so the
//! bound holds from each surface to the other. The fit starts from eight
//! uniform spans and halves each span that misses the target. The carrier's
//! normal is proved nonzero. Each contact is proved inside its support's
//! domain and the closest point of its support to the centre, so the ball
//! touches each support face once and enters neither.
//!
//! The pair search of M10 proves first that the offsets meet in exactly one
//! branch, and that it spans the requested range of `tau`. The face follows
//! that branch. Spans are proved on every core and adopted in order, so the
//! answer and the work do not depend on threads.

use crate::{BSplineSurface, Isometry, PatchPairIntersection, Point3, RationalPatch};
use catallax_numbers::{
    Ball, Decider, EvidenceBound, EvidenceUnit, MPoly, Rational, Real, Refusal, Sign,
};
use std::sync::Arc;

#[path = "spline_fillet/bound.rs"]
mod bound;
#[path = "spline_fillet/carrier.rs"]
mod carrier;
#[path = "spline_fillet/clearance.rs"]
mod clearance;
#[path = "spline_fillet/cut.rs"]
mod cut;
#[path = "spline_fillet/jets.rs"]
mod jets;
#[path = "spline_fillet/order.rs"]
mod order;
#[path = "spline_fillet/surface.rs"]
mod surface;

use jets::{Exact, Jet, Support};

/// The fitted face's target two-sided deviation, in nanometres.
pub const SPLINE_FILLET_TARGET_NM: i64 = 100;
/// The kernel's hard cap; a fit above it refuses with the bound it reached.
const HARD_CAP_NM: i64 = 1_000;
/// Uniform spans of `tau` the fit starts from.
const FIRST_SPANS: usize = 8;
/// Most halvings of one starting span.
const MAX_HALVINGS: u32 = 5;

/// What a fit aims for and accepts: the target and cap of its two-sided
/// bound, in millimetres, and the most halvings of a starting span.
struct Limits {
    target: Rational,
    cap: Rational,
    halvings: u32,
}

impl Limits {
    fn kernel() -> Self {
        let millimetres = |nm: i64| Rational::from_ratio(nm, 1_000_000).expect("nonzero");
        Self {
            target: millimetres(SPLINE_FILLET_TARGET_NM),
            cap: millimetres(HARD_CAP_NM),
            halvings: MAX_HALVINGS,
        }
    }
}

/// A definition the builder cannot fit.
pub const SPLINE_FILLET_DEFINITION: &str = "fillet.spline-definition";
/// The fitted face misses the 1 um cap; the detail names the bound reached.
pub const SPLINE_FILLET_DEVIATION: &str = "fillet.spline-deviation";
/// The spine was not proved a single smooth graph over `tau`.
pub const SPLINE_FILLET_SPINE: &str = "fillet.spline-spine";
/// A support of more than one knot span.
pub const SPLINE_FILLET_SUPPORT_SPANS: &str = "fillet.spline-support-spans";
/// A contact was not proved its support's closest point to the centre.
pub const SPLINE_FILLET_CLEARANCE: &str = "fillet.spline-clearance";
/// The fitted carrier's normal was not proved nonzero.
pub const SPLINE_FILLET_IRREGULAR: &str = "fillet.spline-carrier-irregular";
/// A fitted rail was not proved inside its support's domain.
pub const SPLINE_FILLET_RAIL_DOMAIN: &str = "fillet.spline-rail-domain";
/// A plane's cut of the face was not proved one crossing of every section,
/// ending on both rails at fit nodes.
pub const SPLINE_FILLET_CUT: &str = "fillet.spline-cut";

pub use cut::FilletCut;

use SPLINE_FILLET_SUPPORT_SPANS as SUPPORT_SPANS;

fn definition(detail: impl Into<String>) -> Refusal {
    Refusal::new(SPLINE_FILLET_DEFINITION, detail)
}
fn spine(detail: impl Into<String>) -> Refusal {
    Refusal::new(SPLINE_FILLET_SPINE, detail)
}

/// The exact rolling-ball face of a ball of radius `radius` between two
/// rational patches, over `tau` in `range`, where `tau` is the first patch's
/// parameter `axis`.
///
/// `sides[k]` places the ball's centre at `P + sides[k] r N/|N|` on patch
/// `k`, with `N = P_u x P_v`. With outward normals, a ball inside the
/// material of a convex edge has both sides negative.
#[derive(Clone, Debug)]
pub struct SplineFilletDefinition {
    supports: [RationalPatch; 2],
    radius: Rational,
    sides: [Sign; 2],
    axis: usize,
    range: [Rational; 2],
    borders: Vec<RailBorder>,
}

/// Where a rail crosses one of its support's parameter lines: support
/// `support`'s parameter `parameter` equals `value`. The fit puts a node of
/// `tau` there and gives the rail that parameter exactly, so the rail meets
/// the line at a rational point.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RailBorder {
    /// The support, first or second.
    pub support: usize,
    /// The support's parameter, first or second.
    pub parameter: usize,
    /// The parameter's value on the line.
    pub value: Rational,
}

impl SplineFilletDefinition {
    /// Checks the radius, sides, axis and range; the supports are checked
    /// when the face is fitted.
    pub fn new(
        supports: [RationalPatch; 2],
        radius: Rational,
        sides: [Sign; 2],
        axis: usize,
        range: [Rational; 2],
    ) -> Result<Self, Refusal> {
        if radius.sign() != Sign::Positive {
            return Err(definition("the ball radius must be positive"));
        }
        if sides.contains(&Sign::Zero) {
            return Err(definition("each ball side must be positive or negative"));
        }
        if axis > 1 {
            return Err(definition("tau is one of the first patch's two parameters"));
        }
        let domain = &supports[0].domain()[axis];
        if range[0] >= range[1] || range[0] < domain[0] || range[1] > domain[1] {
            return Err(definition(
                "the range of tau must be a nonempty part of the first patch's domain",
            ));
        }
        Ok(Self {
            supports,
            radius,
            sides,
            axis,
            range,
            borders: Vec::new(),
        })
    }
    /// Adds the lines where the rails must meet their supports' borders.
    /// Each must lie inside the support's domain; a line of `tau` itself
    /// must lie inside the range.
    pub fn with_borders(mut self, borders: Vec<RailBorder>) -> Result<Self, Refusal> {
        for border in &borders {
            if border.support > 1 || border.parameter > 1 {
                return Err(definition("a rail border names a support and a parameter"));
            }
            let side = &self.supports[border.support].domain()[border.parameter];
            if border.value < side[0] || border.value > side[1] {
                return Err(definition("a rail border must lie in its support's domain"));
            }
            if border.support == 0
                && border.parameter == self.axis
                && (border.value < self.range[0] || border.value > self.range[1])
            {
                return Err(definition("a border of tau must lie in the range"));
            }
        }
        self.borders = borders;
        Ok(self)
    }
    /// The two supports, first and second.
    pub fn supports(&self) -> &[RationalPatch; 2] {
        &self.supports
    }
    /// The ball's radius.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }
    /// The side of each support the ball's centre lies on.
    pub fn sides(&self) -> [Sign; 2] {
        self.sides
    }
    /// The first support's parameter that is `tau`.
    pub fn axis(&self) -> usize {
        self.axis
    }
    /// The range of `tau` the face covers.
    pub fn range(&self) -> &[Rational; 2] {
        &self.range
    }
    /// The lines where the rails meet their supports' borders.
    pub fn borders(&self) -> &[RailBorder] {
        &self.borders
    }
}

/// One span of a rail: the Bernstein coefficients of the support's two
/// parameters as cubics in `tau` over the span. The rail is the support at
/// these parameters, exactly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RailSpan {
    /// The span of `tau`.
    pub tau: [Rational; 2],
    /// The cubic Bernstein coefficients of the support's first and second
    /// parameter over the span.
    pub parameters: [[Rational; 4]; 2],
}

/// Fit data at one node of `tau`: value and slope of each fitted row.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Node {
    tau: Rational,
    values: [Rational; ROWS],
    slopes: [Rational; ROWS],
}

/// The exact value and slope of each fitted row at one node.
type NodeData = [[Ball; 2]; ROWS];

/// Fitted rows: the first patch's other parameter, the second patch's two
/// parameters, the middle point and the middle weight.
const ROWS: usize = 7;
const MIDDLE: usize = 3;
const WEIGHT: usize = 6;

/// A checked rolling-ball fillet face between two rational spline faces.
#[derive(Clone, Debug)]
pub struct SplineFillet {
    data: Arc<Data>,
}

/// One placed copy: everything the fit proved, shared by every copy, and
/// the carrier and supports moved by the copy's frame.
#[derive(Debug)]
struct Data {
    fit: Arc<Fitted>,
    carrier: BSplineSurface,
    supports: [RationalPatch; 2],
    /// From the frame of the fit to this copy's.
    frame: Isometry,
}

/// What the fit proved, in the frame it was fitted in.
#[derive(Debug)]
struct Fitted {
    definition: SplineFilletDefinition,
    nodes: Vec<Node>,
    estimates: Vec<[f64; 3]>,
    carrier: BSplineSurface,
    rows: Vec<carrier::SpanRows>,
    /// The carrier's homogeneous controls, enclosed for queries.
    net: Vec<[Ball; 4]>,
    rails: [Vec<RailSpan>; 2],
    borders: Vec<Rational>,
    bound: EvidenceBound,
    spine: PatchPairIntersection,
}

impl SplineFillet {
    /// Fits Hermite spans of `tau`, starting from eight uniform spans and
    /// halving each span whose recomputed bound misses the 100 nm target, at
    /// most five times. A fit that meets only the 1 um cap is accepted and
    /// reported by its bound; anything larger refuses with the bound it
    /// reached.
    pub fn fit(definition: SplineFilletDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::fit_within(definition, &Limits::kernel(), dec)
    }

    fn fit_within(
        definition: SplineFilletDefinition,
        limits: &Limits,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let exact = Exact::new(&definition, dec)?;
        let spine = carrier::spine(&definition, dec)?;
        let seed = carrier::seed(&definition, &spine)?;
        let borders = definition
            .borders
            .iter()
            .map(|border| carrier::border(&exact, &definition, border, seed))
            .collect::<Result<Vec<_>, _>>()?;
        let fit = bound::refine(&exact, &definition.range, seed, &borders, limits, dec)?;
        if fit.upper > limits.cap {
            return Err(Refusal::new(
                SPLINE_FILLET_DEVIATION,
                format!(
                    "the fitted face's two-sided bound reaches {} nm against the {} nm cap",
                    nanometres(&fit.upper),
                    nanometres(&limits.cap)
                ),
            ));
        }
        let (carrier, rails, rows) = carrier::build(&exact, &fit.nodes, dec)?;
        carrier::rails_inside(&rails, &definition.supports, dec)?;
        carrier::regular(&rows, dec)?;
        clearance::check(&exact, &fit.nodes, &fit.estimates, &fit.spans, dec)?;
        let net = surface::net(&carrier);
        let supports = definition.supports.clone();
        Ok(Self {
            data: Arc::new(Data {
                carrier: carrier.clone(),
                supports,
                frame: Isometry::identity(),
                fit: Arc::new(Fitted {
                    definition,
                    nodes: fit.nodes,
                    estimates: fit.estimates,
                    carrier,
                    rows,
                    net,
                    rails,
                    borders: borders.into_iter().map(|(tau, _)| tau).collect(),
                    bound: EvidenceBound {
                        upper: fit.upper,
                        unit: EvidenceUnit::Millimetres,
                    },
                    spine,
                }),
            }),
        })
    }

    /// The definition the face was fitted to, in the frame of the fit.
    pub fn definition(&self) -> &SplineFilletDefinition {
        &self.data.fit.definition
    }
    /// The represented face over `[tau0, tau1] x [0, 1]`: `tau` along the
    /// spine, the second parameter from the first rail to the second.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// The two supports, placed with the face.
    pub fn supports(&self) -> &[RationalPatch; 2] {
        &self.data.supports
    }
    /// The proved two-sided bound between the carrier and the exact face.
    pub fn bound(&self) -> &EvidenceBound {
        &self.data.fit.bound
    }
    /// Fitted spans of `tau`.
    pub fn spans(&self) -> usize {
        self.data.fit.nodes.len() - 1
    }
    /// The rail on support `side`, span by span, as that support's
    /// parameters.
    pub fn rail(&self, side: usize) -> &[RailSpan] {
        &self.data.fit.rails[side]
    }
    /// The M10 meeting of the two offsets that the face follows.
    pub fn spine(&self) -> &PatchPairIntersection {
        &self.data.fit.spine
    }
    /// The node of `tau` at each of the definition's rail borders, in order.
    pub fn borders(&self) -> &[Rational] {
        &self.data.fit.borders
    }
    /// The nodes of `tau`, the carrier's knots between its spans.
    pub fn nodes(&self) -> impl Iterator<Item = &Rational> {
        self.data.fit.nodes.iter().map(|node| &node.tau)
    }
    /// Where the plane `normal . x = offset` cuts the face, or `None` when
    /// the whole face lies on one side of it. A cut crosses every curve of
    /// constant second parameter exactly once, and ends on both rails at fit
    /// nodes.
    pub fn plane_cut(
        &self,
        normal: [Rational; 3],
        offset: Rational,
        dec: &mut Decider,
    ) -> Result<Option<FilletCut>, Refusal> {
        let (local_normal, local_offset) = surface::local_plane(&self.data.frame, &normal, &offset)
            .ok_or_else(|| {
                Refusal::new(
                    SPLINE_FILLET_CUT,
                    "a face moved by an irrational frame has no rational cutting plane",
                )
            })?;
        // The proof runs in the fit's frame; the cut reports the caller's
        // plane, which its ends satisfy on the placed face.
        let cut = cut::plane(
            &self.data.fit.rows,
            &self.data.fit.nodes,
            local_normal,
            local_offset,
            dec,
        )?;
        Ok(cut.map(|cut| FilletCut {
            normal,
            offset,
            ..cut
        }))
    }
    /// Runs the two-sided bound again on the retained fit data.
    pub fn remeasure(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let fit = &self.data.fit;
        let exact = Exact::new(&fit.definition, dec)?;
        bound::recompute(&exact, &fit.nodes, &fit.estimates, dec)
    }
}

/// A bound in millimetres as whole nanometres, rounded up.
fn nanometres(upper: &Rational) -> String {
    let scaled = upper * &Rational::from_i64(1_000_000);
    Rational::from_integer(scaled.ceil()).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    /// The M27 loft edge's two bilinear sides, from their corners.
    fn loft_edge(dec: &mut Decider) -> [RationalPatch; 2] {
        let mut patch = |corners: [[i64; 3]; 4]| {
            let controls = corners
                .map(|[x, y, z]| {
                    Point3::new(
                        Rational::from_i64(x).into(),
                        Rational::from_i64(y).into(),
                        Rational::from_i64(z).into(),
                    )
                })
                .to_vec();
            let knots = vec![
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::one(),
            ];
            let carrier = BSplineSurface::new(
                [1, 1],
                [knots.clone(), knots],
                [2, 2],
                controls,
                vec![Real::one(); 4],
                dec,
            )
            .unwrap();
            RationalPatch::new(carrier, dec).unwrap()
        };
        [
            patch([[-1, -1, 0], [0, -1, 2], [1, -1, 0], [1, 0, 2]]),
            patch([[1, -1, 0], [1, 0, 2], [1, 1, 0], [0, 2, 2]]),
        ]
    }

    /// A review's overshoot: the Hermite cubic of `epsilon + t^4` on `[0, h]`
    /// dips to `epsilon - h^4/27`, between nodes whose values are inside.
    /// Built as a rail is built, it is refused when that dip leaves the
    /// domain, and proved inside by halving when only its control polygon
    /// does.
    #[test]
    fn a_rail_that_overshoots_its_domain_is_refused() {
        let q = |n: i64, d: i64| Rational::from_ratio(n, d).unwrap();
        let h = q(1, 2);
        let hermite = |epsilon: Rational| -> [Rational; 4] {
            let values = [epsilon.clone(), &epsilon + &h.pow(4).unwrap()];
            let slopes = [Rational::zero(), q(4, 1) * h.pow(3).unwrap()];
            carrier::cubic(
                [&Rational::zero(), &h],
                [&values[0], &values[1]],
                [&slopes[0], &slopes[1]],
            )
            .try_into()
            .unwrap()
        };
        let domain = [Rational::zero(), Rational::one()];
        let mut dec = Decider::default();
        // Here h^4/27 = 1/432, so a thousandth dips below zero.
        let leaving = hermite(q(1, 1000));
        assert_eq!(
            carrier::cubic_inside(&leaving, &domain, &mut dec).unwrap(),
            Some(carrier::Escape::Leaves)
        );
        // A hundredth stays inside, though its polygon dips to 1/100 - 1/48.
        let inside = hermite(q(1, 100));
        assert_eq!(inside[2].sign(), Sign::Negative);
        assert_eq!(
            carrier::cubic_inside(&inside, &domain, &mut dec).unwrap(),
            None
        );
        // The fit's check refuses a returned rail with the dip, by its code.
        let supports = loft_edge(&mut dec);
        let rail = RailSpan {
            tau: [Rational::zero(), h.clone()],
            parameters: [leaving, [q(1, 2), q(1, 2), q(1, 2), q(1, 2)]],
        };
        let refusal =
            carrier::rails_inside(&[vec![rail], Vec::new()], &supports, &mut dec).unwrap_err();
        assert_eq!(refusal.code, SPLINE_FILLET_RAIL_DOMAIN);
        assert!(refusal.detail.contains("leaves"), "{}", refusal.detail);
    }

    /// With no halvings and a 1 nm cap, the eight starting spans' bound is
    /// the deviation the refusal reports.
    #[test]
    fn a_fit_above_the_cap_refuses_with_its_bound() {
        let mut dec = Decider::default();
        let definition = SplineFilletDefinition::new(
            loft_edge(&mut dec),
            Rational::from_ratio(1, 10).unwrap(),
            [Sign::Negative, Sign::Negative],
            1,
            [Rational::from_ratio(1, 32).unwrap(), Rational::one()],
        )
        .unwrap();
        let nanometre = Rational::from_ratio(1, 1_000_000).unwrap();
        let limits = Limits {
            target: nanometre.clone(),
            cap: nanometre,
            halvings: 0,
        };
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            SplineFillet::fit_within(definition, &limits, &mut dec)
        });
        let refusal = result.expect_err("eight spans miss a 1 nm cap");
        assert_eq!(refusal.code, SPLINE_FILLET_DEVIATION);
        assert_eq!(
            refusal.detail,
            "the fitted face's two-sided bound reaches 463 nm against the 1 nm cap"
        );
        assert!(!work.exhausted);
    }
}
