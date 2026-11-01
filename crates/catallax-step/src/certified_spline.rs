//! Complete finite spline-prism export, with shared rounded geometry.
//!
//! Independent coordinate rounding would break planar incidence and periodic
//! C1 seams. Round one planar Hermite node/derivative table and one sweep instead.
//! All STEP curves, surfaces, caps and added chart vertices share those values.

use std::io::Write;

use catallax_numbers::{
    Answer, Answered, Decider, Integer, Outcome, Rational, Real, Refusal, Sign, StrictSign,
};
use catallax_shapes::{Point3, Surface, Vector3, certify_convex_cubic_span};
use catallax_topology::{BodyId, Model, Sense, SplinePrism};

use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use crate::{StepFaceRegion, StepReport};

/// Output accuracy, independent of modelling sameness and the body's fitting
/// bound. This exporter preserves exact rational directions and terminating
/// parameter knots of the first checked polynomial spline-prism family.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct CertifiedSplineStepOptions {
    /// Product name written to the AP242 file.
    pub product_name: String,
    /// Required uniform two-way file-to-represented-model bound, in millimetres.
    /// A failed bound refuses; it is never enlarged to finish the export.
    pub file_tolerance: Rational,
    /// Decimal places for the common rounded origin. Local coefficient grids
    /// scale with their exact integer directions; shared handles and knots can
    /// require additional decimal places to retain exact incidence and C1.
    pub decimal_places: u32,
}

impl Default for CertifiedSplineStepOptions {
    fn default() -> Self {
        Self {
            product_name: "catallax".into(),
            file_tolerance: Rational::from_ratio(1, 10_000_000).unwrap(),
            decimal_places: 15,
        }
    }
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("step.spline-proof", detail)
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}
fn positive(value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(value, dec)? == Sign::Positive {
        Ok(())
    } else {
        Err(bad(
            "strict positive orientation or regularity condition failed",
        ))
    }
}

fn integer_direction(vector: &Vector3) -> Result<Vector3, Refusal> {
    let values = vector
        .rational_components()
        .ok_or_else(|| bad("exact rational direction frame required"))?;
    let mut denominator = Integer::one();
    for value in values {
        let gcd = denominator.gcd(value.denominator());
        denominator = denominator.div_rem_floor(&gcd).unwrap().0 * value.denominator().clone();
    }
    let values = values.map(|value| {
        value.numerator() * &denominator.div_rem_floor(value.denominator()).unwrap().0
    });
    let divisor = values
        .iter()
        .fold(Integer::zero(), |gcd, value| gcd.gcd(value));
    if divisor.is_zero() {
        return Err(bad("export direction vanished"));
    }
    let [x, y, z] = values.map(|value| {
        Real::from_rational(Rational::from_integer(
            value.div_rem_floor(&divisor).unwrap().0,
        ))
    });
    Ok(Vector3::new(x, y, z))
}

struct Grid {
    unit: Rational,
    requested: Rational,
}
impl Grid {
    fn new(places: u32, dec: &Decider) -> Result<Self, Refusal> {
        if places > dec.policy().max_precision_bits / 4 {
            return Err(Refusal::new(
                "step.budget-exhausted",
                "decimal grid exceeds the caller's precision allowance",
            ));
        }
        let unit = Rational::one() / Rational::from_integer(Integer::from_i64(10).pow(places));
        let requested = &unit / &Rational::from_i64(16);
        Ok(Self { unit, requested })
    }
    fn for_direction(places: u32, direction: &Vector3, dec: &Decider) -> Result<Self, Refusal> {
        let bits = direction
            .rational_components()
            .unwrap()
            .into_iter()
            .map(|q| q.numerator().bit_length())
            .max()
            .unwrap();
        // 30103/100000 is an upper bound on log10(2); the extra decade covers
        // the three-coordinate L1 scale. This choice precedes all rounding.
        let extra = bits
            .checked_mul(30103)
            .and_then(|b| u32::try_from(b.div_ceil(100000) + 1).ok())
            .ok_or_else(|| bad("integer direction scale exceeds the output grid"))?;
        Self::new(
            places
                .checked_add(extra)
                .ok_or_else(|| bad("decimal grid overflow"))?,
            dec,
        )
    }
    fn round(&self, value: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
        dec.check_deadline()?;
        let midpoint = if let Some(value) = value.rational_ref() {
            value.clone()
        } else {
            let ball = dec.enclose(value, &self.requested)?;
            (&ball.lo() + &ball.hi()).mul_pow2(-1)
        };
        let scaled = &midpoint / &self.unit;
        let integer = (&scaled.abs() + &Rational::from_ratio(1, 2).unwrap()).floor();
        let integer = if scaled.sign() == Sign::Negative {
            -integer
        } else {
            integer
        };
        Ok(Real::from_rational(
            &Rational::from_integer(integer) * &self.unit,
        ))
    }
    fn point(&self, point: &Point3, dec: &mut Decider) -> Result<Point3, Refusal> {
        Ok(Point3::new(
            self.round(&point.x, dec)?,
            self.round(&point.y, dec)?,
            self.round(&point.z, dec)?,
        ))
    }
}

struct Node {
    point: [Real; 2],
    derivative: [Real; 2],
}
struct Rounded {
    controls: Vec<[Point3; 2]>,
    centers: [Point3; 2],
    axis: Vector3,
    reference: Vector3,
    bound: Rational,
}

fn difference_bound(
    a: &Point3,
    b: &Point3,
    bits: u32,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let mut result = Rational::zero();
    for (a, b) in [&a.x, &a.y, &a.z].into_iter().zip([&b.x, &b.y, &b.z]) {
        let difference = a - b;
        let upper = if let Some(value) = difference.rational_ref() {
            value.abs()
        } else {
            let ball = dec.enclose_budgeted_at(&difference, bits)?;
            ball.lo().abs().max(ball.hi().abs())
        };
        result = result + upper;
    }
    Ok(result)
}
fn same_point(a: &Point3, b: &Point3) -> bool {
    matches!((a.rational_components(), b.rational_components()), (Some(a), Some(b)) if a == b)
}
fn same_vector(a: &Vector3, b: &Vector3) -> bool {
    matches!((a.rational_components(), b.rational_components()), (Some(a), Some(b)) if a == b)
}

fn round_geometry(
    view: &SplinePrism,
    options: &CertifiedSplineStepOptions,
    dec: &mut Decider,
) -> Result<Rounded, Refusal> {
    let spline = view.spline();
    let carrier = spline.carrier();
    let Surface::EllipticCylinder(original) = spline.definition().original() else {
        unreachable!("checked spline-prism source");
    };
    let axis = integer_direction(&original.axis)?;
    let reference = integer_direction(&original.major_dir)?;
    let minor = integer_direction(&axis.cross(&reference))?;
    let u2 = reference.norm_squared();
    let v2 = minor.norm_squared();
    let n2 = axis.norm_squared();
    let center_grid = Grid::new(options.decimal_places, dec)?;
    let u_grid = Grid::for_direction(options.decimal_places, &reference, dec)?;
    let v_grid = Grid::for_direction(options.decimal_places, &minor, dec)?;
    let n_grid = Grid::for_direction(options.decimal_places, &axis, dec)?;
    let [lo, hi] = spline.axial();
    let model_center = original.origin.plus(&spline.axis().scaled(lo));
    let model_sweep = spline.axis().scaled(&(hi - lo));
    let model_upper_center = model_center.plus(&model_sweep);
    let center = center_grid.point(&model_center, dec)?;
    let height = n_grid.round(&(&model_sweep.dot(&axis) / &n2), dec)?;
    positive(&height, dec)?;
    let sweep = axis.scaled(&height);
    let upper_center = center.plus(&sweep);
    let widths: Vec<_> = carrier.spans(0).map(|span| span[1] - span[0]).collect();
    let source = carrier.controls();
    let mut nodes = Vec::with_capacity(widths.len());
    for (index, width) in widths.iter().enumerate() {
        let point = &source[6 * index];
        let relative = point.minus(&model_center);
        let derivative = source[6 * index + 2].minus(point);
        let width = Real::from_rational(width.clone());
        nodes.push(Node {
            point: [
                u_grid.round(&(&relative.dot(&reference) / &u2), dec)?,
                v_grid.round(&(&relative.dot(&minor) / &v2), dec)?,
            ],
            derivative: [
                u_grid.round(&(&derivative.dot(&reference) / &(&width * &u2)), dec)?,
                v_grid.round(&(&derivative.dot(&minor) / &(&width * &v2)), dec)?,
            ],
        });
    }
    let world = |point: [&Real; 2]| {
        center
            .plus(&reference.scaled(point[0]))
            .plus(&minor.scaled(point[1]))
    };
    let mut controls = Vec::with_capacity(carrier.control_count()[0]);
    let first = world(nodes[0].point.each_ref());
    let upper = first.plus(&sweep);
    controls.push([first, upper]);
    for (index, width) in widths.iter().enumerate() {
        let next = (index + 1) % nodes.len();
        let width = Real::from_rational(width.clone());
        let a = std::array::from_fn::<_, 2, _>(|i| {
            &nodes[index].point[i] + &(&nodes[index].derivative[i] * &width)
        });
        let b = std::array::from_fn::<_, 2, _>(|i| {
            &nodes[next].point[i] - &(&nodes[next].derivative[i] * &width)
        });
        for point in [
            world(a.each_ref()),
            world(b.each_ref()),
            world(nodes[next].point.each_ref()),
        ] {
            let upper = point.plus(&sweep);
            controls.push([point, upper]);
        }
    }
    if controls.len() != carrier.control_count()[0] || controls.len() * 2 != source.len() {
        return Err(bad(
            "rounded control cover does not match the whole source surface",
        ));
    }
    let bits = options
        .decimal_places
        .saturating_add(16)
        .saturating_mul(4)
        .max(128)
        .min(dec.policy().max_precision_bits);
    let mut bound = difference_bound(&model_center, &center, bits, dec)?.max(difference_bound(
        &model_upper_center,
        &upper_center,
        bits,
        dec,
    )?);
    for (file, model) in controls.iter().flatten().zip(source) {
        bound = bound.max(difference_bound(model, file, bits, dec)?);
    }
    if bound > options.file_tolerance {
        return Err(Refusal::new(
            "step.tolerance-not-reached",
            format!(
                "shared rounded geometry reaches {bound} mm, exceeding the fixed {} mm file tolerance",
                options.file_tolerance
            ),
        ));
    }
    // The same knot cover and a convex-combination correspondence bound the
    // entire surfaces. Check planar incidence, periodic C1, positive curvature
    // and a nonzero-tangent homotopy to the original single-winding profile.
    let mut previous_tangent: Option<Vector3> = None;
    let mut first_tangent = None;
    for (index, width) in widths.iter().enumerate() {
        let base = 3 * index;
        let points = [
            &controls[base][0],
            &controls[base + 1][0],
            &controls[base + 2][0],
            &controls[base + 3][0],
        ];
        let old = [
            &source[6 * index],
            &source[6 * index + 2],
            &source[6 * index + 4],
            &source[6 * index + 6],
        ];
        let halfplane = old[1].minus(old[0]).plus(&old[3].minus(old[2]));
        for pair in old.windows(2) {
            positive(&pair[1].minus(pair[0]).dot(&halfplane), dec)?;
        }
        certify_convex_cubic_span(points, &halfplane, spline.axis(), dec)?;
        let inverse_width = Real::from_rational(Rational::one() / width.clone());
        let start = points[1].minus(points[0]).scaled(&inverse_width);
        let end = points[3].minus(points[2]).scaled(&inverse_width);
        if previous_tangent
            .as_ref()
            .is_some_and(|previous| !same_vector(previous, &start))
        {
            return Err(bad("rounded spline lost an interior C1 seam"));
        }
        if first_tangent.is_none() {
            first_tangent = Some(start);
        }
        previous_tangent = Some(end);
    }
    if !same_point(&controls[0][0], &controls.last().unwrap()[0])
        || !same_vector(
            first_tangent.as_ref().unwrap(),
            previous_tangent.as_ref().unwrap(),
        )
    {
        return Err(bad(
            "rounded spline lost its periodic position or tangent seam",
        ));
    }
    for pair in &controls {
        if sign(&pair[0].minus(&center).dot(&axis), dec)? != Sign::Zero
            || sign(&pair[1].minus(&upper_center).dot(&axis), dec)? != Sign::Zero
            || !same_vector(&pair[1].minus(&pair[0]), &sweep)
        {
            return Err(bad("rounded spline detached a cap or sweep boundary"));
        }
    }
    let signed = match spline.definition().side() {
        catallax_shapes::OffsetSide::Inward => -spline.definition().distance(),
        catallax_shapes::OffsetSide::Outward => spline.definition().distance().clone(),
    };
    let radius =
        &(&original.minor_sq.sqrt() + &signed) - &Real::from_rational(spline.bound().upper.clone());
    let radius_lower = dec.enclose_budgeted_at(&radius, bits)?.lo();
    if &bound * &Rational::from_i64(2) >= radius_lower {
        return Err(bad(
            "rounded cap centre lacks a strict interior correspondence margin",
        ));
    }
    // Both closed convex profiles contain their corresponding centres: the
    // homotopy stays clear by radius_lower-2*bound. Radial cap fans and the
    // common sweep therefore cover both entire solids with this same bound.
    Ok(Rounded {
        controls,
        centers: [center, upper_center],
        axis,
        reference,
        bound,
    })
}

fn numeral(value: &Real) -> Result<String, Refusal> {
    let rational = value
        .rational_ref()
        .ok_or_else(|| bad("written geometry is not exact rational data"))?;
    let decimal = rational
        .to_decimal_exact()
        .ok_or_else(|| bad("shared geometry requires a nonterminating STEP numeral"))?;
    Ok(crate::number::decimal(&decimal))
}
fn coordinates(graph: &mut Graph, name: &'static str, values: [&Real; 3]) -> Result<Id, Refusal> {
    Ok(graph.add(
        name,
        vec![
            text("''"),
            text(format!(
                "({},{},{})",
                numeral(values[0])?,
                numeral(values[1])?,
                numeral(values[2])?
            )),
        ],
    ))
}
fn point(graph: &mut Graph, point: &Point3) -> Result<Id, Refusal> {
    coordinates(graph, "CARTESIAN_POINT", [&point.x, &point.y, &point.z])
}
fn direction(graph: &mut Graph, direction: &Vector3) -> Result<Id, Refusal> {
    coordinates(
        graph,
        "DIRECTION",
        [&direction.x, &direction.y, &direction.z],
    )
}
fn knot_fields(knots: &[Rational]) -> Result<(Field, Field), Refusal> {
    let mut multiplicities = Vec::new();
    let mut values = Vec::new();
    let mut start = 0;
    while start < knots.len() {
        let mut end = start + 1;
        while end < knots.len() && knots[end] == knots[start] {
            end += 1;
        }
        multiplicities.push(text((end - start).to_string()));
        values.push(text(parameter(&knots[start])?));
        start = end;
    }
    Ok((Field::List(multiplicities), Field::List(values)))
}

/// The closed side is written as two open charts, split at its retained u=2
/// knot. A checked carrier ends a span at every quarter, so this slices existing
/// Bezier rows: nothing is refitted and the two halves reassemble exactly.
struct Charts {
    /// Index of the control row at the cut.
    cut: usize,
    /// Expanded U knots of each chart, clamped at the cut.
    knots: [Vec<Rational>; 2],
    /// Exact U domain of each chart.
    domains: [[Rational; 2]; 2],
}
fn open_charts(knots: &[Rational]) -> Result<Charts, Refusal> {
    let cut = Rational::from_i64(2);
    let first = knots
        .iter()
        .position(|knot| *knot == cut)
        .ok_or_else(|| bad("spline side lacks its retained u=2 knot"))?;
    let end = first + knots[first..].iter().take_while(|k| **k == cut).count();
    if first == 0 || end - first != 3 || end == knots.len() {
        return Err(bad("u=2 is not an interior Bezier span boundary"));
    }
    let lower = knots.first().unwrap().clone();
    let upper = knots.last().unwrap().clone();
    Ok(Charts {
        cut: first - 1,
        knots: [
            [&knots[..end], &knots[first..=first]].concat(),
            [&knots[first..=first], &knots[first..]].concat(),
        ],
        domains: [[lower, cut.clone()], [cut, upper]],
    })
}
fn parameter(value: &Rational) -> Result<String, Refusal> {
    let decimal = value.to_decimal_exact().ok_or_else(|| {
        bad("parameter knots must have an exact terminating decimal representation")
    })?;
    Ok(crate::number::decimal(&decimal))
}
fn pcurve_line(
    graph: &mut Graph,
    surface: Id,
    context: Id,
    origin: [String; 2],
    direction: &'static str,
) -> Id {
    let point = graph.add(
        "CARTESIAN_POINT",
        vec![text("''"), text(format!("({},{})", origin[0], origin[1]))],
    );
    let direction = graph.add("DIRECTION", vec![text("''"), text(direction)]);
    let vector = graph.add("VECTOR", vec![text("''"), direction.into(), text("1.")]);
    let line = graph.add("LINE", vec![text("''"), point.into(), vector.into()]);
    let representation = graph.add(
        "DEFINITIONAL_REPRESENTATION",
        vec![text("''"), refs([line]), context.into()],
    );
    graph.add(
        "PCURVE",
        vec![text("''"), surface.into(), representation.into()],
    )
}
fn edge_curve(graph: &mut Graph, ends: [Id; 2], curve: Id, pcurves: Vec<Id>) -> Id {
    let wrapper = graph.add(
        "SURFACE_CURVE",
        vec![text("''"), curve.into(), refs(pcurves), text(".CURVE_3D.")],
    );
    graph.add(
        "EDGE_CURVE",
        vec![
            text("''"),
            ends[0].into(),
            ends[1].into(),
            wrapper.into(),
            boolean(true),
        ],
    )
}
fn advanced_face(
    graph: &mut Graph,
    name: &str,
    uses: &[(Id, bool)],
    surface: Id,
    same: bool,
) -> Field {
    let oriented: Vec<_> = uses
        .iter()
        .map(|(edge, forward)| {
            graph.add(
                "ORIENTED_EDGE",
                vec![
                    text("''"),
                    text("*"),
                    text("*"),
                    (*edge).into(),
                    boolean(*forward),
                ],
            )
        })
        .collect();
    let loop_ = graph.add("EDGE_LOOP", vec![text("''"), refs(oriented)]);
    let bound = graph.add(
        "FACE_OUTER_BOUND",
        vec![text("''"), loop_.into(), boolean(true)],
    );
    graph
        .add(
            "ADVANCED_FACE",
            vec![
                text(string(name)),
                refs([bound]),
                surface.into(),
                boolean(same),
            ],
        )
        .into()
}

/// Writes the first checked finite polynomial spline-prism family as genuine
/// `B_SPLINE_SURFACE_WITH_KNOTS` and shared `B_SPLINE_CURVE_WITH_KNOTS` entities.
/// The one model side becomes two file faces with shared boundaries, recorded
/// in `StepReport::face_map`. This is export topology; the point sets and the
/// model stay as they are.
/// The returned `max_deviation` bounds both directions over the complete solid,
/// relative to the represented model. Compose it separately with that model's
/// definition certificate. It is not the modelling tolerance or a refit bound.
/// All geometry and topology proofs finish before any bytes are written.
pub fn write_spline_prism_step<W: Write>(
    model: &Model,
    body: BodyId,
    options: &CertifiedSplineStepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    if options.file_tolerance.sign() == Sign::Negative {
        return Err(bad("negative file tolerance"));
    }
    let view = SplinePrism::check(model, body, dec)?;
    let rounded = round_geometry(&view, options, dec)?;
    let charts = open_charts(view.spline().carrier().knots()[0])?;
    let mut graph = Graph::default();
    let mut controls: Vec<[Id; 2]> = Vec::with_capacity(rounded.controls.len());
    for (index, pair) in rounded.controls.iter().enumerate() {
        if index + 1 == rounded.controls.len() {
            controls.push(controls[0]);
        } else {
            controls.push([point(&mut graph, &pair[0])?, point(&mut graph, &pair[1])?]);
        }
    }
    let rows = [&controls[..=charts.cut], &controls[charts.cut..]];
    let mut sides = Vec::with_capacity(2);
    let mut arcs = Vec::with_capacity(2);
    for (rows, knots) in rows.into_iter().zip(&charts.knots) {
        let (multiplicities, knots) = knot_fields(knots)?;
        sides.push(graph.add(
            "B_SPLINE_SURFACE_WITH_KNOTS",
            vec![
                text("''"),
                text("3"),
                text("1"),
                Field::List(rows.iter().map(|pair| refs(*pair)).collect()),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                boolean(false),
                multiplicities.clone(),
                text("(2,2)"),
                knots.clone(),
                text("(0.,1.)"),
                text(".UNSPECIFIED."),
            ],
        ));
        arcs.push([0, 1].map(|row| {
            graph.add(
                "B_SPLINE_CURVE_WITH_KNOTS",
                vec![
                    text("''"),
                    text("3"),
                    refs(rows.iter().map(|pair| pair[row])),
                    text(".UNSPECIFIED."),
                    boolean(false),
                    boolean(false),
                    multiplicities.clone(),
                    knots.clone(),
                    text(".UNSPECIFIED."),
                ],
            )
        }));
    }
    let reference = direction(&mut graph, &rounded.reference)?;
    let axis = direction(&mut graph, &rounded.axis)?;
    let faces = view.faces();
    let mut plane = |index| -> Result<_, Refusal> {
        let center = &rounded.centers[index];
        let origin = point(&mut graph, center)?;
        let same = model.face(faces[index]).sense == Sense::Same;
        let normal_positive = if index == 0 { !same } else { same };
        let normal = if normal_positive {
            axis
        } else {
            direction(&mut graph, &rounded.axis.scaled(&Real::from_i64(-1)))?
        };
        let placement = graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), origin.into(), normal.into(), reference.into()],
        );
        Ok(graph.add("PLANE", vec![text("''"), placement.into()]))
    };
    let caps = [plane(0)?, plane(1)?];
    // Both charts of one U parameter value carry its corner rows: index 0 at
    // the periodic u=0/4 join and index 1 at the cut.
    let corners = [controls[0], controls[charts.cut]];
    let vertices =
        corners.map(|pair| pair.map(|row| graph.add("VERTEX_POINT", vec![text("''"), row.into()])));
    let context = graph.add(
        "GEOMETRIC_REPRESENTATION_CONTEXT",
        vec![text("''"), text("''"), text("2")],
    );
    let mut rings = Vec::with_capacity(2);
    for (index, surface) in sides.iter().enumerate() {
        let mut edges = [None, None];
        for row in 0..2 {
            let line = pcurve_line(
                &mut graph,
                *surface,
                context,
                ["0.".to_owned(), format!("{row}.")],
                "(1.,0.)",
            );
            edges[row] = Some(edge_curve(
                &mut graph,
                [vertices[index][row], vertices[1 - index][row]],
                arcs[index][row],
                vec![line],
            ));
        }
        rings.push(edges.map(Option::unwrap));
    }
    let join = |index: usize, graph: &mut Graph| -> Result<Id, Refusal> {
        let curve = graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("1"),
                refs(corners[index]),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                text("(2,2)"),
                text("(0.,1.)"),
                text(".UNSPECIFIED."),
            ],
        );
        let parameters = if index == 0 {
            [&charts.domains[0][0], &charts.domains[1][1]]
        } else {
            [&charts.domains[0][1], &charts.domains[1][0]]
        };
        let mut lines = Vec::with_capacity(2);
        for (surface, value) in sides.iter().zip(parameters) {
            lines.push(pcurve_line(
                graph,
                *surface,
                context,
                [parameter(value)?, "0.".to_owned()],
                "(0.,1.)",
            ));
        }
        Ok(edge_curve(graph, vertices[index], curve, lines))
    };
    let joins = [join(0, &mut graph)?, join(1, &mut graph)?];
    // `SplinePrism::check` proves the side has the surface's sense, travelling
    // the lower ring with increasing u and the upper ring against it. Each cap
    // uses its ring the opposite way, so one shell has one material side.
    let mut face_entities: Vec<Field> = Vec::with_capacity(4);
    let mut face_map = Vec::with_capacity(4);
    for (index, surface) in caps.into_iter().enumerate() {
        let face = faces[index];
        let uses = if index == 0 {
            [(rings[1][0], false), (rings[0][0], false)]
        } else {
            [(rings[0][1], true), (rings[1][1], true)]
        };
        let name = format!("face_{}", face.index());
        let same = model.face(face).sense == Sense::Same;
        face_entities.push(advanced_face(&mut graph, &name, &uses, surface, same));
        face_map.push(StepFaceRegion {
            face,
            file_region: name,
            file_domain: None,
        });
    }
    let side = faces[2];
    let same = model.face(side).sense == Sense::Same;
    for (index, domain) in charts.domains.into_iter().enumerate() {
        let uses = [
            (rings[index][0], true),
            (joins[1 - index], true),
            (rings[index][1], false),
            (joins[index], false),
        ];
        let name = format!("face_{}_patch_{index}", side.index());
        face_entities.push(advanced_face(&mut graph, &name, &uses, sides[index], same));
        face_map.push(StepFaceRegion {
            face: side,
            file_region: name,
            file_domain: Some([domain, [Rational::zero(), Rational::one()]]),
        });
    }
    let shell = graph.add("CLOSED_SHELL", vec![text("''"), Field::List(face_entities)]);
    let solid = graph.add("MANIFOLD_SOLID_BREP", vec![text("''"), shell.into()]);
    let uncertainty = crate::number::decimal(&options.file_tolerance.to_decimal_bounds(17).1);
    let mut file = graph.write(&[solid], |roots| {
        crate::export::header(&options.product_name, roots, &uncertainty)
    });
    file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    dec.check_deadline()?;
    out.write_all(file.as_bytes())
        .map_err(|error| Refusal::new("step.io-error", error.to_string()))?;
    Ok(Answer::exact(StepReport {
        max_deviation: rounded.bound,
        added_vertices: 4,
        added_edges: 4,
        solids: 1,
        faces: 4,
        edges: 6,
        vertices: 4,
        face_map,
        ..StepReport::default()
    }))
}
