//! Finite NURBS definitions as an AP242 geometric surface set, not sewn bodies.

use std::io::Write;
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::spline_tools::RationalSplineSurface;
use crate::graph::{Graph, Field, text, refs, boolean, string};

/// Decimal-file proof for a set of exact spline definitions.
#[derive(Clone, Debug)]
pub struct SplineSurfaceStepReport {
    /// Number of input surface definitions.
    pub surfaces: usize,
    /// Number of emitted Bezier surfaces, one per nonempty tensor knot span.
    pub patches: usize,
    /// Requested maximum distance in millimetres, fixed before rounding.
    pub tolerance_mm: Rational,
    /// Whole-surface L1 rounding bound, hence also a Euclidean distance bound.
    /// Every patch uses [0,1] knots exactly, so no rounded parameter knots enter it.
    pub max_deviation_mm: Rational,
    /// True only when every emitted control and weight equals its source exactly.
    pub exact: bool,
}

fn bad(detail: &str) -> Refusal { Refusal::new("step.spline-surfaces", detail) }

/// Writes an AP242 geometric surface set. It does not claim surface regularity,
/// sewing or solid topology. Each source knot rectangle becomes a Bezier surface
/// with exact [0,1] knots, avoiding an unbounded rounded-knot displacement.
/// Decimal controls and weights have a whole-patch positive-denominator bound.
/// All validation and rounding finish before any bytes reach `out`.
pub fn write_spline_surfaces_step<W: Write>(surfaces: &[RationalSplineSurface],
    tolerance_mm: &Rational, significant_digits: u32, dec: &mut Decider, out: &mut W)
    -> Result<SplineSurfaceStepReport, Refusal> {
    if surfaces.is_empty() || tolerance_mm.sign() != Sign::Positive || !(2..=1000).contains(&significant_digits) {
        return Err(bad("surface set must be nonempty, tolerance positive and decimal precision between 2 and 1000"));
    }
    let mut graph = Graph::default();
    let mut ids = Vec::new();
    let mut maximum = Rational::zero();
    let mut exact = true;
    for surface in surfaces {
        let surface = surface.clone().into_clamped(dec)?;
        let mut breaks = surface.knots().map(<[Rational]>::to_vec);
        for knots in &mut breaks { knots.dedup(); }
        for u in breaks[0].windows(2) { for v in breaks[1].windows(2) {
            dec.check_deadline()?;
            let patch = surface.extracted([[&u[0], &u[1]], [&v[0], &v[1]]], dec)?;
            let columns = patch.control_count()[1];
            let mut points = Vec::with_capacity(patch.controls().len());
            let mut weights = Vec::with_capacity(patch.weights().len());
            let origin = &patch.controls()[0];
            let mut radius = Rational::zero();
            let mut point_error = Rational::zero();
            let mut weight_error = Rational::zero();
            let mut min_weight: Option<Rational> = None;
            for (point, weight) in patch.controls().iter().zip(patch.weights()) {
                dec.charge_work(1)?;
                let numbers = point.each_ref().map(|q| crate::number::rational(q, significant_digits));
                let error = numbers.iter().fold(Rational::zero(), |s,n| &s + &n.error);
                point_error = point_error.max(error);
                radius = radius.max(point.iter().zip(origin).fold(Rational::zero(), |s,(q,o)| s + (q-o).abs()));
                let id = graph.add("CARTESIAN_POINT", vec![text("''"), Field::List(numbers.map(|n| text(n.text)).to_vec())]);
                points.push(Field::from(id));
                let number = crate::number::rational(weight, significant_digits);
                if number.value.sign() != Sign::Positive { return Err(bad("rounded spline weight is nonpositive")); }
                weight_error = weight_error.max(number.error);
                min_weight = Some(min_weight.map_or(number.value.clone(), |w| w.min(number.value)));
                weights.push(text(number.text));
            }
            exact &= point_error.is_zero() && weight_error.is_zero();
            let bound = point_error + Rational::from_i64(2) * radius * weight_error / min_weight.unwrap();
            if bound > *tolerance_mm { return Err(bad("decimal spline control bound exceeds requested tolerance")); }
            maximum = maximum.max(bound);
            let grid = |fields: Vec<Field>| Field::List(fields.chunks(columns).map(|row| Field::List(row.to_vec())).collect());
            let degrees = patch.degrees();
            let multiplicities = degrees.map(|d| Field::List(vec![text((d+1).to_string()), text((d+1).to_string())]));
            let knots = || Field::List(vec![text("0."), text("1.")]);
            ids.push(graph.add_complex(vec![
                ("BOUNDED_SURFACE", vec![]),
                ("B_SPLINE_SURFACE", vec![text(degrees[0].to_string()), text(degrees[1].to_string()), grid(points), text(".UNSPECIFIED."), boolean(false), boolean(false), boolean(false)]),
                ("B_SPLINE_SURFACE_WITH_KNOTS", vec![multiplicities[0].clone(),multiplicities[1].clone(),knots(),knots(),text(".UNSPECIFIED.")]),
                ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
                ("RATIONAL_B_SPLINE_SURFACE", vec![grid(weights)]),
                ("REPRESENTATION_ITEM", vec![text("''")]),
                ("SURFACE", vec![]),
            ]));
        }}
    }
    let name = if exact { "catallax exact spline definitions" } else { "catallax approximate spline definitions" };
    let set = graph.add("GEOMETRIC_SET", vec![text(string(name)), refs(ids.iter().copied())]);
    let mut file = graph.write(&[set], |roots| format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('{}'),'2;1');\nFILE_NAME('surfaces.step','2026-11-17T00:00:00',(''),(''),'catallax','catallax','');\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF'));\nENDSEC;\nDATA;\n#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1);\n#3=PRODUCT_CONTEXT('',#1,'mechanical');\n#4=PRODUCT('surfaces','surfaces','',(#3));\n#5=PRODUCT_DEFINITION_FORMATION('','',#4);\n#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n#7=PRODUCT_DEFINITION('design','',#5,#6);\n#8=PRODUCT_DEFINITION_SHAPE('','',#7);\n#9=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#10=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n#11=(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());\n#12=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({}),#9,'distance_accuracy_value','whole surface rounding bound');\n#13=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#12)) GLOBAL_UNIT_ASSIGNED_CONTEXT((#9,#10,#11)) REPRESENTATION_CONTEXT('',''));\n#14=GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION('{}',(#{}),#13);\n#15=SHAPE_DEFINITION_REPRESENTATION(#8,#14);\n",
        name, crate::number::rational(tolerance_mm, significant_digits).text, name, roots[0]));
    file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    out.write_all(file.as_bytes()).map_err(|e| Refusal::new("step.io-error",e.to_string()))?;
    Ok(SplineSurfaceStepReport { surfaces: surfaces.len(), patches: ids.len(), tolerance_mm: tolerance_mm.clone(), max_deviation_mm: maximum, exact })
}
