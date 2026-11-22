//! Ordered exact spline wires, independent of face loops and material shells.

use crate::{Body, BodyId, BodyKind, Coedge, Edge, Model, Sense, Vertex, VertexId, Wire};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, MarkLog, Rational, Refusal, operation_budget_active,
    with_operation_budget,
};
use catallax_shapes::{Curve, GeomKey, PointDef, RPoint3, spline_tools::RationalSplineCurve};

fn point(model: &mut Model, value: &[Rational; 3]) -> VertexId {
    let point = model.geom_mut().add_point(PointDef::Explicit(
        RPoint3::new(value[0].clone(), value[1].clone(), value[2].clone()).to_real(),
    ));
    model.add_vertex(Vertex {
        point,
        tolerance: None,
    })
}

fn comparison_work(value: &[Rational; 3], dec: &mut Decider) -> Result<(), Refusal> {
    for number in value {
        let bits = number
            .numerator()
            .bit_length()
            .saturating_add(number.denominator().bit_length());
        dec.charge_work(
            bits.saturating_add(1)
                .saturating_mul(bits.saturating_add(1))
                .saturating_mul(4),
        )?;
    }
    Ok(())
}

/// Construct one open or closed directed path. Consecutive endpoints must agree
/// exactly; no endpoint is moved or snapped. Closure shares the first vertex.
/// The fresh model is published only after every piece and incidence is checked.
/// Crossings are permitted: a wire asserts no material or face regularity.
pub fn spline_wire(
    pieces: &[RationalSplineCurve],
    dec: &mut Decider,
) -> Result<(Model, BodyId), Refusal> {
    if !operation_budget_active() {
        return with_operation_budget(crate::numerical::DEFAULT_WORK, || spline_wire(pieces, dec))
            .0;
    }
    dec.check_deadline()?;
    if pieces.is_empty() {
        return Err(Refusal::new(
            "topology.empty-wire",
            "a wire needs at least one curve",
        ));
    }
    let mut model = Model::new();
    dec.charge_work((pieces.len() as u64).saturating_mul(16))?;
    let mut coedges = Vec::with_capacity(pieces.len());
    let mut first = None;
    let mut previous = None;
    for piece in pieces {
        let domain = piece.domain();
        let clamped;
        let piece = if piece.knots().first() == Some(domain[0])
            && piece.knots().last() == Some(domain[1])
        {
            piece
        } else {
            clamped = piece.extracted(domain, dec)?;
            &clamped
        };
        let a = piece.point_at(domain[0], dec)?;
        let b = piece.point_at(domain[1], dec)?;
        comparison_work(&a, dec)?;
        comparison_work(&b, dec)?;
        let start = if let Some((value, vertex)) = &previous {
            if *value != a {
                return Err(Refusal::new(
                    "topology.disconnected-wire",
                    "consecutive curve endpoints differ",
                ));
            }
            *vertex
        } else {
            let vertex = point(&mut model, &a);
            first = Some((a, vertex));
            vertex
        };
        let (first_value, first_vertex) =
            first.as_ref().expect("first piece establishes the start");
        let end = if &b == first_value {
            *first_vertex
        } else {
            point(&mut model, &b)
        };
        let curve = piece.to_curve(dec)?;
        let curve = model.geom_mut().add_curve(Curve::BSpline(curve));
        for vertex in [start, end] {
            let point = model.vertex(vertex).point;
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Point(point),
                b: GeomKey::Curve(curve),
                at: None,
                origin: FactOrigin::Construction("topology.spline-wire-endpoint"),
            });
        }
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
            tolerance: None,
        });
        coedges.push(model.add_coedge(Coedge {
            edge,
            sense: Sense::Same,
        }));
        previous = Some((b, end));
    }
    let body = model.add_body(Body {
        kind: BodyKind::Wire,
        lumps: Vec::new(),
        wires: vec![Wire { coedges }],
        marks: MarkLog::new(),
    });
    if let Some(fault) = model.check(body).into_iter().next() {
        return Err(Refusal::new(fault.code, fault.detail));
    }
    dec.check_deadline()?;
    Ok((model, body))
}
