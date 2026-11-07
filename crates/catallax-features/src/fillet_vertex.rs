//! Exact spherical polygons at fully selected planar vertices.

use super::*;
use catallax_shapes::Sphere;

/// All contacts must share one ball, and each selected pair must be an exposed
/// edge of their positive cone. These predicates prove the whole spherical
/// polygon, not just its trim vertices. Partial selections need end closures.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    old: VertexId,
    active: &[usize],
    runs: &mut [Run],
    incident: &BTreeSet<FaceId>,
    uses: &Uses,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Origin, Surface, bool), Refusal> {
    if active.len() != incident.len()
        || uses.keys().any(|edge| {
            let e = model.edge(*edge);
            (e.start == Some(old) || e.end == Some(old))
                && active.iter().all(|i| runs[*i].edge != *edge)
        })
    {
        return Err(bad(
            "a spherical polygon needs every incident edge selected",
        ));
    }
    let at = model.geom().position(model.vertex(old).point);
    let faces: Vec<_> = incident.iter().copied().collect();
    let convex = runs[active[0]].definition.convex;
    let radius = runs[active[0]].definition.radius.clone();
    let mut directions = Vec::with_capacity(faces.len());
    for face in &faces {
        let p = plane(model, *face)?;
        if sign(&at.minus(&p.origin).dot(&p.normal), dec)? != Sign::Zero {
            return Err(bad("a corner vertex is not on every incident plane"));
        }
        let n = unit(&outward(model, *face, &at, dec)?);
        directions.push(if convex { n } else { n.scaled(&-Real::one()) });
    }
    let first = &runs[active[0]].definition;
    let Surface::Cylinder(axis) = &first.surface else {
        return Err(bad("a planar corner strip must reduce to a cylinder"));
    };
    let mut centre = None;
    for (face, normal) in faces.iter().zip(&directions) {
        let divisor = axis.axis.dot(normal);
        if sign(&divisor, dec)? != Sign::Zero {
            let p = plane(model, *face)?;
            let depth = -(&radius + &axis.origin.minus(&p.origin).dot(normal));
            centre = Some(axis.origin.plus(&axis.axis.scaled(&(depth / divisor))));
            break;
        }
    }
    let centre = centre.ok_or_else(|| bad("corner supports do not fix one ball centre"))?;
    let sum = directions
        .iter()
        .fold(Vector3::zero(), |sum, n| sum.plus(n));
    for (face, n) in faces.iter().zip(&directions) {
        let p = plane(model, *face)?;
        if sign(&(&centre.minus(&p.origin).dot(n) + &radius), dec)? != Sign::Zero {
            return Err(bad("the corner ball does not contact every incident plane"));
        }
        if sign(&sum.dot(n), dec)? != Sign::Positive {
            return Err(bad("corner contact normals need one open hemisphere"));
        }
    }
    let mut pairs = Vec::with_capacity(active.len());
    let mut degree = vec![0; faces.len()];
    for i in active {
        let run = &runs[*i];
        if run.definition.convex != convex
            || sign(&(&run.definition.radius - &radius), dec)? != Sign::Zero
        {
            return Err(bad(
                "a spherical polygon needs one radius and material side",
            ));
        }
        let Surface::Cylinder(cylinder) = &run.definition.surface else {
            return Err(bad("a planar corner strip must reduce to a cylinder"));
        };
        if sign(
            &centre
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
        {
            return Err(bad("corner strip axes do not meet at one ball centre"));
        }
        let Curve::Line(edge) = &run.definition.edge else {
            return Err(bad("a planar corner needs straight input edges"));
        };
        let input = model.edge(run.edge);
        let far = if input.start == Some(old) {
            input.end
        } else {
            input.start
        }
        .ok_or_else(|| bad("a corner strip needs a finite other endpoint"))?;
        let far = model.geom().position(model.vertex(far).point);
        let length = far.minus(&at).dot(&edge.dir);
        let station = centre.minus(&at).dot(&edge.dir);
        let remaining = &length - &station;
        if sign(&(&station * &length), dec)? != Sign::Positive
            || sign(&(&remaining * &length), dec)? != Sign::Positive
        {
            return Err(bad("the corner ball centre lies beyond an input edge end"));
        }
        let pair = run
            .faces
            .map(|f| faces.iter().position(|face| *face == f).unwrap());
        let [a, b] = pair;
        degree[a] += 1;
        degree[b] += 1;
        let boundary = directions[a].cross(&directions[b]);
        let mut side = None;
        for (k, n) in directions.iter().enumerate() {
            if k == a || k == b {
                continue;
            }
            let sign = sign(&boundary.dot(n), dec)?;
            if sign == Sign::Zero || side.is_some_and(|previous| previous != sign) {
                return Err(bad(
                    "a selected pair is not an exposed edge of the contact cone",
                ));
            }
            side = Some(sign);
        }
        pairs.push(pair);
    }
    if degree.iter().any(|degree| *degree != 2) {
        return Err(bad(
            "corner contacts do not form a closed spherical polygon",
        ));
    }
    let reference = directions[1].minus(&directions[0].scaled(&directions[1].dot(&directions[0])));
    let surface = Surface::Sphere(Sphere::new(
        centre.clone(),
        radius.clone(),
        directions[0].clone(),
        reference,
        dec,
    )?);
    for face in model.faces_of_body(body) {
        if incident.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if !separated(&surface, &other, dec)? {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-incident face is not proved clear of the corner sphere",
            ));
        }
    }
    let feet: Vec<_> = directions
        .iter()
        .map(|n| vertex(model, centre.plus(&n.scaled(&radius))))
        .collect();
    for (i, [a, b]) in active.iter().zip(pairs) {
        let run = &mut runs[*i];
        run.ends.insert((old, run.faces[0]), feet[a]);
        run.ends.insert((old, run.faces[1]), feet[b]);
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: centre.clone(),
            normal: directions[a].cross(&directions[b]),
            radius_sq: radius.square(),
            ref_dir: directions[a].clone(),
        }));
        let cap = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(feet[a]),
            end: Some(feet[b]),
        });
        run.caps.insert(old, vec![cap]);
    }
    let name = model
        .vertex_name(body, old)
        .ok_or_else(|| bad("a corner vertex needs an unambiguous origin name"))?;
    Ok((
        Origin {
            step: step.clone(),
            role: Role::Corner { vertex: name },
        },
        surface,
        convex,
    ))
}
