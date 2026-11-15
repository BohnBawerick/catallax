//! Corner occurrences in one face must not hide disconnected vertex fans.
mod common;
use common::*;

fn sheet_with_caps(pinched: bool) -> (Model, BodyId) {
    let mut m = Model::new();
    let points = [
        p(0, 0, 0),
        p(20, 0, 0),
        p(20, 20, 0),
        p(0, 20, 0),
        p(2, 2, 0),
        p(3, 2, 0),
        p(2, 3, 0),
        p(5, 2, 0),
        p(6, 2, 0),
        p(5, 3, 0),
        p(10, 10, 0),
        p(11, 10, 0),
        p(10, 11, 0),
        p(13, 10, 0),
        p(14, 10, 0),
        p(13, 11, 0),
    ];
    let vertices = points
        .iter()
        .map(|p| vertex(&mut m, p.clone()))
        .collect::<Vec<_>>();
    let mut polygons = [
        vec![0, 1, 2, 3],
        vec![4, 5, 6],
        vec![7, 8, 9],
        vec![10, 11, 12],
        vec![13, 14, 15],
    ];
    if pinched {
        polygons[2][0] = 4;
        polygons[4][0] = 10;
    }
    let mut loops = Vec::new();
    let mut caps = Vec::new();
    for (index, polygon) in polygons.iter().enumerate() {
        let edges = polygon
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let b = polygon[(i + 1) % polygon.len()];
                let curve = Curve::Line(Line {
                    origin: points[*a].clone().into(),
                    dir: RVector3::new(
                        &points[b].x - &points[*a].x,
                        &points[b].y - &points[*a].y,
                        &points[b].z - &points[*a].z,
                    )
                    .into(),
                });
                edge(&mut m, curve, Some(vertices[*a]), Some(vertices[b]))
            })
            .collect::<Vec<_>>();
        loops.push(cycle(
            &mut m,
            &edges.iter().map(|e| (*e, Sense::Same)).collect::<Vec<_>>(),
        ));
        if index > 0 {
            let boundary = cycle(
                &mut m,
                &edges
                    .iter()
                    .rev()
                    .map(|e| (*e, Sense::Reversed))
                    .collect::<Vec<_>>(),
            );
            caps.push(face(
                &mut m,
                plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)),
                Sense::Same,
                vec![boundary],
            ));
        }
    }
    caps.push(face(
        &mut m,
        plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)),
        Sense::Same,
        loops,
    ));
    let shell = m.add_shell(Shell { faces: caps });
    let lump = m.add_lump(Lump {
        shells: vec![shell],
    });
    let body = m.add_body(Body {
        kind: BodyKind::Surface,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    (m, body)
}

#[test]
fn separate_corner_cycles_in_the_same_face_refuse() {
    let (model, body) = sheet_with_caps(false);
    assert!(model.check(body).is_empty());
    let (model, body) = sheet_with_caps(true);
    let faults = model.check(body);
    assert!(
        faults
            .iter()
            .any(|fault| fault.code == "topology.vertex-fan"),
        "{faults:?}"
    );
}
