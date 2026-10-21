use super::*;

// Rational paths on a carrier. A plane path is linear; cylinder and sphere
// paths are quaternion rotations, so the implicit stays exactly zero.
fn carrier_candidate(surface: &Surface, p: &RPoint3, t: &Rational, alternate: bool) -> RPoint3 {
    match surface {
        Surface::Plane(s) => {
            let direction = if alternate {
                s.normal.cross(&s.ref_dir)
            } else {
                s.ref_dir.clone()
            };
            p.plus(&direction.scaled(t))
        }
        Surface::Cylinder(s) => {
            let d = p.minus(&s.origin);
            let rotation =
                Rotation::from_quaternion(r(1), &s.axis.x * t, &s.axis.y * t, &s.axis.z * t)
                    .unwrap();
            s.origin.plus(&rotation.apply(&d))
        }
        Surface::Sphere(s) => {
            let d = p.minus(&s.center);
            let axis = if d.cross(&s.axis).is_zero() {
                &s.ref_dir
            } else {
                &s.axis
            };
            let rotation =
                Rotation::from_quaternion(r(1), &axis.x * t, &axis.y * t, &axis.z * t).unwrap();
            s.center.plus(&rotation.apply(&d))
        }
        _ => panic!("unsupported contact carrier"),
    }
}

pub fn touches_another_carrier(m: &Model, f: FaceId, faces: &[FaceId], p: &RPoint3) -> bool {
    let surface = m.face(f).surface;
    let (_, n) = implicit_normal(m.geom().surface(surface), p);
    let mut dec = Decider::default();
    faces.iter().any(|&other| {
        let s = m.face(other).surface;
        if matches!(
            m.geom().same_surface(surface, s, &mut dec),
            Outcome::ProvenCoincident(_)
        ) {
            return false;
        }
        let (value, other_n) = implicit_normal(m.geom().surface(s), p);
        value.is_zero() && n.cross(&other_n).is_zero()
    })
}

pub fn contact_witness(
    m: &Model,
    f: FaceId,
    source: &Primitive,
    other: &Primitive,
    p: &RPoint3,
) -> RPoint3 {
    contact_witness_where(m, f, source, other, p, |_| true)
}

fn contact_witness_where(
    m: &Model,
    f: FaceId,
    source: &Primitive,
    other: &Primitive,
    p: &RPoint3,
    accept: impl Fn(&RPoint3) -> bool,
) -> RPoint3 {
    let surface = m.geom().surface(m.face(f).surface);
    assert_eq!(source.locate(p), Location::On);
    assert_eq!(
        source.constraints(p).iter().filter(|x| x.is_zero()).count(),
        1,
        "a contact fallback cannot excuse a face_point on a primitive edge"
    );
    // Only contact-only faces call this helper. They have the entire source
    // primitive face as their relative interior. A candidate on that carrier,
    // strictly away from all its other constraints, is therefore in the piece.
    // The other primitive's boundary is excluded, as the architect ruled.
    let mut t = r(1);
    loop {
        for alternate in [false, true] {
            let q = carrier_candidate(surface, p, &t, alternate);
            if source.locate(&q) != Location::On || other.locate(&q) == Location::On || !accept(&q)
            {
                continue;
            }
            if source
                .constraints(&q)
                .iter()
                .filter(|x| x.is_zero())
                .count()
                != 1
            {
                continue;
            }
            assert_eq!(implicit_normal(surface, &q).0, r(0));
            return q;
        }
        // The input is strictly inside this primitive face. Its open relative
        // neighbourhood contains a non-contact point along one of these paths.
        // No arbitrary minimum distance or search cutoff is imposed.
        t = t / r(2);
    }
}

pub fn orientation_witness(
    m: &Model,
    b: BodyId,
    f: FaceId,
    operands: [(BodyId, &Primitive); 2],
    p: RPoint3,
) -> RPoint3 {
    let faces = m.faces_of_body(b);
    if !touches_another_carrier(m, f, &faces, &p) {
        return p;
    }
    let surface = m.face(f).surface;
    let mut dec = Decider::default();
    for (index, (body, primitive)) in operands.into_iter().enumerate() {
        let opposite = operands[1 - index].1;
        if m.faces_of_body(body).iter().any(|&g| {
            matches!(
                m.geom().same_surface(surface, m.face(g).surface, &mut dec),
                Outcome::ProvenCoincident(_)
            )
        }) {
            let q = contact_witness_where(m, f, primitive, opposite, &p, |q| {
                !touches_another_carrier(m, f, &faces, q)
            });
            assert!(!touches_another_carrier(m, f, &faces, &q));
            return q;
        }
    }
    panic!("no source face for contact witness");
}
