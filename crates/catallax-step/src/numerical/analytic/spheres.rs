//! Regular spherical charts use longitude and latitude in radians, as in OCCT.

use super::*;
use catallax_shapes::Sphere;

pub(super) fn uv(
    sphere: &Sphere,
    p: &[Rational; 3],
    old: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let [x, y, z] = local(&sphere.center, &sphere.axis, &sphere.ref_dir, p)?;
    let radial = (x.square() + y.square()).sqrt().ok_or_else(|| {
        bad(
            "rounding-topology",
            "written sphere radial coordinate is unresolved",
        )
    })?;
    if radial.lo().sign() != Sign::Positive {
        return Err(bad(
            "rounding-topology",
            "written sphere chart reaches a latitude pole",
        ));
    }
    Ok([
        angle(&x, &y, &old[0], dec)?,
        angle(&radial, &z, &old[1], dec)?,
    ])
}

pub(super) fn rounded(
    w: &mut Writer<'_>,
    key: &str,
    sphere: &Sphere,
) -> Result<(Id, Surface), Refusal> {
    let WrittenFrame {
        entity: frame,
        origin,
        normal,
        reference,
    } = placement(w, &sphere.center, &sphere.axis, &sphere.ref_dir)?;
    let radius = sphere
        .radius
        .as_rational()
        .ok_or_else(|| bad("unsupported", "sphere export requires a rational radius"))?;
    let value = w.scalar(&radius)?;
    let entity = w.graph.add(
        "SPHERICAL_SURFACE",
        vec![text(string(key)), frame.into(), value],
    );
    let rounded = Sphere::new(
        point_from(&origin),
        written(&radius, w.options)?.into(),
        vector_from(&normal),
        vector_from(&reference),
        w.dec,
    )?;
    Ok((entity, Surface::Sphere(rounded)))
}
