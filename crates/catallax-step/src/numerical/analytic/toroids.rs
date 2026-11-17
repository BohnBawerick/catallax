use super::*;
use catallax_shapes::Torus;

pub(super) fn uv(
    torus: &Torus,
    p: &[Rational; 3],
    old: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let [x, y, z] = local(&torus.center, &torus.axis, &torus.ref_dir, p)?;
    let radial = (x.square() + y.square())
        .sqrt()
        .ok_or_else(|| bad("unsupported", "torus radial coordinate is unresolved"))?;
    let major = dec.enclose_budgeted_at(&torus.major, BITS)?;
    Ok([
        angle(&x, &y, &old[0], dec)?,
        angle(&(radial - major), &z, &old[1], dec)?,
    ])
}

pub(super) fn rounded(
    w: &mut Writer<'_>,
    key: &str,
    torus: &Torus,
) -> Result<(Id, Surface), Refusal> {
    let WrittenFrame {
        entity: frame,
        origin,
        normal,
        reference,
    } = placement(w, &torus.center, &torus.axis, &torus.ref_dir)?;
    let major = torus.major.as_rational().ok_or_else(|| {
        bad(
            "unsupported",
            "torus export requires a rational major radius",
        )
    })?;
    let minor = torus.minor.as_rational().ok_or_else(|| {
        bad(
            "unsupported",
            "torus export requires a rational minor radius",
        )
    })?;
    let major_value = w.scalar(&major)?;
    let minor_value = w.scalar(&minor)?;
    let entity = w.graph.add(
        "TOROIDAL_SURFACE",
        vec![text(string(key)), frame.into(), major_value, minor_value],
    );
    let rounded = Torus::new(
        point_from(&origin),
        vector_from(&normal),
        written(&major, w.options)?.into(),
        written(&minor, w.options)?.into(),
        vector_from(&reference),
        w.dec,
    )?;
    Ok((entity, Surface::Torus(rounded)))
}
