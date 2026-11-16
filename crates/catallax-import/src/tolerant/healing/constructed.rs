//! Native face definitions enter the existing import sewer without STEP I/O.

use super::*;
use crate::tolerant::ToleranceOptions;
use catallax_numbers::Sign;

pub(crate) fn sew(
    name: &str,
    faces: &[crate::Face],
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<HealedFile, Error> {
    if name.is_empty() || faces.is_empty() {
        return Err(failure(
            "surface.sew-input",
            "a constructed surface set needs a name and faces",
        ));
    }
    if options.initial_mm.sign() == Sign::Negative
        || options.maximum_mm.sign() != Sign::Positive
        || options.initial_mm > options.maximum_mm
    {
        return Err(failure(
            "surface.sew-tolerance",
            "initial tolerance must lie between zero and the positive sewing cap",
        ));
    }
    dec.check_deadline()?;
    dec.charge_work(name.len() as u64 + faces.len() as u64)?;
    for face in faces {
        dec.charge_work(face.name.len() as u64 + face.loops.len() as u64)?;
        if matches!(face.surface, crate::Surface::Retained(_)) {
            return Err(failure(
                "surface.sew-carrier",
                "constructed faces need native carriers, not opaque imported entities",
            ));
        }
        charge_surface(&face.surface, dec)?;
        for boundary in &face.loops {
            dec.charge_work(1 + boundary.edges.len() as u64)?;
            if let Some(vertex) = &boundary.vertex {
                dec.charge_work(vertex.name.len() as u64)?;
                charge_vector(&vertex.point, dec)?;
            }
            for edge in &boundary.edges {
                dec.charge_work(
                    edge.edge.len() as u64
                        + edge.start.name.len() as u64
                        + edge.end.name.len() as u64,
                )?;
                if matches!(edge.curve, crate::Curve::Retained(_)) {
                    return Err(failure(
                        "surface.sew-carrier",
                        "constructed edges need native carriers, not opaque imported entities",
                    ));
                }
                charge_curve(&edge.curve, dec)?;
                charge_vector(&edge.start.point, dec)?;
                charge_vector(&edge.end.point, dec)?;
            }
        }
    }
    // A diagnostic identity for these immutable in-memory definitions, not a
    // STEP digest or a persistence codec accepted by tolerant::reopen.
    let fingerprint = {
        let encoding = format!("catallax.constructed-surfaces/1\n{name:?}\n{faces:?}");
        dec.charge_work(encoding.len() as u64)?;
        crate::sha256_hex(encoding.as_bytes())
    };
    let mut body = TolerantBody::new(name, false, &options);
    body.shell(name, false, faces, dec)?;
    let source = TolerantFile {
        source_sha256: fingerprint,
        mm_per_unit: Rational::one(),
        options,
        bodies: vec![body],
    };
    heal(&source, dec)
}
