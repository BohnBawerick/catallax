//! Queries on the unique transverse cap cut of a spline fillet.

use crate::*;
use catallax_shapes::{Curve, FilletCutTrace};

pub(crate) fn of_curve(
    geom: &GeomStore,
    curve: CurveId,
    dec: &mut Decider,
) -> Result<Option<FilletCutTrace>, Refusal> {
    let Curve::Intersection(branch) = geom.curve(curve) else {
        return Ok(None);
    };
    let surfaces = branch.surfaces.map(|id| geom.surface(id));
    let (fillet, plane) = match surfaces {
        [
            Surface::BSpline(SplineSupport::SplineFillet(fillet)),
            Surface::Plane(plane),
        ]
        | [
            Surface::Plane(plane),
            Surface::BSpline(SplineSupport::SplineFillet(fillet)),
        ] => (fillet, plane),
        _ => return Ok(None),
    };
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(|v| v.as_rational());
    let [Some(x), Some(y), Some(z)] = normal else {
        return Err(Refusal::new(
            "intersect.fillet-cut-plane",
            "a cap cut requires a rational plane",
        ));
    };
    let origin = Vector3::new(
        plane.origin.x.clone(),
        plane.origin.y.clone(),
        plane.origin.z.clone(),
    );
    let offset = plane.normal.dot(&origin);
    let offset = offset.as_rational().ok_or_else(|| {
        Refusal::new(
            "intersect.fillet-cut-plane",
            "a cap cut requires a rational plane offset",
        )
    })?;
    let trace = fillet.cut_trace([x, y, z], offset, dec)?.ok_or_else(|| {
        Refusal::new(
            "intersect.fillet-cut-plane",
            "the plane misses the finite fillet carrier",
        )
    })?;
    // The unique graph supplies the component; its witness must still belong
    // to both supports, rather than merely name them.
    let witness = geom.position(branch.witness);
    if fillet.side_at(&witness, dec)? != catallax_numbers::Sign::Zero
        || math::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?
            != catallax_numbers::Sign::Zero
    {
        return Err(Refusal::new(
            "intersect.fillet-cut-witness",
            "cap-cut witness is not on both supports",
        ));
    }
    Ok(Some(trace))
}

pub(crate) fn trace(
    geom: &GeomStore,
    curve: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Result<Option<CurveTrace>, Refusal> {
    let Some(cut) = of_curve(geom, curve, dec)? else {
        return Ok(None);
    };
    let (Some(from), Some(to)) = (from, to) else {
        return Err(Refusal::new(
            "intersect.fillet-cut-arc",
            "a cap cut needs both rail endpoints",
        ));
    };
    let same = |p: PointId, side: usize, dec: &mut Decider| match dec
        .sign(&geom.position(p).minus(&cut.endpoint(side)).norm_squared())
    {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(proof)
            if proof.proof != catallax_numbers::ZeroProof::SeparationBound =>
        {
            Ok(true)
        }
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "intersect.fillet-cut-arc",
            "cap-cut endpoints need exact identity",
        )),
        Outcome::Refused(error) => Err(error),
        Outcome::NotExact(error) => Err(error.into()),
    };
    let forward = if same(from, 0, dec)? && same(to, 1, dec)? {
        true
    } else if same(from, 1, dec)? && same(to, 0, dec)? {
        false
    } else {
        return Err(Refusal::new(
            "intersect.fillet-cut-arc",
            "a cap-cut arc must join its two exact rail endpoints",
        ));
    };
    let mut piece = TracePiece::new(std::array::from_fn(|_| {
        chart::E::constant(catallax_numbers::Real::zero())
    }));
    piece.fillet_cut = Some((Arc::new(cut), forward));
    Ok(Some(CurveTrace {
        pieces: vec![piece],
    }))
}
