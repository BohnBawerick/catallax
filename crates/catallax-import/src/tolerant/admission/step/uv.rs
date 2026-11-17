//! Rebuild sampled affine UV relations in the written 3D curve's parameter units.
//! Retained nonlinear pcurves are not replaced by an unmeasured approximation.
use super::super::{NativeTrim, SampledCheck, TolerantImport};
use super::*;
use crate::geometry_source::{Curve, Geometry};

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (0..3).fold(0.0, |distance, axis| distance.hypot(a[axis] - b[axis]))
}
fn number(value: f64) -> Result<Param, Error> {
    let value = Rational::from_f64_exact(value).ok_or_else(|| {
        invalid(
            "import.native-sampling-range",
            "nonfinite native pcurve coordinate",
        )
    })?;
    decimal(&value)
}

fn parameter_kind(
    builder: &Builder<'_>,
    mut entity: u64,
    dec: &mut Decider,
) -> Result<Option<u8>, Error> {
    for _ in 0..32 {
        match builder.parsed.geometry(entity, dec)? {
            Geometry::Curve(Curve::Line(..)) => return Ok(Some(0)),
            Geometry::Curve(Curve::Circle(..) | Curve::Ellipse(..)) => return Ok(Some(1)),
            Geometry::Curve(Curve::BSpline(_)) => return Ok(Some(2)),
            Geometry::Curve(
                Curve::SurfaceCurve(parent, ..)
                | Curve::Replica(parent, _)
                | Curve::Offset3d(parent, ..),
            ) => entity = parent,
            _ => return Ok(None),
        }
    }
    Err(invalid(
        "import.tolerant-broken-reference",
        "cyclic pcurve basis wrappers",
    ))
}

pub(super) struct Writer<'a> {
    import: &'a TolerantImport,
    context: Option<u64>,
    checks: Vec<SampledCheck>,
}

impl<'a> Writer<'a> {
    pub(super) fn new(import: &'a TolerantImport) -> Self {
        Self {
            import,
            context: None,
            checks: Vec::new(),
        }
    }
    pub(super) fn into_checks(self) -> Vec<SampledCheck> {
        self.checks
    }

    fn affine(
        &mut self,
        builder: &mut Builder<'_>,
        edge: &super::super::super::TolerantEdge,
        original: u64,
        surface_id: u64,
        trim: &NativeTrim,
        dec: &mut Decider,
    ) -> Result<Option<u64>, Error> {
        let import = self.import;
        let Some(curve) = import.native_curve(original) else {
            return Ok(None);
        };
        let Some(surface) = import.native_surface(surface_id) else {
            return Ok(None);
        };
        let Some(kind) = parameter_kind(builder, original, dec)? else {
            return Ok(None);
        };
        let first = trim.samples.first().expect("admitted native trim");
        let last = trim.samples.last().expect("admitted native trim");
        let t0 = first.curve.parameters[0];
        let t1 = last.curve.parameters[0];
        let width = t1 - t0;
        if width == 0.0 || !width.is_finite() {
            return Ok(None);
        }
        let start = curve.sample(first.curve, dec)?;
        let end = curve.sample(last.curve, dec)?;
        let scale = match kind {
            0 => distance(&start.point_mm, &end.point_mm) / width.abs(),
            1 => import.source_unit_factors()[1],
            _ => 1.0,
        };
        if scale <= 0.0 || !scale.is_finite() {
            return Ok(None);
        }
        let mut maximum: f64 = 0.0;
        for (index, sample) in trim.samples.iter().enumerate() {
            dec.charge_work(1)?;
            let fraction = (sample.curve.parameters[0] - t0) / width;
            let uv = std::array::from_fn(|axis| {
                first.surface.parameters[axis]
                    + fraction * (last.surface.parameters[axis] - first.surface.parameters[axis])
            });
            let Ok(location) = surface.location(uv) else {
                return Ok(None);
            };
            let actual = match surface.sample(location, dec) {
                Ok(sample) => sample,
                Err(error) if super::super::super::native::control(&error) => return Err(error),
                Err(_) => return Ok(None),
            };
            let source = if index == 0 {
                start
            } else if index + 1 == trim.samples.len() {
                end
            } else {
                curve.sample(sample.curve, dec)?
            };
            let residual = distance(&source.point_mm, &actual.point_mm)
                + source.uncertainty_mm
                + actual.uncertainty_mm;
            maximum = maximum.max(residual + 16.0 * f64::EPSILON * residual.max(1.0));
            let measured = Rational::from_f64_exact(maximum).ok_or_else(|| {
                invalid("import.native-sampling-range", "nonfinite pcurve residual")
            })?;
            if measured > *edge.tolerance_mm() {
                return Ok(None);
            }
        }
        let endpoints = if t0 < t1 {
            [first, last]
        } else {
            [last, first]
        };
        let mut controls = Vec::with_capacity(2);
        for sample in endpoints {
            let coordinates = sample
                .surface
                .parameters
                .into_iter()
                .map(number)
                .collect::<Result<_, _>>()?;
            let point = builder.add("CARTESIAN_POINT", vec![text(""), Param::List(coordinates)])?;
            controls.push(point);
        }
        let knots = if t0 < t1 {
            [t0 * scale, t1 * scale]
        } else {
            [t1 * scale, t0 * scale]
        };
        let spline = builder.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text(""),
                Param::Number("1".into()),
                references(controls),
                Param::Enum("UNSPECIFIED".into()),
                boolean(false),
                boolean(false),
                Param::List(vec![Param::Number("2".into()), Param::Number("2".into())]),
                Param::List(knots.into_iter().map(number).collect::<Result<_, _>>()?),
                Param::Enum("UNSPECIFIED".into()),
            ],
        )?;
        let context = match self.context {
            Some(context) => context,
            None => {
                let created = builder.complex(vec![
                    (
                        "GEOMETRIC_REPRESENTATION_CONTEXT".into(),
                        vec![Param::Number("2".into())],
                    ),
                    ("PARAMETRIC_REPRESENTATION_CONTEXT".into(), vec![]),
                    ("REPRESENTATION_CONTEXT".into(), vec![text(""), text("")]),
                ])?;
                self.context = Some(created);
                created
            }
        };
        let representation = builder.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text(""), references([spline]), reference_to(context)],
        )?;
        let pcurve = builder.add(
            "PCURVE",
            vec![
                text(""),
                reference_to(surface_id),
                reference_to(representation),
            ],
        )?;
        self.checks.push(SampledCheck {
            body: trim.body.clone(),
            scope: format!("{}:{}", trim.face, edge.name),
            check: "export.affine-native-uv-sampled-correspondence",
            samples: trim.samples.len(),
            maximum_residual_mm: maximum,
            converged: trim.samples.len(),
        });
        Ok(Some(pcurve))
    }

    pub(super) fn edge_geometry(
        &mut self,
        builder: &mut Builder<'_>,
        edge: &super::super::super::TolerantEdge,
        original: u64,
        trims: &[&NativeTrim],
        dec: &mut Decider,
    ) -> Result<u64, Error> {
        // STEP's associated_geometry has at most two entries. Do not invent a
        // per-face mapping for a compound's higher-valence shared source edge.
        if trims.len() > 2 {
            return Ok(original);
        }
        let mut generated = false;
        let mut associated = Vec::with_capacity(trims.len());
        let mut surfaces = Vec::with_capacity(trims.len());
        let mut all_pcurves = true;
        for trim in trims {
            let Some(surface) = builder.source_surface(&trim.face, dec)? else {
                continue;
            };
            surfaces.push(surface);
            let curve = self.affine(builder, edge, original, surface, trim, dec)?;
            associated.push(match curve {
                Some(curve) => {
                    generated = true;
                    curve
                }
                None => match trim.source_pcurve {
                    Some(curve) => curve,
                    None => {
                        all_pcurves = false;
                        surface
                    }
                },
            });
        }
        if !generated {
            return Ok(original);
        }
        let mut basis = original;
        for depth in 0..32 {
            match builder.parsed.geometry(basis, dec)? {
                Geometry::Curve(Curve::SurfaceCurve(parent, ..)) => basis = parent,
                _ => break,
            }
            if depth == 31 {
                return Err(invalid(
                    "import.tolerant-broken-reference",
                    "cyclic surface-curve basis",
                ));
            }
        }
        let seam = all_pcurves && surfaces.len() == 2 && surfaces[0] == surfaces[1];
        builder.add(
            if seam { "SEAM_CURVE" } else { "SURFACE_CURVE" },
            vec![
                text("native sampled UV"),
                reference_to(basis),
                references(associated),
                Param::Enum("CURVE_3D".into()),
            ],
        )
    }
}
