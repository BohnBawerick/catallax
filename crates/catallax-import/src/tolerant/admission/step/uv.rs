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
    bases: BTreeMap<u64, u64>,
}

impl<'a> Writer<'a> {
    pub(super) fn new(import: &'a TolerantImport) -> Self {
        Self {
            import,
            context: None,
            checks: Vec::new(),
            bases: BTreeMap::new(),
        }
    }
    pub(super) fn into_checks(self) -> Vec<SampledCheck> {
        self.checks
    }

    fn spline(
        &mut self,
        builder: &mut Builder<'_>,
        surface_id: u64,
        points: impl ExactSizeIterator<Item = ([f64; 2], f64)>,
    ) -> Result<u64, Error> {
        let count = points.len();
        let mut controls = Vec::with_capacity(count);
        let mut knots = Vec::with_capacity(count);
        let mut multiplicities = Vec::with_capacity(count);
        let mut previous = None;
        for (index, (uv, parameter)) in points.enumerate() {
            if !parameter.is_finite() || previous.is_some_and(|value| parameter <= value) {
                return Err(invalid(
                    "import.native-trace-parameter-order",
                    "reconstructed trace knots must increase",
                ));
            }
            previous = Some(parameter);
            let coordinates = uv.into_iter().map(number).collect::<Result<Vec<_>, _>>()?;
            controls
                .push(builder.add("CARTESIAN_POINT", vec![text(""), Param::List(coordinates)])?);
            knots.push(number(parameter)?);
            multiplicities.push(Param::Number(
                if index == 0 || index + 1 == count {
                    "2"
                } else {
                    "1"
                }
                .into(),
            ));
        }
        let spline = builder.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text(""),
                Param::Number("1".into()),
                references(controls),
                Param::Enum("UNSPECIFIED".into()),
                Param::Enum("U".into()),
                Param::Enum("U".into()),
                Param::List(multiplicities),
                Param::List(knots),
                Param::Enum("UNSPECIFIED".into()),
            ],
        )?;
        self.wrap(builder, surface_id, spline)
    }

    fn wrap(
        &mut self,
        builder: &mut Builder<'_>,
        surface_id: u64,
        spline: u64,
    ) -> Result<u64, Error> {
        let context = match self.context {
            Some(context) => context,
            None => {
                let context = builder.complex(vec![
                    (
                        "GEOMETRIC_REPRESENTATION_CONTEXT".into(),
                        vec![Param::Number("2".into())],
                    ),
                    ("PARAMETRIC_REPRESENTATION_CONTEXT".into(), vec![]),
                    ("REPRESENTATION_CONTEXT".into(), vec![text(""), text("")]),
                ])?;
                self.context = Some(context);
                context
            }
        };
        let representation = builder.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text(""), references([spline]), reference_to(context)],
        )?;
        builder.add(
            "PCURVE",
            vec![
                text("reconstructed native trace"),
                reference_to(surface_id),
                reference_to(representation),
            ],
        )
    }

    fn conic(
        &mut self,
        builder: &mut Builder<'_>,
        surface: u64,
        centre: [f64; 2],
        axes: [[f64; 2]; 2],
        radii: [f64; 2],
    ) -> Result<u64, Error> {
        let coordinates = |builder: &mut Builder<'_>, kind, values: [f64; 2]| {
            builder.add(
                kind,
                vec![
                    text(""),
                    Param::List(values.into_iter().map(number).collect::<Result<_, _>>()?),
                ],
            )
        };
        let zero = coordinates(builder, "CARTESIAN_POINT", [0.0; 2])?;
        let unit_x = coordinates(builder, "DIRECTION", [1.0, 0.0])?;
        let placement = builder.add(
            "AXIS2_PLACEMENT_2D",
            vec![text(""), reference_to(zero), reference_to(unit_x)],
        )?;
        let ellipse = builder.add(
            "ELLIPSE",
            vec![
                text(""),
                reference_to(placement),
                number(radii[0])?,
                number(radii[1])?,
            ],
        )?;
        let origin = coordinates(builder, "CARTESIAN_POINT", centre)?;
        let x = coordinates(builder, "DIRECTION", axes[0])?;
        let y = coordinates(builder, "DIRECTION", axes[1])?;
        let transform = builder.add(
            "CARTESIAN_TRANSFORMATION_OPERATOR_2D",
            vec![
                text(""),
                reference_to(x),
                reference_to(y),
                reference_to(origin),
                number(1.0)?,
            ],
        )?;
        let curve = builder.add(
            "CURVE_REPLICA",
            vec![text(""), reference_to(ellipse), reference_to(transform)],
        )?;
        self.wrap(builder, surface, curve)
    }

    fn projected(
        &mut self,
        builder: &mut Builder<'_>,
        surface: u64,
        trace: &super::super::ProjectedTrace,
    ) -> Result<u64, Error> {
        let (source_curve, uv_controls) = match trace {
            super::super::ProjectedTrace::Spline {
                source_curve,
                controls,
            } => (source_curve, controls),
            super::super::ProjectedTrace::Conic {
                centre,
                axes,
                radii,
            } => {
                return self.conic(builder, surface, *centre, *axes, *radii);
            }
        };
        let mut controls = Vec::with_capacity(uv_controls.len());
        for control in uv_controls {
            controls.push(builder.add(
                "CARTESIAN_POINT",
                vec![text(""),
                Param::List(control.iter().map(|value| number(*value)).collect::<Result<_, _>>()?)],
            )?);
        }
        let mut entity = builder
            .parsed
            .source_graph()
            .get(source_curve)
            .ok_or_else(|| {
                invalid(
                    "import.tolerant-broken-reference",
                    "projected spline source is missing",
                )
            })?
            .clone();
        let complex = entity.parts.len() > 1;
        let mut replaced = 0;
        for (kind, parameters) in &mut entity.parts {
            for parameter in parameters.iter_mut() {
                if let Param::List(values) = parameter
                    && values.len() == controls.len()
                    && values.iter().all(|value| matches!(value, Param::Ref(_)))
                {
                    *values = controls
                        .iter()
                        .copied()
                        .map(reference_to)
                        .collect::<Vec<_>>();
                    replaced += 1;
                }
            }
            if kind == "B_SPLINE_CURVE" {
                let count = parameters.len();
                parameters[count - 2] = Param::Enum("U".into());
                parameters[count - 1] = Param::Enum("U".into());
            } else if !complex && kind == "B_SPLINE_CURVE_WITH_KNOTS" {
                parameters[4] = Param::Enum("U".into());
                parameters[5] = Param::Enum("U".into());
            }
        }
        if replaced != 1 {
            return Err(invalid(
                "import.native-projected-controls",
                "one source B-spline control array is required",
            ));
        }
        let curve = builder.complex(entity.parts)?;
        self.wrap(builder, surface, curve)
    }

    fn reconstructed(
        &mut self,
        builder: &mut Builder<'_>,
        original: u64,
        surface_id: u64,
        trim: &NativeTrim,
        dec: &mut Decider,
    ) -> Result<u64, Error> {
        let first = trim.samples.first().expect("finite trace");
        let last = trim.samples.last().expect("finite trace");
        let width = last.curve.parameters[0] - first.curve.parameters[0];
        let kind = if builder.changes.contains_key(&original) {
            Some(2)
        } else {
            parameter_kind(builder, original, dec)?
        };
        let scale = match kind {
            Some(0) => {
                let curve = self.import.native_curve(original).ok_or_else(|| {
                    invalid(
                        "import.native-trace-parameter-map",
                        "native line carrier is missing",
                    )
                })?;
                let a = curve.sample(first.curve, dec)?;
                let b = curve.sample(last.curve, dec)?;
                distance(&a.point_mm, &b.point_mm) / width.abs()
            }
            Some(1) => self.import.source_unit_factors()[1],
            Some(2) => 1.0,
            _ => {
                return Err(invalid(
                    "import.native-trace-parameter-map",
                    "source composite parameter conversion is not established",
                ));
            }
        };
        let reverse = width < 0.0;
        self.spline(
            builder,
            surface_id,
            (0..trim.samples.len()).map(|index| {
                let sample = &trim.samples[if reverse {
                    trim.samples.len() - 1 - index
                } else {
                    index
                }];
                (
                    sample.surface.parameters,
                    sample.curve.parameters[0] * scale,
                )
            }),
        )
    }

    fn basis_curve(
        &mut self,
        builder: &mut Builder<'_>,
        original: u64,
        depth: usize,
        dec: &mut Decider,
    ) -> Result<u64, Error> {
        if let Some(&basis) = self.bases.get(&original) {
            return Ok(basis);
        }
        if depth == 128 {
            return Err(invalid(
                "import.tolerant-broken-reference",
                "cyclic or excessive curve wrappers",
            ));
        }
        if !builder.parsed.source_graph().contains_key(&original) {
            return Ok(original);
        }
        let mut replacements = Vec::new();
        let basis = match builder.parsed.geometry(original, dec)? {
            Geometry::Curve(Curve::SurfaceCurve(parent, ..)) => {
                self.basis_curve(builder, parent, depth + 1, dec)?
            }
            Geometry::Curve(
                Curve::Replica(parent, _)
                | Curve::Offset3d(parent, ..)
                | Curve::Trimmed(parent, ..),
            )
            | Geometry::CompositeSegment(_, _, parent) => {
                let basis = self.basis_curve(builder, parent, depth + 1, dec)?;
                if basis != parent {
                    replacements.push((parent, basis));
                }
                original
            }
            Geometry::Curve(Curve::Composite(segments, _)) => {
                for segment in segments {
                    let basis = self.basis_curve(builder, segment, depth + 1, dec)?;
                    if basis != segment {
                        replacements.push((segment, basis));
                    }
                }
                original
            }
            _ => original,
        };
        let basis = if replacements.is_empty() {
            basis
        } else {
            fn rebind(value: &mut Param, replacements: &[(u64, u64)]) {
                match value {
                    Param::Ref(id) => {
                        if let Some((_, replacement)) =
                            replacements.iter().find(|(source, _)| *source == *id)
                        {
                            *id = *replacement;
                        }
                    }
                    Param::List(values) | Param::Typed(_, values) => {
                        for value in values {
                            rebind(value, replacements);
                        }
                    }
                    _ => {}
                }
            }
            let mut value = builder.parsed.source_graph()[&original].clone();
            for (_, parameters) in &mut value.parts {
                for parameter in parameters {
                    rebind(parameter, &replacements);
                }
            }
            builder.complex(value.parts)?
        };
        self.bases.insert(original, basis);
        Ok(basis)
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
        let pcurve = self.spline(
            builder,
            surface_id,
            endpoints.into_iter().map(|sample| {
                (
                    sample.surface.parameters,
                    sample.curve.parameters[0] * scale,
                )
            }),
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
        // Each emitted topology edge belongs to one body and only its uses.
        if trims.len() > 2 {
            return Err(invalid(
                "import.tolerant-non-manifold",
                "one body's edge has more than two pcurve associations",
            ));
        }
        let mut associated = Vec::with_capacity(trims.len());
        let mut surfaces = Vec::with_capacity(trims.len());
        let mut all_pcurves = true;
        for trim in trims {
            let Some(surface) = builder.source_surface(&trim.face, dec)? else {
                continue;
            };
            surfaces.push(surface);
            let curve = if let Some(projected) = &trim.projected {
                Some(self.projected(builder, surface, projected)?)
            } else if trim.reconstructed {
                Some(self.reconstructed(builder, original, surface, trim, dec)?)
            } else {
                self.affine(builder, edge, original, surface, trim, dec)?
            };
            associated.push(match curve {
                Some(curve) => curve,
                None => match trim.source_pcurve {
                    Some(curve) => curve,
                    None => {
                        all_pcurves = false;
                        surface
                    }
                },
            });
        }
        let basis = self.basis_curve(builder, original, 0, dec)?;
        if associated.is_empty() {
            return Ok(basis);
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
