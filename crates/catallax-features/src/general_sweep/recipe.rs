//! Fitted recipe construction. Only the numerical admission callback publishes a body.

use super::*;
use catallax_recipe::{ElementId, Operation, Origin, PathSegment, Recipe, Role, Step, StepId};
use catallax_shapes::PointDef;
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Shell, Vertex};

/// A private fitted model awaiting numerical embedding and admission.
/// Its geometry is a rational approximation, not an exact answer to the recipe.
pub struct RecipeCandidate {
    model: Model,
    body: BodyId,
    recipe: Recipe,
    options: NumericalOptions,
    checks: BTreeMap<StepId, SweepChecks>,
}

/// An ordinary saved recipe with an admitted numerical sweep result.
pub struct RecipeSwept {
    /// Admitted boundary with its numerical tolerance and checks.
    pub body: NumericalBody,
    /// Original authored definition for save and parameter rebuild.
    pub recipe: Recipe,
    /// Pointwise section-law and profile fitting checks per sweep.
    pub checks: BTreeMap<StepId, SweepChecks>,
    /// Original operation origins of every retained fitted face.
    pub face_origins: BTreeMap<String, Vec<Origin>>,
}

fn require_budget() -> Result<(), Refusal> {
    if catallax_numbers::operation_budget_active() {
        Ok(())
    } else {
        Err(Refusal::new(
            "sweep.operation-budget-required",
            "prepare_recipe and admit require one caller-owned work allowance covering both stages",
        ))
    }
}

impl RecipeCandidate {
    /// Embeds the fitted geometry and admits only a numerical body. This follows
    /// the existing numerical Boolean selection callback convention. No exact
    /// evaluation checkpoint is returned or made available by this producer.
    /// The caller must keep the same operation-budget scope active across
    /// `prepare_recipe` and this admission callback; neither stage resets it.
    pub fn admit<T>(
        self,
        admit: impl FnOnce(
            Model,
            BodyId,
            &NumericalOptions,
            &mut Decider,
        ) -> Result<(NumericalBody, T), Refusal>,
        dec: &mut Decider,
    ) -> Result<(RecipeSwept, T), Refusal> {
        require_budget()?;
        let face_origins = self
            .model
            .faces_of_body(self.body)
            .into_iter()
            .map(|id| (format!("{id:?}"), self.model.face_origins(id).to_vec()))
            .collect();
        let (body, consumer) = admit(self.model, self.body, &self.options, dec)?;
        Ok((
            RecipeSwept {
                body,
                recipe: self.recipe,
                checks: self.checks,
                face_origins,
            },
            consumer,
        ))
    }
}

/// Whether the exact evaluator must not claim authority for this recipe.
pub fn requires_numerical(recipe: &Recipe) -> bool {
    recipe.steps.iter().any(|step| match &step.op {
        Operation::Path { segments, .. } => segments
            .iter()
            .any(|s| matches!(s, PathSegment::Helix { .. })),
        Operation::Sweep { guide, twist, .. } => {
            guide.is_some() || twist.as_ref().is_some_and(|n| !n.value().is_zero())
        }
        _ => false,
    })
}

pub(crate) fn path_needs_numerical(recipe: &Recipe, path: &StepId) -> bool {
    if recipe.step(path).is_some_and(|step| {
        matches!(&step.op,Operation::Path {segments,..}
        if segments.iter().any(|s|matches!(s,PathSegment::Helix {..})))
    }) {
        return true;
    }
    let mut numerical = false;
    for step in &recipe.steps {
        if let Operation::Sweep {
            path: spine,
            guide,
            twist,
            ..
        } = &step.op
        {
            if spine == path {
                if guide.is_none() && twist.as_ref().is_none_or(|n| n.value().is_zero()) {
                    return false;
                }
                numerical = true;
            }
            numerical |= guide.as_ref() == Some(path);
        }
    }
    numerical
}

/// Constructs the fitted recipe privately. The caller must supply the numerical
/// embedding/admission consumer before obtaining an editable numerical body.
/// Requires a caller-owned operation budget that remains active through `admit`.
pub fn prepare_recipe(
    recipe: &Recipe,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<RecipeCandidate, crate::BuildError> {
    require_budget().map_err(|refusal| crate::BuildError::Geometry {
        step: recipe.result.clone(),
        refusal,
    })?;
    let mut context = Context {
        tolerance: options.tolerance_mm.clone(),
        checks: BTreeMap::new(),
    };
    let (model, body) = crate::evaluate::numerical_sweep_impl(recipe, &mut context, dec)?;
    Ok(RecipeCandidate {
        model,
        body,
        recipe: recipe.clone(),
        options,
        checks: context.checks,
    })
}

pub(crate) struct Context {
    pub tolerance: Rational,
    pub checks: BTreeMap<StepId, SweepChecks>,
}
impl Context {
    pub fn sweep(
        &mut self,
        model: &mut Model,
        recipe: &Recipe,
        step: &Step,
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational), Refusal> {
        let Operation::Sweep {
            profile,
            path,
            orientation,
            guide,
            twist,
        } = &step.op
        else {
            return Err(bad("the numerical sweep context needs a sweep step"));
        };
        let sketch = recipe
            .steps
            .iter()
            .find(|s| &s.id == profile)
            .ok_or_else(|| bad("missing sweep profile"))?;
        let Operation::Sketch { plane, loops } = &sketch.op else {
            return Err(bad("sweep needs a sketch"));
        };
        let (id, motion) = recipe_path(recipe, path)?;
        let guide = guide
            .as_ref()
            .map(|id| recipe_path(recipe, id).map(|(_, motion)| motion))
            .transpose()?;
        let definition = SweepRecipe {
            name: step.id.to_string(),
            loops: loops.clone(),
            plane: Some(plane.clone()),
            path: motion,
            orientation: *orientation,
            guide,
            reference: [1, 0, 0].map(Into::into),
            twist_degrees: twist.clone().unwrap_or_else(|| 0.into()),
            end_scale: 1.into(),
            tolerance_mm: self.tolerance.clone().into(),
            max_sections: 4096,
        };
        let (body, radius, checks) = append(model, &definition, &step.id, Some(&id), dec)?;
        self.checks.insert(step.id.clone(), checks);
        Ok((body, radius))
    }
}

fn append(
    model: &mut Model,
    recipe: &SweepRecipe,
    step: &StepId,
    path: Option<&ElementId>,
    dec: &mut Decider,
) -> Result<(BodyId, Rational, SweepChecks), Refusal> {
    let (definition, checks, roles) = super::definition(recipe, dec)?;
    let mut radius = q(0);
    for face in &definition.faces {
        if let Surface::BSpline(support) = &face.surface {
            for point in support.carrier().controls() {
                let mut bound = q(0);
                for value in [&point.x, &point.y, &point.z] {
                    let ball = dec.enclose_budgeted_at(value, 128)?;
                    bound = bound + ball.lo().abs().max(ball.hi().abs());
                }
                radius = radius.max(bound);
            }
        }
    }
    let vertices: Vec<_> = definition
        .vertices
        .iter()
        .map(|v| {
            let [x, y, z] = v.position.clone();
            let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
                x.into(),
                y.into(),
                z.into(),
            )));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let edges: Vec<_> = definition
        .edges
        .iter()
        .map(|e| {
            let curve = model.geom_mut().add_curve(e.curve.clone());
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[e.vertices[0]]),
                end: Some(vertices[e.vertices[1]]),
            })
        })
        .collect();
    let mut faces = Vec::new();
    for (face, mut role) in definition.faces.into_iter().zip(roles) {
        if let Role::Side { path: along, .. } = &mut role {
            *along = path.cloned();
        }
        let paths = face
            .loops
            .into_iter()
            .map(|lp| lp.into_iter().map(|t| (edges[t.edge], t.sense)).collect())
            .collect();
        // These incidences describe only the fitted rational definition: sides
        // share native rows and cap rows use one exact rational affine plane.
        faces.push(crate::extrude::face(
            model,
            face.surface,
            face.sense,
            paths,
            None,
            vec![Origin {
                step: step.clone(),
                role,
            }],
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: Default::default(),
    });
    Ok((body, radius, checks))
}

fn recipe_path(recipe: &Recipe, path: &StepId) -> Result<(ElementId, SweepPath), Refusal> {
    let rail = recipe.step(path).ok_or_else(|| bad("missing sweep path"))?;
    let Operation::Path { start, segments } = &rail.op else {
        return Err(bad("a sweep input must be a path"));
    };
    let [segment] = segments.as_slice() else {
        return Err(bad(
            "a numerical sweep currently needs one complete authored segment per path",
        ));
    };
    let motion = match segment {
        PathSegment::Line { to, .. } => SweepPath::Line {
            start: start.clone(),
            end: to.clone(),
        },
        PathSegment::Helix {
            centre,
            axis,
            angle,
            rise,
            ..
        } => SweepPath::Helix {
            start: start.clone(),
            centre: centre.clone(),
            axis: axis.clone(),
            angle: angle.clone(),
            rise: rise.clone(),
        },
        PathSegment::Spline {
            degree,
            control_points,
            weights,
            knots,
            ..
        } => SweepPath::Spline {
            degree: *degree,
            control_points: control_points.clone(),
            weights: weights.clone(),
            knots: knots.clone(),
        },
        PathSegment::Arc {
            centre, to, axis, ..
        } => {
            let point = |p: &RecipePoint| {
                catallax_shapes::RPoint3::new(
                    p.0[0].value().clone(),
                    p.0[1].value().clone(),
                    p.0[2].value().clone(),
                )
            };
            let a = point(start).minus(&point(centre));
            let b = point(to).minus(&point(centre));
            let normal = catallax_shapes::RVector3::new(
                axis.0[0].value().clone(),
                axis.0[1].value().clone(),
                axis.0[2].value().clone(),
            );
            if normal.is_zero()
                || a.is_zero()
                || a == b
                || a.norm_squared() != b.norm_squared()
                || !a.dot(&normal).is_zero()
                || !b.dot(&normal).is_zero()
            {
                return Err(bad(
                    "an arc needs distinct endpoints on one nonzero perpendicular circle",
                ));
            }
            let av = [
                a.x.to_f64_nearest(),
                a.y.to_f64_nearest(),
                a.z.to_f64_nearest(),
            ];
            let bv = [
                b.x.to_f64_nearest(),
                b.y.to_f64_nearest(),
                b.z.to_f64_nearest(),
            ];
            let n = unit(axis.0.each_ref().map(number))?;
            let mut angle = libm::atan2(dot(n, cross(av, bv)), dot(av, bv));
            if angle <= 0.0 {
                angle += std::f64::consts::TAU;
            }
            SweepPath::Helix {
                start: start.clone(),
                centre: centre.clone(),
                axis: axis.clone(),
                angle: rational(angle.to_degrees())?.into(),
                rise: 0.into(),
            }
        }
    };
    Ok((segment.id().clone(), motion))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plain_sweep_can_share_the_numerical_sweeps_path() {
        let mut input: serde_json::Value =
            serde_json::from_str(include_str!("../../tests/fixtures/straight-sweep.json")).unwrap();
        input["steps"][1]["start"] = serde_json::json!(["7", "-3", "11"]);
        input["steps"][1]["segments"][0]["to"] = serde_json::json!(["7", "-3", "14"]);
        let mut twisted = input["steps"][2].clone();
        twisted["id"] = "twisted".into();
        twisted["twist"] = "90".into();
        input["steps"].as_array_mut().unwrap().insert(2, twisted);
        let recipe: Recipe = serde_json::from_value(input).unwrap();
        let (result, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let candidate = prepare_recipe(
                &recipe,
                NumericalOptions {
                    tolerance_mm: Rational::from_ratio(1, 100).unwrap(),
                },
                &mut dec,
            )
            .unwrap();
            assert_eq!(candidate.model.body(candidate.body).kind, BodyKind::Solid);
            let mut caps = Vec::new();
            for face in candidate.model.faces_of_body(candidate.body) {
                let Surface::Plane(plane) = candidate
                    .model
                    .geom()
                    .surface(candidate.model.face(face).surface)
                else {
                    panic!("the ordinary consumer must still produce a planar prism");
                };
                if dec.sign(&plane.normal.x).is_coincident()
                    && dec.sign(&plane.normal.y).is_coincident()
                {
                    caps.push(plane.origin.z.as_rational().unwrap());
                }
            }
            caps.sort();
            assert_eq!(caps, vec![q(11), q(14)]);
            Ok::<_, Refusal>(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }

    #[test]
    fn recipe_preparation_and_admission_require_a_caller_work_scope() {
        let recipe: Recipe =
            serde_json::from_str(include_str!("../../tests/fixtures/straight-sweep.json")).unwrap();
        let options = || NumericalOptions {
            tolerance_mm: Rational::from_ratio(1, 100).unwrap(),
        };
        let error = prepare_recipe(&recipe, options(), &mut Decider::default())
            .err()
            .expect("unbudgeted construction must refuse");
        let crate::BuildError::Geometry { refusal, .. } = error else {
            panic!("expected a work-scope refusal")
        };
        assert_eq!(refusal.code, "sweep.operation-budget-required");

        let (prepared, work) = catallax_numbers::with_operation_budget(20_000_000, || {
            prepare_recipe(&recipe, options(), &mut Decider::default())
        });
        assert!(!work.exhausted);
        let candidate = prepared.unwrap();
        let error = candidate
            .admit::<()>(
                |_, _, _, _| panic!("an admission callback must not run outside a work scope"),
                &mut Decider::default(),
            )
            .err()
            .expect("unbudgeted admission must refuse");
        assert_eq!(error.code, "sweep.operation-budget-required");

        let (result, work) = catallax_numbers::with_operation_budget(1, || {
            prepare_recipe(&recipe, options(), &mut Decider::default())
        });
        let crate::BuildError::Geometry { refusal, .. } = result.err().unwrap() else {
            panic!("expected caller work exhaustion")
        };
        assert_eq!(refusal.code, "numbers.operation-budget");
        assert!(work.exhausted && work.used <= 1);
    }
}
