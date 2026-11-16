//! Fitted recipe construction. Only the numerical admission callback publishes a body.

use super::*;
use catallax_recipe::{
    ElementId, Operation, Origin, PathSegment, Recipe, Role, StepId, SweepOrientation,
};
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

impl RecipeCandidate {
    /// Embeds the fitted geometry and admits only a numerical body. This follows
    /// the existing numerical Boolean selection callback convention. No exact
    /// evaluation checkpoint is returned or made available by this producer.
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
        Operation::Path {segments,..} => segments.iter().any(|s|matches!(s,PathSegment::Helix {..})),
        Operation::Sweep {guide,twist,..} => guide.is_some() || twist.as_ref().is_some_and(|n|!n.value().is_zero()),
        _=>false,
    })
}

pub(crate) fn path_needs_numerical(recipe:&Recipe,path:&StepId)->bool {
    if recipe.step(path).is_some_and(|step|matches!(&step.op,Operation::Path {segments,..}
        if segments.iter().any(|s|matches!(s,PathSegment::Helix {..})))) {return true;}
    recipe.steps.iter().any(|step|matches!(&step.op,Operation::Sweep {path:spine,guide,twist,..}
        if (guide.is_some() || twist.as_ref().is_some_and(|n|!n.value().is_zero()))
            && (spine==path || guide.as_ref()==Some(path))))
}

/// Constructs the fitted recipe privately. The caller must supply the numerical
/// embedding/admission consumer before obtaining an editable numerical body.
pub fn prepare_recipe(
    recipe: &Recipe,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<RecipeCandidate, crate::BuildError> {
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
        profile: &StepId,
        path: &StepId,
        step: &StepId,
        orientation: SweepOrientation,
        guide: Option<&StepId>,
        twist: Option<&Number>,
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational), Refusal> {
        let sketch = recipe
            .steps
            .iter()
            .find(|s| &s.id == profile)
            .ok_or_else(|| bad("missing sweep profile"))?;
        let Operation::Sketch { plane, loops } = &sketch.op else {
            return Err(bad("sweep needs a sketch"));
        };
        let (id,motion)=recipe_path(recipe,path)?;
        let guide=guide.map(|id|recipe_path(recipe,id).map(|(_,motion)|motion)).transpose()?;
        let definition = SweepRecipe {
            name: step.to_string(),
            loops: loops.clone(),
            plane: Some(plane.clone()),
            path: motion,
            orientation,
            guide,
            reference: [1, 0, 0].map(Into::into),
            twist_degrees: twist.cloned().unwrap_or_else(||0.into()),
            end_scale: 1.into(),
            tolerance_mm: self.tolerance.clone().into(),
            max_sections: 4096,
        };
        let (body, radius, checks) = append(model, &definition, step, Some(&id), dec)?;
        self.checks.insert(step.clone(), checks);
        Ok((body, radius))
    }

    pub fn stock(
        &mut self,
        model: &mut Model,
        operation: &Operation,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<(BodyId, Rational)>, Refusal> {
        let Operation::Cylinder {
            radius,
            height,
            pose,
        } = operation
        else {
            return Ok(None);
        };
        let placement = crate::profile::placement(pose, dec)?;
        let end = placement.apply_point(
            &catallax_shapes::RPoint3::new(q(0), q(0), height.value().clone()).to_real(),
        );
        let coordinate = |x: &Real| {
            x.as_rational()
                .map(Into::into)
                .ok_or_else(|| bad("stock pose is not rational"))
        };
        let end = RecipePoint([
            coordinate(&end.x)?,
            coordinate(&end.y)?,
            coordinate(&end.z)?,
        ]);
        let definition = SweepRecipe {
            name: step.to_string(),
            loops: vec![Loop(vec![catallax_recipe::Segment::Circle {
                id: ElementId::new("wall").expect("fixed identity"),
                centre: catallax_recipe::Point2([0.into(), 0.into()]),
                radius: radius.clone(),
            }])],
            plane: Some(pose.clone()),
            path: SweepPath::Line {
                start: RecipePoint(pose.translation.0.clone()),
                end,
            },
            orientation: SweepOrientation::Fixed,
            guide: None,
            reference: [1, 0, 0].map(Into::into),
            twist_degrees: 0.into(),
            end_scale: 1.into(),
            tolerance_mm: self.tolerance.clone().into(),
            max_sections: 1,
        };
        let (body, radius, checks) = append(model, &definition, step, None, dec)?;
        self.checks.insert(step.clone(), checks);
        Ok(Some((body, radius)))
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
        marks: Default::default(),
    });
    Ok((body, radius, checks))
}

fn recipe_path(recipe:&Recipe,path:&StepId)->Result<(ElementId,SweepPath),Refusal> {
    let rail=recipe.step(path).ok_or_else(||bad("missing sweep path"))?;
    let Operation::Path {start,segments}=&rail.op else {return Err(bad("a sweep input must be a path"));};
    let [segment]=segments.as_slice() else {
        return Err(bad("a numerical sweep currently needs one complete authored segment per path"));
    };
    let motion=match segment {
        PathSegment::Line {to,..}=>SweepPath::Line {start:start.clone(),end:to.clone()},
        PathSegment::Helix {centre,axis,angle,rise,..}=>SweepPath::Helix {
            start:start.clone(),centre:centre.clone(),axis:axis.clone(),angle:angle.clone(),rise:rise.clone()},
        PathSegment::Spline {degree,control_points,weights,knots,..}=>SweepPath::Spline {
            degree:*degree,control_points:control_points.clone(),weights:weights.clone(),knots:knots.clone()},
        PathSegment::Arc {centre,to,axis,..}=>{
            let point=|p:&RecipePoint|catallax_shapes::RPoint3::new(
                p.0[0].value().clone(),p.0[1].value().clone(),p.0[2].value().clone());
            let a=point(start).minus(&point(centre));
            let b=point(to).minus(&point(centre));
            let normal=catallax_shapes::RVector3::new(axis.0[0].value().clone(),
                axis.0[1].value().clone(),axis.0[2].value().clone());
            if normal.is_zero() || a.is_zero() || a==b || a.norm_squared()!=b.norm_squared()
                || !a.dot(&normal).is_zero() || !b.dot(&normal).is_zero() {
                return Err(bad("an arc needs distinct endpoints on one nonzero perpendicular circle"));
            }
            let av=[a.x.to_f64_nearest(),a.y.to_f64_nearest(),a.z.to_f64_nearest()];
            let bv=[b.x.to_f64_nearest(),b.y.to_f64_nearest(),b.z.to_f64_nearest()];
            let n=unit(axis.0.each_ref().map(number))?;
            let mut angle=libm::atan2(dot(n,cross(av,bv)),dot(av,bv));
            if angle<=0.0 {angle+=std::f64::consts::TAU;}
            SweepPath::Helix {start:start.clone(),centre:centre.clone(),axis:axis.clone(),
                angle:rational(angle.to_degrees())?.into(),rise:0.into()}
        }
    };
    Ok((segment.id().clone(),motion))
}
