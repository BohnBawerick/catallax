//! Numerical fillet retention without an exact checkpoint or model conversion.

use super::*;
use catallax_recipe::{EdgeSelection, FilletCorner, FilletOverflow, Origin, Step};
use catallax_topology::numerical::{NumericalBody, NumericalOptions, tolerance};

/// Rule T1's floor for own fitted geometry built from an exact body: the
/// floor in the frame of the body's vertices.
fn exact_floor(model: &Model, body: BodyId, dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut points = Vec::new();
    for vertex in model.vertices_of_body(body) {
        dec.charge_work(1)?;
        let p = model.geom().position(model.vertex(vertex).point);
        points.push([
            dec.enclose_budgeted_at(&p.x, 128)?.mid(),
            dec.enclose_budgeted_at(&p.y, 128)?.mid(),
            dec.enclose_budgeted_at(&p.z, 128)?.mid(),
        ]);
    }
    Ok(tolerance::BodyFrame::of_points(&points).floor_mm())
}

/// Rule T2's tightening near small features: the smallest distance between
/// two entities of `source` that share no topology, `L`, tightens the fit
/// target to `L/100` where that is under the target.
fn near_small_features(
    source: &NumericalBody,
    eps: &Rational,
    target: Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let reach = tolerance::fitted_feature_reach_mm(eps);
    let smallest = tolerance::features_within(source, &reach, dec)?
        .into_iter()
        .map(|pair| pair.distance_mm)
        .min();
    Ok(match smallest {
        Some(size) => tolerance::own_fitted(eps, Some(&size))
            .target_mm
            .min(target),
        None => target,
    })
}

/// An exact continuation checkpoint or a separately admitted numerical result.
#[derive(Clone, Debug)]
pub enum EvaluationOutcome {
    /// Every retained body and continuation state remains on the exact path.
    Exact(EvaluationCheckpoint),
    /// Scoped numerical checks, never an exact continuation checkpoint.
    Numerical(Arc<NumericalFilletEvaluation>),
}

/// Immutable authored controls, numerical geometry and canonical face origins.
#[derive(Clone, Debug)]
pub struct NumericalFilletEvaluation {
    recipe: Arc<Recipe>,
    body: NumericalBody,
    origins: BTreeMap<String, Vec<Origin>>,
}
impl NumericalFilletEvaluation {
    /// Numerical geometry, native trims, declared tolerance and edit history.
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }
    /// Authored controls retained for replay and parameter edits.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }
    /// Canonical origins keyed by the numerical definition's stable face keys.
    pub fn face_origins(&self) -> &BTreeMap<String, Vec<Origin>> {
        &self.origins
    }
    /// Consumes the result into its admitted geometry, authored controls and origins.
    pub fn into_parts(self) -> (NumericalBody, Arc<Recipe>, BTreeMap<String, Vec<Origin>>) {
        (self.body, self.recipe, self.origins)
    }
}

#[derive(Default)]
pub(super) struct Context {
    bodies: BTreeMap<StepId, Arc<NumericalFilletEvaluation>>,
}
impl Context {
    /// Refuses a step that reads a body with tolerant entities and is about
    /// to be built by the exact operations.
    ///
    /// Such a body never reaches the exact intersector (Checkpoint 2
    /// decisions, section 4.5). An operation goes on from one only through a
    /// construction of its own that decides with the filtered predicate. The
    /// hole has one ([`Self::hole`]). Every other operation refuses here by
    /// name until it has one.
    pub(super) fn exact_inputs(&self, step: &Step) -> Result<(), Refusal> {
        if self.bodies.is_empty() {
            return Ok(());
        }
        if let Some((input, _)) = step
            .op
            .inputs()
            .into_iter()
            .find(|(id, _)| self.bodies.contains_key(*id))
        {
            return Err(Refusal::new(
                "features.tolerant-input",
                format!(
                    "operation {} reads body {input}, which has tolerant entities, and has no construction on such a body",
                    step.op.name()
                ),
            ));
        }
        Ok(())
    }

    /// Whether `step` built a body with tolerant entities.
    pub(super) fn has(&self, step: &StepId) -> bool {
        self.bodies.contains_key(step)
    }

    /// Drills a hole in a body with tolerant entities (Checkpoint 2
    /// decisions, section 4.5: "a hole after a numerical fillet works").
    ///
    /// Each placement is a through hole of one radius, written into the
    /// body's definition and admitted again
    /// ([`crate::hole::through_tolerant`]). A blind hole, a head and a tool
    /// that is not square to the faces it crosses refuse by name.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn hole(
        &mut self,
        target: &StepId,
        at: &[catallax_recipe::HolePlacement],
        diameter: &Rational,
        depth: &catallax_recipe::HoleDepth,
        tip: &catallax_recipe::HoleTip,
        head: Option<&catallax_recipe::HoleHead>,
        step: &StepId,
        recipe: &Recipe,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        crate::hole::validate(diameter, depth, tip, head)?;
        if !matches!(depth, catallax_recipe::HoleDepth::Through) || head.is_some() {
            return Err(Refusal::new(
                "hole.numerical-tool",
                "only a through hole with no head has a construction on a body with tolerant entities",
            ));
        }
        let source = Arc::clone(&self.bodies[target]);
        let mut body = source.body.clone();
        let mut origins = source.origins.clone();
        // The tool's bands are cut at a bound of the body about the origin.
        // A through hole is read along its whole axis, so the bound only
        // has to be a size: the largest vertex coordinate sum.
        let mut bound = Rational::zero();
        for vertex in &body.definition().vertices {
            let mut reach = Rational::zero();
            for value in &vertex.position {
                let size = value.abs();
                dec.reserve_exact_arithmetic(&reach, &size, false)?;
                reach = reach + size;
            }
            dec.reserve_uncounted_arithmetic(&bound, &reach, true)?;
            bound = bound.max(reach);
        }
        for placement in at {
            let tool =
                crate::hole::HoleTool::new(placement, diameter, depth, tip, head, &bound, dec)?;
            let key = format!("{step}/{}", placement.id);
            let drilled = crate::hole::through_tolerant::through(&body, &key, &tool, dec)?;
            for (wall, [low, high]) in drilled.walls {
                // The tool's bands meet at the placement point: the entry
                // above it and the wall below.
                let mut parts = Vec::new();
                if low.sign() == catallax_numbers::Sign::Negative {
                    parts.push(catallax_recipe::HolePart::Entry);
                }
                if high.sign() == catallax_numbers::Sign::Positive {
                    parts.push(catallax_recipe::HolePart::Wall);
                }
                origins.insert(
                    wall,
                    parts
                        .into_iter()
                        .map(|part| Origin {
                            step: step.clone(),
                            role: Role::Hole {
                                at: placement.id.clone(),
                                part,
                            },
                        })
                        .collect(),
                );
            }
            body = drilled.body;
        }
        self.bodies.insert(
            step.clone(),
            Arc::new(NumericalFilletEvaluation {
                recipe: Arc::new(recipe.clone()),
                body,
                origins,
            }),
        );
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn fillet(
        &mut self,
        model: &Model,
        target: &BodyValue,
        selection: &EdgeSelection,
        radius: &Rational,
        overflow: Option<&FilletOverflow>,
        corner: Option<FilletCorner>,
        earlier: &[crate::fillet::FilletDefinition],
        step: &StepId,
        recipe: &Recipe,
        exact_refusal: Refusal,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        // A radius that is not positive has no rolling ball to construct.
        if radius.sign() != catallax_numbers::Sign::Positive {
            return Err(exact_refusal);
        }
        let result = (|| {
            if !rolls(overflow, corner) {
                return Err(Refusal::new(
                    "features.numerical-fillet-closure",
                    "the rolling-ball construction cannot honour explicit neighbour limits or a non-cliff corner closure",
                ));
            }
            let selected = crate::fillet_build::select(
                model,
                target.body,
                selection,
                recipe.format != catallax_recipe::FORMAT,
                earlier,
                dec,
            )?;
            if selected.is_empty() {
                return Err(Refusal::new(
                    "features.numerical-fillet-selection",
                    "the selection contains no edges to round",
                ));
            }
            let faults = model.check(target.body);
            if !faults.is_empty() {
                return Err(Refusal::new(
                    "features.numerical-fillet-source",
                    format!("source topology is invalid: {faults:?}"),
                ));
            }
            // Rule T2 for own fitted geometry: the fitters aim at the
            // target, and every body admitted here takes the cap as its
            // class cap. The floor is rule T1's in the frame of the exact
            // body's vertices.
            let eps = exact_floor(model, target.body, dec)?;
            let fitted = tolerance::own_fitted(&eps, None);
            let cap = NumericalOptions {
                tolerance_mm: fitted.cap_mm.clone(),
                vertex_tolerance_mm: None,
            };
            let tolerance = fitted.target_mm.clone();
            // One edge between two retained spline faces keeps the
            // construction that was fitted and measured for it.
            if let [edge] = selected.as_slice()
                && spline_faces(model, target.body, *edge)
            {
                let name = model.edge_name(target.body, *edge).ok_or_else(|| {
                    Refusal::new(
                        "features.numerical-fillet-name",
                        format!("selected edge {edge:?} has no canonical recipe name"),
                    )
                })?;
                let (definition, trims) = crate::numerical_embedding::definition_from_exact_model(
                    model,
                    target.body,
                    recipe.name.clone(),
                    &tolerance,
                    dec,
                )?;
                let source = catallax_topology::numerical::build_with_native_trims(
                    definition,
                    trims,
                    cap.clone(),
                    dec,
                )?;
                let tolerance = near_small_features(&source, &eps, tolerance, dec)?;
                let edge_key = format!("{edge:?}");
                let mut origins: BTreeMap<_, _> = model
                    .faces_of_body(target.body)
                    .into_iter()
                    .map(|face| (format!("{face:?}"), model.face_origins(face).to_vec()))
                    .collect();
                origins.insert(
                    format!("{step}/blend/{edge_key}"),
                    vec![Origin {
                        step: step.clone(),
                        role: Role::Blend { edge: name },
                    }],
                );
                let body = crate::fillet::numerical::build(
                    source,
                    &edge_key,
                    step.as_str(),
                    radius.clone(),
                    tolerance,
                    dec,
                )?
                .body()
                .clone();
                let origins = body
                    .definition()
                    .faces
                    .iter()
                    .map(|face| {
                        origins
                            .remove(&face.key)
                            .map(|origins| (face.key.clone(), origins))
                            .ok_or_else(|| {
                                Refusal::new(
                                    "features.numerical-fillet-name",
                                    format!("output face {} lost its canonical origins", face.key),
                                )
                            })
                    })
                    .collect::<Result<_, _>>()?;
                return Ok(NumericalFilletEvaluation {
                    recipe: Arc::new(recipe.clone()),
                    body,
                    origins,
                });
            }
            // A radius no ball can roll at is refused before anything is built.
            let grouped = crate::fillet_build::chains::group(
                model,
                target.body,
                &selected,
                radius,
                &tolerance,
                dec,
            )?;
            crate::fillet_build::fillet_limit::check_rolling_contacts(
                model,
                target.body,
                &selected,
                grouped.iter().flat_map(|chain| chain.edges.iter()),
                radius,
                dec,
            )?;
            if let Some((body, origins)) = crate::fillet::numerical::native::build(
                model,
                target.body,
                &selected,
                radius,
                step,
                recipe.name.clone(),
                &tolerance,
                &cap,
                dec,
            )? {
                return Ok(NumericalFilletEvaluation {
                    recipe: Arc::new(recipe.clone()),
                    body,
                    origins,
                });
            }
            let chains = crate::fillet::numerical::rolling::construct(
                model,
                target.body,
                &selected,
                radius,
                &tolerance,
                dec,
            )?;
            let (definition, trims) = crate::numerical_embedding::definition_from_exact_model(
                model,
                target.body,
                recipe.name.clone(),
                &tolerance,
                dec,
            )?;
            let source =
                catallax_topology::numerical::build_with_native_trims(definition, trims, cap, dec)?;
            let mut source_origins: BTreeMap<_, _> = model
                .faces_of_body(target.body)
                .into_iter()
                .map(|face| (format!("{face:?}"), model.face_origins(face).to_vec()))
                .collect();
            for strip in chains.iter().flat_map(|chain| &chain.strips) {
                let edge = strip.source.edge;
                let name = model.edge_name(target.body, edge).ok_or_else(|| {
                    Refusal::new(
                        "features.numerical-fillet-name",
                        format!("rolled edge {edge:?} has no canonical recipe name"),
                    )
                })?;
                source_origins.insert(
                    format!("{step}/blend/{edge:?}"),
                    vec![Origin {
                        step: step.clone(),
                        role: Role::Blend { edge: name },
                    }],
                );
            }
            let body = if chains.len() == 1 && chains[0].strips.len() == 1 && !chains[0].closed {
                crate::fillet::numerical::body::cap_strip(
                    &source,
                    &chains[0].strips[0],
                    step.as_str(),
                    &tolerance,
                    dec,
                )?
            } else {
                let corners = crate::fillet::numerical::vertex::construct(
                    model, &chains, radius, &tolerance, dec,
                )?;
                for corner in &corners {
                    let vertex =
                        model
                            .vertex_name(target.body, corner.source)
                            .ok_or_else(|| {
                                Refusal::new(
                                    "features.numerical-fillet-name",
                                    format!(
                                        "corner {:?} has no canonical recipe name",
                                        corner.source
                                    ),
                                )
                            })?;
                    source_origins.insert(
                        format!("{step}/corner/{:?}", corner.source),
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Corner { vertex },
                        }],
                    );
                }
                crate::fillet::numerical::vertex_body::build(
                    &source,
                    model,
                    &chains,
                    &corners,
                    step.as_str(),
                    &tolerance,
                    dec,
                )?
            };
            let origins = body
                .definition()
                .faces
                .iter()
                .map(|face| {
                    source_origins
                        .remove(&face.key)
                        .map(|origins| (face.key.clone(), origins))
                        .ok_or_else(|| {
                            Refusal::new(
                                "features.numerical-fillet-name",
                                format!("output face {} lost its canonical origins", face.key),
                            )
                        })
                })
                .collect::<Result<_, _>>()?;
            Ok(NumericalFilletEvaluation {
                recipe: Arc::new(recipe.clone()),
                body,
                origins,
            })
        })()
        .map_err(|refusal: Refusal| both_refusals(exact_refusal, refusal))?;
        self.bodies.insert(step.clone(), Arc::new(result));
        Ok(())
    }

    /// Bends or twists an exact body (RECIPE.md section 4.20). Both build
    /// only as fitted geometry: the exact body is embedded, every face and
    /// edge is refitted through the warp, and admission measures every gap.
    /// The tolerances are the rolling ball's: rule T2's target for the fits
    /// and its cap for admission, over rule T1's floor.
    pub(super) fn deform(
        &mut self,
        model: &Model,
        target: &BodyValue,
        warp: &crate::deform::Warp,
        step: &StepId,
        recipe: &Recipe,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let faults = model.check(target.body);
        if !faults.is_empty() {
            return Err(Refusal::new(
                "deform.source",
                format!("source topology is invalid: {faults:?}"),
            ));
        }
        let eps = exact_floor(model, target.body, dec)?;
        let fitted = tolerance::own_fitted(&eps, None);
        let (source, trims) = crate::numerical_embedding::definition_from_exact_model(
            model,
            target.body,
            recipe.name.clone(),
            &fitted.target_mm,
            dec,
        )?;
        let (definition, trims, sampled) = crate::deform::warp(
            &source,
            &trims,
            warp,
            fitted.target_mm.to_f64_nearest(),
            dec,
        )?;
        // T1's floor is set by the frame the operation runs in; a warp moves
        // the body, so the floor covers the warped vertices as well.
        let moved: Vec<_> = definition
            .vertices
            .iter()
            .map(|v| v.position.clone())
            .collect();
        let eps = eps.max(tolerance::BodyFrame::of_points(&moved).floor_mm());
        let cap = NumericalOptions {
            tolerance_mm: tolerance::own_fitted(&eps, None).cap_mm,
            vertex_tolerance_mm: None,
        };
        let body = catallax_topology::numerical::build_with_sampled_carriers(
            definition, trims, sampled, cap, dec,
        )?;
        let origins: BTreeMap<_, _> = model
            .faces_of_body(target.body)
            .into_iter()
            .map(|face| (format!("{face:?}"), model.face_origins(face).to_vec()))
            .collect();
        let origins = body
            .definition()
            .faces
            .iter()
            .map(|face| {
                origins
                    .get(&face.key)
                    .map(|origins| (face.key.clone(), origins.clone()))
                    .ok_or_else(|| {
                        Refusal::new(
                            "deform.name",
                            format!("output face {} lost its canonical origins", face.key),
                        )
                    })
            })
            .collect::<Result<_, _>>()?;
        self.bodies.insert(
            step.clone(),
            Arc::new(NumericalFilletEvaluation {
                recipe: Arc::new(recipe.clone()),
                body,
                origins,
            }),
        );
        Ok(())
    }
}

/// Whether both faces of an edge are retained rational spline atlases.
fn spline_faces(model: &Model, body: BodyId, edge: catallax_topology::EdgeId) -> bool {
    let mut faces = 0;
    for face in model.faces_of_body(body) {
        let uses_edge = model.face(face).loops.iter().any(|boundary| {
            matches!(model.loop_(*boundary), catallax_topology::Loop::Coedges(uses)
                if uses.iter().any(|coedge| model.coedge(*coedge).edge == edge))
        });
        if !uses_edge {
            continue;
        }
        if !matches!(
            model.geom().surface(model.face(face).surface),
            catallax_shapes::Surface::BSpline(catallax_shapes::SplineSupport::RationalAtlas(_))
        ) {
            return false;
        }
        faces += 1;
    }
    faces == 2
}

/// Whether the rolling ball builds what these options ask for.
///
/// Rolling over neighbours and a cliff corner are what it does, and what the
/// exact builders do when the options are left out, so naming either changes
/// nothing. Explicit station limits on the neighbours and a smooth corner are
/// constructions it does not have.
fn rolls(overflow: Option<&FilletOverflow>, corner: Option<FilletCorner>) -> bool {
    !matches!(overflow, Some(FilletOverflow::KeepNeighbours(_)))
        && !matches!(corner, Some(FilletCorner::Smooth))
}

/// Reports a fillet neither construction built.
///
/// The exact builders' refusal is the one callers already read, so its code
/// stays the outward one and its text comes first, with the rolling ball's
/// cause after it. A radius limit is the exception, because it is the refusal
/// the user can act on: one the rolling ball proved is reported as such, and
/// one only the exact builders proved keeps its code, and either goes last so
/// that a terminal `; radius limit` stays terminal.
fn both_refusals(mut exact: Refusal, mut numerical: Refusal) -> Refusal {
    const TOO_LARGE: &str = "fillet.radius-too-large";
    if numerical.code == TOO_LARGE {
        numerical.detail = format!(
            "exact path refused {}: {}; numerical continuation refused {}: {}",
            exact.code, exact.detail, numerical.code, numerical.detail
        );
        return numerical;
    }
    exact.detail = if exact.code == TOO_LARGE {
        format!(
            "numerical continuation refused {}: {}; exact path refused {}: {}",
            numerical.code, numerical.detail, exact.code, exact.detail
        )
    } else {
        format!(
            "{}; numerical continuation refused {}: {}",
            exact.detail, numerical.code, numerical.detail
        )
    };
    exact
}

/// Evaluates the matching exact prefix, then each fillet.
///
/// A fillet the exact builders finish stays exact and keeps its continuation
/// checkpoint. Any refusal of theirs, a spent allowance included, returns what
/// it charged and hands the fillet to the rolling-ball construction, which
/// keeps the request allowance it found. A numerical result never enters an
/// exact checkpoint.
pub fn evaluate_incremental_outcome(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<EvaluationOutcome, BuildError> {
    let (state, start) = resume_state(recipe, prefix, dec);
    let mut context = Context::default();
    let state = evaluate_impl(
        recipe,
        files,
        None,
        None,
        state,
        start,
        Some(&mut context),
        dec,
    )?;
    if let Some(result) = context.bodies.remove(&recipe.result) {
        return Ok(EvaluationOutcome::Numerical(result));
    }
    if !context.bodies.is_empty() {
        return Err(BuildError::Geometry {
            step: recipe.result.clone(),
            refusal: Refusal::new(
                "features.numerical-mixed-checkpoint",
                "an exact result with numerical intermediate branches cannot retain an exact-only checkpoint",
            ),
        });
    }
    Ok(EvaluationOutcome::Exact(checkpoint(state, recipe, dec)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_exact_refusal_stays_outward_unless_a_radius_limit_is_proved() {
        let limit = "fillet.radius-too-large";
        let refusal = |code: &'static str, detail: &str| Refusal::new(code, detail);
        // Neither proved a limit: the exact code, its text first.
        let both = both_refusals(
            refusal("shapes.rolling-envelope", "no envelope"),
            refusal("fillet.numerical-spine", "no walk"),
        );
        assert_eq!(both.code, "shapes.rolling-envelope");
        assert_eq!(
            both.detail,
            "no envelope; numerical continuation refused fillet.numerical-spine: no walk"
        );
        // The rolling ball proved the limit: its code, the limit terminal.
        let rolled = both_refusals(
            refusal("fillet.unsupported-topology", "needs one closed ring"),
            refusal(limit, "reaches face #57; radius limit 4 mm"),
        );
        assert_eq!(rolled.code, limit);
        assert!(rolled.detail.starts_with(
            "exact path refused fillet.unsupported-topology: needs one closed ring; numerical continuation refused"
        ));
        assert!(rolled.detail.ends_with("; radius limit 4 mm"));
        // Only the exact builders proved one: their code, the limit terminal.
        let exact = both_refusals(
            refusal(limit, "contact; radius limit 5 mm"),
            refusal("fillet.numerical-chain", "no chain"),
        );
        assert_eq!(exact.code, limit);
        assert!(
            exact
                .detail
                .starts_with("numerical continuation refused fillet.numerical-chain: no chain;")
        );
        assert!(exact.detail.ends_with("; radius limit 5 mm"));
        assert_eq!(
            crate::fillet::radius_limit(&exact),
            Some(Rational::from_i64(5))
        );
        assert_eq!(
            crate::fillet::radius_limit(&rolled),
            Some(Rational::from_i64(4))
        );
    }

    #[test]
    fn naming_the_default_overflow_or_corner_changes_nothing() {
        let roll = FilletOverflow::RollOver;
        let keep = FilletOverflow::KeepNeighbours(Vec::new());
        for overflow in [None, Some(&roll)] {
            for corner in [None, Some(FilletCorner::Cliff)] {
                assert!(rolls(overflow, corner), "{overflow:?} {corner:?}");
            }
            assert!(!rolls(overflow, Some(FilletCorner::Smooth)), "{overflow:?}");
        }
        for corner in [None, Some(FilletCorner::Cliff), Some(FilletCorner::Smooth)] {
            assert!(!rolls(Some(&keep), corner), "{corner:?}");
        }
    }
}
