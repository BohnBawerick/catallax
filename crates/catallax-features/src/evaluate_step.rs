//! One B-rep recipe step and its effects on the build.
//!
//! [`run`] builds one step from the values and records that earlier steps
//! published. Its [`Outcome`] holds the value the step made, what it appends
//! to the records carried from step to step, and which records it read. A
//! refusal comes back as the step's outward refusal. The step never changes a
//! published record in place: the evaluator publishes each outcome in step
//! order and keeps the checkpoints.
//!
//! [`Effects::of`] declares, before a step runs, which records it may read or
//! change. The work-span profiler orders steps by it, and an outcome that
//! reaches past the declaration fails a debug assertion in the evaluator.
//! Recipe input IDs name the step values a step reads; the records are what
//! those IDs cannot name. Reads of the model are tracked by the model itself
//! (`Model::watch`), and a step's model additions can be replayed onto a
//! grown model (`Model::replay`); reuse of a step outcome builds on both.

use super::*;
use crate::{chamfer::ChamferDefinition, fillet::FilletDefinition};
use catallax_recipe::{Origin, Step};
use std::cell::Cell;

/// A record the build carries from step to step, besides the model and the
/// value each step makes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Record {
    /// Generating definitions of fillets and their placed copies.
    Fillets,
    /// Faces that fillets and drafts removed.
    Consumed,
    /// Generating definitions of chamfers and their placed copies.
    Chamfers,
    /// Field-only steps carried as their region envelope.
    StandIns,
}

/// A set of records.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub(super) struct Records(u8);

impl Records {
    const fn of(records: &[Record]) -> Self {
        let mut set = 0;
        let mut i = 0;
        while i < records.len() {
            set |= 1 << records[i] as u8;
            i += 1;
        }
        Self(set)
    }
    fn with(self, record: Record) -> Self {
        Self(self.0 | Self::of(&[record]).0)
    }
    fn covers(self, other: Self) -> bool {
        other.0 & !self.0 == 0
    }
    #[cfg(any(test, feature = "work-span"))]
    fn meets(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

/// The feature registries: the records of feature history that the build
/// keeps in step order and that the profiler treats as one shared output.
#[cfg(any(test, feature = "work-span"))]
const FEATURE_REGISTRIES: Records =
    Records::of(&[Record::Fillets, Record::Consumed, Record::Chamfers]);

/// The records a step reads and the records it appends to.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub(super) struct Effects {
    pub(super) reads: Records,
    pub(super) writes: Records,
}

impl Effects {
    /// The records `op` may read or append to, declared before it runs.
    pub(super) fn of(op: &Operation) -> Self {
        use Record::{Chamfers, Consumed, Fillets, StandIns};
        let (reads, writes): (&[Record], &[Record]) = match op {
            Operation::Draft { .. } => (&[], &[Consumed]),
            // A fillet selection can name the strips of earlier fillets.
            Operation::Fillet { .. } | Operation::FilletVariable { .. } => {
                (&[Fillets], &[Fillets, Consumed])
            }
            Operation::Chamfer { .. } => (&[], &[Chamfers]),
            // A placed copy carries the fillet and chamfer records of its source.
            Operation::Pattern { .. } | Operation::Mirror { .. } => {
                (&[Fillets, Chamfers], &[Fillets, Chamfers])
            }
            Operation::Lattice { .. } | Operation::Gyroid { .. } => (&[], &[StandIns]),
            _ => (&[], &[]),
        };
        Self {
            reads: Records::of(reads),
            writes: Records::of(writes),
        }
    }

    /// Whether the step reads or extends the feature registries, and so
    /// waits for the last step that changed them.
    #[cfg(any(test, feature = "work-span"))]
    pub(super) fn feature_registries(self) -> bool {
        Records(self.reads.0 | self.writes.0).meets(FEATURE_REGISTRIES)
    }

    /// Whether every read and write of `actual` is declared here.
    pub(super) fn covers(self, actual: Self) -> bool {
        self.reads.covers(actual.reads) && self.writes.covers(actual.writes)
    }
}

/// Read access to the published records. It notes each record read.
pub(super) struct Published<'a> {
    fillets: &'a [FilletDefinition],
    chamfers: &'a [ChamferDefinition],
    read: Cell<Records>,
}

impl<'a> Published<'a> {
    pub(super) fn new(fillets: &'a [FilletDefinition], chamfers: &'a [ChamferDefinition]) -> Self {
        Self {
            fillets,
            chamfers,
            read: Cell::default(),
        }
    }
    fn note(&self, record: Record) {
        self.read.set(self.read.get().with(record));
    }
    fn fillets(&self) -> &'a [FilletDefinition] {
        self.note(Record::Fillets);
        self.fillets
    }
    fn chamfers(&self) -> &'a [ChamferDefinition] {
        self.note(Record::Chamfers);
        self.chamfers
    }
}

/// What earlier steps published: their values by step ID, and the records.
pub(super) struct Inputs<'a> {
    pub(super) profiles: &'a OrdMap<StepId, Profile>,
    pub(super) paths: &'a OrdMap<StepId, crate::sweep::Path>,
    pub(super) bodies: &'a OrdMap<StepId, BodyValue>,
    pub(super) meshes: &'a OrdMap<StepId, crate::mesh_ops::MeshValue>,
    pub(super) records: Published<'a>,
}

/// The optional constructions of one evaluation.
pub(super) struct Modes<'a, 'c> {
    /// Certified offsets and the booleans that read them.
    pub(super) certification: Option<&'a mut crate::certified::CertificationContext<'c>>,
    /// Numerical sweeps along helices and guides.
    pub(super) numerical_sweep: Option<&'a mut crate::general_sweep::Context>,
    /// The numerical continuation of fillets, holes and warps.
    pub(super) numerical: Option<&'a mut numerical_fillet::Context>,
}

/// The value one step makes.
pub(super) enum Value {
    /// Reference data, or a path that only the numerical sweep reads.
    Nothing,
    Profile(Profile),
    Path(crate::sweep::Path),
    Body(BodyValue),
    /// A mesh step's body and its exact triangles.
    Mesh(BodyValue, crate::mesh_ops::MeshValue),
    /// An import's body, its repair record and the placed lumps of each
    /// source solid.
    Import(
        BodyValue,
        catallax_import::RepairRecord,
        Vec<crate::import_edit::ImportedSolid>,
    ),
    /// A body with tolerant entities. The numerical continuation holds it;
    /// it never enters an exact checkpoint.
    Numerical,
}

impl Value {
    fn body(&self) -> Option<BodyId> {
        match self {
            Self::Body(value) | Self::Mesh(value, _) | Self::Import(value, ..) => Some(value.body),
            Self::Nothing | Self::Profile(_) | Self::Path(_) | Self::Numerical => None,
        }
    }
}

/// What one step appends to the records, in the order it made them.
#[derive(Default)]
pub(super) struct Changes {
    pub(super) fillets: Vec<FilletDefinition>,
    pub(super) consumed: Vec<crate::ConsumedFace>,
    pub(super) chamfers: Vec<ChamferDefinition>,
    pub(super) stand_ins: Vec<StandIn>,
}

impl Changes {
    fn written(&self) -> Records {
        let mut written = Records::default();
        for (record, empty) in [
            (Record::Fillets, self.fillets.is_empty()),
            (Record::Consumed, self.consumed.is_empty()),
            (Record::Chamfers, self.chamfers.is_empty()),
            (Record::StandIns, self.stand_ins.is_empty()),
        ] {
            if !empty {
                written = written.with(record);
            }
        }
        written
    }
}

/// One step's value, its record changes and the records it touched.
pub(super) struct Outcome {
    pub(super) value: Value,
    pub(super) changes: Changes,
    pub(super) effects: Effects,
}

/// Builds one step. A refusal is the step's outward refusal: a fillet or
/// distance feature that found no contact reports an unknown radius limit.
pub(super) fn run(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
    modes: &mut Modes<'_, '_>,
    step: &Step,
    inputs: &Inputs<'_>,
    model: &mut Model,
    dec: &mut Decider,
) -> Result<Outcome, Refusal> {
    let span = catallax_numbers::OperationSpan::begin("brep", &step.id, step.op.name(), dec);
    let mut changes = Changes::default();
    let result = build(recipe, files, modes, step, inputs, &mut changes, model, dec);
    span.finish(
        dec,
        result.as_ref().err(),
        result
            .as_ref()
            .ok()
            .and_then(Value::body)
            .map(|id| &model.body(id).marks),
    );
    let value = result.map_err(|refusal| outward(&step.op, refusal))?;
    Ok(Outcome {
        value,
        effects: Effects {
            reads: inputs.records.read.get(),
            writes: changes.written(),
        },
        changes,
    })
}

fn outward(op: &Operation, refusal: Refusal) -> Refusal {
    let limit_feature = matches!(
        op,
        Operation::Fillet { .. }
            | Operation::FilletVariable { .. }
            | Operation::FilletG2 { .. }
            | Operation::Shell { .. }
            | Operation::Thicken { .. }
    );
    let contact_refusal = refusal.code.starts_with("fillet.")
        || refusal.code.starts_with("shapes.rolling-")
        || refusal.code.starts_with("distance.")
        || refusal.code.starts_with("offset.");
    if limit_feature && contact_refusal && refusal.undecided().is_none() {
        crate::refusal_limit::unknown(refusal)
    } else {
        refusal
    }
}

#[allow(clippy::too_many_arguments)]
fn build(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
    modes: &mut Modes<'_, '_>,
    step: &Step,
    inputs: &Inputs<'_>,
    changes: &mut Changes,
    model: &mut Model,
    dec: &mut Decider,
) -> Result<Value, Refusal> {
    match &step.op {
        Operation::Sketch { plane, loops } => {
            return Profile::new(plane, loops, dec).map(Value::Profile);
        }
        Operation::Path { start, segments } => {
            return if modes.numerical_sweep.is_some()
                && crate::general_sweep::path_needs_numerical(recipe, &step.id)
            {
                Ok(Value::Nothing)
            } else {
                crate::sweep::Path::new(start, segments).map(Value::Path)
            };
        }
        // Reference data, RECIPE.md section 4.19.1: no body, and no body
        // operation may read it, so the B-rep has nothing to build.
        Operation::ScanImport { .. } => return Ok(Value::Nothing),
        _ => {}
    }
    let Inputs {
        profiles,
        paths,
        bodies,
        meshes,
        records,
    } = inputs;
    let Modes {
        certification,
        numerical_sweep,
        numerical,
    } = modes;
    if let Some(context) = numerical.as_deref_mut() {
        // A hole in a body with tolerant entities has a construction
        // of its own. Nothing else that reads one does yet.
        if let Operation::Hole {
            target,
            at,
            diameter,
            depth,
            tip,
            head,
        } = &step.op
            && context.has(target)
        {
            context.hole(
                target,
                at,
                diameter.value(),
                depth,
                tip,
                head.as_ref(),
                &step.id,
                recipe,
                dec,
            )?;
            return Ok(Value::Numerical);
        }
        context.exact_inputs(step)?;
    }
    if let Some((body, radius)) = crate::stock::build(model, &step.op, &step.id, dec)? {
        return Ok(Value::Body(BodyValue { body, radius }));
    }
    Ok(Value::Body(match &step.op {
        Operation::Sketch { .. } => return Ok(Value::Nothing),
        Operation::ImportStep {
            source_solid,
            file,
            pose,
            decisions,
        }
        | Operation::ImportIges {
            source_solid,
            file,
            pose,
            decisions,
        } => {
            let format = if matches!(step.op, Operation::ImportIges { .. }) {
                catallax_import::SourceFormat::Iges
            } else {
                catallax_import::SourceFormat::Step
            };
            let decisions = crate::import_edit::import_decisions(&file.sha256, decisions)?;
            let mut imported = catallax_import::apply_source(
                model,
                format,
                files[&file.sha256],
                &file.name,
                &Default::default(),
                &decisions,
                Some(&step.id),
            )
            .map_err(import_refusal)?;
            if let Some(entity) = source_solid {
                imported.bodies.retain(|(key, _)| key == entity);
                imported.left_open.retain(|key| key == entity);
                if imported.bodies.is_empty() && imported.left_open.is_empty() {
                    return Err(Refusal::new(
                        "import.source-solid",
                        format!("{entity} is not a source solid in {}", file.sha256),
                    ));
                }
            }
            if !imported.left_open.is_empty() {
                return Err(Refusal::new(
                    "import.left-open",
                    format!(
                        "{} not imported until the gap is resolved. The original file was never changed",
                        imported.left_open.join(", ")
                    ),
                ));
            }
            let raw = crate::import_edit::collect_solids(model, &imported.bodies, dec)?;
            let body = if pose == &catallax_recipe::Pose::identity() {
                raw
            } else {
                let placement = crate::profile::placement(pose, dec)?;
                catallax_topology::transformed_body(
                    model,
                    raw,
                    &Isometry::from(&placement),
                    Clone::clone,
                )
            };
            let bounds = catallax_measure::bounding_box(model, body, &Default::default(), dec)?
                .value
                .ok_or_else(|| {
                    Refusal::new("import.not-solid", "the imported body has no bounds")
                })?;
            let radius = (0..3).fold(Rational::zero(), |sum, i| {
                sum + bounds.min[i].abs().max(bounds.max[i].abs())
            });
            let mut first_lump = 0;
            let membership = imported
                .bodies
                .iter()
                .map(|(entity, raw)| {
                    let count = model.body(*raw).lumps.len();
                    let lumps = model.body(body).lumps[first_lump..first_lump + count].to_vec();
                    first_lump += count;
                    crate::import_edit::ImportedSolid {
                        entity: entity.clone(),
                        lumps,
                    }
                })
                .collect();
            return Ok(Value::Import(
                BodyValue { body, radius },
                imported.record,
                membership,
            ));
        }
        Operation::ImportMesh { file, pose, scale } => {
            let mesh =
                catallax_import::mesh::ClosedMesh::read(files[&file.sha256], scale.value(), dec)?;
            let placement = crate::profile::placement(pose, dec)?;
            let placement = Isometry::from(&placement);
            let (body, radius) =
                crate::mesh_import::build(model, &mesh, &placement, &step.id, dec)?;
            return Ok(Value::Mesh(
                BodyValue { body, radius },
                crate::mesh_ops::MeshValue::imported(mesh, placement, &step.id),
            ));
        }
        Operation::MeshBoolean { kind, target, tool } => {
            let (body, radius, mesh) =
                crate::mesh_ops::boolean(model, kind, &meshes[target], &meshes[tool], dec)?;
            return Ok(Value::Mesh(BodyValue { body, radius }, mesh));
        }
        Operation::MeshDecimate { body, triangles } => {
            let (body, radius, mesh) =
                crate::mesh_ops::decimate(model, &meshes[body], *triangles, dec)?;
            return Ok(Value::Mesh(BodyValue { body, radius }, mesh));
        }
        Operation::Extrude {
            profile,
            distance,
            direction,
        } => {
            let profile = &profiles[profile];
            let direction = direction.as_ref().map(vector);
            let body = crate::extrude::extrude_with_decider(
                model,
                profile,
                &step.id,
                distance.value(),
                direction.as_ref(),
                dec,
            )?;
            BodyValue {
                body,
                radius: profile_radius(profile) + distance.value().abs(),
            }
        }
        Operation::ExtrudeTapered {
            profile,
            distance,
            angle,
        } => {
            let (body, radius) = crate::tapered_extrude::with_radius(
                model,
                &profiles[profile],
                &step.id,
                distance.value(),
                angle.value(),
                dec,
            )?;
            BodyValue { body, radius }
        }
        Operation::ExtrudeOnFace {
            profile,
            body,
            face,
            normal_direction,
            x_direction,
            offset,
            distance,
            placement,
        } => {
            let profile = &profiles[profile];
            crate::check_face_profile(&profile.pose, dec)?;
            let support = &bodies[body];
            let planes = model
                .faces_of_body(support.body)
                .into_iter()
                .filter(|id| {
                    model
                        .face_origins(*id)
                        .iter()
                        .any(|origin| face.matches(origin))
                })
                .map(|id| {
                    let native = model.face(id);
                    let catallax_shapes::Surface::Plane(plane) =
                        model.geom().surface(native.surface)
                    else {
                        return Err(Refusal::new(
                            "extrude.nonplanar-face",
                            "the sketch support is not planar",
                        ));
                    };
                    Ok(plane.clone())
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            let map = crate::face_placement(
                &planes,
                &vector(normal_direction),
                &vector(x_direction),
                offset.value(),
                placement,
                dec,
            )?;
            let map = Isometry::from(&map);
            let local = crate::extrude::extrude_with_decider(
                model,
                profile,
                &step.id,
                distance.value(),
                None,
                dec,
            )?;
            let body = catallax_topology::transformed_body(model, local, &map, Clone::clone);
            BodyValue {
                body,
                radius: support.radius.clone()
                    + profile_radius(profile)
                    + offset.value().abs()
                    + distance.value().abs()
                    + placement
                        .translation
                        .0
                        .iter()
                        .fold(Rational::zero(), |sum, n| sum + n.value().abs()),
            }
        }
        Operation::Loft { sections, shape } => {
            let (body, radius) =
                crate::loft::build(model, recipe, sections, *shape, &step.id, profiles, dec)?;
            BodyValue { body, radius }
        }
        Operation::Sweep {
            profile,
            path,
            orientation,
            guide,
            twist,
        } => {
            let helical = recipe
                .steps
                .iter()
                .find(|s| &s.id == path)
                .is_some_and(|s| {
                    matches!(&s.op,Operation::Path {segments,..}
                    if segments.iter().any(|s|matches!(s,catallax_recipe::PathSegment::Helix {..})))
                });
            if helical || guide.is_some() || twist.as_ref().is_some_and(|n| !n.value().is_zero()) {
                let context = numerical_sweep.as_deref_mut().ok_or_else(|| {
                    Refusal::new(
                        "sweep.requires-numerical-evaluation",
                        "helical sweeps require numerical admission under a stated tolerance",
                    )
                })?;
                let (body, radius) = context.sweep(model, recipe, step, dec)?;
                BodyValue { body, radius }
            } else {
                let profile = &profiles[profile];
                let path = &paths[path];
                let body = path.build(model, profile, &step.id, *orientation, dec)?;
                BodyValue {
                    body,
                    radius: profile_radius(profile) + path.extent(),
                }
            }
        }
        Operation::Revolve {
            profile,
            axis,
            angle,
        } => {
            let profile = &profiles[profile];
            let body = revolve(model, profile, &step.id, axis, angle.value(), dec)?;
            let origin = RPoint3::new(
                axis.point.0[0].value().clone(),
                axis.point.0[1].value().clone(),
                axis.point.0[2].value().clone(),
            );
            BodyValue {
                body,
                radius: profile_radius(profile) + Rational::from_i64(2) * norm_upper(&origin),
            }
        }
        Operation::ExtrudeThrough {
            target,
            profile,
            direction,
            side,
        } => {
            let target = &bodies[target];
            let body = through(
                model,
                &profiles[profile],
                &step.id,
                target,
                direction.as_ref(),
                *side,
                dec,
            )?;
            BodyValue {
                body,
                radius: target.radius.clone(),
            }
        }
        Operation::ExtrudeFrom {
            target,
            profile,
            direction,
            start,
            distance,
        } => {
            profiles[profile].require_circular("extrude-from.ellipse-profile")?;
            let stock = &bodies[target];
            let (body, radius) = crate::end_condition::face::start(
                model,
                recipe,
                (target, stock.body, &stock.radius),
                (profile, &profiles[profile]),
                direction.as_ref(),
                start,
                distance.value(),
                &step.id,
                dec,
            )?;
            BodyValue { body, radius }
        }
        Operation::ExtrudeTo {
            target,
            profile,
            direction,
            ..
        }
        | Operation::ExtrudeNext {
            target,
            profile,
            direction,
        } => {
            use crate::end_condition::{Capped, curved::AuthoredCap, face, toroidal::AnnularCap};
            let limit = match &step.op {
                Operation::ExtrudeTo { limit, .. } => Some(limit),
                _ => None,
            };
            let profile_value = &profiles[profile];
            let stock = &bodies[target];
            profile_value.require_circular(if limit.is_some() {
                "extrude-to.ellipse-profile"
            } else {
                "extrude-next.ellipse-profile"
            })?;
            let result = (|| -> Result<BodyValue, Refusal> {
                match recipe
                    .step(target)
                    .map(|s| &s.op)
                    .filter(|_| !profile_value.is_empty())
                {
                    Some(Operation::Torus { .. }) => {
                        let cap = AnnularCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            limit,
                            dec,
                        )?;
                        Ok(BodyValue {
                            body: cap.build(model, stock.body, profile_value, &step.id, dec)?,
                            radius: cap.radius().clone(),
                        })
                    }
                    Some(
                        Operation::Sphere { .. }
                        | Operation::Cylinder { .. }
                        | Operation::Cone { .. },
                    ) => {
                        let cap = AuthoredCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            limit,
                            dec,
                        )?;
                        Ok(BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius().clone(),
                        })
                    }
                    _ => {
                        let cap = Capped::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            limit,
                            dec,
                        )?;
                        Ok(BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius,
                        })
                    }
                }
            })();
            match result {
                Err(error)
                    if recipe.format == catallax_recipe::FORMAT_V6 && face::may_retry(&error) =>
                {
                    let selected = limit
                        .map(|reference| face::select(model, stock.body, reference, dec))
                        .transpose()?;
                    let (body, radius) = face::limit(
                        model,
                        recipe,
                        (target, stock.body, &stock.radius),
                        (profile, profile_value),
                        direction.as_ref(),
                        selected,
                        &step.id,
                        dec,
                    )?;
                    BodyValue { body, radius }
                }
                result => result?,
            }
        }
        Operation::Hole {
            target,
            at,
            diameter,
            depth,
            tip,
            head,
        } => {
            crate::hole::validate(diameter.value(), depth, tip, head.as_ref())?;
            let target = &bodies[target];
            let mut body = target.body;
            let tools = at
                .iter()
                .map(|placement| {
                    crate::hole::HoleTool::new(
                        placement,
                        diameter.value(),
                        depth,
                        tip,
                        head.as_ref(),
                        &target.radius,
                        dec,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            if model.body(body).kind != BodyKind::Empty {
                let mut disjoint = true;
                for i in 0..tools.len() {
                    for other in &tools[..i] {
                        if !tools[i].disjoint(other, dec)? {
                            disjoint = false;
                            break;
                        }
                    }
                    if !disjoint {
                        break;
                    }
                }
                if disjoint && !tools.is_empty() {
                    let mut lumps = Vec::new();
                    for (tool, placement) in tools.iter().zip(at) {
                        let tool = tool.build(model, &step.id, &placement.id);
                        lumps.extend_from_slice(&model.body(tool).lumps);
                    }
                    let tool = model.add_body(catallax_topology::Body {
                        kind: BodyKind::Solid,
                        lumps,
                        wires: Vec::new(),
                        marks: Default::default(),
                    });
                    body = boolean(model, body, tool, BoolOp::Subtract, dec)
                        .map_err(|e| e.refusal)?
                        .body;
                } else {
                    for (tool, placement) in tools.iter().zip(at) {
                        let tool = tool.build(model, &step.id, &placement.id);
                        body = boolean(model, body, tool, BoolOp::Subtract, dec)
                            .map_err(|e| e.refusal)?
                            .body;
                    }
                }
            }
            BodyValue {
                body,
                radius: target.radius.clone(),
            }
        }
        Operation::Boolean { kind, target, tool } => {
            let a = &bodies[target];
            let b = &bodies[tool];
            require_material(model, &[a.body, b.body])?;
            if let Some(context) = certification.as_deref_mut()
                && (context.has_certificate(target) || context.has_certificate(tool))
            {
                let (body, radius) = context.boolean(
                    model,
                    &step.id,
                    kind,
                    [(target, a.body), (tool, b.body)],
                    dec,
                )?;
                return Ok(Value::Body(BodyValue { body, radius }));
            }
            let (op, radius) = match kind {
                BooleanKind::Union => (BoolOp::Unite, a.radius.clone().max(b.radius.clone())),
                BooleanKind::Subtract => (BoolOp::Subtract, a.radius.clone()),
                BooleanKind::Intersect => {
                    (BoolOp::Intersect, a.radius.clone().min(b.radius.clone()))
                }
            };
            // Recipes allow repeated operands. Apply exact set identities
            // before the two-distinct-operands intersection graph.
            let body = if a.body == b.body {
                if *kind == BooleanKind::Subtract {
                    // An empty result has no fitted face.
                    model.add_body(catallax_topology::Body {
                        kind: BodyKind::Empty,
                        lumps: Vec::new(),
                        wires: Vec::new(),
                        marks: catallax_topology::without_fit_mark(&model.body(a.body).marks),
                    })
                } else {
                    a.body
                }
            } else {
                // Constructing the part merges a sub-tolerance gap between
                // the operands once, here, rather than leaving a sliver.
                catallax_boolean::boolean_merging_gaps(model, a.body, b.body, op, dec)
                    .map_err(|e| e.refusal)?
                    .body
            };
            BodyValue { body, radius }
        }
        Operation::Draft {
            body,
            faces,
            pull,
            neutral,
            angle,
            keep_regions,
        } => {
            let definition = crate::draft::DraftDefinition::new(
                vector(pull).to_real(),
                crate::draft::neutral_plane(neutral),
                angle.value(),
                dec,
            )?;
            let target = &bodies[body];
            let (body, history, radius) = crate::draft_build::build(
                model,
                target.body,
                faces,
                &definition,
                keep_regions.as_deref(),
                &step.id,
                dec,
            )?;
            changes.consumed.extend(history);
            BodyValue { body, radius }
        }
        Operation::Fillet {
            body,
            edges,
            radius,
            overflow,
            corner,
        } => {
            let target = &bodies[body];
            let mut exact = |model: &mut Model, dec: &mut Decider| {
                if let Some(catallax_recipe::FilletOverflow::KeepNeighbours(limits)) = overflow {
                    crate::fillet_build::fillet_setback::build(
                        model,
                        target.body,
                        edges,
                        &radius.value().clone().into(),
                        limits,
                        &step.id,
                        dec,
                    )
                } else {
                    crate::fillet_build::build(
                        model,
                        target.body,
                        edges,
                        &radius.value().clone().into(),
                        recipe.format != catallax_recipe::FORMAT,
                        corner.unwrap_or(catallax_recipe::FilletCorner::Cliff),
                        &step.id,
                        records.fillets(),
                        &mut changes.consumed,
                        dec,
                    )
                }
            };
            let (body, definitions) = match numerical.as_deref_mut() {
                None => exact(model, dec)?,
                // The exact builders are an optional proof of the
                // same fillet: a body they finish is the exact result,
                // charged as it always was. A refusal of theirs is no
                // verdict on the fillet and returns what it charged,
                // so the rolling-ball construction keeps the whole
                // allowance it found. The attempt runs under the cap
                // of its tier, which the builders set once they have
                // resolved their selection (`predictor.rs`): reaching
                // it is a refusal of the attempt like any other.
                Some(context) => {
                    let work = catallax_numbers::current_operation_work;
                    // Everything counted for the request so far:
                    // refunded work only grows, so this never falls.
                    let counted = || {
                        let work = work();
                        work.used + work.refunded + work.uncharged
                    };
                    let (used, refunded, before) = (work().used, work().refunded, counted());
                    match catallax_numbers::capped_attempt(|| exact(model, dec)) {
                        Ok(result) => {
                            crate::predictor::record(
                                &step.id,
                                crate::predictor::Path::Exact,
                                work().used - used,
                                counted() - before,
                            );
                            result
                        }
                        Err(refusal) => {
                            let attempt = work().refunded - refunded;
                            let fitted = context.fillet(
                                model,
                                target,
                                edges,
                                radius.value(),
                                overflow.as_ref(),
                                *corner,
                                records.fillets(),
                                &step.id,
                                recipe,
                                refusal,
                                dec,
                            );
                            crate::predictor::record(
                                &step.id,
                                if fitted.is_ok() {
                                    crate::predictor::Path::Fitted
                                } else {
                                    crate::predictor::Path::Refused
                                },
                                attempt,
                                counted() - before,
                            );
                            fitted?;
                            return Ok(Value::Numerical);
                        }
                    }
                }
            };
            changes.fillets.extend(definitions);
            BodyValue {
                body,
                radius: &target.radius + &(radius.value() * &Rational::from_i64(2)),
            }
        }
        Operation::FilletVariable { body, edges } => {
            let target = &bodies[body];
            let mut maximum = Rational::zero();
            for law in edges {
                let catallax_recipe::RadiusLaw::Linear { start, end } = &law.radius else {
                    return Err(Refusal::new(
                        "fillet.unsupported-variable-contact",
                        "the radius law has no construction definition",
                    ));
                };
                maximum = maximum.max(start.value().clone()).max(end.value().clone());
            }
            let (body, definitions) = crate::fillet_build::build_variable(
                model,
                target.body,
                edges,
                &step.id,
                records.fillets(),
                &mut changes.consumed,
                dec,
            )?;
            changes.fillets.extend(definitions);
            BodyValue {
                body,
                radius: &target.radius + &(&maximum * &Rational::from_i64(2)),
            }
        }
        Operation::FilletG2 {
            body,
            edges,
            radius,
            fullness,
        } => {
            let target = &bodies[body];
            let fullness = fullness
                .as_ref()
                .map_or_else(Rational::one, |f| f.value().clone());
            let body = crate::fillet_build::build_g2(
                model,
                target.body,
                edges,
                &radius.value().clone().into(),
                &fullness.into(),
                &step.id,
                dec,
            )?;
            BodyValue {
                body,
                radius: &target.radius + &(radius.value() * &Rational::from_i64(2)),
            }
        }
        Operation::Chamfer {
            body,
            edges,
            distance,
            other_distance,
        } => {
            let target = &bodies[body];
            let other = other_distance.as_ref().unwrap_or(distance);
            let distances = [
                distance.value().clone().into(),
                other.value().clone().into(),
            ];
            let (body, definitions) = crate::chamfer_build::build(
                model,
                target.body,
                edges,
                &distances,
                other_distance.is_some(),
                &step.id,
                dec,
            )?;
            changes.chamfers.extend(definitions);
            BodyValue {
                body,
                radius: &(&target.radius + distance.value()) + other.value(),
            }
        }
        Operation::Pattern { body, layout } => {
            // Validate even count-one and empty patterns before set identities.
            crate::pattern_placement(layout, 0, dec)?;
            let source = &bodies[body];
            if layout.count() > 1 {
                require_material(model, &[source.body])?;
            }
            let mut body = source.body;
            let mut radius = source.radius.clone();
            let copied = if layout.count() > 1 {
                CopiedRecords::on(model, source.body, records)
            } else {
                CopiedRecords::default()
            };
            let source_bounds = if layout.count() > 1 {
                catallax_measure::bounding_box(model, source.body, &Default::default(), dec)?.value
            } else {
                None
            };
            let mut placed_bounds = source_bounds.iter().cloned().collect::<Vec<_>>();
            for instance in 1..layout.count() {
                let map = crate::pattern_placement(layout, instance, dec)?;
                let shift = [&map.translation.x, &map.translation.y, &map.translation.z]
                    .into_iter()
                    .try_fold(Rational::zero(), |sum, v| {
                        let b = dec.enclose_at(v, 128)?;
                        Ok::<_, Refusal>(sum + b.lo().abs().max(b.hi().abs()))
                    })?;
                radius = radius.max(&source.radius + &shift);
                let copy =
                    catallax_topology::transformed_body(model, source.body, &map, |origin| {
                        crate::copy_origin(&step.id, instance, origin)
                    });
                copied.place(records, &map, &step.id, instance, changes)?;
                let copy_bounds = source_bounds
                    .as_ref()
                    .map(|bounds| transformed_bounds(bounds, &map, dec))
                    .transpose()?;
                body = if copy_bounds
                    .as_ref()
                    .is_some_and(|copy| placed_bounds.iter().all(|prior| separated(prior, copy)))
                {
                    unite_separated(model, body, copy)
                } else {
                    boolean(model, body, copy, BoolOp::Unite, dec)
                        .map_err(|e| e.refusal)?
                        .body
                };
                if let Some(bounds) = copy_bounds {
                    placed_bounds.push(bounds);
                }
            }
            BodyValue { body, radius }
        }
        Operation::Mirror {
            body,
            plane,
            keep_original,
        } => {
            let map = crate::mirror_placement(plane)?;
            let source = &bodies[body];
            if *keep_original {
                require_material(model, &[source.body])?;
            }
            let copied = CopiedRecords::on(model, source.body, records);
            let source_bounds = if *keep_original {
                catallax_measure::bounding_box(model, source.body, &Default::default(), dec)?.value
            } else {
                None
            };
            let shift = [&map.translation.x, &map.translation.y, &map.translation.z]
                .into_iter()
                .try_fold(Rational::zero(), |sum, v| {
                    let b = dec.enclose_at(v, 128)?;
                    Ok::<_, Refusal>(sum + b.lo().abs().max(b.hi().abs()))
                })?;
            let copy = catallax_topology::transformed_body(model, source.body, &map, |origin| {
                crate::copy_origin(&step.id, 1, origin)
            });
            copied.place(records, &map, &step.id, 1, changes)?;
            let body = if *keep_original {
                let copy_bounds = source_bounds
                    .as_ref()
                    .map(|bounds| transformed_bounds(bounds, &map, dec))
                    .transpose()?;
                if source_bounds
                    .as_ref()
                    .zip(copy_bounds.as_ref())
                    .is_some_and(|(a, b)| separated(a, b))
                {
                    unite_separated(model, source.body, copy)
                } else {
                    boolean(model, source.body, copy, BoolOp::Unite, dec)
                        .map_err(|e| e.refusal)?
                        .body
                }
            } else {
                copy
            };
            BodyValue {
                body,
                radius: &source.radius + &shift,
            }
        }
        Operation::Lattice {
            region,
            cell,
            size,
            radius,
            ..
        } => {
            // The lattice is field-only. The B-rep carries the PRD section 7
            // stand-in: the region is the outer envelope, and the field kernel
            // adds bodies matching mass, centre of gravity and inertia.
            crate::lattice::validate(cell, size.value(), radius)?;
            let source = &bodies[region];
            changes.stand_ins.push(StandIn {
                step: step.id.clone(),
                region: region.clone(),
            });
            BodyValue {
                body: source.body,
                radius: source.radius.clone(),
            }
        }
        Operation::Gyroid {
            region,
            period,
            level,
            ..
        } => {
            // The gyroid is field-only too, carried by the same region envelope.
            crate::lattice::validate_gyroid(period.value(), level)?;
            let source = &bodies[region];
            changes.stand_ins.push(StandIn {
                step: step.id.clone(),
                region: region.clone(),
            });
            BodyValue {
                body: source.body,
                radius: source.radius.clone(),
            }
        }
        Operation::Offset { body, distance } if certification.is_some() => {
            let input = &bodies[body];
            let (body, radius) = certification
                .as_deref_mut()
                .expect("certified branch")
                .offset(
                    model,
                    profiles,
                    &step.id,
                    (body, input.body, &input.radius),
                    distance.value(),
                    dec,
                )?;
            BodyValue { body, radius }
        }
        Operation::Twist { body, .. } | Operation::Bend { body, .. } => {
            let warp = crate::deform::Warp::of(&step.op)?.expect("a bend or twist step");
            let context = numerical.as_deref_mut().ok_or_else(|| {
                Refusal::new(
                    "deform.requires-numerical-evaluation",
                    "bend and twist build only as numerical geometry; build through the numerical outcome",
                )
            })?;
            context.deform(model, &bodies[body], &warp, &step.id, recipe, dec)?;
            return Ok(Value::Numerical);
        }
        Operation::Offset { body, distance } => {
            let input = &bodies[body];
            let result = crate::distance_build::offset(
                model,
                input.body,
                &distance.value().clone().into(),
                &step.id,
                dec,
            )?;
            BodyValue {
                body: crate::distance_build::inherit_marks(model, input.body, result),
                radius: &input.radius + &distance.value().abs(),
            }
        }
        Operation::Shell {
            body,
            open,
            thickness,
        } => {
            let input = &bodies[body];
            let result = crate::distance_build::shell(
                model,
                input.body,
                open,
                &thickness.value().clone().into(),
                &step.id,
                dec,
            )?;
            BodyValue {
                body: crate::distance_build::inherit_marks(model, input.body, result),
                radius: input.radius.clone(),
            }
        }
        Operation::Thicken {
            body,
            faces,
            thickness,
            side,
        } => {
            let input = &bodies[body];
            let result = crate::distance_build::thicken(
                model,
                input.body,
                faces,
                &thickness.value().clone().into(),
                *side,
                &step.id,
                recipe.format != catallax_recipe::FORMAT,
                dec,
            )?;
            BodyValue {
                body: crate::distance_build::inherit_marks(model, input.body, result),
                radius: &input.radius + &thickness.value().abs(),
            }
        }
        Operation::ShellGraded {
            body,
            open,
            thickness,
        } => {
            use crate::graded_wall::Thickness;
            let input = &bodies[body];
            let result = match Thickness::new(thickness)? {
                Thickness::Constant(thickness) => crate::distance_build::shell(
                    model,
                    input.body,
                    open,
                    &thickness.into(),
                    &step.id,
                    dec,
                )?,
                Thickness::Graded(grade) => crate::distance_build::graded_shell(
                    model, input.body, open, &grade, &step.id, dec,
                )?,
            };
            BodyValue {
                body: crate::distance_build::inherit_marks(model, input.body, result),
                radius: input.radius.clone(),
            }
        }
        Operation::ThickenGraded {
            body,
            faces,
            thickness,
            side,
        } => {
            use crate::graded_wall::Thickness;
            let input = &bodies[body];
            let (result, thickness) = match Thickness::new(thickness)? {
                Thickness::Constant(thickness) => (
                    crate::distance_build::thicken(
                        model,
                        input.body,
                        faces,
                        &thickness.clone().into(),
                        *side,
                        &step.id,
                        true,
                        dec,
                    )?,
                    thickness,
                ),
                Thickness::Graded(grade) => (
                    crate::distance_build::graded_thicken(
                        model, input.body, faces, &grade, *side, &step.id, dec,
                    )?,
                    grade.maximum().clone(),
                ),
            };
            BodyValue {
                body: crate::distance_build::inherit_marks(model, input.body, result),
                radius: &input.radius + &thickness,
            }
        }
        _ => {
            return Err(Refusal::new(
                "features.unsupported-operation",
                format!(
                    "B-rep recipe operation {} is not implemented",
                    step.op.name()
                ),
            ));
        }
    }))
}

/// The fillet and chamfer records on a source body, by index into the
/// published records. Pattern and mirror give every placed copy of the
/// body a placed copy of each, named for that copy.
#[derive(Default)]
struct CopiedRecords {
    fillets: Vec<usize>,
    chamfers: Vec<usize>,
}

impl CopiedRecords {
    fn on(model: &Model, body: BodyId, records: &Published<'_>) -> Self {
        Self {
            fillets: on_body(model, body, records.fillets()),
            chamfers: on_body(model, body, records.chamfers()),
        }
    }

    /// Appends the copies for `instance` of `step`: fillets, then chamfers.
    fn place(
        &self,
        records: &Published<'_>,
        map: &Isometry,
        step: &StepId,
        instance: u32,
        changes: &mut Changes,
    ) -> Result<(), Refusal> {
        place(
            records.fillets(),
            &self.fillets,
            map,
            step,
            instance,
            &mut changes.fillets,
        )?;
        place(
            records.chamfers(),
            &self.chamfers,
            map,
            step,
            instance,
            &mut changes.chamfers,
        )
    }
}

/// A feature definition that a placed copy carries.
trait Placed: Clone {
    fn origin(&self) -> &Origin;
    fn origin_mut(&mut self) -> &mut Origin;
    fn transform(&mut self, map: &Isometry) -> Result<(), Refusal>;
}

impl Placed for FilletDefinition {
    fn origin(&self) -> &Origin {
        &self.origin
    }
    fn origin_mut(&mut self) -> &mut Origin {
        &mut self.origin
    }
    fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        FilletDefinition::transform(self, map)
    }
}

impl Placed for ChamferDefinition {
    fn origin(&self) -> &Origin {
        &self.origin
    }
    fn origin_mut(&mut self) -> &mut Origin {
        &mut self.origin
    }
    fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        ChamferDefinition::transform(self, map)
    }
}

fn on_body<T: Placed>(model: &Model, body: BodyId, records: &[T]) -> Vec<usize> {
    if records.is_empty() {
        return Vec::new();
    }
    let faces = model.faces_of_body(body);
    records
        .iter()
        .enumerate()
        .filter_map(|(index, record)| {
            faces
                .iter()
                .any(|face| model.face_origins(*face).contains(record.origin()))
                .then_some(index)
        })
        .collect()
}

fn place<T: Placed>(
    records: &[T],
    indices: &[usize],
    map: &Isometry,
    step: &StepId,
    instance: u32,
    out: &mut Vec<T>,
) -> Result<(), Refusal> {
    for &index in indices {
        let mut copy = records[index].clone();
        copy.transform(map)?;
        *copy.origin_mut() = crate::copy_origin(step, instance, copy.origin());
        out.push(copy);
    }
    Ok(())
}

fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    RVector3::new(
        v.0[0].value().clone(),
        v.0[1].value().clone(),
        v.0[2].value().clone(),
    )
}

fn transformed_bounds(
    bounds: &BoundingBox,
    map: &Isometry,
    dec: &mut Decider,
) -> Result<BoundingBox, Refusal> {
    let local = std::array::from_fn::<_, 3, _>(|i| {
        Ball::from_rational(&bounds.min[i], 128).hull(&Ball::from_rational(&bounds.max[i], 128))
    });
    let mut enclose = |value: &Real| {
        if let Some(value) = value.rational_ref() {
            Ok(Ball::from_rational(value, 128))
        } else {
            dec.enclose_budgeted_at(value, 128)
        }
    };
    let translation = [&map.translation.x, &map.translation.y, &map.translation.z];
    let mut min = std::array::from_fn(|_| Rational::zero());
    let mut max = min.clone();
    for i in 0..3 {
        let mut coordinate = enclose(translation[i])?;
        for (j, local) in local.iter().enumerate() {
            coordinate = &coordinate + &(enclose(&map.matrix()[i][j])? * local.clone());
        }
        min[i] = coordinate.lo();
        max[i] = coordinate.hi();
    }
    Ok(BoundingBox { min, max })
}

fn separated(a: &BoundingBox, b: &BoundingBox) -> bool {
    (0..3).any(|i| a.max[i] < b.min[i] || b.max[i] < a.min[i])
}

fn require_material(model: &Model, bodies: &[BodyId]) -> Result<(), Refusal> {
    if bodies
        .iter()
        .any(|body| model.body(*body).kind == BodyKind::Surface)
    {
        return Err(Refusal::new(
            "boolean.not-solid",
            "a surface body has no solid material",
        ));
    }
    Ok(())
}

// Strictly separated enclosing boxes prove that regularised union is just the
// disjoint collection of the existing immutable lumps. No support intersection
// or replacement geometry is needed, and contact never takes this path.
fn unite_separated(model: &mut Model, a: BodyId, b: BodyId) -> BodyId {
    let mut result = model.body(a).clone();
    result.lumps.extend_from_slice(&model.body(b).lumps);
    result.marks.extend(&model.body(b).marks);
    result.kind = if result.lumps.is_empty() {
        BodyKind::Empty
    } else {
        BodyKind::Solid
    };
    model.add_body(result)
}

fn through(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    target: &BodyValue,
    direction: Option<&catallax_recipe::Vector3>,
    side: Option<ThroughSide>,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let direction = direction.map(vector).unwrap_or_else(|| normal.clone());
    if direction.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion direction must not be zero",
        ));
    }
    let dot = direction.dot(&normal);
    if dot.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "extrusion direction lies in the sketch plane",
        ));
    }
    if model.body(target.body).kind == BodyKind::Empty || profile.is_empty() {
        return Ok(target.body);
    }
    let origin = crate::profile::rational_point(&profile.pose, &RPoint3::origin());
    let normal_l1 = normal.x.abs() + normal.y.abs() + normal.z.abs();
    // For every target point p, |(p-o).n| <= (radius + |o|_1)|n|_1.
    // The strict extra unit leaves both construction caps outside the target.
    // Rescaling direction changes t inversely, not the resulting through-all set.
    let extent =
        ((&target.radius + &norm_upper(&origin)) * normal_l1) / dot.abs() + Rational::one();
    let (low, high, caps) = match side.unwrap_or(ThroughSide::Both) {
        ThroughSide::Both => (-&extent, extent.clone(), [None, None]),
        ThroughSide::Forward => (Rational::zero(), extent.clone(), [Some(Role::Start), None]),
        ThroughSide::Backward => (-&extent, Rational::zero(), [None, Some(Role::Start)]),
    };
    let mut shifted = profile.clone();
    shifted.pose.translation = shifted
        .pose
        .translation
        .plus(&direction.scaled(&low).to_real());
    let tool = crate::extrude::extrude_vector(
        model,
        &shifted,
        step,
        &direction.scaled(&(high - low)),
        caps,
        None,
        dec,
    )?;
    Ok(boolean(model, target.body, tool, BoolOp::Subtract, dec)
        .map_err(|e| e.refusal)?
        .body)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn effects(step: serde_json::Value) -> Effects {
        Effects::of(&serde_json::from_value::<Step>(step).unwrap().op)
    }

    // Pattern and mirror read and extend both records they copy. The
    // profiler orders exactly the steps that touch the feature registries.
    #[test]
    fn copies_declare_the_records_they_carry() {
        let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
        let end = json!({"step":"b","role":"end"});
        let wall = json!({"step":"b","role":{"side":{"segment":"wall"}}});
        let edge = json!({"between":[end, wall]});
        let copied = Records::of(&[Record::Fillets, Record::Chamfers]);
        for step in [
            json!({"id":"p","op":"pattern","body":"b","layout":{"linear":{"vector":["3","0","0"],"count":2}}}),
            json!({"id":"m","op":"mirror","body":"b","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":true}),
        ] {
            let declared = effects(step);
            assert_eq!(declared.reads, copied);
            assert_eq!(declared.writes, copied);
        }
        for (step, registries) in [
            (
                json!({"id":"r","op":"fillet","body":"b","radius":"2","edges":[edge]}),
                true,
            ),
            (
                json!({"id":"r","op":"fillet_variable","body":"b","edges":[{"edge":edge,"from":{"faces":[end,wall]},"radius":{"kind":"linear","start":"1","end":"2"}}]}),
                true,
            ),
            (
                json!({"id":"c","op":"chamfer","body":"b","distance":"2","edges":[edge]}),
                true,
            ),
            (
                json!({"id":"d","op":"draft","body":"b","faces":[{"face":end}],"pull":["0","0","1"],"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"5"}),
                true,
            ),
            (
                json!({"id":"l","op":"lattice","region":"b","cell":"cubic","size":"2","frame":pose,"radius":{"constant":"1/4"}}),
                false,
            ),
            (
                json!({"id":"u","op":"boolean","kind":"subtract","target":"a","tool":"b"}),
                false,
            ),
            (
                json!({"id":"b","op":"block","size":["1","1","1"],"pose":pose}),
                false,
            ),
        ] {
            assert_eq!(
                effects(step.clone()).feature_registries(),
                registries,
                "{step}"
            );
        }
    }
}
