//! Read the AP242 occurrence relationship used by the STEP assembly writer
//! and OpenCASCADE XDE. Parts stay in definition coordinates.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_step::assembly::{Assembly, Colour, Occurrence, Product, ProductContent};

use crate::{
    ImportRefusal, StepFile, StepPmi,
    part21::{Graph, Param, parse},
    pmi::{assembly_name as string, read_assembly_pmi},
    recognize::{number, read_graph, source_length_scale, source_placement},
};

/// A product's recognized source geometry and metadata. The original STEP
/// entity keys tie PMI to its faces and can be passed to the part builder.
#[derive(Clone, Debug)]
pub struct Part {
    /// Geometry in product-local coordinates.
    pub geometry: StepFile,
    /// Uniform surface colour inherited by occurrences of this part.
    pub colour: Option<Colour>,
    /// Semantic PMI owned by this product.
    pub pmi: StepPmi,
}

#[derive(Clone, Copy)]
struct Representation {
    shape: u64,
    representation: u64,
}

/// Reads one rooted assembly. Rigid placement matrices are recovered from
/// STEP frames and converted into millimetres. Unsupported shapes or layouts
/// refuse instead of dropping geometry or returning a flattened assembly.
pub fn read_step_assembly(
    bytes: &[u8],
    dec: &mut Decider,
) -> Result<Assembly<Part>, ImportRefusal> {
    let graph = parse(bytes)?;
    dec.charge_work(graph.len() as u64)
        .map_err(crate::recognize::source_error)?;
    let scale = source_length_scale(&graph, dec)?;

    let mut representations = BTreeMap::new();
    let mut shape_definitions = BTreeMap::new();
    let mut uses = BTreeMap::new();
    let mut contexts = BTreeMap::new();
    let mut styles = BTreeMap::<u64, Vec<u64>>::new();
    let mut leaf_relationships = BTreeMap::<u64, Vec<u64>>::new();
    for (&id, entry) in &graph {
        match entry.type_name() {
            "SHAPE_DEFINITION_REPRESENTATION" => {
                let fields = fields(&graph, id, "SHAPE_DEFINITION_REPRESENTATION")?;
                if matches!(fields.get(1), Some(Param::Ref(representation))
                    if graph.get(representation).is_some_and(|e|
                        e.part("SHAPE_REPRESENTATION_WITH_PARAMETERS").is_some()))
                {
                    continue;
                }
                let shape = reference(fields, 0)?;
                // A property's shape definition is a PMI carrier, such as a
                // datum target's parameters, never a product's geometry. The
                // PMI reader reports a broken one with its item.
                if graph.get(&shape).is_some_and(|e| {
                    e.part("PROPERTY_DEFINITION").is_some()
                        && e.part("PRODUCT_DEFINITION_SHAPE").is_none()
                }) {
                    continue;
                }
                let shape_fields = self::fields(&graph, shape, "PRODUCT_DEFINITION_SHAPE")?;
                let definition = reference(shape_fields, 2)?;
                if !graph
                    .get(&definition)
                    .is_some_and(|e| e.part("PRODUCT_DEFINITION").is_some())
                {
                    continue;
                }
                let representation = reference(fields, 1)?;
                if representations
                    .insert(
                        definition,
                        Representation {
                            shape,
                            representation,
                        },
                    )
                    .is_some()
                {
                    return Err(invalid(format!(
                        "product #{definition} has multiple representations"
                    )));
                }
            }
            "PRODUCT_DEFINITION_SHAPE" => {
                let fields = fields(&graph, id, "PRODUCT_DEFINITION_SHAPE")?;
                shape_definitions.insert(id, reference(fields, 2)?);
            }
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE" => {
                uses.insert(id, fields(&graph, id, "NEXT_ASSEMBLY_USAGE_OCCURRENCE")?);
            }
            "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION" => {
                let fields = fields(&graph, id, "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")?;
                let shape = reference(fields, 1)?;
                if contexts.insert(shape, reference(fields, 0)?).is_some() {
                    return Err(invalid(format!(
                        "shape #{shape} has multiple occurrence transforms"
                    )));
                }
            }
            "STYLED_ITEM" => {
                let fields = fields(&graph, id, "STYLED_ITEM")?;
                styles
                    .entry(reference(fields, 2)?)
                    .or_default()
                    .extend(references(fields, 1)?);
            }
            _ => {
                if entry.part("SHAPE_REPRESENTATION_RELATIONSHIP").is_some()
                    && entry
                        .part("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
                        .is_none()
                {
                    let relation = entry
                        .part("REPRESENTATION_RELATIONSHIP")
                        .unwrap_or(&entry.parts[0].1);
                    let a = reference(relation, 2)?;
                    let b = reference(relation, 3)?;
                    leaf_relationships.entry(a).or_default().push(b);
                    leaf_relationships.entry(b).or_default().push(a);
                }
            }
        }
    }
    let mut transforms = BTreeMap::new();
    for (shape, relationship) in contexts {
        let usage = *shape_definitions
            .get(&shape)
            .ok_or_else(|| invalid(format!("missing shape #{shape}")))?;
        if transforms.insert(usage, relationship).is_some() {
            return Err(invalid(format!("usage #{usage} has multiple transforms")));
        }
    }

    let indices: BTreeMap<_, _> = representations
        .keys()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect();
    let representation_ids: BTreeSet<_> = representations
        .values()
        .flat_map(|rep| {
            std::iter::once(rep.representation).chain(
                leaf_relationships
                    .get(&rep.representation)
                    .into_iter()
                    .flatten()
                    .copied(),
            )
        })
        .collect();
    let mut representation_items = BTreeMap::new();
    for id in representation_ids {
        let (frames, geometry): (BTreeSet<_>, BTreeSet<_>) = items(&graph, id)?
            .into_iter()
            .partition(|item| is_frame(&graph, *item));
        representation_items.insert(id, (frames, geometry));
    }
    let mut child_definitions = BTreeSet::new();
    let mut children = BTreeMap::<u64, Vec<Occurrence>>::new();
    for (id, args) in uses {
        let parent = reference(args, 3)?;
        let child = reference(args, 4)?;
        let locate = |definition| {
            representations
                .get(&definition)
                .ok_or_else(|| invalid(format!("missing representation for product #{definition}")))
        };
        let parent_rep = locate(parent)?.representation;
        let child_rep = locate(child)?.representation;
        let relationship = *transforms
            .get(&id)
            .ok_or_else(|| invalid(format!("occurrence #{id} has no transform")))?;
        let rel = fields(&graph, relationship, "REPRESENTATION_RELATIONSHIP")?;
        let a = reference(rel, 2)?;
        let b = reference(rel, 3)?;
        if (a, b) != (child_rep, parent_rep) && (a, b) != (parent_rep, child_rep) {
            return Err(invalid(format!(
                "occurrence #{id} relates the wrong representations"
            )));
        }
        let transform = reference(
            fields(
                &graph,
                relationship,
                "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
            )?,
            0,
        )?;
        let frames = fields(&graph, transform, "ITEM_DEFINED_TRANSFORMATION")?;
        let first = reference(frames, 2)?;
        let second = reference(frames, 3)?;
        if !representation_items[&a].0.contains(&first)
            || !representation_items[&b].0.contains(&second)
        {
            return Err(invalid(format!(
                "occurrence #{id} uses frames outside its representations"
            )));
        }
        let pose = |frame, dec: &mut Decider| {
            let placement = source_placement(&graph, frame, &scale, dec)?;
            let vector = |[x, y, z]: crate::Vec3| RVector3::new(x, y, z);
            let rotation = Rotation::from_frame(
                &vector(placement.axis).to_real(),
                &vector(placement.ref_dir).to_real(),
                dec,
            )
            .map_err(crate::recognize::source_error)?;
            Ok::<_, ImportRefusal>(Pose::from_rational(rotation, vector(placement.origin)))
        };
        let motion = pose(first, dec)?.inverse().then(&pose(second, dec)?);
        children.entry(parent).or_default().push(Occurrence {
            name: string(args.get(1))
                .ok_or_else(|| invalid(format!("occurrence #{id} has no name")))?,
            product: indices[&child],
            pose: if a == child_rep {
                motion
            } else {
                motion.inverse()
            },
        });
        child_definitions.insert(child);
    }
    let roots: Vec<_> = representations
        .keys()
        .filter(|id| !child_definitions.contains(id))
        .collect();
    let [root] = roots.as_slice() else {
        return Err(invalid("expected a single assembly root"));
    };
    let root = indices[root];
    let leaf_shapes = representations
        .iter()
        .filter(|(definition, _)| !children.contains_key(definition))
        .map(|(&definition, rep)| (rep.shape, definition))
        .collect();
    let mut pmi = read_assembly_pmi(&graph, &leaf_shapes)?;
    let mut products = Vec::new();
    for (&id, rep) in &representations {
        let definition = fields(&graph, id, "PRODUCT_DEFINITION")?;
        let formation = reference(definition, 2)?;
        let formation = graph
            .get(&formation)
            .ok_or_else(|| invalid("missing product formation"))?;
        let args = formation
            .part("PRODUCT_DEFINITION_FORMATION")
            .or_else(|| formation.part("PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE"))
            .ok_or_else(|| invalid("unsupported product formation"))?;
        let product = reference(args, 2)?;
        let name = string(fields(&graph, product, "PRODUCT")?.get(1))
            .ok_or_else(|| invalid("missing product name"))?;
        let mut related = vec![rep.representation];
        related.extend(
            leaf_relationships
                .get(&rep.representation)
                .into_iter()
                .flatten()
                .copied(),
        );
        let content = if let Some(children) = children.remove(&id) {
            for representation in related {
                if !representation_items[&representation].1.is_empty() {
                    return Err(invalid(format!(
                        "assembly #{id} also contains direct geometry"
                    )));
                }
            }
            ProductContent::Assembly(children)
        } else {
            let mut local = Graph::new();
            let geometry_items = related
                .iter()
                .flat_map(|representation| representation_items[representation].1.iter().copied())
                .collect::<BTreeSet<_>>();
            for &item in &geometry_items {
                let entity = graph
                    .get(&item)
                    .ok_or_else(|| invalid(format!("missing representation item #{item}")))?;
                if entity.type_name() != "MANIFOLD_SOLID_BREP" {
                    return Err(ImportRefusal::Unsupported {
                        entity: format!("#{item}"),
                        type_name: entity.type_name().into(),
                    });
                }
            }
            let mut pending = related;
            while let Some(entity) = pending.pop() {
                if local.contains_key(&entity) {
                    continue;
                }
                let entry = graph
                    .get(&entity)
                    .ok_or_else(|| invalid(format!("missing entity #{entity}")))?;
                for (_, params) in &entry.parts {
                    collect_refs(params, &mut pending);
                }
                local.insert(entity, entry.clone());
            }
            let geometry = read_graph(&local, dec)?;
            if geometry.solids.is_empty()
                || geometry.solids.len() != geometry_items.len()
                || !geometry.surface_models.is_empty()
            {
                return Err(invalid(format!(
                    "part #{id} has unsupported representation items"
                )));
            }
            let colour = part_colour(&graph, &styles, &geometry_items)?;
            let pmi = pmi.remove(&rep.shape).expect("indexed leaf PMI");
            let faces: BTreeSet<_> = geometry
                .named_faces()
                .map(|(_, face)| face.name.as_str())
                .collect();
            for item in &pmi.pmi.items {
                if item
                    .kind
                    .faces()
                    .into_iter()
                    .any(|face| !faces.contains(face.as_str()))
                {
                    return Err(invalid(format!(
                        "PMI item {:?} references a face outside product #{id}",
                        item.id
                    )));
                }
            }
            ProductContent::Part(Part {
                geometry,
                colour,
                pmi,
            })
        };
        products.push(Product { name, content });
    }
    let assembly = Assembly { root, products };
    assembly.validate().map_err(|r| invalid(r.detail))?;
    Ok(assembly)
}

fn invalid(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::Malformed {
        detail: detail.into(),
    }
}

fn fields<'a>(graph: &'a Graph, id: u64, kind: &str) -> Result<&'a [Param], ImportRefusal> {
    graph
        .get(&id)
        .and_then(|e| e.part(kind))
        .ok_or_else(|| invalid(format!("expected {kind} at #{id}")))
}

fn reference(fields: &[Param], index: usize) -> Result<u64, ImportRefusal> {
    match fields.get(index) {
        Some(Param::Ref(id)) => Ok(*id),
        _ => Err(invalid(format!(
            "expected reference at parameter {}",
            index + 1
        ))),
    }
}

fn references(fields: &[Param], index: usize) -> Result<Vec<u64>, ImportRefusal> {
    let Some(Param::List(list)) = fields.get(index) else {
        return Err(invalid(format!(
            "expected a list at parameter {}",
            index + 1
        )));
    };
    (0..list.len()).map(|i| reference(list, i)).collect()
}

fn collect_refs(fields: &[Param], out: &mut Vec<u64>) {
    for field in fields {
        match field {
            Param::Ref(id) => out.push(*id),
            Param::List(list) | Param::Typed(_, list) => collect_refs(list, out),
            _ => {}
        }
    }
}

fn items(graph: &Graph, representation: u64) -> Result<Vec<u64>, ImportRefusal> {
    let entity = graph
        .get(&representation)
        .ok_or_else(|| invalid(format!("missing representation #{representation}")))?;
    if !matches!(
        entity.type_name(),
        "SHAPE_REPRESENTATION" | "ADVANCED_BREP_SHAPE_REPRESENTATION"
    ) {
        return Err(invalid(format!(
            "unsupported representation #{representation}"
        )));
    }
    references(&entity.parts[0].1, 1)
}

fn is_frame(graph: &Graph, id: u64) -> bool {
    graph
        .get(&id)
        .is_some_and(|entity| entity.part("AXIS2_PLACEMENT_3D").is_some())
}

fn part_colour(
    graph: &Graph,
    styles: &BTreeMap<u64, Vec<u64>>,
    solids: &BTreeSet<u64>,
) -> Result<Option<Colour>, ImportRefusal> {
    let mut colours = BTreeSet::new();
    for solid in solids {
        let mut side_colours = [None, None];
        let mut has_one_sided_style = false;
        for assignment in styles.get(solid).into_iter().flatten() {
            let fields = fields(graph, *assignment, "PRESENTATION_STYLE_ASSIGNMENT")?;
            for style in references(fields, 0)? {
                let Some(usage) = graph
                    .get(&style)
                    .and_then(|entry| entry.part("SURFACE_STYLE_USAGE"))
                else {
                    continue;
                };
                let sides: &[usize] = match usage.first() {
                    Some(Param::Enum(side)) if side == "BOTH" => &[0, 1],
                    Some(Param::Enum(side)) if side == "POSITIVE" => &[0],
                    Some(Param::Enum(side)) if side == "NEGATIVE" => &[1],
                    _ => return Err(invalid("unknown surface style side")),
                };
                has_one_sided_style |= sides.len() == 1;
                let side = reference(usage, 1)?;
                for fill in references(self::fields(graph, side, "SURFACE_SIDE_STYLE")?, 1)? {
                    let Some(fill) = graph
                        .get(&fill)
                        .and_then(|entry| entry.part("SURFACE_STYLE_FILL_AREA"))
                    else {
                        continue;
                    };
                    let area = reference(fill, 0)?;
                    for paint in references(self::fields(graph, area, "FILL_AREA_STYLE")?, 1)? {
                        let Some(paint) = graph
                            .get(&paint)
                            .and_then(|entry| entry.part("FILL_AREA_STYLE_COLOUR"))
                        else {
                            continue;
                        };
                        let rgb = reference(paint, 1)?;
                        let args = self::fields(graph, rgb, "COLOUR_RGB")?;
                        let read = |i| {
                            let q = args
                                .get(i)
                                .ok_or_else(|| invalid("expected an RGB channel"))
                                .and_then(number)?;
                            if q < Rational::zero() || q > Rational::one() {
                                return Err(invalid("RGB channels must be between zero and one"));
                            }
                            Ok(q)
                        };
                        let current = [read(1)?, read(2)?, read(3)?];
                        for &side in sides {
                            let colour = &mut side_colours[side];
                            if colour.as_ref().is_some_and(|previous| previous != &current) {
                                return Err(invalid("one side of a solid has conflicting colours"));
                            }
                            *colour = Some(current.clone());
                        }
                    }
                }
            }
        }
        let [positive, negative] = side_colours;
        if positive != negative || (has_one_sided_style && positive.is_none()) {
            return Err(invalid("a part colour must resolve equally on both sides"));
        }
        colours.insert(positive);
    }
    if colours.len() > 1 {
        return Err(invalid("bodies of one product have different colours"));
    }
    Ok(colours.pop_first().flatten())
}
