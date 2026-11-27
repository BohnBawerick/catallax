//! AP242 product definitions and rigid occurrences.

use std::collections::BTreeSet;
use std::io::Write;

use catallax_numbers::{Decider, MarkLog, Rational, Refusal};
use catallax_shapes::{Point3, Pose, Vector3};
use catallax_topology::{BodyId, BodyKind, FaceId, Model};

use crate::document::{self, MILLI, Preamble, Product as ProductHandles};
use crate::graph::{Field, Graph, Id, Layout, refs, string, text};
use crate::{StepOptions, StepReport, pmi::Pmi};

/// A product and its children or part data. Indices are local to the document.
#[derive(Clone, Debug)]
pub struct Product<T> {
    /// Product definition name, shared by all its occurrences.
    pub name: String,
    /// A part definition or a group of placed occurrences.
    pub content: ProductContent<T>,
}

/// The contents of a product definition.
#[derive(Clone, Debug)]
pub enum ProductContent<T> {
    /// Geometry and metadata in the product's local coordinates.
    Part(T),
    /// Ordered placements of child definitions relative to this product.
    Assembly(Vec<Occurrence>),
}

/// An instance of a product in an assembly.
#[derive(Clone, Debug)]
pub struct Occurrence {
    /// The occurrence name, independent of its product name.
    pub name: String,
    /// Index of the child product definition.
    pub product: usize,
    /// Motion from child coordinates into parent coordinates.
    pub pose: Pose,
}

/// A rooted product graph. Several occurrences can share the same definition.
#[derive(Clone, Debug)]
pub struct Assembly<T> {
    /// Index of the root assembly.
    pub root: usize,
    /// Local definitions, each reachable from the root.
    pub products: Vec<Product<T>>,
}

impl<T> Assembly<T> {
    /// Checks references, reachability and cycles without flattening instances.
    pub fn validate(&self) -> Result<(), Refusal> {
        if !matches!(
            self.products.get(self.root).map(|p| &p.content),
            Some(ProductContent::Assembly(_))
        ) {
            return Err(invalid("the root must be an assembly"));
        }
        let mut active = BTreeSet::new();
        let mut visited = BTreeSet::new();
        let mut pending = vec![(self.root, false)];
        while let Some((index, leaving)) = pending.pop() {
            if leaving {
                active.remove(&index);
                visited.insert(index);
                continue;
            }
            if visited.contains(&index) {
                continue;
            }
            if !active.insert(index) {
                return Err(invalid("cycle in product definitions"));
            }
            let product = self
                .products
                .get(index)
                .ok_or_else(|| invalid("missing product definition"))?;
            pending.push((index, true));
            if let ProductContent::Assembly(children) = &product.content {
                if children.is_empty() {
                    return Err(invalid("an assembly has no children"));
                }
                pending.extend(children.iter().rev().map(|child| (child.product, false)));
            }
        }
        if visited.len() != self.products.len() {
            return Err(invalid("a product is unreachable from the root"));
        }
        Ok(())
    }
}

fn invalid(detail: &str) -> Refusal {
    Refusal::new("step.assembly-invalid", detail)
}

/// Surface colour for a part definition, in STEP's normalized sRGB channels.
pub type Colour = [Rational; 3];

/// A kernel part definition with optional appearance and semantic PMI.
pub struct Part<'a> {
    /// Model containing the part. Different definitions may use different models.
    pub model: &'a Model,
    /// Solid in product-local coordinates.
    pub body: BodyId,
    /// Colour inherited by all occurrences of this definition.
    pub colour: Option<Colour>,
    /// Manufacturing information attached to faces in this model.
    pub pmi: Option<&'a Pmi<FaceId>>,
}

/// Reports the underlying part exports in product order.
pub struct AssemblyReport {
    /// Part reports paired with their definition indices, one per part.
    pub parts: Vec<(usize, StepReport)>,
    /// Rounding bounds for the occurrence frames in product and child order.
    pub placements: Vec<PlacementRounding>,
    /// Marks issued by the part writers.
    pub marks: MarkLog,
}

/// Bounds for decimal rounding of one local placement.
pub struct PlacementRounding {
    /// Index of the assembly that owns the occurrence.
    pub parent: usize,
    /// Index in that assembly's occurrence list.
    pub child: usize,
    /// Bound in millimetres on displacement of the local origin.
    pub origin_error: Rational,
    /// Bound on the sum of the three basis-vector errors after normalization.
    pub basis_error: Rational,
}

/// Writes a rooted assembly as one AP242 document, using the ordinary part
/// writers for each unique definition. Validation completes before writing.
pub fn write_step_assembly<W: Write>(
    assembly: &Assembly<Part<'_>>,
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Result<AssemblyReport, Refusal> {
    dec.charge_work(assembly.products.len() as u64)?;
    assembly.validate()?;
    let name = &assembly.products[assembly.root].name;
    let mut graph = Graph::default();
    let mut layout = Layout::default();
    let root = Preamble::add(
        &mut graph,
        name,
        "SHAPE_REPRESENTATION",
        &[],
        None,
        MILLI,
        true,
    );
    layout.prefix(&root.slots);
    let mut products = Vec::with_capacity(assembly.products.len());
    let mut reports = Vec::new();
    let mut placements = Vec::new();
    let mut marks = MarkLog::new();
    for (index, product) in assembly.products.iter().enumerate() {
        let entry = match &product.content {
            ProductContent::Assembly(_) if index == assembly.root => {
                (root.references.product, None)
            }
            ProductContent::Assembly(_) => {
                let (product, entities) = root.subassembly(&mut graph, &product.name);
                layout.reserve(&entities);
                (product, Some(entities))
            }
            ProductContent::Part(part) => {
                let body = part.model.body(part.body);
                if body.kind != BodyKind::Solid {
                    return Err(invalid("a part definition must be a solid"));
                }
                if body
                    .lumps
                    .iter()
                    .any(|&lump| part.model.lump(lump).shells.len() > 1)
                {
                    return Err(Refusal::new(
                        "step.assembly-unsupported-part",
                        format!(
                            "part {:?} has a void shell outside the assembly reader's solid subset",
                            product.name
                        ),
                    ));
                }
                let bodies = [part.body];
                let mut settings = options.clone();
                settings.product_name.clone_from(&product.name);
                let answer = crate::build_part(part.model, &bodies, &settings, part.pmi, dec)?;
                if settings.seams == crate::Seams::Required {
                    // Keep part-local seam diagnostics and their refusal order.
                    crate::seams::check(&answer.value.file(&settings.product_name))?;
                }
                dec.charge_work(answer.value.layout.part_span() as u64)?;
                let (product, context, items, report) =
                    answer.value.append(&mut graph, &mut layout);
                if let Some(colour) = &part.colour {
                    add_colour(&mut graph, &mut layout, colour, context, items)?;
                }
                reports.push((index, report));
                marks.extend(&answer.marks);
                (product, None)
            }
        };
        products.push(entry);
    }

    for (index, product) in assembly.products.iter().enumerate() {
        let ProductContent::Assembly(children) = &product.content else {
            continue;
        };
        let (parent, entities) = &products[index];
        let mut items = Vec::with_capacity(children.len());
        dec.charge_work(children.len() as u64)?;
        for (child_index, child) in children.iter().enumerate() {
            let target = products[child.product].0;
            let (frame, origin_error, basis_error) = add_frame(
                &mut graph,
                &mut layout,
                &child.pose,
                options.significant_digits,
                dec,
            )?;
            placements.push(PlacementRounding {
                parent: index,
                child: child_index,
                origin_error,
                basis_error,
            });
            items.push(frame);
            layout.extend(add_occurrence(
                &mut graph,
                &child.name,
                *parent,
                target,
                frame,
            ));
        }
        parent.add_frames(&mut graph, items);
        if let Some(entities) = entities {
            layout.emit(entities);
        }
    }

    let complete_body_bound = reports
        .iter()
        .all(|(_, report)| report.definition_deviation.is_some());
    let description: &[_] = if complete_body_bound {
        &["catallax assembly export"]
    } else {
        &[
            "catallax assembly export",
            "at least one part lacks a complete deviation bound; complete_body_bound=unestablished",
        ]
    };
    graph.unique_draughting_model_names();
    let file = document::render(&graph, &layout, name, description);
    out.write_all(file.as_bytes())
        .map_err(|error| Refusal::new("step.io-error", error.to_string()))?;
    Ok(AssemblyReport {
        parts: reports,
        placements,
        marks,
    })
}

fn add_occurrence(
    graph: &mut Graph,
    name: &str,
    parent: ProductHandles,
    child: ProductHandles,
    frame: Id,
) -> [Id; 5] {
    const OCCURRENCE: &str = "NEXT_ASSEMBLY_USAGE_OCCURRENCE";
    let occurrence = graph.add(OCCURRENCE, Vec::new());
    graph.replace(
        occurrence,
        OCCURRENCE,
        vec![
            Field::QuotedNumber(occurrence),
            text(string(name)),
            text("''"),
            parent.definition.into(),
            child.definition.into(),
            text("$"),
        ],
    );
    let shape = graph.add(
        "PRODUCT_DEFINITION_SHAPE",
        vec![text("''"), text("''"), occurrence.into()],
    );
    // Map the child frame to the occurrence frame in the parent representation.
    let transform = graph.add(
        "ITEM_DEFINED_TRANSFORMATION",
        vec![text("''"), text("''"), child.frame.into(), frame.into()],
    );
    let relationship = graph.add_complex(vec![
        (
            "REPRESENTATION_RELATIONSHIP",
            vec![
                text("''"),
                text("''"),
                child.representation.into(),
                parent.representation.into(),
            ],
        ),
        (
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
            vec![transform.into()],
        ),
        ("SHAPE_REPRESENTATION_RELATIONSHIP", vec![]),
    ]);
    let dependent = graph.add(
        "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION",
        vec![relationship.into(), shape.into()],
    );
    [occurrence, shape, transform, relationship, dependent]
}

fn add_frame(
    graph: &mut Graph,
    layout: &mut Layout,
    pose: &Pose,
    digits: u32,
    dec: &mut Decider,
) -> Result<(Id, Rational, Rational), Refusal> {
    let matrix = pose.rotation.matrix();
    let basis = |column: usize| {
        Vector3::new(
            matrix[0][column].clone(),
            matrix[1][column].clone(),
            matrix[2][column].clone(),
        )
    };
    let origin = crate::number::point(&Point3::origin().plus(&pose.translation), digits, dec)?;
    let z = crate::number::direction(&basis(2), digits, dec)?;
    let x = crate::number::direction(&basis(0), digits, dec)?;
    let triple = |point: &crate::number::Point| {
        let [a, b, c] = point.xyz.each_ref();
        text(format!("({},{},{})", a.text, b.text, c.text))
    };
    let p = graph.add("CARTESIAN_POINT", vec![text("''"), triple(&origin)]);
    let axis = graph.add("DIRECTION", vec![text("''"), triple(&z)]);
    let reference = graph.add("DIRECTION", vec![text("''"), triple(&x)]);
    let frame = graph.add(
        "AXIS2_PLACEMENT_3D",
        vec![text("''"), p.into(), axis.into(), reference.into()],
    );
    layout.extend([p, axis, reference, frame]);
    Ok((frame, origin.error(), crate::number::frame_error(&x, &z)))
}

fn add_colour(
    graph: &mut Graph,
    layout: &mut Layout,
    colour: &Colour,
    context: Id,
    targets: Vec<Id>,
) -> Result<(), Refusal> {
    if colour
        .iter()
        .any(|value| value < &Rational::zero() || value > &Rational::one())
    {
        return Err(invalid("RGB channels must be between zero and one"));
    }
    if targets.is_empty() {
        return Err(invalid("part has no supported B-rep representation items"));
    }
    let [r, g, b] = colour
        .each_ref()
        .map(|c| crate::number::rational(c, 17).text);
    let mut add = |name, fields| {
        let id = graph.add(name, fields);
        layout.push(id);
        id
    };
    let rgb = add("COLOUR_RGB", vec![text("''"), text(r), text(g), text(b)]);
    let colour = add("FILL_AREA_STYLE_COLOUR", vec![text("''"), rgb.into()]);
    let fill = add("FILL_AREA_STYLE", vec![text("''"), refs([colour])]);
    let area = add("SURFACE_STYLE_FILL_AREA", vec![fill.into()]);
    let side = add("SURFACE_SIDE_STYLE", vec![text("''"), refs([area])]);
    let usage = add("SURFACE_STYLE_USAGE", vec![text(".BOTH."), side.into()]);
    let assignment = add("PRESENTATION_STYLE_ASSIGNMENT", vec![refs([usage])]);
    let styles = targets
        .into_iter()
        .map(|target| {
            add(
                "STYLED_ITEM",
                vec![text("''"), refs([assignment]), target.into()],
            )
        })
        .collect::<Vec<_>>();
    add(
        "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION",
        vec![text("''"), refs(styles), context.into()],
    );
    Ok(())
}
