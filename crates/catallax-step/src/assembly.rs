//! AP242 product definitions and rigid occurrences.

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::io::Write;

use catallax_numbers::{Decider, MarkLog, Rational, Refusal};
use catallax_shapes::{Point3, Pose, Vector3};
use catallax_topology::{BodyId, BodyKind, FaceId, Model};

use crate::graph::string;
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

#[derive(Clone, Copy)]
struct ProductIds {
    definition: usize,
    representation: usize,
    frame: usize,
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
    let mut step = String::new();
    let mut next = 21; // The assembly header reserves the original STEP prefix.
    let mut products = Vec::with_capacity(assembly.products.len());
    let mut reports = Vec::new();
    let mut placements = Vec::new();
    let mut marks = MarkLog::new();
    for (index, product) in assembly.products.iter().enumerate() {
        let ids = match &product.content {
            ProductContent::Assembly(_) if index == assembly.root => ProductIds {
                definition: 7,
                representation: 16,
                frame: 17,
            },
            ProductContent::Assembly(_) => {
                let id = next;
                next += 7;
                ProductIds {
                    definition: id + 2,
                    representation: id + 5,
                    frame: 17,
                }
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
                let mut file = Vec::new();
                let mut settings = options.clone();
                settings.product_name.clone_from(&product.name);
                let answer = match part.pmi {
                    Some(pmi) => crate::write_step_with_pmi(
                        part.model, part.body, &settings, pmi, dec, &mut file,
                    )?,
                    None => crate::write_step(part.model, &[part.body], &settings, dec, &mut file)?,
                };
                let file = std::str::from_utf8(&file).expect("the STEP writer emits ASCII");
                let base = next - 1;
                let entities = append_part(file, base, &mut step)?;
                dec.charge_work(entities as u64)?;
                next += entities;
                if let Some(colour) = &part.colour {
                    write_colour(colour, file, base, &mut next, &mut step)?;
                }
                reports.push((index, answer.value));
                marks.extend(&answer.marks);
                ProductIds {
                    definition: base + 7,
                    representation: base + 16,
                    frame: base + 17,
                }
            }
        };
        products.push(ids);
    }

    let mut root_items = Vec::new();
    for (index, product) in assembly.products.iter().enumerate() {
        let ProductContent::Assembly(children) = &product.content else {
            continue;
        };
        let parent = products[index];
        let mut items = Vec::with_capacity(children.len());
        dec.charge_work(children.len() as u64)?;
        for (child_index, child) in children.iter().enumerate() {
            let target = products[child.product];
            let (frame, origin_error, basis_error) = write_frame(
                &child.pose,
                options.significant_digits,
                dec,
                &mut next,
                &mut step,
            )?;
            placements.push(PlacementRounding {
                parent: index,
                child: child_index,
                origin_error,
                basis_error,
            });
            items.push(frame);
            let id = next;
            next += 5;
            let name = string(&child.name);
            // OCCT's ProSTEP assembly relationship maps the child frame to the
            // occurrence frame listed in the parent representation.
            writeln!(
                step,
                "#{id}=NEXT_ASSEMBLY_USAGE_OCCURRENCE('{id}',{name},'',#{},#{},$);",
                parent.definition, target.definition
            )
            .unwrap();
            writeln!(step, "#{}=PRODUCT_DEFINITION_SHAPE('','',#{id});", id + 1).unwrap();
            writeln!(
                step,
                "#{}=ITEM_DEFINED_TRANSFORMATION('','',#{},#{frame});",
                id + 2,
                target.frame
            )
            .unwrap();
            writeln!(step, "#{}=(REPRESENTATION_RELATIONSHIP('','',#{},#{})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{})SHAPE_REPRESENTATION_RELATIONSHIP());",
                     id + 3, target.representation, parent.representation, id + 2).unwrap();
            writeln!(
                step,
                "#{}=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{},#{});",
                id + 4,
                id + 3,
                id + 1
            )
            .unwrap();
        }
        if index == assembly.root {
            root_items = items;
        } else {
            let id = parent.definition - 2;
            let name = string(&product.name);
            let items = items.iter().fold(String::from("#17"), |mut list, id| {
                write!(list, ",#{id}").unwrap();
                list
            });
            writeln!(step, "#{id}=PRODUCT({name},{name},'',(#3));").unwrap();
            writeln!(
                step,
                "#{}=PRODUCT_DEFINITION_FORMATION('','',#{id});",
                id + 1
            )
            .unwrap();
            writeln!(
                step,
                "#{}=PRODUCT_DEFINITION('design','',#{},#6);",
                id + 2,
                id + 1
            )
            .unwrap();
            writeln!(
                step,
                "#{}=PRODUCT_DEFINITION_SHAPE('','',#{});",
                id + 3,
                id + 2
            )
            .unwrap();
            writeln!(
                step,
                "#{}=SHAPE_DEFINITION_REPRESENTATION(#{},#{});",
                id + 4,
                id + 3,
                id + 5
            )
            .unwrap();
            writeln!(
                step,
                "#{}=SHAPE_REPRESENTATION({name},({items}),#14);",
                id + 5
            )
            .unwrap();
            writeln!(
                step,
                "#{}=PRODUCT_RELATED_PRODUCT_CATEGORY('assembly',$,(#{id}));",
                id + 6
            )
            .unwrap();
        }
    }

    let complete_body_bound = reports
        .iter()
        .all(|(_, report)| report.definition_deviation.is_some());
    let header = crate::export::assembly_header(
        &assembly.products[assembly.root].name,
        &root_items,
        complete_body_bound,
    );
    out.write_all(header.as_bytes())
        .and_then(|_| out.write_all(step.as_bytes()))
        .and_then(|_| out.write_all(b"ENDSEC;\nEND-ISO-10303-21;\n"))
        .map_err(|error| Refusal::new("step.io-error", error.to_string()))?;
    Ok(AssemblyReport {
        parts: reports,
        placements,
        marks,
    })
}

/// Part writers all emit one Part 21 DATA section with symbolic references
/// already numbered. Relocate those numbers without touching quoted strings.
fn append_part(file: &str, base: usize, output: &mut String) -> Result<usize, Refusal> {
    let data = file
        .split_once("\nDATA;\n")
        .and_then(|(_, data)| data.strip_suffix("ENDSEC;\nEND-ISO-10303-21;\n"))
        .ok_or_else(|| invalid("unexpected part writer envelope"))?;
    let mut quoted = false;
    let mut highest = 0;
    let mut at = 0;
    for (index, byte) in data.bytes().enumerate() {
        if byte == b'\'' {
            quoted = !quoted;
        }
        if byte != b'#' || quoted {
            continue;
        }
        let end = index
            + 1
            + data[index + 1..]
                .bytes()
                .take_while(u8::is_ascii_digit)
                .count();
        let entity: usize = data[index + 1..end]
            .parse()
            .map_err(|_| invalid("invalid part entity number"))?;
        output.push_str(&data[at..index]);
        write!(output, "#{}", base + entity).unwrap();
        highest = highest.max(entity);
        at = end;
    }
    output.push_str(&data[at..]);
    Ok(highest)
}

fn write_frame(
    pose: &Pose,
    digits: u32,
    dec: &mut Decider,
    next: &mut usize,
    output: &mut String,
) -> Result<(usize, Rational, Rational), Refusal> {
    let id = *next;
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
        format!("{},{},{}", a.text, b.text, c.text)
    };
    writeln!(output, "#{id}=CARTESIAN_POINT('',({}));", triple(&origin)).unwrap();
    writeln!(output, "#{}=DIRECTION('',({}));", id + 1, triple(&z)).unwrap();
    writeln!(output, "#{}=DIRECTION('',({}));", id + 2, triple(&x)).unwrap();
    writeln!(
        output,
        "#{}=AXIS2_PLACEMENT_3D('',#{id},#{},#{});",
        id + 3,
        id + 1,
        id + 2
    )
    .unwrap();
    *next += 4;
    Ok((id + 3, origin.error(), crate::number::frame_error(&x, &z)))
}

fn write_colour(
    colour: &Colour,
    part: &str,
    base: usize,
    next: &mut usize,
    output: &mut String,
) -> Result<(), Refusal> {
    if colour
        .iter()
        .any(|value| value < &Rational::zero() || value > &Rational::one())
    {
        return Err(invalid("RGB channels must be between zero and one"));
    }
    // One style can target every solid lump in the part's B-rep representation.
    let line = part
        .lines()
        .find(|line| line.starts_with("#16="))
        .and_then(|line| line.rsplit_once(",(#17,"))
        .and_then(|(_, items)| items.strip_suffix("),#14);"))
        .ok_or_else(|| invalid("part has no supported B-rep representation items"))?;
    let id = *next;
    *next += 6;
    let [r, g, b] = colour
        .each_ref()
        .map(|c| crate::number::rational(c, 17).text);
    writeln!(output, "#{id}=COLOUR_RGB('',{r},{g},{b});").unwrap();
    writeln!(output, "#{}=FILL_AREA_STYLE_COLOUR('',#{id});", id + 1).unwrap();
    writeln!(output, "#{}=FILL_AREA_STYLE('',(#{}));", id + 2, id + 1).unwrap();
    writeln!(output, "#{}=SURFACE_STYLE_FILL_AREA(#{});", id + 3, id + 2).unwrap();
    writeln!(output, "#{}=SURFACE_SIDE_STYLE('',(#{}));", id + 4, id + 3).unwrap();
    writeln!(
        output,
        "#{}=SURFACE_STYLE_USAGE(.BOTH.,#{});",
        id + 5,
        id + 4
    )
    .unwrap();
    let assignment = *next;
    *next += 1;
    writeln!(
        output,
        "#{assignment}=PRESENTATION_STYLE_ASSIGNMENT((#{}));",
        id + 5
    )
    .unwrap();
    let mut styles = String::new();
    for item in line.split(',') {
        let target: usize = item
            .strip_prefix('#')
            .and_then(|n| n.parse().ok())
            .ok_or_else(|| invalid("unexpected B-rep representation item"))?;
        if !styles.is_empty() {
            styles.push(',');
        }
        let id = *next;
        *next += 1;
        write!(styles, "#{id}").unwrap();
        writeln!(
            output,
            "#{id}=STYLED_ITEM('',(#{assignment}),#{});",
            base + target
        )
        .unwrap();
    }
    let id = *next;
    *next += 1;
    writeln!(
        output,
        "#{id}=MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('',({styles}),#{});",
        base + 14
    )
    .unwrap();
    Ok(())
}
