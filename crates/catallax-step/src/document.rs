//! Shared AP242 product and representation declarations.

use catallax_numbers::Refusal;

use crate::StepReport;
use crate::graph::{Field, Graph, Id, Layout, refs, string, text};

pub(crate) const MILLI: &str = ".MILLI.";
pub(crate) const BREP: &str = "ADVANCED_BREP_SHAPE_REPRESENTATION";
pub(crate) const FILE_END: &str = "ENDSEC;\nEND-ISO-10303-21;\n";

#[derive(Clone, Copy)]
pub(crate) struct Product {
    pub(crate) definition: Id,
    pub(crate) shape: Id,
    pub(crate) representation: Id,
    pub(crate) frame: Id,
}

impl Product {
    pub(crate) fn add_frames(self, graph: &mut Graph, frames: impl IntoIterator<Item = Id>) {
        graph.extend_list(self.representation, 1, frames);
    }
}

#[derive(Clone, Copy)]
pub(crate) struct PartReferences {
    pub(crate) product: Product,
    pub(crate) context: Id,
    pub(crate) length_unit: Id,
}

pub(crate) struct Preamble {
    pub(crate) references: PartReferences,
    pub(crate) slots: [Option<Id>; 20],
    product_context: Id,
    definition_context: Id,
}

impl Preamble {
    pub(crate) fn add(
        graph: &mut Graph,
        name: &str,
        representation: &'static str,
        items: &[Id],
        uncertainty: Option<&str>,
        length_prefix: &'static str,
        assembly: bool,
    ) -> Self {
        let app = graph.add(
            "APPLICATION_CONTEXT",
            vec![text("'managed model based 3d engineering'")],
        );
        let protocol = graph.add(
            "APPLICATION_PROTOCOL_DEFINITION",
            vec![
                text("'international standard'"),
                text("'ap242_managed_model_based_3d_engineering'"),
                text("2014"),
                app.into(),
            ],
        );
        let product_context = graph.add(
            "PRODUCT_CONTEXT",
            vec![text("''"), app.into(), text("'mechanical'")],
        );
        let definition_context = graph.add(
            "PRODUCT_DEFINITION_CONTEXT",
            vec![text("'part definition'"), app.into(), text("'design'")],
        );
        let length_unit = graph.add_complex(vec![
            ("LENGTH_UNIT", vec![]),
            ("NAMED_UNIT", vec![text("*")]),
            ("SI_UNIT", vec![text(length_prefix), text(".METRE.")]),
        ]);
        let angle = graph.add_complex(vec![
            ("NAMED_UNIT", vec![text("*")]),
            ("PLANE_ANGLE_UNIT", vec![]),
            ("SI_UNIT", vec![text("$"), text(".RADIAN.")]),
        ]);
        let solid_angle = graph.add_complex(vec![
            ("NAMED_UNIT", vec![text("*")]),
            ("SI_UNIT", vec![text("$"), text(".STERADIAN.")]),
            ("SOLID_ANGLE_UNIT", vec![]),
        ]);
        let uncertainty = uncertainty.map(|value| {
            graph.add(
                "UNCERTAINTY_MEASURE_WITH_UNIT",
                vec![
                    Field::Typed("LENGTH_MEASURE", vec![text(value.to_owned())]),
                    length_unit.into(),
                    text("'distance_accuracy_value'"),
                    text("'maximum deviation from the exact model'"),
                ],
            )
        });
        let mut contexts = vec![("GEOMETRIC_REPRESENTATION_CONTEXT", vec![text("3")])];
        if let Some(uncertainty) = uncertainty {
            contexts.push((
                "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
                vec![refs([uncertainty])],
            ));
        }
        contexts.extend([
            (
                "GLOBAL_UNIT_ASSIGNED_CONTEXT",
                vec![refs([length_unit, angle, solid_angle])],
            ),
            ("REPRESENTATION_CONTEXT", vec![text("''"), text("''")]),
        ]);
        let context = graph.add_complex(contexts);
        let origin = graph.add("CARTESIAN_POINT", vec![text("''"), text("(0.,0.,0.)")]);
        let z = graph.add("DIRECTION", vec![text("''"), text("(0.,0.,1.)")]);
        let x = graph.add("DIRECTION", vec![text("''"), text("(1.,0.,0.)")]);
        let frame = graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), origin.into(), z.into(), x.into()],
        );
        let (product, [entity, formation, definition, shape, link, rep, category]) = add_product(
            graph,
            name,
            [product_context, definition_context, context, frame],
            representation,
            items,
            assembly,
        );
        let slots = [
            Some(app),
            Some(protocol),
            Some(product_context),
            Some(entity),
            Some(formation),
            Some(definition_context),
            Some(definition),
            Some(shape),
            Some(link),
            Some(length_unit),
            Some(angle),
            Some(solid_angle),
            uncertainty,
            Some(context),
            Some(category),
            Some(rep),
            Some(frame),
            Some(origin),
            Some(z),
            Some(x),
        ];
        Self {
            references: PartReferences {
                product,
                context,
                length_unit,
            },
            slots,
            product_context,
            definition_context,
        }
    }

    pub(crate) fn subassembly(&self, graph: &mut Graph, name: &str) -> (Product, [Id; 7]) {
        add_product(
            graph,
            name,
            [
                self.product_context,
                self.definition_context,
                self.references.context,
                self.references.product.frame,
            ],
            "SHAPE_REPRESENTATION",
            &[],
            true,
        )
    }

    pub(crate) fn layout(&self, graph: &Graph, items: &[Id], annotations: &[Id]) -> Layout {
        let mut layout = Layout::default();
        layout.prefix(&self.slots);
        graph.depth_first(
            &mut layout,
            &self.slots,
            items.iter().chain(annotations).copied(),
        );
        layout
    }
}

fn add_product(
    graph: &mut Graph,
    name: &str,
    [product_context, definition_context, context, frame]: [Id; 4],
    representation: &'static str,
    items: &[Id],
    assembly: bool,
) -> (Product, [Id; 7]) {
    let name = text(string(name));
    let entity = graph.add(
        "PRODUCT",
        vec![
            name.clone(),
            name.clone(),
            text("''"),
            refs([product_context]),
        ],
    );
    let formation = graph.add(
        "PRODUCT_DEFINITION_FORMATION",
        vec![text("''"), text("''"), entity.into()],
    );
    let definition = graph.add(
        "PRODUCT_DEFINITION",
        vec![
            text("'design'"),
            text("''"),
            formation.into(),
            definition_context.into(),
        ],
    );
    let shape = graph.add(
        "PRODUCT_DEFINITION_SHAPE",
        vec![text("''"), text("''"), definition.into()],
    );
    let rep = graph.add(
        representation,
        vec![
            name,
            refs(std::iter::once(frame).chain(items.iter().copied())),
            context.into(),
        ],
    );
    let link = graph.add(
        "SHAPE_DEFINITION_REPRESENTATION",
        vec![shape.into(), rep.into()],
    );
    let category = graph.add(
        "PRODUCT_RELATED_PRODUCT_CATEGORY",
        vec![
            text(if assembly { "'assembly'" } else { "'part'" }),
            text("$"),
            refs([entity]),
        ],
    );
    (
        Product {
            definition,
            shape,
            representation: rep,
            frame,
        },
        [entity, formation, definition, shape, link, rep, category],
    )
}

pub(crate) fn render(graph: &Graph, layout: &Layout, name: &str, description: &[&str]) -> String {
    let description = description
        .iter()
        .map(|entry| string(entry))
        .collect::<Vec<_>>()
        .join(",");
    let filename = string(&format!("{name}.step"));
    let header = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(({description}),'2;1');\nFILE_NAME({filename},'1970-01-01T00:00:00',(''),(''),'catallax','catallax','');\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF {{ 1 0 10303 442 1 1 4 }}'));\nENDSEC;\nDATA;\n"
    );
    let mut file = graph.write_layout(layout, header);
    file.push_str(FILE_END);
    file
}

pub(crate) struct PartDocument {
    pub(crate) graph: Graph,
    pub(crate) layout: Layout,
    pub(crate) references: PartReferences,
    pub(crate) items: Vec<Id>,
    pub(crate) description: &'static [&'static str],
    pub(crate) report: StepReport,
}

impl PartDocument {
    pub(crate) fn file(&self, name: &str) -> String {
        render(&self.graph, &self.layout, name, self.description)
    }

    pub(crate) fn write(
        self,
        name: &str,
        seams: crate::Seams,
        out: &mut impl std::io::Write,
    ) -> Result<StepReport, Refusal> {
        let file = self.file(name);
        if seams == crate::Seams::Required {
            crate::seams::check(&file)?;
        }
        out.write_all(file.as_bytes())
            .map_err(|error| Refusal::new("step.io-error", error.to_string()))?;
        Ok(self.report)
    }

    pub(crate) fn append(
        self,
        graph: &mut Graph,
        layout: &mut Layout,
    ) -> (Product, Id, Vec<Id>, StepReport) {
        let offset = graph.append(self.graph);
        layout.append(self.layout, offset);
        let product = self.references.product;
        let product = Product {
            definition: offset.apply(product.definition),
            shape: offset.apply(product.shape),
            representation: offset.apply(product.representation),
            frame: offset.apply(product.frame),
        };
        let mut items = self.items;
        for item in &mut items {
            *item = offset.apply(*item);
        }
        (
            product,
            offset.apply(self.references.context),
            items,
            self.report,
        )
    }
}
