//! The final lattice stage of a finished body.
//!
//! A [`Stage`] is built from the body's recipe, the lattice parameters and the
//! user's output settings. It adds one lattice step that reads the body's result
//! and makes that step the result. The stage owns no cache: editing the body means
//! calling [`Stage::new`] again with the edited recipe, so nothing computed from
//! the earlier body can be reused.

use crate::{Plan, printing};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{
    Grade, LatticeCell, Number, Operation, Pose, Recipe, RecipeError, Step, StepId,
};

/// The lattice the user chose. These survive returning before the lattice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LatticeSpec {
    /// The cell.
    pub cell: LatticeCell,
    /// Edge length of the cubic cell.
    pub size: Number,
    /// Places the cell grid.
    pub frame: Pose,
    /// Beam radius.
    pub radius: Grade,
}

/// The user's output settings. Only these two relax for a lattice export.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Settings {
    /// Two-sided distance limit of the output in millimetres.
    pub accuracy: Rational,
    /// Smallest feature in millimetres the output is promised to keep.
    pub minimum_detail: Rational,
    /// Maximum octree cells visited across all refinement attempts.
    pub max_cells: usize,
    /// Maximum mesh vertices.
    pub max_vertices: usize,
}
impl Settings {
    /// The print controls these settings select.
    pub fn print_options(&self) -> printing::Options {
        printing::Options {
            tolerance: self.accuracy.clone(),
            max_cells: self.max_cells,
            max_vertices: self.max_vertices,
        }
    }
}

/// What the export checked about its mesh. `Stage::export` returns an [`Output`] only
/// when every flag is true.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Verdict {
    /// Every edge has exactly two facets of opposite sense and every vertex link is a
    /// cycle.
    pub closed: bool,
    /// The surface is a closed 2-manifold on a conforming tetrahedral grid whose cells
    /// each hold at most one patch, so facets meet only along shared edges.
    pub non_self_intersecting: bool,
    /// The mesh has exactly one edge-connected component.
    pub connected: bool,
    /// The number of components found.
    pub components: usize,
    /// Proven two-sided distance bound of the output, in millimetres.
    pub max_deviation: Rational,
}

/// A certified lattice mesh with its verdict.
#[derive(Clone, Debug)]
pub struct Output {
    /// The mesh; write it with `Mesh::write_stl` or `Mesh::write_obj`.
    pub mesh: printing::Mesh,
    /// The checks that passed.
    pub verdict: Verdict,
}

/// A body with its lattice stage compiled.
#[derive(Clone, Debug)]
pub struct Stage {
    body: Recipe,
    spec: LatticeSpec,
    settings: Settings,
    recipe: Recipe,
    lattice: StepId,
    plan: Plan,
}

impl Stage {
    /// Builds the stage over `body`, whose result is the finished body.
    ///
    /// Refusals: `lattice.cad-after-lattice` when the body already contains a
    /// lattice, `lattice.settings` for non-positive accuracy, detail or limits,
    /// `lattice.impossible` for impossible lattice parameters, and the field's own
    /// codes when the body cannot be compiled.
    pub fn new(
        body: &Recipe,
        spec: LatticeSpec,
        settings: Settings,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if body
            .steps
            .iter()
            .any(|s| matches!(s.op, Operation::Lattice { .. }))
        {
            return Err(Refusal::new(
                "lattice.cad-after-lattice",
                "the body already contains a lattice; return before it to edit",
            ));
        }
        if settings.accuracy.sign() != Sign::Positive
            || settings.minimum_detail.sign() != Sign::Positive
            || settings.max_cells == 0
            || settings.max_vertices < 4
        {
            return Err(Refusal::new(
                "lattice.settings",
                "positive accuracy, minimum detail and cell limit, and at least four vertices required",
            ));
        }
        catallax_features::lattice::validate(&spec.cell, spec.size.value(), &spec.radius)?;
        let lattice = (0..)
            .map(|n| match n {
                0 => "lattice".to_string(),
                n => format!("lattice-{n}"),
            })
            .find_map(|name| {
                let id = StepId::new(&name).ok()?;
                body.step(&id).is_none().then_some(id)
            })
            .expect("an unused step id exists");
        let mut recipe = body.clone();
        recipe.steps.push(Step {
            id: lattice.clone(),
            op: Operation::Lattice {
                region: body.result.clone(),
                cell: spec.cell.clone(),
                size: spec.size.clone(),
                frame: spec.frame.clone(),
                radius: spec.radius.clone(),
            },
        });
        recipe.result = lattice.clone();
        let plan = Plan::compile_checked(&recipe, dec)?;
        Ok(Self {
            body: body.clone(),
            spec,
            settings,
            recipe,
            lattice,
            plan,
        })
    }

    /// The body recipe plus the lattice step, which is its result.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }
    /// The id of the lattice step.
    pub fn lattice_step(&self) -> &StepId {
        &self.lattice
    }
    /// The compiled field of the staged recipe.
    pub fn plan(&self) -> &Plan {
        &self.plan
    }
    /// The output settings.
    pub fn settings(&self) -> &Settings {
        &self.settings
    }

    /// Builds the certified mesh at the settings' accuracy. The mesh is closed, free of
    /// self-intersection and connected, or the call refuses; nothing partial is
    /// returned. Refusals: `field.print-disconnected`, `field.print-unresolved` (also
    /// when the minimum detail forbids the refinement the boundary needs),
    /// `field.print-topology`, `field.print-budget`, `field.print-grid`,
    /// `field.print-vertices`, and the decider's deadline and cancellation outcomes.
    pub fn export(&self, dec: &mut Decider) -> Result<Output, Refusal> {
        let floor = &self.settings.minimum_detail * &Rational::from_ratio(1, 4).unwrap();
        let mesh = printing::mesh_plan_floor(
            &self.plan,
            &self.settings.print_options(),
            Some(&floor),
            dec,
        )?;
        let components = mesh.component_count();
        if components != 1 {
            return Err(Refusal::new(
                "field.print-disconnected",
                format!("the print output has {components} connected components, expected one"),
            ));
        }
        let verdict = Verdict {
            closed: true,
            non_self_intersecting: true,
            connected: true,
            components,
            max_deviation: mesh.max_deviation.clone(),
        };
        Ok(Output { mesh, verdict })
    }

    /// Leaves the stage: the unchanged body recipe, with the lattice parameters and
    /// settings to put back after the body is edited.
    pub fn return_before_lattice(self) -> (Recipe, LatticeSpec, Settings) {
        (self.body, self.spec, self.settings)
    }
}

/// Appends nothing and reports whether a CAD step may follow `recipe`'s result.
/// A lattice anywhere in a history makes every later step a CAD-after-lattice.
pub fn check_no_cad_after_lattice(recipe: &Recipe) -> Result<(), Refusal> {
    recipe.check_lattice_terminal().map_err(|e| match e {
        RecipeError::LatticeNotTerminal { step, lattice } => Refusal::new(
            "lattice.cad-after-lattice",
            format!("step {step} reads lattice {lattice}; return before the lattice to edit"),
        ),
        other => Refusal::new("field.recipe", other.to_string()),
    })
}
