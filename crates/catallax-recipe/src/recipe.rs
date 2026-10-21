//! The recipe, its steps, its canonical text and its errors.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::ids::{ElementId, StepId};
use crate::operation::{Operation, ValueKind};
use crate::text::Canonical;

/// The format identifier every recipe carries.
pub const FORMAT: &str = "catallax-recipe/2";

/// A part: ordered steps and the step whose body is the part.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Recipe {
    /// Always [`FORMAT`].
    pub format: String,
    /// The part's name, for people.
    pub name: String,
    /// What the part is, for people.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    /// The steps, in order. A step reads only steps before it.
    pub steps: Vec<Step>,
    /// The step whose body is the part.
    pub result: StepId,
}

/// One operation with its id.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct Step {
    /// Unique in the recipe.
    pub id: StepId,
    /// The operation.
    #[serde(flatten)]
    pub op: Operation,
}

#[derive(Deserialize)]
struct FormatOnly {
    format: Option<String>,
}

impl Recipe {
    /// A recipe of the current format.
    pub fn new(name: impl Into<String>, steps: Vec<Step>, result: StepId) -> Self {
        Recipe {
            format: FORMAT.into(),
            name: name.into(),
            description: String::new(),
            steps,
            result,
        }
    }

    /// Reads a recipe and checks its structure. Accepts any JSON layout, key order and
    /// spelling of numbers; rejects unknown keys and every JSON number where an exact
    /// number belongs.
    pub fn parse(text: &str) -> Result<Recipe, RecipeError> {
        let parse = |e: serde_json::Error| RecipeError::Parse {
            message: e.to_string(),
        };
        let found = serde_json::from_str::<FormatOnly>(text)
            .map_err(parse)?
            .format;
        if found.as_deref() != Some(FORMAT) {
            return Err(RecipeError::Format {
                found: found.unwrap_or_default(),
            });
        }
        let recipe: Recipe = serde_json::from_str(text).map_err(parse)?;
        recipe.check()?;
        Ok(recipe)
    }

    /// The canonical text: the same recipe always gives the same bytes. Keys in schema
    /// order, absent optional keys left out, numbers in canonical form, edge and vertex
    /// names in canonical order, the layout of `text.rs`, and one final newline.
    pub fn to_text(&self) -> String {
        let mut out = Vec::new();
        let mut serializer = serde_json::Serializer::with_formatter(&mut out, Canonical::new());
        self.serialize(&mut serializer)
            .expect("a recipe always serialises");
        out.push(b'\n');
        String::from_utf8(out).expect("serde_json writes UTF-8")
    }

    /// Checks everything that needs no geometry: the format, unique step ids, inputs of
    /// the right kind from earlier steps, unique element ids, and names that can exist.
    pub fn check(&self) -> Result<(), RecipeError> {
        crate::check::check(self)
    }

    /// The step with this id.
    pub fn step(&self, id: &StepId) -> Option<&Step> {
        self.steps.iter().find(|s| &s.id == id)
    }
}

/// Why a recipe is not valid. These are recipe errors, not geometric refusals.
#[derive(Clone, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum RecipeError {
    /// The text is not a recipe of this format: bad JSON, a missing or unknown key, a
    /// JSON number where an exact number belongs, a bad id or a malformed name.
    Parse {
        /// What the reader found, with its line and column.
        message: String,
    },
    /// The format is not [`FORMAT`].
    Format {
        /// The format the text names, empty when it names none.
        found: String,
    },
    /// Two steps have the same id.
    DuplicateStep {
        /// The id.
        step: StepId,
    },
    /// A step reads a step that does not come before it.
    UnknownInput {
        /// The reading step.
        step: StepId,
        /// The id it reads.
        input: StepId,
    },
    /// A step reads a value of the wrong kind.
    WrongInput {
        /// The reading step.
        step: StepId,
        /// The step it reads.
        input: StepId,
        /// The kind it needs.
        expected: ValueKind,
        /// The kind the input makes.
        found: ValueKind,
    },
    /// Two elements of one step have the same id.
    DuplicateElement {
        /// The step.
        step: StepId,
        /// The id.
        element: ElementId,
    },
    /// A step names an element that another step does not have.
    UnknownElement {
        /// The naming step.
        step: StepId,
        /// The step that should have the element.
        of: StepId,
        /// The element id.
        element: ElementId,
    },
    /// A loft has fewer than two sections.
    LoftSections {
        /// The loft step.
        step: StepId,
    },
    /// A loft section does not name exactly one start segment in each loop.
    LoftStarts {
        /// The loft step.
        step: StepId,
        /// The section's sketch step.
        profile: StepId,
    },
    /// A name mentions a step that does not come before the step it is used in.
    UnknownOrigin {
        /// The step using the name.
        step: StepId,
        /// The step the name mentions.
        origin: StepId,
    },
    /// A reference matches the faces of a step that makes no faces.
    NoFaces {
        /// The step using the reference.
        step: StepId,
        /// The step the reference mentions.
        origin: StepId,
    },
    /// A name gives a step a role its operation never makes.
    WrongRole {
        /// The step using the name.
        step: StepId,
        /// The step the name mentions.
        origin: StepId,
        /// The role, such as `"side"` or `"hole.floor"`.
        role: String,
    },
    /// A name mentions a copy that its pattern or mirror does not make.
    NoSuchCopy {
        /// The step using the name.
        step: StepId,
        /// The pattern or mirror step.
        origin: StepId,
        /// The copy.
        instance: u32,
    },
    /// A two-distance chamfer selects every edge, so no face is first.
    AsymmetricAll {
        /// The chamfer step.
        step: StepId,
    },
    /// A file's SHA-256 is not 64 lowercase hexadecimal digits.
    FileHash {
        /// The import step.
        step: StepId,
    },
    /// The result is not a step that makes a body.
    Result {
        /// The result id.
        result: StepId,
    },
}

impl RecipeError {
    /// The stable code, `recipe.<reason>`.
    pub fn code(&self) -> &'static str {
        match self {
            RecipeError::Parse { .. } => "recipe.parse",
            RecipeError::Format { .. } => "recipe.format",
            RecipeError::DuplicateStep { .. } => "recipe.duplicate-step",
            RecipeError::UnknownInput { .. } => "recipe.unknown-input",
            RecipeError::WrongInput { .. } => "recipe.wrong-input",
            RecipeError::DuplicateElement { .. } => "recipe.duplicate-element",
            RecipeError::UnknownElement { .. } => "recipe.unknown-element",
            RecipeError::LoftSections { .. } => "recipe.loft-sections",
            RecipeError::LoftStarts { .. } => "recipe.loft-starts",
            RecipeError::UnknownOrigin { .. } => "recipe.unknown-origin",
            RecipeError::NoFaces { .. } => "recipe.no-faces",
            RecipeError::WrongRole { .. } => "recipe.wrong-role",
            RecipeError::NoSuchCopy { .. } => "recipe.no-such-copy",
            RecipeError::AsymmetricAll { .. } => "recipe.asymmetric-all",
            RecipeError::FileHash { .. } => "recipe.file-hash",
            RecipeError::Result { .. } => "recipe.result",
        }
    }
}

impl fmt::Display for RecipeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: ", self.code())?;
        match self {
            RecipeError::Parse { message } => write!(f, "{message}"),
            RecipeError::Format { found } => write!(f, "format is {found:?}, expected {FORMAT:?}"),
            RecipeError::DuplicateStep { step } => write!(f, "two steps have the id {step}"),
            RecipeError::UnknownInput { step, input } => {
                write!(f, "step {step} reads {input}, which is not an earlier step")
            }
            RecipeError::WrongInput {
                step,
                input,
                expected,
                found,
            } => write!(
                f,
                "step {step} needs a {} from {input}, which makes a {}",
                expected.name(),
                found.name()
            ),
            RecipeError::DuplicateElement { step, element } => {
                write!(f, "step {step} has two elements with the id {element}")
            }
            RecipeError::UnknownElement { step, of, element } => {
                write!(
                    f,
                    "step {step} names element {element}, which {of} does not have"
                )
            }
            RecipeError::LoftSections { step } => {
                write!(f, "loft {step} has fewer than two sections")
            }
            RecipeError::LoftStarts { step, profile } => write!(
                f,
                "loft {step}: section {profile} must name one start segment in each loop"
            ),
            RecipeError::UnknownOrigin { step, origin } => {
                write!(
                    f,
                    "step {step} names a face of {origin}, which is not an earlier step"
                )
            }
            RecipeError::NoFaces { step, origin } => {
                write!(
                    f,
                    "step {step} names faces of {origin}, which makes no faces"
                )
            }
            RecipeError::WrongRole { step, origin, role } => {
                write!(
                    f,
                    "step {step} names a {role} face of {origin}, which makes none"
                )
            }
            RecipeError::NoSuchCopy {
                step,
                origin,
                instance,
            } => write!(
                f,
                "step {step} names copy {instance} of {origin}, which makes none"
            ),
            RecipeError::AsymmetricAll { step } => write!(
                f,
                "chamfer {step} has two distances, so it must name its edges"
            ),
            RecipeError::FileHash { step } => {
                write!(
                    f,
                    "step {step}: sha256 must be 64 lowercase hexadecimal digits"
                )
            }
            RecipeError::Result { result } => {
                write!(f, "result {result} is not a step that makes a body")
            }
        }
    }
}

impl std::error::Error for RecipeError {}
