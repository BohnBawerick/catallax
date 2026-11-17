//! Immutable application retention, separate from geometry construction.

use crate::{EvaluatedApplicationBody, NumericalBody};
use catallax_recipe::{Origin, Recipe};
use std::{collections::BTreeMap, sync::Arc};

/// Exact evaluation or separately admitted numerical geometry.
///
/// A numerical result has no exact model or reusable exact checkpoint.
#[derive(Clone, Debug)]
pub enum EvaluatedApplicationOutcome {
    /// Checked exact geometry and its continuation checkpoint.
    Exact(EvaluatedApplicationBody),
    /// Numerical geometry with declared tolerance and scoped evidence.
    Numerical(Arc<NumericalApplicationBody>),
}

/// Retained numerical geometry, authored controls and canonical face origins.
///
/// This type adds no qualification to its body. Numerical consumers use
/// [`Self::body`]; exact exporters and checkpoints cannot consume this wrapper.
#[derive(Debug)]
pub struct NumericalApplicationBody {
    body: NumericalBody,
    recipe: Arc<Recipe>,
    face_origins: BTreeMap<String, Vec<Origin>>,
}

impl NumericalApplicationBody {
    /// Retains an already admitted body without rebuilding or discarding history.
    ///
    /// The producer supplies its source recipe and complete origins keyed by
    /// `NumericalFace.key`, not exact-model picking handles. This constructor
    /// does not construct geometry or verify that provenance against the recipe.
    pub fn new(
        body: NumericalBody,
        recipe: Arc<Recipe>,
        face_origins: BTreeMap<String, Vec<Origin>>,
    ) -> Self {
        Self {
            body,
            recipe,
            face_origins,
        }
    }

    /// Geometry, native trims, declared tolerance, checks and repair history.
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }

    /// Authored controls retained for replay, never an exact checkpoint.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }

    /// Canonical origins keyed by the numerical definition's stable face keys.
    pub fn face_origins(&self) -> &BTreeMap<String, Vec<Origin>> {
        &self.face_origins
    }
}
