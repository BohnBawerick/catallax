//! Saved views: a named camera and the items it shows (CAx-IF PMI
//! recommended practices 4.1, section 9.4).

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::Rational;

/// A saved view of the part's manufacturing information.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SavedView {
    /// The view's name, unique within the part.
    pub name: String,
    /// Where the view is seen from.
    pub camera: Camera,
    /// The names of the items whose callouts the view shows, by
    /// [`super::Item::id`]. A name is never rebound to another item.
    pub items: Vec<String>,
}

/// A camera as `camera_model_d3` and `view_volume` state it, in part
/// millimetres (section 9.4.2.6).
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Camera {
    /// The eye: the projection point.
    pub eye: [Rational; 3],
    /// The direction of view, from the eye towards the target. Not a length.
    pub direction: [Rational; 3],
    /// The screen's up direction, square to `direction`. Not a length.
    pub up: [Rational; 3],
    /// The distance from the eye to the view plane through the target.
    pub distance: Rational,
    /// The width and height of the view window on that plane.
    pub window: [Rational; 2],
    /// A central (perspective) rather than a parallel projection.
    pub perspective: bool,
}

/// Unique view names in order: the first of each authored name keeps it, and
/// a repeat takes `name (n)` with the lowest free `n` from 2 that no authored
/// or earlier name holds. Each base keeps its next suffix.
#[doc(hidden)]
pub fn unique_view_names<'n>(
    bases: impl Iterator<Item = (&'n str, &'n str)> + Clone,
) -> (Vec<String>, usize) {
    let mut reserved: BTreeSet<String> = bases.clone().map(|(name, _)| name.to_owned()).collect();
    let mut seen = BTreeSet::new();
    let mut next = BTreeMap::<&str, usize>::new();
    let mut names = Vec::new();
    let mut tried = 0;
    for (base, suffix_base) in bases {
        if seen.insert(base) {
            names.push(base.to_owned());
            continue;
        }
        let suffix = next.entry(suffix_base).or_insert(2);
        let name = loop {
            let name = format!("{suffix_base} ({suffix})");
            *suffix += 1;
            tried += 1;
            if reserved.insert(name.clone()) {
                break name;
            }
        };
        names.push(name);
    }
    (names, tried)
}
