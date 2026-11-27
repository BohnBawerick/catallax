//! Saved views: a named camera and the items it shows (CAx-IF PMI
//! recommended practices 4.1, section 9.4).

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
