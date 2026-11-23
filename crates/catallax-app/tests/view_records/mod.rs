use catallax_app::views::{HiddenLineView, stroke_drawing};
use catallax_numbers::{Answer, Decider};
use catallax_shapes::RVector3;

/// Emit the ordered drawing curves and their marks for a same-platform
/// comparison of the baseline and candidate test binaries.
pub fn emit(answer: &Answer<HiddenLineView>, direction: &RVector3) {
    if std::env::var_os("CATALLAX_RECORD_VIEWS").is_none() {
        return;
    }
    let mut dec = Decider::default();
    let drawn: Vec<_> = answer
        .value
        .strokes
        .iter()
        .map(|stroke| {
            (
                stroke.source,
                stroke.visible,
                stroke.smooth,
                stroke.seam,
                stroke_drawing(&answer.value, stroke, &mut dec),
            )
        })
        .collect();
    let test = std::thread::current();
    eprintln!(
        "VIEW_RECORD {} {direction:?} {:?} {drawn:?}",
        test.name().unwrap(),
        answer.marks,
    );
}
