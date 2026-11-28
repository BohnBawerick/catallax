//! Prints the manufacturing information of STEP files: what was read, by
//! kind, and every item that was left out with the reason.
//!
//! ```sh
//! cargo run -p catallax-import --example step_pmi -- part.step
//! ```

use catallax_import::read_step_pmi;
use catallax_step::pmi::ItemKind;

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: step_pmi <file.step>...");
        std::process::exit(2);
    }
    for path in paths {
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) => {
                println!("{path}: {error}");
                continue;
            }
        };
        let read = match read_step_pmi(&bytes) {
            Ok(read) => read,
            Err(refusal) => {
                println!("{path}: not read: {refusal:?}");
                continue;
            }
        };
        let count = |wanted: fn(&ItemKind<String>) -> bool| {
            read.pmi
                .items
                .iter()
                .filter(|item| wanted(&item.kind))
                .count()
        };
        println!(
            "{path}: {} datums, {} sizes, {} distances, {} basic angles, {} tolerances, {} textures, {} skipped",
            count(|kind| matches!(kind, ItemKind::Datum(_))),
            count(|kind| matches!(
                kind,
                ItemKind::Size(_)
                    | ItemKind::Width { .. }
                    | ItemKind::Chamfer { .. }
                    | ItemKind::EdgeBreak { .. }
            )),
            count(|kind| matches!(kind, ItemKind::Location(_))),
            count(|kind| matches!(kind, ItemKind::BasicAngle { .. })),
            count(|kind| matches!(kind, ItemKind::Tolerance(_))),
            count(|kind| matches!(kind, ItemKind::Texture(_))),
            read.skipped.len(),
        );
        let part = &read.pmi.part;
        for (name, value) in [
            ("general tolerance", &part.general_tolerance),
            ("material", &part.material),
        ] {
            if let Some(value) = value {
                println!("  {name}: {value}");
            }
        }
        for note in &part.notes {
            println!("  note: {note}");
        }
        for item in &read.pmi.items {
            println!("  {}: {:?}", item.id, item.kind);
        }
        for view in &read.pmi.views {
            println!(
                "  view {}: {} items, eye {:?}, looking {:?}",
                view.name,
                view.items.len(),
                view.camera
                    .eye
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>(),
                view.camera
                    .direction
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
            );
        }
        println!("  {} graphical-only callouts", read.graphics.len());
        for graphic in &read.graphics {
            println!(
                "  graphic {} {:?} ({}): {:?}, presents {:?}, in {:?}",
                graphic.entity,
                graphic.name,
                graphic.kind,
                graphic.text,
                graphic.presents,
                graphic.views
            );
        }
        for skipped in &read.skipped {
            println!(
                "  skipped {} {}: {}",
                skipped.entity, skipped.kind, skipped.reason
            );
        }
    }
}
