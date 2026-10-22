//! Build a SolidDraft recipe and export its named display OBJ.
use catallax_app::build_soliddraft;
use catallax_mesh::MeshOptions;
use catallax_numbers::Decider;
use std::{env, fs, io};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err(io::Error::other("usage: build SOLIDDRAFT_RECIPE.json DISPLAY.obj").into());
    }
    let json = fs::read_to_string(&args[0])?;
    let result = build_soliddraft(&json, &MeshOptions::default(), &mut Decider::default())?;
    result
        .body
        .display
        .value
        .write_obj(&mut fs::File::create(&args[1])?)?;
    println!(
        "{} faces, {} edges, {} triangles, {} display marks",
        result.body.names.faces.len(),
        result.body.names.edges.len(),
        result.body.display.value.triangles.len(),
        result.body.display.marks.len()
    );
    for (face, origins) in &result.body.names.faces {
        println!("face_{} {}", face.index(), serde_json::to_string(origins)?);
    }
    Ok(())
}
