//! Bounded M23 reconstruction and unproved two-way fidelity demonstration.
#[path = "reconstruction/fixture.rs"]
mod fixture;
use catallax_field::reconstruction::{Sampling, filled_envelope_control, reconstruct};
use catallax_numbers::{Decider, Rational};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 2 && (args[0] == "--write-fixture" || args[0] == "--write-cone-fixture") {
        std::fs::write(
            &args[1],
            fixture::fixture(if args[0] == "--write-cone-fixture" {
                3
            } else {
                2
            }),
        )?;
        return Ok(());
    }
    if args.len() > 1 {
        return Err("usage: reconstruct_mesh [mesh.obj] | --write-fixture output.obj | --write-cone-fixture output.obj".into());
    }
    let inputs = if let Some(path) = args.first() {
        vec![("supplied mesh", std::fs::read(path)?)]
    } else {
        vec![
            ("cylinder control", fixture::fixture(2)),
            ("cone control", fixture::fixture(3)),
        ]
    };
    for (name, bytes) in inputs {
        let mut dec = Decider::default();
        let part = reconstruct(&bytes, &Rational::one(), &mut dec)?;
        println!(
            "{name}: source={} regions={:?} centre={:?} radii={:?} B-rep_faces={}",
            part.source_sha256(),
            part.region_triangles(),
            part.centre(),
            part.radii(),
            part.brep().model.faces_of_body(part.brep().result).len()
        );
        let sampling = Sampling {
            triangle_steps: 4,
            surface_steps: 8,
            angles: 128,
        };
        println!("UNPROVED sampled fidelity mm: {:?}", part.fidelity());
        println!(
            "UNPROVED adverse filled-envelope control mm: {:?}",
            filled_envelope_control(&bytes, &part, sampling, &mut dec)?
        );
        println!("inferred definition: {}", part.recipe().to_text());
    }
    Ok(())
}
