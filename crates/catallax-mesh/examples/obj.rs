//! `cargo run --release -p catallax-mesh --example obj -- sphere sphere.obj 0.01`

use std::{
    env,
    fs::File,
    io::{BufWriter, Write},
};

use catallax_boolean::{BoolOp, boolean};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{self as topology, Model};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if !(2..=3).contains(&args.len()) {
        return Err("usage: obj <block|cylinder|cone|frustum|sphere|torus|drill|lens|cross> <file.obj> [chord-mm]".into());
    }
    let mut model = Model::new();
    let mut dec = Decider::default();
    let pose = Pose::identity();
    let shifted = Pose::new(
        Rotation::identity(),
        (RVector3::new(1.into(), 0.into(), 0.into())).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let body = match args[0].as_str() {
        "block" => topology::block(&mut model, [2.into(), 3.into(), 5.into()], &pose)?,
        "cylinder" => topology::cylinder(&mut model, 2.into(), 5.into(), &pose)?,
        "cone" => topology::cone(&mut model, 3.into(), 0.into(), 4.into(), &pose)?,
        "frustum" => topology::cone(&mut model, 3.into(), 1.into(), 4.into(), &pose)?,
        "sphere" => topology::sphere(&mut model, 2.into(), &pose)?,
        "torus" => topology::torus(&mut model, 3.into(), 1.into(), &pose)?,
        "lens" => {
            let a = topology::sphere(&mut model, 2.into(), &pose)?;
            let b = topology::sphere(&mut model, 2.into(), &shifted)?;
            boolean(&mut model, a, b, BoolOp::Intersect, &mut dec)
                .map_err(|e| e.refusal)?
                .body
        }
        "cross" => {
            let pose = Pose::new(
                Rotation::identity(),
                (RVector3::new(0.into(), 0.into(), (-3).into())).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap();
            let a = topology::cylinder(&mut model, 1.into(), 6.into(), &pose)?;
            let pose = Pose::new(
                Rotation::identity(),
                (RVector3::new(2.into(), 0.into(), 0.into())).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap();
            let b = topology::sphere(&mut model, 2.into(), &pose)?;
            boolean(&mut model, a, b, BoolOp::Intersect, &mut dec)
                .map_err(|e| e.refusal)?
                .body
        }
        "drill" => {
            let a = topology::block(&mut model, [6.into(), 6.into(), 4.into()], &pose)?;
            let pose = Pose::new(
                Rotation::identity(),
                (RVector3::new(3.into(), 3.into(), (-1).into())).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap();
            let b = topology::cylinder(&mut model, 1.into(), 6.into(), &pose)?;
            boolean(&mut model, a, b, BoolOp::Subtract, &mut dec)
                .map_err(|e| e.refusal)?
                .body
        }
        _ => return Err("unknown example shape".into()),
    };
    let options = MeshOptions {
        chord_tolerance: args
            .get(2)
            .map(|s| Rational::parse(s))
            .transpose()?
            .unwrap_or_else(|| MeshOptions::default().chord_tolerance),
        ..MeshOptions::default()
    };
    let answer = mesh_body(&model, body, &options, &mut dec)?;
    let mut out = BufWriter::new(File::create(&args[1])?);
    answer.value.write_obj(&mut out)?;
    out.flush()?;
    println!(
        "{} vertices, {} triangles, {} faces, {} shared edges",
        answer.value.positions.len(),
        answer.value.triangles.len(),
        answer.value.faces.len(),
        answer.value.edges.len()
    );
    for face in &answer.value.faces {
        println!(
            "face {}: proven chord error <= {} mm",
            face.face.index(),
            face.max_chord_error.to_f64_nearest()
        );
    }
    Ok(())
}
