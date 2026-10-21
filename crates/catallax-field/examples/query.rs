//! Query a recipe at one exact point: `query part.json x y z`.

use catallax_field::{Batch, Plan, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 4 {
        return Err("usage: query part.json x y z".into());
    }
    let recipe = Recipe::parse(&std::fs::read_to_string(&args[0])?)?;
    let plan = Plan::compile(&recipe)?;
    let point = RPoint3::new(
        Rational::parse(&args[1])?,
        Rational::parse(&args[2])?,
        Rational::parse(&args[3])?,
    );
    let mut storage = Workspace::default();
    println!("proven: {:?}", plan.proven(&point, &mut storage));
    println!(
        "fast: {:?}",
        plan.fast(
            [
                point.x.to_f64_nearest(),
                point.y.to_f64_nearest(),
                point.z.to_f64_nearest()
            ],
            &mut storage
        )
    );
    for face in plan.faces() {
        let id = plan.select(&face.origin).expect("registered face");
        println!(
            "{}: {:?}",
            face.origin.key(),
            plan.face_query(id, &point, &mut storage)?
        );
    }
    let mut batch = Batch::default();
    assert_eq!(
        plan.proven_batch(std::slice::from_ref(&point), &mut batch)[0],
        plan.proven(&point, &mut storage)
    );
    println!("batch: stable scalar/batch result");
    Ok(())
}
