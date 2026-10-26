use super::{Options, qualification};
use catallax_features::{BuildError, Built};
use catallax_field::printing::Mesh;
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::Recipe;
use std::collections::BTreeMap;

fn build_error(error: BuildError) -> Refusal {
    match error {
        BuildError::Geometry { step, refusal } => {
            Refusal::new(refusal.code, format!("step {step}: {}", refusal.detail))
        }
        BuildError::Recipe(error) => Refusal::new("cli.invalid-recipe", error.to_string()),
    }
}
pub(super) struct Output {
    pub mesh: Mesh,
    pub model_bound: Rational,
    pub tessellation_bound: Rational,
    pub rounding_bound: Rational,
}
pub(super) fn mesh(
    recipe: &Recipe,
    options: &Options,
    certified: bool,
    dec: &mut Decider,
) -> Result<Output, Refusal> {
    if super::native_field(recipe) {
        return Err(Refusal::new(
            "print.brep-field-definition",
            "a field envelope is not a B-rep print definition; request the field route",
        ));
    }
    let certificate = if certified {
        Some(catallax_features::evaluate_certified(recipe, dec).map_err(build_error)?)
    } else {
        None
    };
    let exact;
    let built: &Built = if let Some(c) = &certificate {
        c.verify(dec)?;
        c.built()
    } else {
        exact = catallax_features::evaluate(recipe, dec).map_err(build_error)?;
        &exact
    };
    if !built.model.check(built.result).is_empty() {
        return Err(Refusal::new(
            "print.brep-topology",
            "the model topology is not proven closed",
        ));
    }
    let model_bound = if let Some(c) = &certificate {
        c.result_certificate()
            .record()
            .uniform_bound
            .as_ref()
            .ok_or_else(|| {
                Refusal::new(
                    "print.brep-model-bound",
                    "no complete model-to-definition bound",
                )
            })?
            .upper
            .clone()
    } else {
        if !built.model.body(built.result).marks.is_empty() {
            return Err(Refusal::new(
                "print.brep-model-bound",
                "marked geometry requires its complete geometry-bound certificate, not a mark interpreted as a distance",
            ));
        }
        Rational::zero()
    };
    if model_bound >= options.tolerance {
        return Err(Refusal::new(
            "print.brep-accuracy",
            "model-to-definition bound consumes the requested print accuracy",
        ));
    }
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions {
            chord_tolerance: (&options.tolerance - &model_bound).mul_pow2(-1),
            max_vertices: options.max_vertices,
            ..MeshOptions::default()
        },
        dec,
    )?
    .value;
    let face_bound = mesh
        .faces
        .iter()
        .map(|face| &face.max_chord_error)
        .max()
        .cloned()
        .unwrap_or_else(Rational::zero);
    let trim_bound = mesh
        .edges
        .iter()
        .map(|edge| &edge.max_chord_error)
        .max()
        .cloned()
        .unwrap_or_else(Rational::zero);
    // Include the finite-domain replacement, not only distance to the infinite
    // surface carrier. Adding the two certified bounds is conservative where
    // their endpoint-rounding contributions overlap.
    let tessellation_bound = &face_bound + &trim_bound;
    let mut rounding_bound = Rational::zero();
    let mut positions = Vec::with_capacity(mesh.positions.len());
    for point in &mesh.positions {
        dec.check_deadline()?;
        let rounded = point.map(|x| x as f32);
        if !rounded.iter().all(|x| x.is_finite()) {
            return Err(Refusal::new(
                "print.brep-encoding-range",
                "coordinate exceeds finite binary32 range",
            ));
        }
        // L1 dominates Euclidean displacement. Barycentric interpolation gives
        // the same two-sided bound for every point of each corresponding facet.
        let error = (0..3).fold(Rational::zero(), |sum, i| {
            &sum + &(&Rational::from_f64_exact(point[i]).unwrap()
                - &Rational::from_f64_exact(f64::from(rounded[i])).unwrap())
                .abs()
        });
        rounding_bound = rounding_bound.max(error);
        positions.push(rounded);
    }
    let bound = &(&model_bound + &tessellation_bound) + &rounding_bound;
    if bound > options.tolerance {
        return Err(Refusal::new(
            "print.brep-accuracy",
            "model, finite-trim tessellation and binary32 rounding exceed requested accuracy",
        ));
    }
    let triangles: Vec<_> = mesh.triangles.iter().map(|t| t.vertices).collect();
    qualification::topology(&triangles, positions.len(), dec)?;
    qualification::motion(&mesh.positions, &positions, &triangles, dec)?;
    // Preserve the producer's material sense. Voids must have negative signed
    // volume; never independently flip a connected shell to positive volume.
    let mut face_shell = BTreeMap::new();
    let mut shells = Vec::new();
    let mut owners = Vec::new();
    let mut shell_facets = Vec::<Vec<[usize; 3]>>::new();
    for (owner, lump) in built.model.body(built.result).lumps.iter().enumerate() {
        for (i, shell) in built.model.lump(*lump).shells.iter().enumerate() {
            let index = shells.len();
            shells.push((i == 0, Rational::zero()));
            owners.push((owner, i == 0));
            shell_facets.push(Vec::new());
            for face in &built.model.shell(*shell).faces {
                face_shell.insert(*face, index);
            }
        }
    }
    let q: Vec<_> = positions
        .iter()
        .map(|p| p.map(|x| Rational::from_f64_exact(f64::from(x)).unwrap()))
        .collect();
    for t in &mesh.triangles {
        dec.check_deadline()?;
        let [a, b, c] = t.vertices.map(|i| &q[i]);
        let volume = (0..3).fold(Rational::zero(), |sum, i| {
            &sum + &(&a[i]
                * &(&(&b[(i + 1) % 3] * &c[(i + 2) % 3]) - &(&b[(i + 2) % 3] * &c[(i + 1) % 3])))
        });
        shell_facets[face_shell[&t.face]].push(t.vertices);
        let entry = &mut shells[face_shell[&t.face]];
        entry.1 = &entry.1 + &volume;
    }
    if shells.iter().any(|(outer, v)| {
        if *outer {
            v <= &Rational::zero()
        } else {
            v >= &Rational::zero()
        }
    }) {
        return Err(Refusal::new(
            "print.brep-orientation",
            "encoded outer/cavity shell orientation disagrees with material sense",
        ));
    }
    super::nesting::check(&q, &shell_facets, &owners, dec)?;
    Ok(Output {
        mesh: Mesh {
            positions,
            triangles,
            max_deviation: bound,
            cells: 0,
        },
        model_bound,
        tessellation_bound,
        rounding_bound,
    })
}
