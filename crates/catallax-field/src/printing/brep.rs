//! Qualified print meshes from an already-built B-rep. No recipe replay or measurement.

use super::{Mesh, Options};
use catallax_features::{Built, CertifiedBuilt};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, Refusal};
use std::collections::BTreeMap;

mod nesting;
mod qualification;

/// Geometry and contributions to its complete file-distance bound.
#[derive(Clone, Debug)]
pub struct Output {
    /// Qualified binary32 geometry, with material-oriented cavity shells.
    mesh: Mesh,
    /// Complete model-to-definition bound.
    model_bound: Rational,
    /// Finite-face and finite-trim tessellation bound.
    tessellation_bound: Rational,
    /// Exact binary64-to-binary32 displacement bound.
    rounding_bound: Rational,
    options: Options,
    field_route: bool,
}

/// Qualifies the exact supplied body, without rebuilding it or using its display.
/// Marked geometry requires [`mesh_certified`], never an assumed zero model bound.
pub fn mesh(built: &Built, options: &Options, dec: &mut Decider) -> Result<Output, Refusal> {
    check_body(built, options)?;
    if !built.model.body(built.result).marks.is_empty() {
        return Err(Refusal::new(
            "print.brep-model-bound",
            "marked geometry requires its complete geometry-bound certificate, not a mark interpreted as a distance",
        ));
    }
    qualify(built, Rational::zero(), options, dec)
}

/// Qualifies the supplied certified body, verifying its geometry-bound certificate.
pub fn mesh_certified(
    certified: &CertifiedBuilt,
    options: &Options,
    dec: &mut Decider,
) -> Result<Output, Refusal> {
    certified.verify(dec)?;
    check_body(certified.built(), options)?;
    let model_bound = certified
        .result_certificate()
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
        .clone();
    qualify(certified.built(), model_bound, options, dec)
}

fn check_body(built: &Built, options: &Options) -> Result<(), Refusal> {
    if options.tolerance <= Rational::zero() || options.max_vertices < 4 {
        return Err(Refusal::new(
            "print.invalid-options",
            "positive accuracy and at least four vertices required",
        ));
    }
    if !built.model.check(built.result).is_empty() {
        return Err(Refusal::new(
            "print.brep-topology",
            "the model topology is not proven closed",
        ));
    }
    Ok(())
}

fn qualify(
    built: &Built,
    model_bound: Rational,
    options: &Options,
    dec: &mut Decider,
) -> Result<Output, Refusal> {
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
    nesting::check(&q, &shell_facets, &owners, dec)?;
    Ok(Output {
        options: options.clone(),
        field_route: false,
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

fn upward(bound: &Rational) -> Result<(f64, Rational), Refusal> {
    let mut value = bound.to_f64_nearest();
    if !value.is_finite() {
        return Err(Refusal::new(
            "print.metadata-range",
            "certificate bound exceeds binary64 metadata range",
        ));
    }
    if Rational::from_f64_exact(value).unwrap() < *bound {
        value = value.next_up();
    }
    Rational::from_f64_exact(value)
        .map(|q| (value, q))
        .ok_or_else(|| {
            Refusal::new(
                "print.metadata-range",
                "certificate bound exceeds binary64 metadata range",
            )
        })
}
fn header(output: &mut Output) -> Result<[u8; 80], Refusal> {
    let options = &output.options;
    let field_route = output.field_route;
    let accuracy = options.tolerance.to_string();
    let (n, d) = accuracy.split_once('/').unwrap_or((&accuracy, "1"));
    let n: i64 = n.parse().map_err(|_| Refusal::new("print.metadata-range","accuracy numerator does not fit STL metadata; use an exactly representable requested accuracy"))?;
    let d: u64 = d.parse().map_err(|_| {
        Refusal::new(
            "print.metadata-range",
            "accuracy denominator does not fit STL metadata; requested accuracy was not truncated",
        )
    })?;
    let (u, model) = upward(&output.model_bound)?;
    let (v, tessellation) = upward(&output.tessellation_bound)?;
    let (r, rounding) = upward(&output.rounding_bound)?;
    let (w, total) = upward(&(&(&model + &tessellation) + &rounding))?;
    if total > options.tolerance {
        return Err(Refusal::new(
            "print.metadata-accuracy",
            "lossless conservative header bounds exceed requested accuracy",
        ));
    }
    output.model_bound = model;
    output.tessellation_bound = tessellation;
    output.rounding_bound = rounding;
    output.mesh.max_deviation = total;
    let mut bytes = [0u8; 80];
    bytes[..8].copy_from_slice(b"CXPRINT1");
    bytes[8] = if !field_route { 1 } else { 2 };
    bytes[9] = 1; // millimetres
    bytes[10] = 1; // binary32 little endian
    bytes[11] = if field_route { 1 } else { 0 }; // cell cap applicability
    bytes[12..20].copy_from_slice(&n.to_le_bytes());
    bytes[20..28].copy_from_slice(&d.to_le_bytes());
    bytes[28..36].copy_from_slice(&(options.max_cells as u64).to_le_bytes());
    bytes[36..44].copy_from_slice(&(options.max_vertices as u64).to_le_bytes());
    for (offset, value) in [(44, u), (52, v), (60, r), (68, w)] {
        bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    }
    Ok(bytes)
}

impl Output {
    /// The selected route. A refusal never switches this route.
    pub fn route(&self) -> &'static str {
        if self.field_route { "field" } else { "brep" }
    }
    /// Creates the metadata-bearing STL in memory. No destination is touched on refusal.
    /// Returns the rounded certificate with the bytes so both describe the same encoding.
    pub fn stl(mut self) -> Result<(Vec<u8>, Self), Refusal> {
        let header = header(&mut self)?;
        let mut bytes = Vec::new();
        self.mesh
            .write_stl(&mut bytes)
            .map_err(|e| Refusal::new("print.encoding", e.to_string()))?;
        bytes[..80].copy_from_slice(&header);
        Ok((bytes, self))
    }
}

impl Output {
    /// Exact controls used for this proof, retained independently of caller edits.
    pub fn options(&self) -> &Options {
        &self.options
    }
    /// Read-only certified binary32 geometry; no display mesh is substituted.
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }
    /// Complete model-to-definition contribution.
    pub fn model_bound(&self) -> &Rational {
        &self.model_bound
    }
    /// Finite-face and trim contribution.
    pub fn tessellation_bound(&self) -> &Rational {
        &self.tessellation_bound
    }
    /// Binary32 encoding contribution.
    pub fn rounding_bound(&self) -> &Rational {
        &self.rounding_bound
    }
    /// The existing print header; all contributions are conservatively rounded together.
    pub fn header(&mut self) -> Result<[u8; 80], Refusal> {
        header(self)
    }
}

/// Qualifies the native-field route with the same retained output controls.
/// This never retries a failed B-rep request as a field request.
pub fn field(plan: &crate::Plan, options: &Options, dec: &mut Decider) -> Result<Output, Refusal> {
    let mesh = super::mesh_plan(plan, options, dec)?;
    let tessellation_bound = mesh.max_deviation.clone();
    Ok(Output {
        mesh,
        model_bound: Rational::zero(),
        tessellation_bound,
        rounding_bound: Rational::zero(),
        options: options.clone(),
        field_route: true,
    })
}
