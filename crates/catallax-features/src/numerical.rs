//! Opt-in planar numerical edits. Saved definitions replay all checks, never evidence claims.
use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::{Curve, Line, Plane, Point3, Surface, Vector3};
use catallax_topology::{
    Sense,
    numerical::{
        self, NumericalBody, NumericalDefinition, NumericalEdge, NumericalFace, NumericalOptions,
        NumericalTrim, NumericalVertex,
    },
};
use serde_json::{Value, json};

/// A source and ordered edits, separate from exact recipes and checkpoints.
#[derive(Clone, Debug)]
pub struct NumericalRecipe {
    /// Immutable input, before explicit sewing.
    pub source: NumericalDefinition,
    /// Explicit `[keeper, removed]` endpoint pairs in the source table.
    pub sew: Vec<[usize; 2]>,
    /// Source admission tolerance.
    pub options: NumericalOptions,
    /// Edits in evaluation order.
    pub edits: Vec<NumericalEdit>,
}

/// Controls belong to the edit that uses them; changing them rechecks the whole body.
#[derive(Clone, Debug)]
pub struct NumericalEdit {
    /// Unique nonempty edit identity, used for cut origins.
    pub id: String,
    /// Current tolerance fixed before this edit.
    pub options: NumericalOptions,
    /// The requested edit.
    pub operation: NumericalOperation,
}

/// Finite operations supported by the first planar numerical continuation.
#[derive(Clone, Debug)]
pub enum NumericalOperation {
    /// Change the total size on one coordinate axis, anchored at its lower extent.
    Resize { axis: usize, length_mm: Rational },
    /// Remove a half-space from a two-layer straight prism. The cap retains its
    /// source origins and adds this edit's cap origin. No general Boolean claim.
    Cut {
        axis: usize,
        at_mm: Rational,
        keep_lower: bool,
    },
}

/// A checked numerical result and the source/edit definitions needed to replay it.
#[derive(Clone, Debug)]
pub struct NumericalBuilt {
    body: NumericalBody,
    recipe: NumericalRecipe,
    origins: BTreeMap<String, BTreeSet<String>>,
}
impl NumericalBuilt {
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }
    pub fn recipe(&self) -> &NumericalRecipe {
        &self.recipe
    }
    /// Source and cut identities for each retained face key.
    pub fn face_origins(&self) -> &BTreeMap<String, BTreeSet<String>> {
        &self.origins
    }
    /// Transactional edit. A refusal returns no partial result and leaves this one intact.
    pub fn edit(&self, edit: NumericalEdit, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut recipe = self.recipe.clone();
        recipe.edits.push(edit);
        crate::evaluate_numerical(&recipe, dec)
    }
}
fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("features.numerical-definition", detail)
}
fn unsupported() -> Refusal {
    Refusal::new(
        "features.numerical-unsupported",
        "edit/replay requires rational lines and unit coordinate-axis plane frames; cuts require a two-layer straight prism",
    )
}
fn coordinates(p: &Point3) -> Result<[Rational; 3], Refusal> {
    [&p.x, &p.y, &p.z]
        .map(|v| v.as_rational().ok_or_else(unsupported))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unsupported())
}
fn vector(v: &Vector3) -> Result<[Rational; 3], Refusal> {
    coordinates(&Point3::new(v.x.clone(), v.y.clone(), v.z.clone()))
}
fn point(p: [Rational; 3]) -> Point3 {
    let [x, y, z] = p;
    Point3::new(x.into(), y.into(), z.into())
}
fn direction(p: [Rational; 3]) -> Vector3 {
    let [x, y, z] = p;
    Vector3::new(x.into(), y.into(), z.into())
}
fn axis_of(v: &[Rational; 3]) -> Result<usize, Refusal> {
    let nonzero: Vec<_> = (0..3).filter(|i| v[*i] != Rational::zero()).collect();
    if nonzero.len() != 1 || v[nonzero[0]].abs() != Rational::one() {
        return Err(unsupported());
    }
    Ok(nonzero[0])
}
fn frame(plane: &Plane) -> Result<([Rational; 3], [Rational; 3]), Refusal> {
    let n = vector(&plane.normal)?;
    let x = vector(&plane.ref_dir)?;
    if axis_of(&n)? == axis_of(&x)? {
        return Err(unsupported());
    }
    let y = [
        &n[1] * &x[2] - &n[2] * &x[1],
        &n[2] * &x[0] - &n[0] * &x[2],
        &n[0] * &x[1] - &n[1] * &x[0],
    ];
    Ok((x, y))
}
fn extents(d: &NumericalDefinition, axis: usize) -> Result<(Rational, Rational), Refusal> {
    if axis >= 3 {
        return Err(bad("axis must be 0, 1 or 2"));
    }
    let min = d
        .vertices
        .iter()
        .map(|v| &v.position[axis])
        .min()
        .ok_or_else(unsupported)?
        .clone();
    let max = d
        .vertices
        .iter()
        .map(|v| &v.position[axis])
        .max()
        .ok_or_else(unsupported)?
        .clone();
    if min == max {
        return Err(unsupported());
    }
    Ok((min, max))
}
fn prism(
    d: &NumericalDefinition,
    axis: usize,
    lo: &Rational,
    hi: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for v in &d.vertices {
        dec.charge_work(1)?;
        if &v.position[axis] != lo && &v.position[axis] != hi {
            return Err(unsupported());
        }
    }
    for edge in &d.edges {
        dec.charge_work(1)?;
        let a = &d.vertices[edge.vertices[0]].position;
        let b = &d.vertices[edge.vertices[1]].position;
        if a[axis] != b[axis] && (0..3).any(|i| i != axis && a[i] != b[i]) {
            return Err(unsupported());
        }
    }
    // Exactly two planar caps. Side supports are parallel to the extrusion axis.
    let mut caps = [0, 0];
    for f in &d.faces {
        dec.charge_work(1)?;
        let Surface::Plane(p) = &f.surface else {
            return Err(unsupported());
        };
        let n = vector(&p.normal)?;
        if axis_of(&n)? == axis {
            let z = coordinates(&p.origin)?[axis].clone();
            if &z == lo {
                caps[0] += 1;
            } else if &z == hi {
                caps[1] += 1;
            } else {
                return Err(unsupported());
            }
        }
    }
    if caps != [1, 1] {
        return Err(unsupported());
    }
    Ok(())
}
fn deform(
    d: &NumericalDefinition,
    axis: usize,
    lo: &Rational,
    hi: &Rational,
    next_lo: &Rational,
    next_hi: &Rational,
    dec: &mut Decider,
) -> Result<NumericalDefinition, Refusal> {
    let scale = (next_hi - next_lo) / (hi - lo);
    if scale.sign() != Sign::Positive {
        return Err(bad("dimension must remain positive"));
    }
    let map = |p: &mut [Rational; 3]| {
        let delta = &p[axis] - lo;
        p[axis] = next_lo + &(&delta * &scale);
    };
    let mut d = d.clone();
    for v in &mut d.vertices {
        dec.charge_work(1)?;
        map(&mut v.position);
    }
    for edge in &mut d.edges {
        dec.charge_work(1)?;
        let Curve::Line(line) = &edge.curve else {
            return Err(unsupported());
        };
        let mut p = coordinates(&line.origin)?;
        let mut v = vector(&line.dir)?;
        map(&mut p);
        v[axis] = &v[axis] * &scale;
        edge.curve = Curve::Line(Line::new(point(p), direction(v), dec)?);
    }
    for face in &mut d.faces {
        dec.charge_work(1)?;
        let Surface::Plane(plane) = &face.surface else {
            return Err(unsupported());
        };
        let (x, y) = frame(plane)?;
        let mut p = coordinates(&plane.origin)?;
        map(&mut p);
        face.surface = Surface::Plane(Plane::new(
            point(p),
            plane.normal.clone(),
            plane.ref_dir.clone(),
            dec,
        )?);
        for trim in face.loops.iter_mut().flatten() {
            dec.charge_work(1)?;
            for uv in &mut trim.uv {
                if axis_of(&x)? == axis {
                    uv[0] = &uv[0] * &scale;
                }
                if axis_of(&y)? == axis {
                    uv[1] = &uv[1] * &scale;
                }
            }
        }
    }
    Ok(d)
}

pub(crate) fn evaluate(
    recipe: &NumericalRecipe,
    dec: &mut Decider,
) -> Result<NumericalBuilt, Refusal> {
    let mut ids = BTreeSet::new();
    for edit in &recipe.edits {
        dec.charge_work(1)?;
        if edit.id.trim().is_empty() || !ids.insert(&edit.id) {
            return Err(bad("edit IDs must be nonempty and unique"));
        }
    }
    let mut body = if recipe.sew.is_empty() {
        numerical::build(recipe.source.clone(), recipe.options.clone(), dec)?
    } else {
        numerical::sew(&recipe.source, &recipe.sew, recipe.options.clone(), dec)?
    };
    let mut origins: BTreeMap<String, BTreeSet<String>> = body
        .definition()
        .faces
        .iter()
        .map(|f| {
            (
                f.key.clone(),
                BTreeSet::from([format!("source:{}/{}", recipe.source.name, f.key)]),
            )
        })
        .collect();
    for edit in &recipe.edits {
        let d = body.definition();
        let (axis, next_lo, next_hi) = match &edit.operation {
            NumericalOperation::Resize { axis, length_mm } => {
                let (lo, _) = extents(d, *axis)?;
                (*axis, lo.clone(), &lo + length_mm)
            }
            NumericalOperation::Cut {
                axis,
                at_mm,
                keep_lower,
            } => {
                let (lo, hi) = extents(d, *axis)?;
                prism(d, *axis, &lo, &hi, dec)?;
                if at_mm <= &lo || at_mm >= &hi {
                    return Err(bad("cut must lie strictly between the prism caps"));
                }
                for f in &d.faces {
                    dec.charge_work(1)?;
                    let Surface::Plane(p) = &f.surface else {
                        return Err(unsupported());
                    };
                    if axis_of(&vector(&p.normal)?)? == *axis
                        && coordinates(&p.origin)?[*axis]
                            == if *keep_lower { hi.clone() } else { lo.clone() }
                    {
                        origins
                            .get_mut(&f.key)
                            .expect("retained key")
                            .insert(format!("edit:{}/cap", edit.id));
                    }
                }
                if *keep_lower {
                    (*axis, lo, at_mm.clone())
                } else {
                    (*axis, at_mm.clone(), hi)
                }
            }
        };
        let (lo, hi) = extents(d, axis)?;
        let definition = deform(d, axis, &lo, &hi, &next_lo, &next_hi, dec)?;
        body = body.rebuild(definition, edit.options.clone(), dec)?;
    }
    Ok(NumericalBuilt {
        body,
        recipe: recipe.clone(),
        origins,
    })
}

fn q(value: &Value) -> Result<Rational, Refusal> {
    Rational::parse(
        value
            .as_str()
            .ok_or_else(|| bad("numbers must be strings"))?,
    )
    .map_err(|_| bad("invalid rational"))
}
fn array<const N: usize>(value: &Value) -> Result<[Rational; N], Refusal> {
    let values = value
        .as_array()
        .ok_or_else(|| bad("expected number array"))?;
    if values.len() != N {
        return Err(bad("wrong array length"));
    }
    values
        .iter()
        .map(q)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| bad("wrong array length"))
}
fn list(value: &Value) -> Result<&[Value], Refusal> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| bad("expected array"))
}
fn text(value: &Value) -> Result<String, Refusal> {
    value
        .as_str()
        .map(str::to_owned)
        .ok_or_else(|| bad("expected string"))
}
fn index(value: &Value) -> Result<usize, Refusal> {
    value
        .as_u64()
        .and_then(|i| usize::try_from(i).ok())
        .ok_or_else(|| bad("expected index"))
}
fn sense(value: &Value) -> Result<Sense, Refusal> {
    match value.as_str() {
        Some("same") => Ok(Sense::Same),
        Some("reversed") => Ok(Sense::Reversed),
        _ => Err(bad("invalid sense")),
    }
}
fn sense_text(s: Sense) -> &'static str {
    if s == Sense::Same { "same" } else { "reversed" }
}
fn strings<const N: usize>(a: &[Rational; N]) -> [String; N] {
    a.each_ref().map(ToString::to_string)
}

impl NumericalRecipe {
    /// Save source and edit definitions as `catallax-numerical-recipe/1`.
    /// Evidence is recomputed on replay. Unsupported carriers refuse before output.
    pub fn to_json(&self) -> Result<String, Refusal> {
        let vertices: Vec<_> = self
            .source
            .vertices
            .iter()
            .map(|v| json!({"key":v.key,"position":strings(&v.position)}))
            .collect();
        let mut edges = Vec::new();
        for e in &self.source.edges {
            let Curve::Line(l) = &e.curve else {
                return Err(unsupported());
            };
            edges.push(json!({"key":e.key,"origin":strings(&coordinates(&l.origin)?),"direction":strings(&vector(&l.dir)?),"vertices":e.vertices,"parameters":strings(&e.parameters)}));
        }
        let mut faces = Vec::new();
        for f in &self.source.faces {
            let Surface::Plane(p) = &f.surface else {
                return Err(unsupported());
            };
            frame(p)?;
            let loops: Vec<Vec<_>> = f.loops.iter().map(|b| b.iter().map(|t| json!({"edge":t.edge,"sense":sense_text(t.sense),"uv":t.uv.each_ref().map(strings)})).collect()).collect();
            faces.push(json!({"key":f.key,"origin":strings(&coordinates(&p.origin)?),"normal":strings(&vector(&p.normal)?),"reference":strings(&vector(&p.ref_dir)?),"sense":sense_text(f.sense),"loops":loops}));
        }
        let edits: Vec<_> = self.edits.iter().map(|e| {
            let operation = match &e.operation {
                NumericalOperation::Resize {axis,length_mm} => json!({"kind":"resize","axis":axis,"length_mm":length_mm.to_string()}),
                NumericalOperation::Cut {axis,at_mm,keep_lower} => json!({"kind":"cut","axis":axis,"at_mm":at_mm.to_string(),"keep_lower":keep_lower}),
            };
            json!({"id":e.id,"tolerance_mm":e.options.tolerance_mm.to_string(),"operation":operation})
        }).collect();
        Ok(json!({"format":"catallax-numerical-recipe/1","source":{"name":self.source.name,"vertices":vertices,"edges":edges,"faces":faces},"sew":self.sew,"tolerance_mm":self.options.tolerance_mm.to_string(),"edits":edits}).to_string())
    }
    /// Parse generating data only. Call `evaluate_numerical` to admit it under a
    /// shared request allowance. Saved check claims are never accepted as evidence.
    pub fn parse(input: &str, dec: &mut Decider) -> Result<Self, Refusal> {
        let v: Value = serde_json::from_str(input).map_err(|e| bad(e.to_string()))?;
        if v["format"] != "catallax-numerical-recipe/1" {
            return Err(bad("unsupported numerical recipe format"));
        }
        let source = &v["source"];
        let mut vertices = Vec::new();
        for x in list(&source["vertices"])? {
            dec.charge_work(1)?;
            vertices.push(NumericalVertex {
                key: text(&x["key"])?,
                position: array(&x["position"])?,
            });
        }
        let mut edges = Vec::new();
        for x in list(&source["edges"])? {
            dec.charge_work(1)?;
            let ends = list(&x["vertices"])?;
            if ends.len() != 2 {
                return Err(bad("edge needs two vertices"));
            }
            edges.push(NumericalEdge {
                key: text(&x["key"])?,
                curve: Curve::Line(Line::new(
                    point(array(&x["origin"])?),
                    direction(array(&x["direction"])?),
                    dec,
                )?),
                vertices: [index(&ends[0])?, index(&ends[1])?],
                parameters: array(&x["parameters"])?,
            });
        }
        let mut faces = Vec::new();
        for x in list(&source["faces"])? {
            dec.charge_work(1)?;
            let mut loops = Vec::new();
            for b in list(&x["loops"])? {
                let mut boundary = Vec::new();
                for t in list(b)? {
                    dec.charge_work(1)?;
                    let uv = list(&t["uv"])?;
                    if uv.len() != 2 {
                        return Err(bad("trim needs two endpoints"));
                    }
                    boundary.push(NumericalTrim {
                        edge: index(&t["edge"])?,
                        sense: sense(&t["sense"])?,
                        uv: [array(&uv[0])?, array(&uv[1])?],
                    });
                }
                loops.push(boundary);
            }
            faces.push(NumericalFace {
                key: text(&x["key"])?,
                surface: Surface::Plane(Plane::new(
                    point(array(&x["origin"])?),
                    direction(array(&x["normal"])?),
                    direction(array(&x["reference"])?),
                    dec,
                )?),
                sense: sense(&x["sense"])?,
                loops,
            });
        }
        let mut sew = Vec::new();
        for pair in list(&v["sew"])? {
            dec.charge_work(1)?;
            let pair = list(pair)?;
            if pair.len() != 2 {
                return Err(bad("sewing needs two indices"));
            }
            sew.push([index(&pair[0])?, index(&pair[1])?]);
        }
        let mut edits = Vec::new();
        for e in list(&v["edits"])? {
            dec.charge_work(1)?;
            let op = &e["operation"];
            let axis = index(&op["axis"])?;
            let operation = match op["kind"].as_str() {
                Some("resize") => NumericalOperation::Resize {
                    axis,
                    length_mm: q(&op["length_mm"])?,
                },
                Some("cut") => NumericalOperation::Cut {
                    axis,
                    at_mm: q(&op["at_mm"])?,
                    keep_lower: op["keep_lower"]
                        .as_bool()
                        .ok_or_else(|| bad("expected keep_lower boolean"))?,
                },
                _ => return Err(bad("unsupported numerical edit")),
            };
            edits.push(NumericalEdit {
                id: text(&e["id"])?,
                options: NumericalOptions {
                    tolerance_mm: q(&e["tolerance_mm"])?,
                },
                operation,
            });
        }
        Ok(Self {
            source: NumericalDefinition {
                name: text(&source["name"])?,
                vertices,
                edges,
                faces,
            },
            sew,
            options: NumericalOptions {
                tolerance_mm: q(&v["tolerance_mm"])?,
            },
            edits,
        })
    }
}
