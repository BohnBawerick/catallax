//! V7d-K, fixed by the independent V7d designer: one clamped linear grade along a part
//! axis drives the cubic lattice radius and the network gyroid level. Recipes, settings,
//! expected values and refusal codes are the designer's; the builder may not edit them.
use catallax_field::{
    Location, Plan, Workspace,
    lattice_stage::{
        GyroidSpec, LatticeSpec, Output, Proof, Settings, Stage, default_preview_options,
    },
    mass,
    printing::{BoxState, classify_box},
};
use catallax_numbers::{Decider, OperationWork, Rational, Refusal, with_operation_budget};
use catallax_recipe::{Recipe, RecipeError};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

const APP_WORK: u64 = 19_999_999;
const KERNEL_WORK: u64 = 99_999_999;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

/// Designer section 3, D0.
fn d0() -> Value {
    json!({"format":"catallax-recipe/4","name":"V7d-D0","steps":[
        {"id":"region","op":"block","size":["6","6","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-3","-3","20"]}},
        {"id":"field","op":"lattice","region":"region","cell":"cubic","size":"100",
         "frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
         "radius":{"linear":{"from":["0","0","20"],"to":["0","0","30"],"start":"1","end":"2"}}}],"result":"field"})
}
fn d1() -> Value {
    let mut r = d0();
    r["name"] = json!("V7d-D1");
    r["steps"][0]["size"][2] = json!("12");
    r
}
fn d2() -> Value {
    let mut r = d0();
    r["name"] = json!("V7d-D2");
    r["steps"][1]["radius"]["linear"]["start"] = json!("2");
    r["steps"][1]["radius"]["linear"]["end"] = json!("1");
    r
}
/// Designer section 3, G0.
fn g0() -> Value {
    json!({"format":"catallax-recipe/4","name":"V7d-G0","steps":[
        {"id":"region","op":"block","size":["2","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["-1","-1","-1"]}},
        {"id":"field","op":"gyroid","region":"region","period":"10",
         "level":{"linear":{"from":["0","0","-1"],"to":["0","0","1"],"start":"-1/10","end":"1/10"}},
         "kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"field"})
}
fn g1() -> Value {
    let mut r = g0();
    r["name"] = json!("V7d-G1");
    r["steps"][0]["size"][0] = json!("12/5");
    r["steps"][0]["pose"]["translation"][0] = json!("-6/5");
    r
}
fn g2() -> Value {
    let mut r = g0();
    r["name"] = json!("V7d-G2");
    r["steps"][1]["level"]["linear"]["start"] = json!("1/10");
    r["steps"][1]["level"]["linear"]["end"] = json!("-1/10");
    r
}
fn states() -> [(&'static str, Value); 6] {
    [
        ("D0", d0()),
        ("D1", d1()),
        ("D2", d2()),
        ("G0", g0()),
        ("G1", g1()),
        ("G2", g2()),
    ]
}

/// Cyclic coordinate permutation: new axis `(i + shift) % 3` takes old axis `i`.
fn permute(recipe: &Value, shift: usize) -> Value {
    let mut r = recipe.clone();
    let perm = |v: &Value| {
        let mut out = vec![Value::Null; 3];
        for i in 0..3 {
            out[(i + shift) % 3] = v[i].clone();
        }
        Value::Array(out)
    };
    for step in r["steps"].as_array_mut().unwrap() {
        if step["op"] == "block" {
            step["size"] = perm(&step["size"]);
            step["pose"]["translation"] = perm(&step["pose"]["translation"]);
        } else {
            let key = if step["op"] == "lattice" {
                "radius"
            } else {
                "level"
            };
            let linear = &mut step[key]["linear"];
            if linear.is_object() {
                linear["from"] = perm(&linear["from"]);
                linear["to"] = perm(&linear["to"]);
            }
            step["frame"]["translation"] = perm(&step["frame"]["translation"]);
        }
    }
    r
}

/// Translates a D fixture by `(7, -11, 13)`: region, lattice frame and anchors.
fn translate(recipe: &Value) -> Value {
    let mut r = recipe.clone();
    let offset = [q("7"), q("-11"), q("13")];
    let shift = |v: &Value| {
        Value::Array(
            (0..3)
                .map(|i| json!((&q(v[i].as_str().unwrap()) + &offset[i]).to_string()))
                .collect(),
        )
    };
    for step in r["steps"].as_array_mut().unwrap() {
        if step["op"] == "block" {
            step["pose"]["translation"] = shift(&step["pose"]["translation"]);
        } else {
            step["frame"]["translation"] = shift(&step["frame"]["translation"]);
            let linear = &mut step["radius"]["linear"];
            linear["from"] = shift(&linear["from"]);
            linear["to"] = shift(&linear["to"]);
        }
    }
    r
}

fn parse(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}

/// The body recipe (everything before the terminal field) and the field's parameters.
enum Field {
    Lattice(LatticeSpec),
    Gyroid(GyroidSpec),
}
fn split(value: &Value) -> (Recipe, Field) {
    let mut body = value.clone();
    let field = body["steps"].as_array_mut().unwrap().pop().unwrap();
    body["result"] = field["region"].clone();
    let spec = if field["op"] == "lattice" {
        Field::Lattice(LatticeSpec {
            cell: serde_json::from_value(field["cell"].clone()).unwrap(),
            size: serde_json::from_value(field["size"].clone()).unwrap(),
            frame: serde_json::from_value(field["frame"].clone()).unwrap(),
            radius: serde_json::from_value(field["radius"].clone()).unwrap(),
        })
    } else {
        Field::Gyroid(GyroidSpec {
            period: serde_json::from_value(field["period"].clone()).unwrap(),
            level: serde_json::from_value(field["level"].clone()).unwrap(),
            kind: serde_json::from_value(field["kind"].clone()).unwrap(),
            frame: serde_json::from_value(field["frame"].clone()).unwrap(),
        })
    };
    (parse(&body), spec)
}

/// Designer section 3: the positive native Stage settings.
fn settings() -> Settings {
    Settings {
        accuracy: q("1/10"),
        minimum_detail: q("1/10"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}

/// A staged field of either family, for the calls the tests share.
enum Staged {
    Lattice(Stage<LatticeSpec>),
    Gyroid(Stage<GyroidSpec>),
}
impl Staged {
    fn new(value: &Value, settings: Settings, dec: &mut Decider) -> Result<Self, Refusal> {
        let (body, spec) = split(value);
        Ok(match spec {
            Field::Lattice(spec) => Self::Lattice(Stage::new(&body, spec, settings, dec)?),
            Field::Gyroid(spec) => Self::Gyroid(Stage::new(&body, spec, settings, dec)?),
        })
    }
    fn plan(&self) -> &Plan {
        match self {
            Self::Lattice(s) => s.plan(),
            Self::Gyroid(s) => s.plan(),
        }
    }
    fn recipe(&self) -> &Recipe {
        match self {
            Self::Lattice(s) => s.recipe(),
            Self::Gyroid(s) => s.recipe(),
        }
    }
    fn export(&self, dec: &mut Decider) -> Result<Output, Refusal> {
        match self {
            Self::Lattice(s) => s.export(dec),
            Self::Gyroid(s) => s.export(dec),
        }
    }
    fn preview_is_unproved(&self, dec: &mut Decider) -> bool {
        let options = default_preview_options(self.plan(), dec).unwrap();
        let preview = match self {
            Self::Lattice(s) => s.preview(&options, dec),
            Self::Gyroid(s) => s.preview(&options, dec),
        }
        .unwrap();
        preview.proof == Proof::Unproved && !preview.mesh.triangles.is_empty()
    }
}

fn stage_error(value: &Value, settings: Settings) -> Refusal {
    match Staged::new(value, settings, &mut Decider::default()) {
        Ok(_) => panic!("stage accepted {value}"),
        Err(error) => error,
    }
}

/// One complete app action under the shared allowance: build, unproved preview, export
/// and labelling.
fn app_export(value: &Value, settings: Settings) -> ((Vec<u8>, Output), OperationWork) {
    with_operation_budget(APP_WORK, || {
        let mut dec = Decider::default();
        let stage = Staged::new(value, settings, &mut dec).unwrap();
        assert!(stage.preview_is_unproved(&mut dec));
        stage.export(&mut dec).unwrap().labelled_stl().unwrap()
    })
}

/// Independent reading of a labelled binary STL: header fields, then the topology and
/// placement checks of designer section 6 on the actual emitted coordinates.
fn check_file(bytes: &[u8], accuracy: &str, bounds: [[Rational; 2]; 3]) -> (usize, Rational) {
    assert_eq!(&bytes[..8], b"CXPRINT1");
    let word = |start: usize| u64::from_le_bytes(bytes[start..start + 8].try_into().unwrap());
    let requested = q(accuracy);
    let (numerator, denominator) = (word(12) as i64, word(20) as i64);
    assert_eq!(
        Rational::from_i64(numerator)
            .checked_div(&Rational::from_i64(denominator))
            .unwrap(),
        requested
    );
    assert_eq!((word(28), word(36)), (2_000_000, 1_000_000));
    let bound =
        Rational::from_f64_exact(f64::from_le_bytes(bytes[68..76].try_into().unwrap())).unwrap();
    assert!(bound > q("0") && bound <= requested, "file bound {bound}");
    let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    assert_eq!(bytes.len(), 84 + 50 * count);
    let mut ids = std::collections::BTreeMap::new();
    let mut points: Vec<[f32; 3]> = Vec::new();
    let mut triangles = Vec::new();
    for facet in bytes[84..].as_chunks::<50>().0 {
        let corner = |v: usize| -> [f32; 3] {
            std::array::from_fn(|i| {
                let start = 12 + v * 12 + i * 4;
                f32::from_le_bytes(facet[start..start + 4].try_into().unwrap())
            })
        };
        let t: [usize; 3] = std::array::from_fn(|v| {
            let p = corner(v);
            let key = p.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
            let next = ids.len();
            *ids.entry(key).or_insert_with(|| {
                points.push(p);
                next
            })
        });
        assert!(
            t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
            "degenerate facet"
        );
        triangles.push(t);
    }
    // Every edge has exactly two facets of opposite sense; no facet repeats.
    let mut directed = std::collections::BTreeSet::new();
    let mut faces = std::collections::BTreeSet::new();
    for t in &triangles {
        let mut key = *t;
        key.sort_unstable();
        assert!(faces.insert(key), "duplicate facet");
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            assert!(directed.insert((a, b)), "edge used twice in one sense");
        }
    }
    for &(a, b) in &directed {
        assert!(directed.contains(&(b, a)), "open edge");
    }
    // Every vertex link is one cycle.
    let mut links = vec![Vec::new(); points.len()];
    for t in &triangles {
        for k in 0..3 {
            links[t[k]].push((t[(k + 1) % 3], t[(k + 2) % 3]));
        }
    }
    for link in &links {
        let mut next = std::collections::BTreeMap::new();
        for &(a, b) in link {
            assert!(next.insert(a, b).is_none(), "vertex link is not a cycle");
        }
        let start = link[0].0;
        let mut at = start;
        let mut steps = 0;
        loop {
            at = next[&at];
            steps += 1;
            if at == start {
                break;
            }
        }
        assert_eq!(steps, link.len(), "vertex link has several cycles");
    }
    // One connected component, Euler characteristic 2, nonzero facet areas and outward
    // orientation by a positive enclosed volume.
    let mut seen = vec![false; points.len()];
    let mut pending = vec![0];
    while let Some(v) = pending.pop() {
        if !seen[v] {
            seen[v] = true;
            pending.extend(links[v].iter().map(|&(a, _)| a));
        }
    }
    assert!(seen.iter().all(|&v| v), "mesh has several components");
    let edges = directed.len() / 2;
    assert_eq!(
        points.len() as i64 - edges as i64 + triangles.len() as i64,
        2
    );
    let exact = |p: [f32; 3]| p.map(|x| Rational::from_f64_exact(f64::from(x)).unwrap());
    let mut volume = q("0");
    for t in &triangles {
        let [a, b, c] = t.map(|i| exact(points[i]));
        let u: [Rational; 3] = std::array::from_fn(|i| &b[i] - &a[i]);
        let v: [Rational; 3] = std::array::from_fn(|i| &c[i] - &a[i]);
        let n: [Rational; 3] = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            &(&u[j] * &v[k]) - &(&u[k] * &v[j])
        });
        assert!(n.iter().any(|x| !x.is_zero()), "zero-area facet");
        volume = &volume + &(0..3).fold(q("0"), |s, i| &s + &(&a[i] * &n[i]));
    }
    assert!(volume > q("0"), "facets are not outward");
    for p in &points {
        for (i, x) in exact(*p).iter().enumerate() {
            assert!(*x >= &bounds[i][0] - &bound && *x <= &bounds[i][1] + &bound);
        }
    }
    (triangles.len(), bound)
}

fn tight_bounds(name: &str) -> [[Rational; 2]; 3] {
    match name {
        "D0" | "D2" => [[q("-2"), q("2")], [q("-2"), q("2")], [q("20"), q("30")]],
        "D1" => [[q("-2"), q("2")], [q("-2"), q("2")], [q("20"), q("32")]],
        "G1" => [[q("-6/5"), q("6/5")], [q("-1"), q("1")], [q("-1"), q("1")]],
        _ => [[q("-1"), q("1")], [q("-1"), q("1")], [q("-1"), q("1")]],
    }
}

/// Designer sections 3, 6 and 8: every state exports through the native stage at 1/10
/// mm, D0 and G0 also at 1/2 mm, inside one app allowance, with an unproved preview and
/// a labelled file whose actual bytes pass the independent topology checks.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d exports run in release mode")]
fn six_states_export_through_the_stage_within_one_app_allowance() {
    for (name, value) in states() {
        let accuracies: &[&str] = if name.ends_with('0') {
            &["1/10", "1/2"]
        } else {
            &["1/10"]
        };
        for accuracy in accuracies {
            let ((bytes, output), work) = app_export(
                &value,
                Settings {
                    accuracy: q(accuracy),
                    ..settings()
                },
            );
            assert!(!work.exhausted && work.used < APP_WORK, "{name} {work:?}");
            let verdict = output.verdict();
            assert!(verdict.closed && verdict.non_self_intersecting && verdict.connected);
            assert_eq!(verdict.components, 1);
            assert!(output.mesh().cells <= 2_000_000 && output.mesh().positions.len() <= 1_000_000);
            let (facets, bound) = check_file(&bytes, accuracy, tight_bounds(name));
            assert_eq!(bound, verdict.max_deviation);
            println!(
                "{name} at {accuracy}: facets {facets} cells {} work {}",
                output.mesh().cells,
                work.used
            );
        }
    }
}

/// Designer section 6: cyclic permutations of every state, and D translated by
/// `(7, -11, 13)`, export at 1/10 mm inside the same allowance.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d exports run in release mode")]
fn permuted_and_translated_states_export_at_a_tenth() {
    let permute_bounds = |b: [[Rational; 2]; 3], shift: usize| {
        let mut out = b.clone();
        for i in 0..3 {
            out[(i + shift) % 3] = b[i].clone();
        }
        out
    };
    for (name, value) in states() {
        for shift in [1, 2] {
            let ((bytes, output), work) = app_export(&permute(&value, shift), settings());
            assert!(
                !work.exhausted && work.used < APP_WORK,
                "{name} {shift} {work:?}"
            );
            assert_eq!(output.verdict().components, 1);
            check_file(&bytes, "1/10", permute_bounds(tight_bounds(name), shift));
        }
        if name.starts_with('D') {
            let ((bytes, _), work) = app_export(&translate(&value), settings());
            assert!(!work.exhausted && work.used < APP_WORK, "{name} T {work:?}");
            let offset = [q("7"), q("-11"), q("13")];
            let moved = tight_bounds(name).map(|b| b);
            let moved: [[Rational; 2]; 3] =
                std::array::from_fn(|i| [&moved[i][0] + &offset[i], &moved[i][1] + &offset[i]]);
            check_file(&bytes, "1/10", moved);
        }
    }
}

const PI_LO: &str = "3141592653589793/1000000000000000";
const PI_HI: &str = "3141592653589794/1000000000000000";

fn accuracy() -> mass::Accuracy {
    let width = |absolute: &str| catallax_measure::Accuracy {
        relative: q("0"),
        absolute: q(absolute),
    };
    mass::Accuracy {
        volume: width("1/100"),
        centroid: width("1/1000"),
        inertia: width("1/10"),
    }
}

/// Native moments of one staged state inside one app allowance.
fn moments(value: &Value) -> (mass::Properties, OperationWork) {
    let (properties, work) = with_operation_budget(APP_WORK, || {
        let mut dec = Decider::default();
        let stage = Staged::new(value, settings(), &mut dec).unwrap();
        mass::properties(
            stage.plan(),
            &accuracy(),
            &mass::Options {
                max_cells: 2_000_000,
            },
            &mut dec,
        )
    });
    (properties.unwrap(), work)
}

fn contains(e: &catallax_measure::Enclosure, v: &Rational) -> bool {
    e.lo <= *v && *v <= e.hi
}
fn contains_pi(e: &catallax_measure::Enclosure, a: &Rational) -> bool {
    e.lo <= a * &q(PI_LO) && a * &q(PI_HI) <= e.hi
}
fn overlaps(e: &catallax_measure::Enclosure, [lo, hi]: [&str; 2]) -> bool {
    e.lo <= q(hi) && q(lo) <= e.hi
}
fn widths(p: &mass::Properties) {
    assert!(
        p.volume.width() <= q("1/100"),
        "volume width {}",
        p.volume.width()
    );
    for c in p.centroid.as_ref().unwrap() {
        assert!(c.width() <= q("1/1000"), "centroid width {}", c.width());
    }
    for row in p.inertia.as_ref().unwrap() {
        for v in row {
            assert!(v.width() <= q("1/10"), "inertia width {}", v.width());
        }
    }
}

/// Undoes a cyclic permutation: original axis `i` is reported on axis `(i + shift) % 3`.
fn unpermute(
    p: &mass::Properties,
    shift: usize,
) -> (
    [catallax_measure::Enclosure; 3],
    [[catallax_measure::Enclosure; 3]; 3],
) {
    let c = p.centroid.as_ref().unwrap();
    let t = p.inertia.as_ref().unwrap();
    (
        std::array::from_fn(|i| c[(i + shift) % 3].clone()),
        std::array::from_fn(|i| {
            std::array::from_fn(|j| t[(i + shift) % 3][(j + shift) % 3].clone())
        }),
    )
}

/// Designer section 5, D references: `V/π`, centroid z, `Ixx/π = Iyy/π`, `Izz/π`.
fn d_reference(name: &str) -> [&'static str; 4] {
    match name {
        "D0" => ["70/3", "365/14", "1321/7", "31"],
        "D1" => ["94/3", "2569/94", "48517/141", "47"],
        _ => ["70/3", "335/14", "1321/7", "31"],
    }
}

fn check_d(name: &str, p: &mass::Properties, shift: usize, offset: [Rational; 3]) {
    widths(p);
    let [v, cz, ixx, izz] = d_reference(name).map(q);
    assert!(contains_pi(&p.volume, &v), "{name} volume");
    let (c, t) = unpermute(p, shift);
    let expected = [&offset[0] + &q("0"), offset[1].clone(), &cz + &offset[2]];
    for i in 0..3 {
        assert!(contains(&c[i], &expected[i]), "{name} centroid {i}");
    }
    assert!(
        contains_pi(&t[0][0], &ixx) && contains_pi(&t[1][1], &ixx),
        "{name} Ixx Iyy"
    );
    assert!(contains_pi(&t[2][2], &izz), "{name} Izz");
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        assert!(contains(&t[i][j], &q("0")), "{name} I{i}{j}");
    }
}

/// Designer section 5, G references: exact volume, exact raw second moments, and the
/// independent centroid and tensor intervals (xy, xz, yz off-diagonal order).
type GReference = (&'static str, [&'static str; 3], [[&'static str; 2]; 9]);
fn g_reference(name: &str) -> GReference {
    match name {
        "G0" => (
            "4",
            ["4/3", "4/3", "4/3"],
            [
                ["-7323723/25000000", "-140831983/500000000"],
                ["-293376209/1000000000", "-282357203/1000000000"],
                ["-118964903/500000000", "-113220591/500000000"],
                ["2095945899/1000000000", "2142661869/1000000000"],
                ["2096948019/1000000000", "2144225871/1000000000"],
                ["1979111989/1000000000", "2030425947/1000000000"],
                ["159059699/500000000", "171888487/500000000"],
                ["127560643/500000000", "139402559/500000000"],
                ["63937299/250000000", "279211777/1000000000"],
            ],
        ),
        "G1" => (
            "24/5",
            ["288/125", "8/5", "8/5"],
            [
                ["-10000177/25000000", "-386465873/1000000000"],
                ["-268183199/1000000000", "-128022813/500000000"],
                ["-43057467/200000000", "-203326091/1000000000"],
                ["52645997/20000000", "1343437931/500000000"],
                ["2913499361/1000000000", "1494326311/500000000"],
                ["2790746121/1000000000", "2872406877/1000000000"],
                ["29685869/62500000", "257460427/500000000"],
                ["377177257/1000000000", "206679499/500000000"],
                ["24989163/100000000", "277134941/1000000000"],
            ],
        ),
        _ => (
            "4",
            ["4/3", "4/3", "4/3"],
            [
                ["-25674223/100000000", "-49246889/200000000"],
                ["-256756431/1000000000", "-247014493/1000000000"],
                ["-6191257/20000000", "-18749169/62500000"],
                ["126228411/62500000", "2062633937/1000000000"],
                ["1009841873/500000000", "2064172967/1000000000"],
                ["2139304919/1000000000", "2180076419/1000000000"],
                ["243293907/1000000000", "131840437/500000000"],
                ["295468239/1000000000", "12716457/40000000"],
                ["148202127/500000000", "31792901/100000000"],
            ],
        ),
    }
}

fn check_g(name: &str, p: &mass::Properties, shift: usize) {
    use catallax_measure::Enclosure;
    widths(p);
    let (volume, raw, reference) = g_reference(name);
    assert!(contains(&p.volume, &q(volume)), "{name} volume");
    let (c, t) = unpermute(p, shift);
    let add = |a: &Enclosure, b: &Enclosure| Enclosure {
        lo: &a.lo + &b.lo,
        hi: &a.hi + &b.hi,
    };
    let product = |a: &Enclosure, b: &Enclosure| {
        let v = [&a.lo * &b.lo, &a.lo * &b.hi, &a.hi * &b.lo, &a.hi * &b.hi];
        Enclosure {
            lo: v.iter().min().unwrap().clone(),
            hi: v.iter().max().unwrap().clone(),
        }
    };
    for i in 0..3 {
        // The raw second moment, reconstructed from V, C and I: Σᵢᵢ + V cᵢ².
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let sum = add(&t[j][j], &t[k][k]);
        let sigma = Enclosure {
            lo: (&sum.lo - &t[i][i].hi).mul_pow2(-1),
            hi: (&sum.hi - &t[i][i].lo).mul_pow2(-1),
        };
        let raw_enclosure = add(&sigma, &product(&p.volume, &product(&c[i], &c[i])));
        assert!(
            contains(&raw_enclosure, &q(raw[i])),
            "{name} raw second moment {i}"
        );
        assert!(overlaps(&c[i], reference[i]), "{name} centroid {i}");
        assert!(overlaps(&t[i][i], reference[3 + i]), "{name} I{i}{i}");
    }
    for (n, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        assert!(overlaps(&t[i][j], reference[6 + n]), "{name} I{i}{j}");
    }
}

/// Designer sections 5 and 8: requested native moments of every state meet the widths
/// and contain the independent references, inside one app allowance.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d moments run in release mode")]
fn six_states_measure_within_the_requested_widths() {
    for (name, value) in states() {
        let (properties, work) = moments(&value);
        assert!(!work.exhausted && work.used < APP_WORK, "{name} {work:?}");
        assert!(properties.cells <= 2_000_000);
        if name.starts_with('D') {
            check_d(name, &properties, 0, [q("0"), q("0"), q("0")]);
        } else {
            check_g(name, &properties, 0);
        }
        println!("{name}: cells {} work {}", properties.cells, work.used);
    }
}

/// Designer section 6: permuted states transform C and I by the permutation, and the
/// translated D states move only the centroid.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d moments run in release mode")]
fn permuted_and_translated_moments_transform_with_the_part() {
    for (name, value) in states() {
        for shift in [1, 2] {
            let (properties, work) = moments(&permute(&value, shift));
            assert!(
                !work.exhausted && work.used < APP_WORK,
                "{name} {shift} {work:?}"
            );
            if name.starts_with('D') {
                check_d(name, &properties, shift, [q("0"), q("0"), q("0")]);
            } else {
                check_g(name, &properties, shift);
            }
        }
        if name.starts_with('D') {
            let (properties, work) = moments(&translate(&value));
            assert!(!work.exhausted && work.used < APP_WORK, "{name} T {work:?}");
            check_d(name, &properties, 0, [q("7"), q("-11"), q("13")]);
        }
    }
}

fn location(value: &Value, point: [&str; 3]) -> Location {
    let plan = Plan::compile_checked(&parse(value), &mut Decider::default()).unwrap();
    let [x, y, z] = point.map(q);
    plan.proven_checked(
        &RPoint3::new(x, y, z),
        &mut Workspace::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .location
}

/// Designer section 6: exact point classifications. The radius is read at the query
/// point and anchors stay in part coordinates through the upstream edit.
#[test]
fn exact_points_classify_on_both_sides_of_each_boundary() {
    assert_eq!(location(&d0(), ["3/2", "0", "25"]), Location::Boundary);
    assert_eq!(location(&d0(), ["1501/1000", "0", "25"]), Location::Outside);
    assert_eq!(location(&d0(), ["1499/1000", "0", "25"]), Location::Inside);
    assert_eq!(location(&d0(), ["6/5", "0", "22"]), Location::Boundary);
    assert_eq!(location(&d2(), ["9/5", "0", "22"]), Location::Boundary);
    assert_eq!(location(&d1(), ["2", "0", "31"]), Location::Boundary);
    for value in [g0(), g1(), g2()] {
        assert_eq!(location(&value, ["0", "0", "0"]), Location::Boundary);
    }
    for value in [g0(), g1()] {
        assert_eq!(location(&value, ["0", "0", "1/1000"]), Location::Outside);
        assert_eq!(location(&value, ["0", "0", "-1/1000"]), Location::Inside);
    }
}

/// Designer section 6: grade queries return the exact clamped value at each point.
#[test]
fn grade_queries_return_exact_clamped_values() {
    let level: catallax_recipe::Grade =
        serde_json::from_value(g0()["steps"][1]["level"].clone()).unwrap();
    for (z, expected) in [
        ("-2", "-1/10"),
        ("-1", "-1/10"),
        ("0", "0"),
        ("1", "1/10"),
        ("2", "1/10"),
    ] {
        let p = [q("0"), q("0"), q(z)];
        let range = catallax_field::grade_range(&level, &p, &p, &mut Decider::default()).unwrap();
        assert_eq!(range, [q(expected), q(expected)], "z = {z}");
    }
}

/// Q-clamp: a box across the end anchor plane keeps the exact clamped range, not the
/// value at its centre or an unclamped endpoint. After that successful evaluation, the
/// same evaluation with an isolated zero operation allowance refuses at its first
/// checkpoint (designer ruling on the zero-work rows).
#[test]
fn a_box_across_an_anchor_plane_keeps_the_clamped_range() {
    let level: catallax_recipe::Grade =
        serde_json::from_value(g0()["steps"][1]["level"].clone()).unwrap();
    let lo = [q("0"), q("0"), q("999/1000")];
    let hi = [q("0"), q("0"), q("1001/1000")];
    let range = catallax_field::grade_range(&level, &lo, &hi, &mut Decider::default()).unwrap();
    assert_eq!(range, [q("999/10000"), q("1/10")]);
    let (result, work) = with_operation_budget(0, || {
        catallax_field::grade_range(&level, &lo, &hi, &mut Decider::default())
    });
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted && work.used == 0);
}

/// Q-uncertain: a box holding material and void witnesses is unresolved, never wholly
/// inside or outside. After a separately recorded successful compilation, the proven
/// point query at `(0, 0, 1/1000)` with an isolated zero operation allowance refuses
/// at its first checkpoint (designer ruling on the zero-work rows).
#[test]
fn a_box_holding_both_sides_is_unresolved() {
    let ((plan, state), work) = with_operation_budget(KERNEL_WORK, || {
        let mut dec = Decider::default();
        let plan = Plan::compile_checked(&parse(&g0()), &mut dec).unwrap();
        let state = classify_box(
            &plan,
            &[q("-1/1000"), q("-1/1000"), q("-1/1000")],
            &[q("1/1000"), q("1/1000"), q("1/1000")],
            &mut dec,
        )
        .unwrap();
        (plan, state)
    });
    assert!(!work.exhausted);
    assert_eq!(state, BoxState::Unknown);
    let point = RPoint3::new(q("0"), q("0"), q("1/1000"));
    let answer = plan
        .proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
        .unwrap();
    assert_eq!(answer.location, Location::Outside);
    let (result, work) = with_operation_budget(0, || {
        plan.proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
    });
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted && work.used == 0);
}

fn with_radius(base: Value, change: impl FnOnce(&mut Value)) -> Value {
    let mut r = base;
    change(&mut r["steps"][1]["radius"]);
    r
}
fn with_level(base: Value, change: impl FnOnce(&mut Value)) -> Value {
    let mut r = base;
    change(&mut r["steps"][1]["level"]);
    r
}

/// Designer section 7: the grade-validation rows at native stage validation.
#[test]
fn stage_validation_names_each_grade_refusal() {
    let rows: Vec<(&str, Value, Option<&str>)> = vec![
        (
            "L-zero-span",
            with_radius(d0(), |g| g["linear"]["to"] = json!(["0", "0", "20"])),
            Some("lattice.impossible"),
        ),
        (
            "L-zero-radius",
            with_radius(d0(), |g| g["linear"]["start"] = json!("0")),
            Some("lattice.impossible"),
        ),
        (
            "L-negative-radius",
            with_radius(d0(), |g| g["linear"]["start"] = json!("-1/1000")),
            Some("lattice.impossible"),
        ),
        (
            "L-K-edge",
            with_radius(d0(), |g| g["linear"]["end"] = json!("6")),
            None,
        ),
        (
            "L-K-over",
            with_radius(d0(), |g| g["linear"]["end"] = json!("600001/100000")),
            Some("lattice.grade-range-unsupported"),
        ),
        (
            "L-radius-edge",
            with_radius(d0(), |g| {
                g["linear"]["start"] = json!("24");
                g["linear"]["end"] = json!("25");
            }),
            None,
        ),
        (
            "L-radius-over",
            with_radius(d0(), |g| {
                g["linear"]["start"] = json!("24");
                g["linear"]["end"] = json!("25000001/1000000");
            }),
            Some("lattice.grade-range-unsupported"),
        ),
        (
            "L-axis-other",
            with_radius(d0(), |g| g["linear"]["to"] = json!(["1", "0", "30"])),
            Some("lattice.grade-function-unsupported"),
        ),
        (
            "G-zero-span",
            with_level(g0(), |g| g["linear"]["to"] = json!(["0", "0", "-1"])),
            Some("gyroid.impossible"),
        ),
        (
            "G-zero-period",
            {
                let mut r = g0();
                r["steps"][1]["period"] = json!("0");
                r
            },
            Some("gyroid.impossible"),
        ),
        ("G-K-edge", g0(), None),
        (
            "G-K-over",
            with_level(g0(), |g| g["linear"]["end"] = json!("50001/500000")),
            Some("gyroid.grade-range-unsupported"),
        ),
        (
            "G-level-edge",
            with_level(g0(), |g| {
                g["linear"]["start"] = json!("6/25");
                g["linear"]["end"] = json!("1/4");
            }),
            None,
        ),
        (
            "G-level-over",
            with_level(g0(), |g| {
                g["linear"]["start"] = json!("6/25");
                g["linear"]["end"] = json!("250001/1000000");
            }),
            Some("gyroid.grade-range-unsupported"),
        ),
        (
            "G-axis-other",
            with_level(g0(), |g| g["linear"]["to"] = json!(["1", "0", "1"])),
            Some("gyroid.grade-function-unsupported"),
        ),
        (
            "G-sheet",
            {
                let mut r = g0();
                r["steps"][1]["kind"] = json!("sheet");
                r
            },
            Some("gyroid.unsupported"),
        ),
        (
            "G-frame",
            {
                let mut r = g0();
                r["steps"][1]["frame"]["translation"] = json!(["1", "0", "0"]);
                r
            },
            Some("gyroid.unsupported"),
        ),
        (
            "U-distance-L",
            with_radius(
                d0(),
                |g| *g = json!({"distance":{"faces":[],"near":"1/10","far":"1/5","span":"1"}}),
            ),
            Some("lattice.unsupported-grade"),
        ),
        (
            "U-distance-G",
            with_level(
                g0(),
                |g| *g = json!({"distance":{"faces":[],"near":"1/10","far":"1/5","span":"1"}}),
            ),
            Some("gyroid.unsupported"),
        ),
    ];
    for (row, mut value, expected) in rows {
        value["name"] = json!(row);
        match (
            Staged::new(&value, settings(), &mut Decider::default()),
            expected,
        ) {
            (Ok(_), None) => {}
            (Err(error), Some(code)) => assert_eq!(error.code, code, "{row}: {}", error.detail),
            (Ok(_), Some(code)) => panic!("{row}: accepted, expected {code}"),
            (Err(error), None) => panic!("{row}: refused {}", error.code),
        }
    }
    // The oblique grade stays an ordinary lattice definition for point queries.
    let oblique = with_radius(d0(), |g| g["linear"]["to"] = json!(["1", "0", "30"]));
    assert_eq!(location(&oblique, ["0", "0", "25"]), Location::Inside);
}

/// Designer section 7: equal endpoints are their constants, with no graded limits.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d exports run in release mode")]
fn equal_endpoints_equal_their_constants() {
    let equal_d = with_radius(d0(), |g| {
        g["linear"]["start"] = json!("1");
        g["linear"]["end"] = json!("1");
    });
    let constant_d = with_radius(d0(), |g| *g = json!({"constant":"1"}));
    let equal_g = with_level(g0(), |g| {
        g["linear"]["start"] = json!("0");
        g["linear"]["end"] = json!("0");
    });
    let constant_g = with_level(g0(), |g| *g = json!({"constant":"0"}));
    for (equal, constant) in [(equal_d, constant_d), (equal_g, constant_g)] {
        let export = |value: &Value| {
            let mut dec = Decider::default();
            Staged::new(value, settings(), &mut dec)
                .unwrap()
                .export(&mut dec)
                .map(|output| output.labelled_stl().unwrap().0)
        };
        match (export(&equal), export(&constant)) {
            (Ok(a), Ok(b)) => assert_eq!(a, b),
            (Err(a), Err(b)) => assert_eq!((a.code, a.detail), (b.code, b.detail)),
            _ => panic!("an equal-endpoint grade differs from its constant"),
        }
    }
    // A constant level outside the graded range keeps its constant support.
    let high = with_level(g0(), |g| {
        g["linear"]["start"] = json!("1/2");
        g["linear"]["end"] = json!("1/2");
    });
    assert!(Staged::new(&high, settings(), &mut Decider::default()).is_ok());
}

/// Designer section 7: settings, cell, vertex and moment limits refuse with their own
/// typed codes, never a partial mesh or fabricated bounds.
#[test]
fn settings_and_caps_refuse_with_typed_codes() {
    for change in [
        Settings {
            accuracy: q("0"),
            ..settings()
        },
        Settings {
            minimum_detail: q("0"),
            ..settings()
        },
        Settings {
            max_cells: 0,
            ..settings()
        },
        Settings {
            max_vertices: 3,
            ..settings()
        },
    ] {
        assert_eq!(stage_error(&d0(), change).code, "lattice.settings");
    }
    for value in [d0(), g0()] {
        let mut dec = Decider::default();
        let stage = Staged::new(
            &value,
            Settings {
                max_cells: 1,
                ..settings()
            },
            &mut dec,
        )
        .unwrap();
        assert_eq!(
            stage.export(&mut dec).unwrap_err().code,
            "field.print-budget"
        );
        let stage = Staged::new(
            &value,
            Settings {
                max_vertices: 4,
                accuracy: q("1/1000"),
                ..settings()
            },
            &mut dec,
        )
        .unwrap();
        let code = stage.export(&mut dec).unwrap_err().code;
        assert!(
            ["field.print-budget", "field.print-vertices"].contains(&code),
            "{code}"
        );
        let exact = || catallax_measure::Accuracy {
            relative: q("0"),
            absolute: q("0"),
        };
        let zero = mass::Accuracy {
            volume: exact(),
            centroid: exact(),
            inertia: exact(),
        };
        let code = mass::properties(
            stage.plan(),
            &zero,
            &mass::Options { max_cells: 1 },
            &mut dec,
        )
        .unwrap_err()
        .code;
        assert!(
            ["field.mass-budget", "field.mass-accuracy"].contains(&code),
            "{code}"
        );
    }
    let mut dec = Decider::default();
    let stage = Staged::new(&g0(), settings(), &mut dec).unwrap();
    assert_eq!(
        mass::properties(
            stage.plan(),
            &accuracy(),
            &mass::Options { max_cells: 0 },
            &mut dec
        )
        .unwrap_err()
        .code,
        "field.mass-options"
    );
}

/// The neck constructor of designer section 7: two blocks joined by a 0.02 mm square
/// neck inside the central beam, so region and field coincide.
fn neck() -> Value {
    let block = |id: &str, size: [&str; 3], translation: [&str; 3]| json!({"id":id,"op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":translation}});
    let union = |id: &str, a: &str, b: &str| json!({"id":id,"op":"boolean","kind":"union","target":a,"tool":b});
    json!({"format":"catallax-recipe/4","name":"Neck","steps":[
        block("low", ["1", "1", "4"], ["-1/2", "-1/2", "20"]),
        block("neck", ["1/50", "1/50", "2"], ["-1/100", "-1/100", "24"]),
        union("join", "low", "neck"),
        block("high", ["1", "1", "4"], ["-1/2", "-1/2", "26"]),
        union("region", "join", "high"),
        d0()["steps"][1].clone()],"result":"field"})
}

/// Neck-cap and Neck-detail: a typed resource or coverage refusal, or a single proven
/// component; never a hang or a disconnected claim about the connected solid.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d neck runs in release mode")]
fn the_neck_refuses_by_type_or_keeps_one_component() {
    let neck_settings = |max_cells| Settings {
        accuracy: q("1/1000"),
        minimum_detail: q("1/10"),
        max_cells,
        max_vertices: 1_000_000,
    };
    let (result, _) = with_operation_budget(KERNEL_WORK, || {
        let mut dec = Decider::default();
        Staged::new(&neck(), neck_settings(1), &mut dec)?.export(&mut dec)
    });
    let code = result.unwrap_err().code;
    assert!(
        ["field.print-budget", "field.print-unresolved"].contains(&code),
        "{code}"
    );
    println!("Neck-cap: {code}");
    let (result, work) = with_operation_budget(KERNEL_WORK, || {
        let mut dec = Decider::default();
        Staged::new(&neck(), neck_settings(2_000_000), &mut dec)?.export(&mut dec)
    });
    match result {
        Ok(output) => assert_eq!(output.verdict().components, 1),
        Err(error) => {
            println!(
                "Neck-detail: {} {} work {}",
                error.code, error.detail, work.used
            );
            assert!(
                ["field.print-budget", "field.print-unresolved"].contains(&error.code),
                "{}: {}",
                error.code,
                error.detail
            )
        }
    }
}

/// Disconnected-G: the documented corner piece stays in the definition; the export
/// refuses rather than deleting it, or a resource limit preempts the determination.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d exports run in release mode")]
fn a_disconnected_constant_gyroid_refuses_without_deleting_the_corner() {
    let mut value = g0();
    value["steps"][0]["size"] = json!(["20", "20", "20"]);
    value["steps"][0]["pose"]["translation"] = json!(["-10", "-10", "-5"]);
    value["steps"][1]["level"] = json!({"constant":"1/2"});
    let (result, _) = with_operation_budget(KERNEL_WORK, || {
        let mut dec = Decider::default();
        Staged::new(&value, settings(), &mut dec)?.export(&mut dec)
    });
    let error = result.unwrap_err();
    assert!(
        ["field.print-disconnected", "field.print-budget"].contains(&error.code),
        "{}: {}",
        error.code,
        error.detail
    );
    println!("Disconnected-G: {} {}", error.code, error.detail);
    // A constant level keeps its uniform route below a detail floor of a quarter of
    // the accuracy. With the V7c detail of 1/2 the certified route settles the
    // topology: the corner piece refuses the export.
    let (result, _) = with_operation_budget(KERNEL_WORK, || {
        let mut dec = Decider::default();
        Staged::new(
            &value,
            Settings {
                minimum_detail: q("1/2"),
                ..settings()
            },
            &mut dec,
        )?
        .export(&mut dec)
    });
    assert_eq!(result.unwrap_err().code, "field.print-disconnected");
}

/// After-field: a step that reads the terminal field is refused by the recipe's terminal
/// policy and by the stage over the already terminal body.
#[test]
fn a_step_after_the_terminal_field_is_refused() {
    for base in [d0(), g0()] {
        let mut value = base.clone();
        value["steps"].as_array_mut().unwrap().push(json!({"id":"after","op":"lattice","region":"field","cell":"cubic","size":"100",
            "frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"radius":{"constant":"1"}}));
        value["result"] = json!("after");
        let recipe = parse(&value);
        assert!(matches!(
            recipe.check_lattice_terminal(),
            Err(RecipeError::LatticeNotTerminal { .. })
        ));
        let mut dec = Decider::default();
        let terminal = Staged::new(&base, settings(), &mut dec).unwrap();
        let (_, spec) = split(&d0());
        let Field::Lattice(spec) = spec else {
            unreachable!()
        };
        let error = Stage::new(terminal.recipe(), spec, settings(), &mut dec).unwrap_err();
        assert_eq!(error.code, "lattice.cad-after-lattice");
    }
}

/// Negative-envelope: subtracting the field from its region reverses containment, so
/// the region is not an outer envelope for an explicit stand-in.
#[test]
fn a_subtracted_field_is_not_covered_by_its_region_envelope() {
    let mut value = d0();
    value["steps"].as_array_mut().unwrap().push(
        json!({"id":"negative","op":"boolean","kind":"subtract","target":"region","tool":"field"}),
    );
    value["result"] = json!("negative");
    let plan = Plan::compile_checked(&parse(&value), &mut Decider::default()).unwrap();
    assert_eq!(
        mass::check_region_envelope(&plan).unwrap_err().code,
        "field.stand-in-envelope-coverage"
    );
}

/// The app's route to the same region: a rectangle sketched on the source box's lower
/// face and extruded through its height, then the same terminal field.
fn sketched(value: &Value) -> Value {
    let region = &value["steps"][0];
    let [sx, sy, sz] = [0, 1, 2].map(|i| q(region["size"][i].as_str().unwrap()));
    let [tx, ty, tz] = [0, 1, 2].map(|i| q(region["pose"]["translation"][i].as_str().unwrap()));
    let corner = |x: &Rational, y: &Rational| json!([x.to_string(), y.to_string()]);
    let [x0, x1, y0, y1] = [tx.clone(), &tx + &sx, ty.clone(), &ty + &sy];
    let mut r = value.clone();
    r["steps"][0] = json!({"id":"profile","op":"sketch",
        "plane":{"rotation":["1","0","0","0"],"translation":["0","0",tz.to_string()]},
        "loops":[[{"kind":"line","id":"a","to":corner(&x1, &y0)},{"kind":"line","id":"b","to":corner(&x1, &y1)},
                  {"kind":"line","id":"c","to":corner(&x0, &y1)},{"kind":"line","id":"d","to":corner(&x0, &y0)}]]});
    let steps = r["steps"].as_array_mut().unwrap();
    steps.insert(
        1,
        json!({"id":"region","op":"extrude","profile":"profile","distance":sz.to_string()}),
    );
    r
}

/// Designer sections 3 and 8: the app's sketch-and-extrude region states the same
/// field, so it exports within the allowance and its moments meet the same references.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d moments run in release mode")]
fn a_sketched_region_states_the_same_field() {
    for (name, value) in [("D0", d0()), ("G0", g0())] {
        let value = sketched(&value);
        let ((bytes, output), work) = app_export(&value, settings());
        assert!(!work.exhausted && work.used < APP_WORK, "{name} {work:?}");
        assert_eq!(output.verdict().components, 1);
        check_file(&bytes, "1/10", tight_bounds(name));
        let (properties, work) = moments(&value);
        assert!(!work.exhausted && work.used < APP_WORK, "{name} {work:?}");
        if name == "D0" {
            check_d(name, &properties, 0, [q("0"), q("0"), q("0")]);
        } else {
            check_g(name, &properties, 0);
        }
    }
}
