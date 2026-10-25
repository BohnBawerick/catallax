//! Saved hole/chamfer expansions survive a bounded edit sequence like fresh replay.
use catallax_app::pattern::{FeatureInstance, expand_feature};
use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Operation, Recipe, Role, Step, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn decider() -> Decider {
    Decider::default()
}
#[allow(dead_code)]
fn wall(step: &str) -> Value {
    json!({"step":step,"role":{"hole":{"at":"c","part":"wall"}}})
}
fn source() -> Recipe {
    Recipe::parse(&json!({"format":catallax_recipe::FORMAT,"name":"saved-holes","steps":[
        {"id":"base","op":"block","size":["30","8","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"drill","op":"hole","target":"base","at":[{"id":"c","point":["3","4","4"],"direction":["0","0","-1"]}],
            "diameter":"2","depth":"through","tip":"flat"}
    ],"result":"drill"}).to_string())
    .unwrap()
}
fn graph() -> [StepId; 1] {
    [id("drill")]
}
// Copy `k` sits 5k along x; `previous` is the result it reads.
fn instance(k: i64, name: &str, previous: &str) -> FeatureInstance {
    FeatureInstance {
        placement: Isometry::translation(RVector3::new(
            Rational::from_i64(5 * k),
            0.into(),
            0.into(),
        )),
        steps: BTreeMap::from([(id("drill"), id(&format!("drill{name}")))]),
        inputs: BTreeMap::from([(id("base"), id(previous))]),
        origins: BTreeMap::from([(id("base"), id("base"))]),
    }
}
fn expand(recipe: &Recipe, copies: &[FeatureInstance]) -> Recipe {
    expand_feature(recipe, &graph(), &id("drill"), copies, &mut decider()).unwrap()
}
fn build(recipe: &Recipe) -> Built {
    let saved = Recipe::parse(&recipe.to_text()).unwrap();
    assert_eq!(&saved, recipe);
    let built = evaluate(&saved, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    built
}
fn volume(built: &Built) -> Rational {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    assert!(mass.volume.width() <= Rational::parse("1e-9").unwrap());
    mass.volume.lo
}
// Creating step of every origin on every surviving face, with the face's hole parts.
fn origin_steps(built: &Built) -> BTreeSet<String> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).to_vec())
        .map(|o| o.step.as_str().to_string())
        .collect()
}
// Per face: sorted origin descriptions, so split pieces and merges are compared by name.
fn names_of(built: &Built) -> Vec<String> {
    let mut names = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|f| {
            let mut o = built
                .model
                .face_origins(f)
                .iter()
                .map(|o| format!("{}:{:?}", o.step, o.role))
                .collect::<Vec<_>>();
            o.sort();
            o.join("|")
        })
        .collect::<Vec<_>>();
    names.sort();
    names
}
fn remove_copy(recipe: &Recipe, name: &str, next: Option<(&str, &str)>) -> Recipe {
    let mut r = recipe.clone();
    r.steps.retain(|s| s.id.as_str() != format!("drill{name}"));
    if let Some((reader, previous)) = next {
        let step = r.steps.iter_mut().find(|s| s.id == id(reader)).unwrap();
        let Operation::Hole { target, .. } = &mut step.op else {
            panic!("reader is a hole");
        };
        *target = id(previous);
    }
    r.result = r.steps.last().unwrap().id.clone();
    r.check().unwrap();
    r
}
fn through_slot(after: &str, x: &str, name: &str) -> Vec<Step> {
    let sketch = json!({"id":format!("p{name}"),"op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":[x,"7/2","0"]},"loops":[[
        {"kind":"line","id":"a","to":["2","0"]},{"kind":"line","id":"b","to":["2","1"]},
        {"kind":"line","id":"c","to":["0","1"]},{"kind":"line","id":"d","to":["0","0"]}]]});
    let cut = json!({"id":format!("s{name}"),"op":"extrude_through","target":after,"profile":format!("p{name}")});
    [sketch, cut]
        .map(|s| serde_json::from_value(s).unwrap())
        .into()
}

#[test]
fn deleted_and_reinserted_copies_keep_origins_and_match_fresh_replay() {
    let three = expand(
        &source(),
        &[
            instance(1, "1", "drill"),
            instance(2, "2", "drill1"),
            instance(3, "3", "drill2"),
        ],
    );
    let full = build(&three);
    let steps = origin_steps(&full);
    for name in ["drill", "drill1", "drill2", "drill3"] {
        assert!(steps.contains(name), "{name} is named in the full result");
    }
    let full_volume = volume(&full);

    // Delete the middle copy: the survivors keep their own origins, and no neighbour takes over.
    let gap = remove_copy(&three, "2", Some(("drill3", "drill1")));
    let built = build(&gap);
    let steps = origin_steps(&built);
    assert!(!steps.iter().any(|s| s.ends_with('2')));
    for name in ["drill", "drill", "drill1", "drill3"] {
        assert!(steps.contains(name), "{name} survives the deletion");
    }
    assert!(volume(&built) > full_volume);
    let fresh_gap = expand(
        &source(),
        &[instance(1, "1", "drill"), instance(3, "3", "drill1")],
    );
    assert_eq!(gap.to_text(), fresh_gap.to_text());
    assert_eq!(names_of(&built), names_of(&build(&fresh_gap)));

    // Reinsert a copy under a new id at the old place: it is a new feature, named separately.
    let mut back = gap.clone();
    let inserted = expand_feature(
        &gap,
        &graph(),
        &id("drill"),
        &[instance(2, "2b", "drill1")],
        &mut decider(),
    )
    .unwrap();
    let at = back
        .steps
        .iter()
        .position(|s| s.id == id("drill3"))
        .unwrap();
    back.steps
        .insert(at, inserted.steps.last().unwrap().clone());
    let Operation::Hole { target, .. } = &mut back.steps[at + 1].op else {
        panic!("reader is a hole");
    };
    *target = id("drill2b");
    back.result = id("drill3");
    back.check().unwrap();
    let built = build(&back);
    let steps = origin_steps(&built);
    assert!(steps.contains("drill2b") && !steps.contains("drill2"));
    assert_eq!(volume(&built), full_volume);
    let fresh = expand(
        &source(),
        &[
            instance(1, "1", "drill"),
            instance(2, "2b", "drill1"),
            instance(3, "3", "drill2b"),
        ],
    );
    assert_eq!(back.to_text(), fresh.to_text());
    assert_eq!(names_of(&built), names_of(&build(&fresh)));
}

#[test]
fn later_cut_splits_and_merges_hole_origins_like_fresh_replay() {
    let copies = [instance(1, "1", "drill")];
    let mut edited = expand(&source(), &copies);
    // A slot across the middle of hole 1 cuts its wall into two separate faces.
    let before = build(&edited);
    let slot = through_slot("drill1", "7", "x");
    assert_eq!(slot.len(), 2);
    edited.steps.extend(slot);
    edited.result = id("sx");
    let built = build(&edited);
    let names = names_of(&built);
    let pieces = |built: &Built, step: &str| {
        names_of(built)
            .iter()
            .filter(|n| n.contains(&format!("{step}:Hole")))
            .count()
    };
    assert_eq!(pieces(&before, "drill1"), 1);
    assert_eq!(pieces(&built, "drill1"), 2, "{names:?}");
    assert_eq!(pieces(&built, "drill"), pieces(&before, "drill"));
    // Fresh replay of the same final list gives the same recipe and names.
    let mut fresh = expand(&source(), &copies);
    fresh.steps.extend(through_slot("drill1", "7", "x"));
    fresh.result = id("sx");
    assert_eq!(edited.to_text(), fresh.to_text());
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f).to_vec())
            .all(|o| !matches!(o.role, Role::Copy { .. }))
    );
}

#[test]
fn reflected_copy_keeps_source_origins_and_missing_bindings_stay_atomic() {
    let mirror = Isometry::reflection(
        &RPoint3::new(Rational::from_i64(15), 0.into(), 0.into()),
        &RVector3::new(1.into(), 0.into(), 0.into()),
    )
    .unwrap();
    let mut copy = instance(0, "m", "drill");
    copy.placement = mirror;
    let mirrored = expand(&source(), &[copy.clone()]);
    let built = build(&mirrored);
    let steps = origin_steps(&built);
    assert!(steps.contains("drill") && steps.contains("drillm") && steps.contains("base"));
    let single = build(&source());
    assert!(volume(&built) < volume(&single));

    // Without the explicit input binding for the base the expansion refuses and changes nothing.
    let text = source().to_text();
    let mut missing = copy.clone();
    missing.inputs.clear();
    let err = expand_feature(
        &source(),
        &graph(),
        &id("drill"),
        &[missing],
        &mut decider(),
    )
    .unwrap_err();
    assert_eq!(err.code, "pattern.expansion");
    assert_eq!(source().to_text(), text);
    // A half-turn of 45 degrees has no rational spelling.
    let mut bad = copy;
    bad.placement = Isometry::around_axis(
        &RPoint3::origin(),
        &RVector3::new(0.into(), 0.into(), 1.into()),
        &Rational::from_i64(45),
        &mut decider(),
    )
    .unwrap();
    assert_eq!(
        expand_feature(&source(), &graph(), &id("drill"), &[bad], &mut decider())
            .unwrap_err()
            .code,
        "pattern.expansion-exact-parameter"
    );
}

// A chamfered hole copy on its own plate, with the plate named by an explicit origin binding.
// One face cannot hold two chamfered hole rims yet, so each copy gets its own plate.
fn plate_source() -> Recipe {
    let mut steps = vec![
        json!({"id":"plate","op":"block","size":["8","8","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}),
    ];
    for k in 1..=2 {
        steps.push(json!({"id":format!("plate{k}"),"op":"block","size":["8","8","4"],"pose":{"rotation":["1","0","0","0"],"translation":[(10 * k).to_string(),"0","0"]}}));
    }
    steps.push(json!({"id":"drill","op":"hole","target":"plate","at":[{"id":"c","point":["4","4","4"],"direction":["0","0","-1"]}],"diameter":"2","depth":"through","tip":"flat"}));
    steps.push(
        json!({"id":"bevel","op":"chamfer","body":"drill","distance":"1/4","edges":[{"between":[
        wall("drill"),{"step":"plate","role":"end"}]}]}),
    );
    Recipe::parse(
        &json!({"format":catallax_recipe::FORMAT,"name":"plates","steps":steps,"result":"bevel"})
            .to_string(),
    )
    .unwrap()
}
fn plate_copy(k: i64, name: &str) -> FeatureInstance {
    FeatureInstance {
        placement: Isometry::translation(RVector3::new(
            Rational::from_i64(10 * k),
            0.into(),
            0.into(),
        )),
        steps: BTreeMap::from([
            (id("drill"), id(&format!("drill{name}"))),
            (id("bevel"), id(&format!("bevel{name}"))),
        ]),
        inputs: BTreeMap::from([(id("plate"), id(&format!("plate{k}")))]),
        origins: BTreeMap::from([(id("plate"), id(&format!("plate{k}")))]),
    }
}
fn show(recipe: &Recipe, result: &str) -> Built {
    let mut r = recipe.clone();
    r.result = id(result);
    build(&r)
}
fn bevel_faces(built: &Built) -> usize {
    names_of(built)
        .iter()
        .filter(|n| n.contains("Bevel"))
        .count()
}

#[test]
fn chamfer_copies_keep_explicit_bindings_through_deletion_and_reinsertion() {
    let graph = [id("drill"), id("bevel")];
    let source = plate_source();
    let all = expand_feature(
        &source,
        &graph,
        &id("bevel"),
        &[plate_copy(1, "1"), plate_copy(2, "2")],
        &mut decider(),
    )
    .unwrap();
    let one = show(&all, "bevel1");
    let steps = origin_steps(&one);
    assert!(steps.contains("plate1") && steps.contains("bevel1") && steps.contains("drill1"));
    assert!(!steps.contains("plate2") && !steps.contains("bevel2") && !steps.contains("plate"));
    assert_eq!(bevel_faces(&one), 1);
    let two = show(&all, "bevel2");
    assert!(origin_steps(&two).contains("plate2"));
    assert_eq!(names_of(&one).len(), names_of(&two).len());

    // Delete copy 1: copy 2 is unchanged and nothing takes over the freed plate.
    let mut gap = all.clone();
    gap.steps
        .retain(|s| s.id != id("drill1") && s.id != id("bevel1"));
    gap.check().unwrap();
    let kept = show(&gap, "bevel2");
    assert_eq!(names_of(&kept), names_of(&two));
    let fresh = expand_feature(
        &source,
        &graph,
        &id("bevel"),
        &[plate_copy(2, "2")],
        &mut decider(),
    )
    .unwrap();
    assert_eq!(gap.steps, fresh.steps);

    // Reinsert under a new id on the same plate: a separate feature, same geometry.
    let back = expand_feature(
        &gap,
        &graph,
        &id("bevel"),
        &[plate_copy(1, "1b")],
        &mut decider(),
    )
    .unwrap();
    let again = show(&back, "bevel1b");
    let steps = origin_steps(&again);
    assert!(steps.contains("bevel1b") && !steps.contains("bevel1") && !steps.contains("bevel2"));
    assert_eq!(names_of(&again).len(), names_of(&one).len());
    assert_eq!(volume(&again), volume(&one));
}
