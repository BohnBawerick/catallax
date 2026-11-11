//! Public finite-mouth clearance controls; full semantic worker comparisons.
use catallax_field::{FaceQuery, FaceState, Location, Plan, Query, Workspace};
use catallax_numbers::{Decider, OperationAllowance, Rational, Refusal, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub radius: String,
    pub placement: usize,
    pub result: Result<(Vec<Query>, FaceQuery, Vec<String>), Refusal>,
    pub work: u64,
    pub exhausted: bool,
}
fn recipe(radius: &str, placement: usize) -> Recipe {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../../../../parts/P5-flanged-housing.json")).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    for s in steps {
        let id = s["id"].as_str().unwrap().to_owned();
        if id == "filleted" {
            s["radius"] = radius.into();
        }
        if placement == 0 {
            continue;
        }
        let path = match id.as_str() {
            "section" => Some("plane"),
            "drill" | "sink" => Some("pose"),
            _ => None,
        };
        if let Some(path) = path {
            let a = s[path]["translation"].as_array().unwrap();
            let coords: Vec<Rational> = a.iter().map(|v| q(v.as_str().unwrap())).collect();
            s[path]["translation"] = serde_json::json!(
                place(
                    [coords[0].clone(), coords[1].clone(), coords[2].clone()],
                    placement
                )
                .map(|v| v.to_string())
            );
            if placement == 2 {
                s[path]["rotation"] = if id == "section" {
                    serde_json::json!(["1", "0", "0", "1"])
                } else {
                    serde_json::json!(["1", "1", "1", "1"])
                };
            }
        }
        let axis = if id == "housing" {
            Some(&mut s["axis"])
        } else if id == "hole-pattern" {
            Some(&mut s["layout"]["circular"]["axis"])
        } else {
            None
        };
        if let Some(axis) = axis {
            axis["point"] = serde_json::json!(
                place([q("0"), q("0"), q("0")], placement).map(|v| v.to_string())
            );
            if placement == 2 {
                axis["direction"] = serde_json::json!(["-1", "0", "0"]);
            }
        }
    }
    Recipe::parse(&value.to_string()).unwrap()
}
fn place([x, y, z]: [Rational; 3], placement: usize) -> [Rational; 3] {
    match placement {
        0 => [x, y, z],
        1 => [x + q("7"), y + q("11"), z - q("13")],
        2 => [-y, x, z],
        _ => unreachable!(),
    }
}
fn run_one(
    radius: &str,
    placement: usize,
) -> Result<(Vec<Query>, FaceQuery, Vec<String>), Refusal> {
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe(radius, placement), &mut dec)?;
    let mut storage = Workspace::default();
    let mut answers = Vec::new();
    // In the meridian at direction (3/5,4/5): ring centre (12.5,5.5), radius .5.
    // Boundary is the exact 3-4-5 radius point (12.2,5.1).
    let rows = [
        (["363/50", "51/10", "242/25"], Location::Inside),
        (["183/25", "51/10", "244/25"], Location::Boundary),
        (["369/50", "26/5", "246/25"], Location::Outside),
        (["17", "4", "0"], Location::Outside),
        (["0", "10", "0"], Location::Outside),
        (["10", "10", "0"], Location::Inside),
    ];
    for (xyz, want) in rows {
        let [x, y, z] = place(xyz.map(q), placement);
        let answer = plan.proven_checked(&RPoint3::new(x, y, z), &mut storage, &mut dec)?;
        assert_eq!(
            answer.location, want,
            "radius={radius} placement={placement} xyz={xyz:?}"
        );
        answers.push(answer);
    }
    let origin = plan
        .faces()
        .iter()
        .find(|f| matches!(f.origin.role, Role::Blend { .. }))
        .unwrap()
        .origin
        .clone();
    let [x, y, z] = place([q("183/25"), q("51/10"), q("244/25")], placement);
    let face = plan.face_query_checked(
        plan.select(&origin).unwrap(),
        &RPoint3::new(x, y, z),
        &mut storage,
        &mut dec,
    )?;
    assert_eq!(face.state, FaceState::OnFace);
    Ok((
        answers,
        face,
        plan.faces().iter().map(|f| f.origin.key()).collect(),
    ))
}
pub fn matrix(workers: Option<usize>, job_cap: u64) -> (Vec<Record>, usize, u64) {
    let cases = [
        ("1/2", 0),
        ("1/2", 1),
        ("1/2", 2),
        ("1", 0),
        ("2", 0),
        ("5/2", 0),
        ("7/2", 0),
        ("0", 0),
    ];
    let (result, operation) = with_operation_budget(20_000_000, || {
        let mut tokens = OperationAllowance::reserve(job_cap * cases.len() as u64).unwrap();
        let count = workers
            .unwrap_or_else(|| std::thread::available_parallelism().unwrap().get())
            .min(cases.len());
        let mut partitions: Vec<Vec<_>> = (0..count).map(|_| Vec::new()).collect();
        for (i, (radius, placement)) in cases.into_iter().enumerate() {
            partitions[i % count].push((i, radius, placement, tokens.split_off(job_cap).unwrap()));
        }
        let mut records = std::thread::scope(|scope| {
            let handles: Vec<_> = partitions
                .into_iter()
                .map(|jobs| {
                    scope.spawn(move || {
                        jobs.into_iter()
                            .map(|(index, radius, placement, token)| {
                                let (mut result, work) = token.run(|| run_one(radius, placement));
                                if work.exhausted && result.is_ok() {
                                    result = Err(Refusal::new(
                                        "m18.exhausted",
                                        "discarded result after child allowance exhaustion",
                                    ));
                                }
                                (
                                    index,
                                    Record {
                                        radius: radius.into(),
                                        placement,
                                        result,
                                        work: work.used,
                                        exhausted: work.exhausted,
                                    },
                                )
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            assert_eq!(handles.len(), count);
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        records.sort_by_key(|r| r.0);
        (records.into_iter().map(|r| r.1).collect(), count)
    });
    assert!(!operation.exhausted);
    (result.0, result.1, operation.used)
}
