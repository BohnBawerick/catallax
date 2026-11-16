//! Edit one body from an actual imported STEP file and export before/after B-reps.
//! Usage: direct_edits SOURCE SHA OUT MODE [FACE [numeric arguments...]]
//! MODE: inspect, delete-face, delete-fillet, delete-hole, move, offset, replace.
//! Each read/open, edit, and aggregate export has an explicit 100M work ceiling.
use catallax_import::{
    Surface,
    direct::{DirectSolid, Edit},
    numerical::ImportError,
    sha256_hex, tolerant,
};
use catallax_numbers::{Decider, Rational, current_operation_work, with_operation_budget};
use catallax_recipe::StepId;
use serde_json::json;
use std::path::PathBuf;

fn error(error: ImportError) -> serde_json::Value {
    match error {
        ImportError::Source(r) => serde_json::to_value(r).unwrap(),
        ImportError::Numerical(r) => json!({"code":r.code,"detail":r.detail}),
    }
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let bytes = std::fs::read(&args[0]).unwrap();
    assert_eq!(sha256_hex(&bytes), args[1], "immutable census source SHA");
    let out = PathBuf::from(&args[2]);
    std::fs::create_dir_all(&out).unwrap();
    let save = |name: &str, value: &serde_json::Value| {
        std::fs::write(out.join(name), serde_json::to_vec_pretty(value).unwrap()).unwrap();
        println!("{value}");
    };
    if args[3] == "inspect" {
        let (result, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let source = tolerant::read_step(&bytes, &args[1], Default::default(), &mut dec)?;
            source.definition(&mut dec)
        });
        match result {
            Ok(file) => save(
                "inventory.json",
                &json!({"inspection_only":true,"work":work.used,
                "solids":file.solids.iter().map(|s| json!({"name":s.name,"faces":s.faces.iter().map(|f| json!({
                    "key":f.name,"surface":format!("{:?}",f.surface),"same_sense":f.same_sense,
                    "bounds":f.loops.iter().map(|b|json!({"outer":b.outer,"edges":b.edges.iter().map(|e|json!({
                        "key":e.edge,"vertices":[e.start.name,e.end.name],"curve":format!("{:?}",e.curve)})).collect::<Vec<_>>()
                    })).collect::<Vec<_>>()
                })).collect::<Vec<_>>()})).collect::<Vec<_>>() }),
            ),
            Err(e) => save(
                "inventory.json",
                &json!({"refusal":error(e),"work":work.used}),
            ),
        }
        return;
    }
    let (opened, open_work) = with_operation_budget(100_000_000, || {
        tolerant::opened::open_step(
            &bytes,
            &args[1],
            &StepId::new("direct-edit-census").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let opened = match opened {
        Ok(value) => value,
        Err(e) => {
            save(
                "report.json",
                &json!({"stage":"open","refusal":error(e),"work":open_work.used}),
            );
            return;
        }
    };
    save(
        "open.json",
        &json!({"work":open_work.used,"exhausted":open_work.exhausted,
        "tiers":opened.bodies().iter().map(|b|b.tier()).collect::<Vec<_>>()}),
    );
    let parse = |index: usize| Rational::parse(&args[index]).unwrap();
    let key = args[4].clone();
    let approximate_source = opened.bodies()[0].approximate();
    let mut stage = "prepare";
    let (result, edit_work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let before = DirectSolid::from_opened(std::sync::Arc::new(opened), 0, &mut dec)?;
        let prepare_work = current_operation_work().used;
        stage = "edit";
        save("prepared.json", &json!({"work":prepare_work}));
        let edit = match args[3].as_str() {
            "delete-face" => Edit::DeleteFace { face: key },
            "delete-fillet" => Edit::DeleteFillet { face: key },
            "delete-hole" => Edit::DeleteHole { face: key },
            "move" => Edit::MoveFace {
                face: key,
                translation_mm: [parse(5), parse(6), parse(7)],
            },
            "offset" => Edit::OffsetFace {
                face: key,
                distance_mm: parse(5),
            },
            "replace" => {
                let mut surface = before
                    .faces()
                    .iter()
                    .find(|f| f.name == key)
                    .unwrap()
                    .surface
                    .clone();
                match &mut surface {
                    Surface::Plane(frame) => {
                        frame.origin[0] = &frame.origin[0] + &parse(5);
                        frame.origin[1] = &frame.origin[1] + &parse(6);
                        frame.origin[2] = &frame.origin[2] + &parse(7);
                    }
                    Surface::Cylinder(frame, radius) => {
                        frame.origin[0] = &frame.origin[0] + &parse(5);
                        frame.origin[1] = &frame.origin[1] + &parse(6);
                        frame.origin[2] = &frame.origin[2] + &parse(7);
                        *radius = parse(8);
                    }
                    _ => panic!("observer replacement family"),
                }
                Edit::ReplaceFace {
                    face: key,
                    surface: Box::new(surface),
                }
            }
            _ => panic!("unknown edit"),
        };
        let after = before.apply(edit, Default::default(), Default::default(), &mut dec)?;
        Ok::<_, ImportError>((before, after, prepare_work))
    });
    let (before, after, prepare_work) = match result {
        Ok(value) => value,
        Err(e) => {
            save(
                "report.json",
                &json!({"stage":stage,"refusal":error(e),"work":edit_work.used,"exhausted":edit_work.exhausted}),
            );
            return;
        }
    };
    assert_eq!(before.source_sha256(), after.source_sha256());
    assert_eq!(before.source_faces(), after.source_faces());
    assert!(before.history().is_empty());
    let (exports, export_work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let mut files = Vec::new();
        for (name, solid) in [("before.step", &before), ("after.step", &after)] {
            let bytes = match solid.body() {
                tolerant::opened::BodyView::Certified { model, body } => {
                    let mut bytes = Vec::new();
                    let mut options = catallax_step::StepOptions::default();
                    if approximate_source {
                        options.product_name = "catallax approximate source repair".into();
                    }
                    catallax_step::write_step(model, &[body], &options, &mut dec, &mut bytes)?;
                    bytes
                }
                tolerant::opened::BodyView::Numerical(body) => {
                    let file =
                        catallax_step::numerical::build(body, &Default::default(), &mut dec)?;
                    file.check_written(file.bytes(), &mut dec)?;
                    file.bytes().to_vec()
                }
                _ => unreachable!("the editor retains only admitted material"),
            };
            std::fs::write(out.join(name), &bytes).unwrap();
            assert_eq!(std::fs::read(out.join(name)).unwrap(), bytes);
            files.push(json!({"file":name,"sha256":sha256_hex(&bytes),"bytes":bytes.len()}));
        }
        Ok::<_, catallax_numbers::Refusal>(files)
    });
    let exports = exports.unwrap();
    let tolerant::opened::BodyView::Numerical(material) = after.body() else {
        unreachable!("edits are numerical")
    };
    let admission = after.admission().expect("completed derivative admission");
    save(
        "report.json",
        &json!({"outcome":"exported","authority":"numerical","source_sha256":after.source_sha256(),
        "tolerance_mm":material.options().tolerance_mm.to_string(),"work_limit":100_000_000,
        "open_work":open_work.used,"prepare_work":prepare_work,"edit_total_work":edit_work.used,
        "export_work":export_work.used,"exhausted":edit_work.exhausted||export_work.exhausted,
        "before_faces":before.faces().len(),"after_faces":after.faces().len(),
        "history":format!("{:?}",after.history()),"source_pair_checks":admission.source_pair_checks(),
        "topology_repairs":admission.topology_repairs(),"files":exports}),
    );
}
