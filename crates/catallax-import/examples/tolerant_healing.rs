//! Observe transactional source healing, without assigning a solid tier.
//! Usage: cargo run --release -p catallax-import --example tolerant_healing --
//! SOURCE.step EXPECTED_SHA256 REPORT.json
//! The complete read/heal query shares 20M work. Observer JSON is outside it.
use std::collections::BTreeMap;

use catallax_import::tolerant::{self, healing};
use catallax_numbers::{Decider, with_operation_budget};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("expected SOURCE.step EXPECTED_SHA256 REPORT.json".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let (attempt, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let source = tolerant::read_step(&bytes, &args[1], Default::default(), &mut dec)?;
        let result = healing::attempt(&source, &mut dec);
        Ok::<_, tolerant::Error>((source, result))
    });
    let mut report = serde_json::json!({"format":"catallax-healing-observation/1","source_sha256":args[1],"kernel_work_used":work.used,"work_exhausted":work.exhausted,"solid_validity_checked":false,"body_tiers_assigned":false});
    match attempt {
        Err(error) => {
            report["status"] = "source_read_refused".into();
            report["refusal"] = format!("{error:?}").into();
        }
        Ok((source, result)) => {
            report["source_bodies"]=source.bodies().iter().map(|body|serde_json::json!({"name":body.name(),"declared_solid":body.declared_solid(),"faces":body.shells().iter().map(|s|s.faces.len()).sum::<usize>(),"edges":body.edges().len(),"vertices":body.vertices().len()})).collect::<Vec<_>>().into();
            match result {
                Err(failure) => {
                    report["status"] = "healing_refused_source_retained".into();
                    report["refusal"] = format!("{:?}", failure.error).into();
                    report["failed_stage"] = failure.failed_stage.into();
                    report["failed_body"] = failure.body.into();
                    let measured = failure.preparation_work_used
                        + failure.stages.iter().map(|s| s.work_used).sum::<u64>()
                        + failure.publication_work_used;
                    assert!(measured <= work.used);
                    report["healing_work_used"] = measured.into();
                    report["preparation_work_used"] = failure.preparation_work_used.into();
                    report["publication_work_used"] = failure.publication_work_used.into();
                    report["stages"] = serde_json::to_value(failure.stages)?;
                    report["source_tolerances"] = serde_json::to_value(source.record())?;
                }
                Ok(healed) => {
                    let definition = healed
                        .definition(&mut Decider::default())
                        .map_err(|e| format!("{e:?}"))?;
                    let mut surfaces = BTreeMap::new();
                    let mut curves = BTreeMap::new();
                    let mut vertices = BTreeMap::new();
                    for body in source.bodies() {
                        for shell in body.shells() {
                            for face in &shell.faces {
                                surfaces.insert((body.name(), face.name.as_str()), &face.surface);
                            }
                        }
                        for edge in body.edges() {
                            curves.insert((body.name(), edge.name.as_str()), &edge.curve);
                        }
                        for vertex in body.vertices() {
                            vertices.insert(
                                (body.name(), vertex.source.name.as_str()),
                                &vertex.source.point,
                            );
                        }
                    }
                    for (body, face) in definition.named_faces() {
                        assert_eq!(
                            surfaces[&(body, face.name.as_str())],
                            &face.surface,
                            "source surface changed"
                        );
                        for boundary in &face.loops {
                            for edge in &boundary.edges {
                                assert_eq!(
                                    curves[&(body, edge.edge.as_str())],
                                    &edge.curve,
                                    "retained carrier changed"
                                );
                                for vertex in [&edge.start, &edge.end] {
                                    assert_eq!(
                                        vertices[&(body, vertex.name.as_str())],
                                        &vertex.point,
                                        "retained point changed"
                                    );
                                }
                            }
                        }
                    }
                    report["status"] = "candidate_complete".into();
                    report["retained_source_geometry_verified"] = true.into();
                    report["healing"] = serde_json::to_value(healed.report())?;
                    report["tolerances"] = serde_json::to_value(healed.tolerances())?;
                    let measured = healed.report().preparation_work_used
                        + healed
                            .report()
                            .stages
                            .iter()
                            .map(|s| s.work_used)
                            .sum::<u64>()
                        + healed.report().publication_work_used;
                    assert!(measured <= work.used);
                    report["healing_work_used"] = measured.into();
                }
            }
        }
    }
    std::fs::write(&args[2], serde_json::to_vec_pretty(&report)?)?;
    println!("{report}");
    Ok(())
}
