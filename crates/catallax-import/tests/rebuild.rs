//! Rebuilding an import from its source and persisted record.

mod common;

use catallax_import::{
    Decision, ImportOptions, ImportRefusal, Imported, RepairRecord, SourceFormat, apply,
    apply_source, rebuild, rebuild_source, sha256_hex,
};
use catallax_numbers::Rational;
use catallax_topology::Model;
use common::*;

fn text(origin: [&str; 3], unit: Unit, by: &str) -> String {
    let mut b = BoxFixture::new(["10", "12", "14"]);
    b.origin = origin.map(r);
    b.unit = unit;
    b.shift = Some((1, [r(by), Rational::zero(), Rational::zero()]));
    box_step(b)
}

fn first(text: &str, decisions: &[Decision]) -> Imported {
    apply(
        &mut Model::new(),
        text.as_bytes(),
        "p.step",
        &ImportOptions::default(),
        decisions,
        None,
    )
    .unwrap()
}

fn reopen(record: &RepairRecord) -> RepairRecord {
    serde_json::from_str(&serde_json::to_string(record).unwrap()).unwrap()
}

fn again(text: &str, record: &RepairRecord) -> Result<Imported, ImportRefusal> {
    rebuild(
        &mut Model::new(),
        text.as_bytes(),
        &ImportOptions::default(),
        record,
        None,
    )
}

#[test]
fn an_approved_record_survives_reopen_and_reproduces_the_derivative() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.02");
    let done = first(
        &t,
        &[Decision::Accept {
            id: "g0".into(),
            tolerance_mm: "0.02".into(),
            approver: "jev".into(),
        }],
    );
    let record = reopen(&done.record);
    assert_eq!(record, done.record);
    let rebuilt = again(&t, &record).unwrap();
    assert_eq!(rebuilt.derivative, done.derivative);
    assert_eq!(rebuilt.record, done.record);
    assert_eq!(rebuilt.bodies.len(), 1);
}

#[test]
fn automatic_and_declined_records_rebuild_too() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let done = first(&t, &[]);
    assert_eq!(
        again(&t, &reopen(&done.record)).unwrap().derivative,
        done.derivative
    );
    let declined = first(&t, &[Decision::Decline { id: "g0".into() }]);
    let back = again(&t, &reopen(&declined.record)).unwrap();
    assert!(back.bodies.is_empty() && back.derivative.is_none());
}

#[test]
fn other_bytes_or_a_forged_record_are_refused() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let record = first(&t, &[]).record;
    let other = t.replace("'box'", "'bax'");
    assert!(matches!(
        again(&other, &record).unwrap_err(),
        ImportRefusal::HashMismatch { .. }
    ));
    let mut forged = record.clone();
    forged.derivative.as_mut().unwrap().sha256 = sha256_hex(b"x");
    assert!(matches!(
        again(&t, &forged).unwrap_err(),
        ImportRefusal::HashMismatch { expected, .. } if expected == sha256_hex(b"x")
    ));
    let mut moved = record.clone();
    moved.changes[0].max_movement_mm = "0.001".into();
    assert!(matches!(
        again(&t, &moved).unwrap_err(),
        ImportRefusal::ValidationFailed { .. }
    ));
    let mut format = record;
    format.format = "catallax-import-repair/9".into();
    assert!(matches!(
        again(&t, &format).unwrap_err(),
        ImportRefusal::Malformed { .. }
    ));
}

#[test]
fn moved_scaled_and_threaded_copies_repair_the_same_way() {
    // A translated copy and an inch copy of the same part repair by the same movement.
    let far = text(["100.5", "-30", "7"], Unit::Millimetre, "0.004");
    assert_eq!(first(&far, &[]).record.changes[0].max_movement_mm, "0.004");
    // 0.004 mm is 0.004/25.4 inch, a terminating decimal 0.00015748...: use an exact inch shift.
    let inch = text(["0", "0", "0"], Unit::Inch, "0.0001");
    let done = first(&inch, &[]);
    assert_eq!(done.record.changes[0].max_movement_mm, "0.00254");
    // Four threads give byte-identical derivatives and hashes.
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let runs: Vec<_> = (0..4)
        .map(|_| {
            let t = t.clone();
            std::thread::spawn(move || first(&t, &[]).derivative.unwrap())
        })
        .collect();
    let all: Vec<_> = runs.into_iter().map(|h| h.join().unwrap()).collect();
    assert!(all.windows(2).all(|w| w[0] == w[1]));
}

#[test]
fn ties_at_the_cap_repair_and_just_above_ask() {
    let at = first(&text(["0", "0", "0"], Unit::Millimetre, "0.01"), &[]);
    assert_eq!(at.record.changes[0].max_movement_mm, "0.01");
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.010001");
    assert!(matches!(
        apply(
            &mut Model::new(),
            t.as_bytes(),
            "p.step",
            &ImportOptions::default(),
            &[],
            None
        )
        .unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
}

/// A hand-made IGES unit cube. The x=1 face owns copies of its four vertices
/// and edges, translated by `by` in x, so its four seams are open by exactly that
/// distance; the other five faces share one vertex list and one edge list, as a
/// closed cube's do.
fn gapped_cube_iges(by: &str) -> Vec<u8> {
    const CORNERS: [[&str; 3]; 8] = [
        ["0", "0", "0"],
        ["1", "0", "0"],
        ["1", "1", "0"],
        ["0", "1", "0"],
        ["0", "0", "1"],
        ["1", "0", "1"],
        ["1", "1", "1"],
        ["0", "1", "1"],
    ];
    const FACES: [([usize; 4], [&str; 3]); 6] = [
        ([0, 3, 2, 1], ["0", "0", "-1"]),
        ([4, 5, 6, 7], ["0", "0", "1"]),
        ([0, 1, 5, 4], ["0", "-1", "0"]),
        ([1, 2, 6, 5], ["1", "0", "0"]),
        ([2, 3, 7, 6], ["0", "1", "0"]),
        ([3, 0, 4, 7], ["-1", "0", "0"]),
    ];
    // The x=1 face moves outward along its own normal, so its four seams open
    // without any face crossing another.
    const SHIFTED: usize = 3;

    struct Writer {
        entities: Vec<(u32, u32, String)>,
    }
    impl Writer {
        fn add(&mut self, kind: u32, form: u32, params: String) -> usize {
            self.entities.push((kind, form, params));
            2 * self.entities.len() - 1
        }
        fn bytes(&self) -> Vec<u8> {
            let record = |out: &mut String, text: &str, section: char, n: usize| {
                out.push_str(&format!("{text:<72}{section}{n:>7}\n"));
            };
            let mut out = String::new();
            record(&mut out, "hand-made fixture", 'S', 1);
            let global = "1H,,1H;,4Htest,7Hfixture,4Htest,4Htest,32,38,6,308,15,4Htest,1.0,2,2HMM,1,1.0,15H19700101.000000,1.0E-9,1.0E3,0H,0H,11,0,15H19700101.000000;";
            let global: Vec<char> = global.chars().collect();
            let g = global.chunks(72).count();
            for (i, chunk) in global.chunks(72).enumerate() {
                record(&mut out, &chunk.iter().collect::<String>(), 'G', i + 1);
            }
            let mut directory = String::new();
            let mut parameters = String::new();
            let mut line = 1;
            for (i, (kind, form, params)) in self.entities.iter().enumerate() {
                let sequence = 2 * i + 1;
                let transform = 0usize;
                let text: Vec<char> = format!("{kind},{params};").chars().collect();
                let lines = text.chunks(64).count();
                directory.push_str(&format!(
                    "{kind:>8}{line:>8}{:>8}{:>8}{:>8}{:>8}{transform:>8}{:>8}{:>8}D{sequence:>7}\n",
                    0, 0, 0, 0, 0, "00000000"
                ));
                directory.push_str(&format!(
                    "{kind:>8}{:>8}{:>8}{lines:>8}{form:>8}{:>8}{:>8}{:>8}{:>8}D{:>7}\n",
                    0,
                    0,
                    "",
                    "",
                    "",
                    0,
                    sequence + 1
                ));
                for chunk in text.chunks(64) {
                    let chunk: String = chunk.iter().collect();
                    parameters.push_str(&format!("{chunk:<64} {sequence:>7}P{line:>7}\n"));
                    line += 1;
                }
            }
            out.push_str(&directory);
            out.push_str(&parameters);
            let d = 2 * self.entities.len();
            record(
                &mut out,
                &format!("S{:>7}G{g:>7}D{d:>7}P{:>7}", 1, line - 1),
                'T',
                1,
            );
            out.into_bytes()
        }
    }

    let mut w = Writer {
        entities: Vec::new(),
    };
    let dec = |v: &Rational| v.to_decimal_exact().expect("decimal");
    let coords: Vec<String> = CORNERS
        .iter()
        .map(|c| {
            c.iter()
                .map(|v| dec(&Rational::parse(v).unwrap()))
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    let vertices = w.add(502, 1, format!("8,{}", coords.join(",")));
    let shift = Rational::parse(by).unwrap();
    let shifted = |c: usize| -> String {
        CORNERS[c]
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let p = Rational::parse(v).unwrap();
                dec(&if i == 0 { &p + &shift } else { p })
            })
            .collect::<Vec<_>>()
            .join(",")
    };
    // The five closed faces share their straight edges.
    let mut shared: Vec<(usize, usize, usize)> = Vec::new();
    for (index, (cycle, _)) in FACES.iter().enumerate() {
        if index == SHIFTED {
            continue;
        }
        for k in 0..4 {
            let (a, b) = (cycle[k], cycle[(k + 1) % 4]);
            let key = (a.min(b), a.max(b));
            if !shared.iter().any(|e| (e.0, e.1) == key) {
                let line = w.add(110, 0, format!("{},{}", coords[key.0], coords[key.1]));
                shared.push((key.0, key.1, line));
            }
        }
    }
    let rows: Vec<String> = shared
        .iter()
        .map(|(a, b, line)| format!("{line},{vertices},{},{vertices},{}", a + 1, b + 1))
        .collect();
    let shared_edges = w.add(504, 1, format!("{},{}", shared.len(), rows.join(",")));
    let mut shell = Vec::new();
    for (index, (cycle, normal)) in FACES.iter().enumerate() {
        let at = if index == SHIFTED {
            shifted(cycle[0])
        } else {
            coords[cycle[0]].clone()
        };
        let point = w.add(116, 0, at);
        let direction = w.add(123, 0, format!("{},{},{}", normal[0], normal[1], normal[2]));
        let plane = w.add(190, 0, format!("{point},{direction}"));
        let (list, uses) = if index == SHIFTED {
            let own: Vec<String> = cycle.iter().map(|&c| shifted(c)).collect();
            let vp = w.add(502, 1, format!("4,{}", own.join(",")));
            let mut rows = Vec::new();
            for k in 0..4 {
                let line = w.add(110, 0, format!("{},{}", own[k], own[(k + 1) % 4]));
                rows.push(format!("{line},{vp},{},{vp},{}", k + 1, (k + 1) % 4 + 1));
            }
            let el = w.add(504, 1, format!("4,{}", rows.join(",")));
            (el, (0..4).map(|k| (k, 1)).collect::<Vec<_>>())
        } else {
            let uses = (0..4)
                .map(|k| {
                    let (a, b) = (cycle[k], cycle[(k + 1) % 4]);
                    let index = shared
                        .iter()
                        .position(|e| (e.0, e.1) == (a.min(b), a.max(b)))
                        .expect("shared edge");
                    (index, usize::from(a < b))
                })
                .collect::<Vec<_>>();
            (shared_edges, uses)
        };
        let list_text = uses
            .iter()
            .map(|(i, f)| format!("0,{list},{},{f},0", i + 1))
            .collect::<Vec<_>>()
            .join(",");
        let lp = w.add(508, 1, format!("4,{list_text}"));
        let face = w.add(510, 1, format!("{plane},1,1,{lp}"));
        shell.push(format!("{face},1"));
    }
    let shell = w.add(514, 1, format!("6,{}", shell.join(",")));
    w.add(186, 0, format!("{shell},1,0"));
    w.bytes()
}

#[test]
fn iges_sources_rebuild_from_their_records() {
    let bytes = gapped_cube_iges("0.004");
    let done = apply_source(
        &mut Model::new(),
        SourceFormat::Iges,
        &bytes,
        "gap.igs",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap();
    assert_eq!(done.record.changes.len(), 1);
    assert_eq!(done.record.changes[0].max_movement_mm, "0.004");
    assert!(done.record.derivative.is_some());
    let rebuilt = rebuild_source(
        &mut Model::new(),
        SourceFormat::Iges,
        &bytes,
        &ImportOptions::default(),
        &reopen(&done.record),
        None,
    )
    .unwrap();
    assert_eq!(rebuilt.derivative, done.derivative);
    assert_eq!(rebuilt.record, done.record);
    assert_eq!(rebuilt.bodies.len(), 1);
}
