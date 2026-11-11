//! Worker threads run kernel code ahead of the request's own thread (see
//! docs/RUN-AHEAD.md). Nothing they do may change what that thread is charged,
//! so every cell, lock, static and thread-local in the crates they run is
//! listed here with the reason it is safe. A new one fails this test until it
//! has a line: decide which kind it is first.

use std::{collections::BTreeSet, fs, path::Path};

/// Why state a worker can reach cannot change the writer's work.
#[derive(Clone, Copy, Debug)]
enum Kind {
    /// Its fill spends no decision work, so who fills it changes no count.
    Uncharged,
    /// Its fill is charged; a worker asks `shared_caches()` and leaves it alone.
    Guarded,
    /// A worker has its own: per thread, emptied or replaced for each item.
    PerThread,
    /// Set by the caller or compiled out; no proof is ever stored in it.
    Setting,
}

const CRATES: [&str; 7] = [
    "catallax-numbers",
    "catallax-shapes",
    "catallax-topology",
    "catallax-intersect",
    "catallax-boolean",
    "catallax-recipe",
    "catallax-step",
];

// The run-ahead machinery itself: its tables, pool and thread roles.
const MACHINERY: [&str; 2] = ["catallax-numbers/src/ahead.rs", "catallax-step/src/uv.rs"];

const KNOWN: &[(&str, &str, Kind)] = &[
    // Certificates kept on a root. A worker records its cold proofs instead.
    (
        "catallax-numbers/src/algebraic.rs",
        "enclosures",
        Kind::Guarded,
    ),
    // The test clock; a thread that has one never uses workers.
    ("catallax-numbers/src/decider.rs", "CLOCK", Kind::Setting),
    // The request's counter and proof caches. A worker's are its own.
    (
        "catallax-numbers/src/operation.rs",
        "ACTIVE",
        Kind::PerThread,
    ),
    (
        "catallax-numbers/src/operation.rs",
        "ROOT_PROOFS",
        Kind::PerThread,
    ),
    (
        "catallax-numbers/src/operation.rs",
        "ROOT_LIST_PROOFS",
        Kind::PerThread,
    ),
    (
        "catallax-numbers/src/operation.rs",
        "ENCLOSURE_PROOFS",
        Kind::PerThread,
    ),
    (
        "catallax-numbers/src/operation.rs",
        "PI_PROOFS",
        Kind::PerThread,
    ),
    // The work profile, a development feature that turns workers off.
    ("catallax-numbers/src/profile.rs", "ROWS", Kind::Setting),
    ("catallax-numbers/src/profile.rs", "FRAMES", Kind::Setting),
    ("catallax-numbers/src/profile.rs", "SECTIONS", Kind::Setting),
    // Construction metadata of an expression node.
    ("catallax-numbers/src/real.rs", "coarse", Kind::Uncharged),
    ("catallax-numbers/src/real.rs", "rational", Kind::Uncharged),
    ("catallax-numbers/src/real.rs", "ZERO", Kind::Uncharged),
    ("catallax-numbers/src/real.rs", "ONE", Kind::Uncharged),
    (
        "catallax-numbers/src/real.rs",
        "NEGATIVE_ONE",
        Kind::Uncharged,
    ),
    // The operation log belongs to the thread that asked for it.
    ("catallax-numbers/src/trace.rs", "LOG", Kind::PerThread),
    (
        "catallax-numbers/src/ball/pi_constant.rs",
        "VALUE",
        Kind::Uncharged,
    ),
    (
        "catallax-shapes/src/rational_patch.rs",
        "normal",
        Kind::Uncharged,
    ),
    // Filling the spans is charged: the writer fills one cell, workers the other.
    (
        "catallax-shapes/src/rational_patch.rs",
        "beziers",
        Kind::Guarded,
    ),
    (
        "catallax-shapes/src/rational_patch.rs",
        "ahead_beziers",
        Kind::Guarded,
    ),
    (
        "catallax-shapes/src/rolling_envelope.rs",
        "grid",
        Kind::Guarded,
    ),
    (
        "catallax-shapes/src/rolling_envelope.rs",
        "moments",
        Kind::Guarded,
    ),
    (
        "catallax-shapes/src/rolling_saddle.rs",
        "grid",
        Kind::Guarded,
    ),
    (
        "catallax-shapes/src/rolling_saddle.rs",
        "turns",
        Kind::Uncharged,
    ),
    (
        "catallax-shapes/src/store.rs",
        "EquationCells",
        Kind::Uncharged,
    ),
    // Local to one moment computation, never shared.
    (
        "catallax-shapes/src/rolling_envelope/moments.rs",
        "known",
        Kind::PerThread,
    ),
    ("catallax-intersect/src/chart.rs", "balls", Kind::Guarded),
    // Uncharged for a rational constant; guarded where an enclosure fills it.
    ("catallax-intersect/src/chart.rs", "interval", Kind::Guarded),
    ("catallax-intersect/src/lib.rs", "jets", Kind::Guarded),
    ("catallax-intersect/src/lib.rs", "program", Kind::Uncharged),
    // A plan hit skips sign decisions. Each thread has its own plans, and a
    // recorded computation that reads them is dropped.
    ("catallax-intersect/src/pair.rs", "PLANS", Kind::PerThread),
    // Pure rational arithmetic, never charged.
    ("catallax-intersect/src/solve.rs", "BASES", Kind::Uncharged),
    (
        "catallax-intersect/src/span_hull.rs",
        "controls",
        Kind::Guarded,
    ),
    ("catallax-step/src/fit.rs", "WEIGHTS", Kind::Uncharged),
];

const MARKERS: [&str; 9] = [
    "Mutex<",
    "RwLock<",
    "OnceLock<",
    "LazyLock<",
    "OnceCell<",
    "RefCell<",
    "Cell<",
    "Atomic",
    "static mut ",
];

// The name a line declares: a static, a field or a tuple struct.
fn declared(line: &str) -> Option<&str> {
    let line = line.trim_start();
    let line = line.strip_prefix("pub(crate) ").unwrap_or(line);
    let line = line.strip_prefix("pub ").unwrap_or(line);
    let line = line.strip_prefix("static ").unwrap_or(line);
    let line = line.strip_prefix("mut ").unwrap_or(line);
    let line = line.strip_prefix("struct ").unwrap_or(line);
    let end = line.find([':', '('])?;
    let name = &line[..end];
    (!name.is_empty() && name.chars().all(|c| c.is_alphanumeric() || c == '_')).then_some(name)
}

fn scan(root: &Path, dir: &Path, found: &mut BTreeSet<(String, String)>) {
    let mut entries: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect();
    entries.sort();
    for path in entries {
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        if path.is_dir() {
            if name != "tests" {
                scan(root, &path, found);
            }
            continue;
        }
        if !name.ends_with(".rs") || name == "tests.rs" {
            continue;
        }
        let file = path
            .strip_prefix(root)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        if MACHINERY.contains(&file.as_str()) {
            continue;
        }
        let text = fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        for (index, line) in lines.iter().enumerate() {
            // Unit tests at the end of a file keep what state they like.
            if line.trim() == "#[cfg(test)]"
                && lines
                    .get(index + 1)
                    .is_some_and(|next| next.trim_start().starts_with("mod "))
            {
                break;
            }
            if line.trim_start().starts_with("//") {
                continue;
            }
            if MARKERS.iter().any(|marker| line.contains(marker))
                && let Some(name) = declared(line)
            {
                found.insert((file.clone(), name.to_string()));
            }
        }
    }
}

#[test]
fn every_shared_cell_a_worker_can_reach_is_accounted_for() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let mut found = BTreeSet::new();
    for name in CRATES {
        scan(crates, &crates.join(name).join("src"), &mut found);
    }
    let known: BTreeSet<(String, String)> = KNOWN
        .iter()
        .map(|(file, name, _)| (file.to_string(), name.to_string()))
        .collect();
    let unlisted: Vec<_> = found.difference(&known).collect();
    assert!(
        unlisted.is_empty(),
        "state a run-ahead worker can reach, with no reason listed here: {unlisted:?}. \
         If filling it spends decision work, guard it with catallax_numbers::shared_caches() \
         (docs/RUN-AHEAD.md, section 4), then add it to KNOWN."
    );
    let gone: Vec<_> = known.difference(&found).collect();
    assert!(
        gone.is_empty(),
        "listed here but no longer in the source: {gone:?}"
    );
    // Every kind is in use; a list of one kind would mean the scan is blind.
    for kind in [
        Kind::Uncharged,
        Kind::Guarded,
        Kind::PerThread,
        Kind::Setting,
    ] {
        assert!(
            KNOWN
                .iter()
                .any(|(_, _, other)| std::mem::discriminant(other) == std::mem::discriminant(&kind)),
            "{kind:?}"
        );
    }
}
