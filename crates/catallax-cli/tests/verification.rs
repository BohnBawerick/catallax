//! Exercise the persistent public command through stdin and output artifacts.
use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
    process::{Command, Stdio},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn question() -> Value {
    json!({"format":"exact-brep-question/2","id":"stream-test",
        "operands":{"a":{"shape":"block","size":["1","2","3"]}},
        "steps":[],"result":"a","limits":{"deadline_ms":0}})
}

#[test]
fn stream_ignores_wall_limit_and_resets_operation_budget() {
    let dir = std::env::temp_dir().join(format!(
        "catallax-stream-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&dir).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    for (i, limit) in [100_000_000u64, 0, 100_000_000].into_iter().enumerate() {
        writeln!(
            input,
            "{}",
            json!({"format":"catallax-request/1","source":"question",
            "problem":question(),"max_work":limit,"out":dir.join(i.to_string())})
        )
        .unwrap();
    }
    drop(input);
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let rows: Vec<Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0]["manifest"]["status"], "answered");
    assert_eq!(
        rows[1]["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert!(!dir.join("1/answer.step").exists());
    assert_eq!(rows[0]["work"], rows[2]["work"]);
    let mut a = rows[0]["manifest"].clone();
    let mut b = rows[2]["manifest"].clone();
    a.as_object_mut().unwrap().remove("stats");
    b.as_object_mut().unwrap().remove("stats");
    assert_eq!(a, b);
    assert_eq!(
        fs::read(dir.join("0/answer.step")).unwrap(),
        fs::read(dir.join("2/answer.step")).unwrap()
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn operation_scope_restores_after_panic_and_is_thread_local() {
    use catallax_numbers::{Decider, with_operation_budget};
    let _ = std::panic::catch_unwind(|| with_operation_budget(0, || panic!("test")));
    assert!(Decider::default().check_deadline().is_ok());
    let handles: Vec<_> = [0, 5]
        .into_iter()
        .map(|limit| {
            std::thread::spawn(move || {
                let (result, work) =
                    with_operation_budget(limit, || Decider::default().check_deadline());
                (result.is_ok(), work.used, work.exhausted)
            })
        })
        .collect();
    let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    assert_eq!(results, [(false, 0, true), (true, 1, false)]);
}
