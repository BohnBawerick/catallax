//! Embeds the source commit without a runtime dependency on Git.
use std::process::Command;

fn main() {
    let commit = git(&["rev-parse", "HEAD"]);
    println!("cargo:rustc-env=CATALLAX_COMMIT={commit}");
    for name in ["HEAD", "packed-refs"] {
        println!(
            "cargo:rerun-if-changed={}",
            git(&["rev-parse", "--git-path", name])
        );
    }
    let branch = Command::new("git")
        .args(["symbolic-ref", "-q", "HEAD"])
        .output()
        .expect("git identifies the build branch");
    if branch.status.success() {
        let branch = String::from_utf8(branch.stdout).expect("git branch is UTF-8");
        println!(
            "cargo:rerun-if-changed={}",
            git(&["rev-parse", "--git-path", branch.trim()])
        );
    }
}

fn git(args: &[&str]) -> String {
    let output = Command::new("git")
        .args(args)
        .output()
        .expect("git is required to identify the kernel build");
    assert!(output.status.success(), "cannot identify the kernel commit");
    String::from_utf8(output.stdout)
        .expect("git output is UTF-8")
        .trim()
        .into()
}
