//! Replay an opt-in numerical recipe and publish approximate STEP after its checks.
use catallax_numbers::{Decider, with_operation_budget};
use catallax_step::numerical::{self, NumericalStepOptions};
use serde_json::json;
use std::{
    error::Error,
    fs::{self, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

static STAGE_ID: AtomicU64 = AtomicU64::new(0);
struct StagedFile(PathBuf);
impl Drop for StagedFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}
fn distinct_files(input: &Path, output: &Path) -> io::Result<()> {
    let source_path = fs::canonicalize(input)?;
    let source = fs::metadata(input)?;
    let destination = match fs::metadata(output) {
        Ok(value) => value,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error),
    };
    let same = fs::canonicalize(output)? == source_path;
    #[cfg(unix)]
    let same = {
        use std::os::unix::fs::MetadataExt;
        same || (source.dev() == destination.dev() && source.ino() == destination.ino())
    };
    #[cfg(not(unix))]
    let _ = (source, destination);
    if same {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input recipe and STEP output must be different files",
        ));
    }
    Ok(())
}
fn stage(output: &Path, bytes: &[u8]) -> io::Result<StagedFile> {
    let mut name = std::ffi::OsString::from(".");
    name.push(output.file_name().ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "STEP output needs a file name")
    })?);
    name.push(format!(
        ".{}-{}.pending",
        std::process::id(),
        STAGE_ID.fetch_add(1, Ordering::Relaxed)
    ));
    let path = output.with_file_name(name);
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)?;
    let staged = StagedFile(path);
    file.write_all(bytes)?;
    file.sync_all()?;
    Ok(staged)
}

fn run(input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    distinct_files(input, output)?;
    if let Err(error) = publish(input, output) {
        distinct_files(input, output)?;
        match fs::remove_file(output) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error.into()),
        }
        return Err(error);
    }
    Ok(())
}

fn publish(input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let saved = std::fs::read_to_string(input)?;
    let (result, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let built = catallax_features::replay_numerical(&saved, &mut dec)?;
        let artifact = numerical::build(built.body(), &NumericalStepOptions::default(), &mut dec)?;
        let staged = stage(output, artifact.bytes())?;
        artifact.check_written(&std::fs::read(&staged.0)?, &mut dec)?;
        Ok::<_, Box<dyn Error>>((built, staged))
    });
    let (built, staged) = result?;
    if work.exhausted {
        return Err("shared replay/export allowance exhausted".into());
    }
    distinct_files(input, output)?;
    fs::rename(&staged.0, output)?;
    println!(
        "{}",
        json!({"status":"numerical","exact":false,"step":"approximate",
        "tolerance_mm":built.body().evidence().tolerance_mm().to_string(),
        "checks":built.body().evidence().checks().len(),
        "repairs":built.body().evidence().repairs().len(),
        "history":built.body().history().iter().map(|stage| json!({"tolerance_mm":stage.options.tolerance_mm.to_string(),"checks":stage.evidence.checks().len(),"repairs":stage.evidence.repairs().len()})).collect::<Vec<_>>(),
        "origins":built.face_origins(),"work":work.used})
    );
    Ok(())
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 3 {
        return Err("usage: numerical_replay INPUT.json OUTPUT.step".into());
    }
    let output = Path::new(&args[2]);
    if let Err(error) = run(Path::new(&args[1]), output) {
        eprintln!(
            "{}",
            json!({"status":"refused","exact":false,"step":null,"detail":error.to_string()})
        );
        return Err(error);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_is_preserved_when_output_names_the_recipe() {
        let dir =
            std::env::temp_dir().join(format!("catallax-n04-identity-{}", std::process::id()));
        std::fs::create_dir(&dir).unwrap();
        let input = dir.join("recipe.json");
        let saved = include_str!("../../../evidence/phase1/N04/recipe.json");
        std::fs::write(&input, saved).unwrap();
        let result = run(&input, &input);
        let unchanged = std::fs::read_to_string(&input).unwrap() == saved;
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(
            unchanged,
            "replay publication overwrote its generating recipe"
        );
        assert!(result.is_err(), "same-file publication must refuse");
    }
    #[test]
    fn refused_replay_preserves_source_and_cleans_stale_output() {
        let dir = std::env::temp_dir().join(format!("catallax-n04-refusal-{}", std::process::id()));
        fs::create_dir(&dir).unwrap();
        let input = dir.join("recipe.json");
        let output = dir.join("result.step");
        fs::write(&input, "invalid JSON").unwrap();
        fs::write(&output, "stale").unwrap();
        assert!(run(&input, &output).is_err());
        assert_eq!(fs::read_to_string(&input).unwrap(), "invalid JSON");
        assert!(!output.exists());
        assert!(run(&input, &input).is_err());
        assert_eq!(fs::read_to_string(&input).unwrap(), "invalid JSON");
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
        fs::remove_dir_all(&dir).unwrap();
    }
    #[test]
    fn verified_stage_replaces_output_without_leftover_files() {
        let dir = std::env::temp_dir().join(format!("catallax-n04-publish-{}", std::process::id()));
        fs::create_dir(&dir).unwrap();
        let input = dir.join("recipe.json");
        let output = dir.join("result.step");
        let saved = include_str!("../../../evidence/phase1/N04/recipe.json");
        fs::write(&input, saved).unwrap();
        fs::write(&output, "stale").unwrap();
        run(&input, &output).unwrap();
        assert_eq!(fs::read_to_string(&input).unwrap(), saved);
        assert_eq!(
            fs::read(&output).unwrap(),
            include_bytes!("../../../evidence/phase1/N04/edited.step")
        );
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 2);
        fs::remove_dir_all(&dir).unwrap();
    }
    #[cfg(unix)]
    #[test]
    fn source_aliases_refuse_before_publication_or_cleanup() {
        use std::os::unix::fs::symlink;
        let dir = std::env::temp_dir().join(format!("catallax-n04-alias-{}", std::process::id()));
        fs::create_dir(&dir).unwrap();
        let input = dir.join("recipe.json");
        fs::write(&input, "invalid JSON").unwrap();
        let hard = dir.join("hard.step");
        let link = dir.join("link.step");
        fs::hard_link(&input, &hard).unwrap();
        symlink(&input, &link).unwrap();
        for output in [&hard, &link, &dir.join("./recipe.json")] {
            assert!(run(&input, output).is_err());
            assert_eq!(fs::read_to_string(&input).unwrap(), "invalid JSON");
            assert_eq!(fs::read_to_string(output).unwrap(), "invalid JSON");
        }
        fs::remove_dir_all(&dir).unwrap();
    }
}
