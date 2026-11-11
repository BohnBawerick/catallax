//! The `catallax` and legacy `exact-brep` programs. See `docs/JUDGE-INTERFACE.md`.

use std::{
    ffi::OsString,
    fs,
    io::{self, Write},
    panic::{AssertUnwindSafe, catch_unwind},
    path::Path,
    process::ExitCode,
    time::Instant,
};

use catallax_cli::{Status, failure, parse_question, render_manifest};

fn main() -> ExitCode {
    let started = Instant::now();
    let mut args: Vec<OsString> = std::env::args_os().skip(1).collect();
    // The kernel uses every available core unless a count is forced here or
    // by CATALLAX_THREADS. The count changes no answer, refusal or work count.
    if args.first().is_some_and(|flag| flag == "--threads") {
        let Some(count) = args
            .get(1)
            .and_then(|count| count.to_str())
            .and_then(|count| count.parse::<std::num::NonZeroUsize>().ok())
        else {
            eprintln!("--threads needs a count of at least 1");
            return ExitCode::from(64);
        };
        catallax_cli::set_threads(count);
        args.drain(..2);
    }
    let code = match args.as_slice() {
        [command] if command == "about" => {
            match writeln!(io::stdout().lock(), "{}", catallax_cli::about()) {
                Ok(()) => 0,
                Err(error) => {
                    eprintln!("{error}");
                    74
                }
            }
        }
        [command, file] if command == "check" => check(Path::new(file)),
        [command, file, flag, dir] if command == "answer" && flag == "--out" => {
            answer(Path::new(file), Path::new(dir), started)
        }
        _ => {
            eprintln!(
                "usage: exact-brep [--threads N] answer <question.json> --out <dir>\n       exact-brep check <question.json>\n       exact-brep about"
            );
            64
        }
    };
    ExitCode::from(code)
}

fn check(file: &Path) -> u8 {
    let result = fs::read_to_string(file)
        .map_err(|error| error.to_string())
        .and_then(|text| parse_question(&text).map_err(|error| error.detail));
    match result {
        Ok(_) => 0,
        Err(error) => {
            eprintln!("{error}");
            65
        }
    }
}

fn answer(file: &Path, dir: &Path, started: Instant) -> u8 {
    if let Err(error) = catallax_cli::prepare_output(dir) {
        eprintln!("{error}");
        return 74;
    }
    let mut question_id = None;
    let mut format = catallax_cli::MANIFEST_FORMAT;
    let result = catch_unwind(AssertUnwindSafe(|| match fs::read_to_string(file) {
        Ok(text) => match parse_question(&text) {
            Ok(question) => {
                question_id = Some(question.id.clone());
                format = question.manifest_format();
                catallax_cli::answer_since(&question, dir, started)
            }
            Err(error) => {
                question_id.clone_from(&error.question_id);
                Ok(failure(
                    error.question_id,
                    Status::InvalidQuestion,
                    "cli.invalid-question",
                    error.detail,
                ))
            }
        },
        Err(error) => Ok(failure(
            None,
            Status::InvalidQuestion,
            "cli.invalid-question",
            error.to_string(),
        )),
    }));
    let (manifest, code) = match result {
        Ok(Ok(manifest)) => (manifest, 0),
        Ok(Err(error)) => {
            eprintln!("{error}");
            return 74;
        }
        Err(panic) => {
            if let Err(error) = catallax_cli::prepare_output(dir) {
                eprintln!("{error}");
                return 74;
            }
            let detail = panic
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| panic.downcast_ref::<&str>().copied())
                .unwrap_or("kernel panicked");
            let mut manifest = failure(
                question_id,
                Status::InternalError,
                "cli.internal-error",
                detail,
            );
            manifest.format = format;
            (manifest, 70)
        }
    };
    match fs::write(dir.join("manifest.json"), render_manifest(&manifest)) {
        Ok(()) => code,
        Err(error) => {
            eprintln!("{error}");
            74
        }
    }
}
