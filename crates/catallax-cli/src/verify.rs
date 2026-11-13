//! Persistent verification command. See docs/VERIFY-INTERFACE.md.
use serde_json::json;
use std::io::{self, BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<_> = std::env::args().skip(1).collect();
    // The kernel uses every available core unless a count is forced here or
    // by CATALLAX_THREADS. The count changes no answer, refusal or work count.
    if args.first().is_some_and(|flag| flag == "--threads") {
        let Some(count) = args
            .get(1)
            .and_then(|count| count.parse::<std::num::NonZeroUsize>().ok())
        else {
            eprintln!("--threads needs a count of at least 1");
            std::process::exit(64);
        };
        catallax_cli::set_threads(count);
        args.drain(..2);
    }
    match args.as_slice() {
        [command] if command == "about" => {
            println!(
                "{}",
                json!({"name":"catallax-verify","request":"catallax-request/1",
                "response":"catallax-response/1","manifest":"catallax-manifest/3",
                "certified_recipe":"catallax-certified-recipe/1","certified_manifest":"catallax-manifest/4",
                "exact_ruled_manifest":"catallax-manifest/5",
                "operation_log":"catallax-operation-log/1","question":"catallax-question/1",
                "field_print":"catallax-print/1","iges":"catallax-iges/1","stand_in":"catallax-stand-in/1",
                "seamed_step":"catallax-step/1",
                "sources":["question","recipe","field","print","iges","compare","name-edit"],"default_work":catallax_cli::verification::DEFAULT_WORK,
                "kernel":catallax_cli::about()})
            );
        }
        [command, file, flag, out] if command == "replay" && flag == "--out" => {
            let question: catallax_cli::verification::ExportedQuestion =
                serde_json::from_slice(&std::fs::read(file)?)?;
            let request = question.request(out.into())?;
            println!("{}", catallax_cli::verification::answer(&request)?);
        }
        [command, file, flag, limit] if command == "steps" && flag == "--max-work" => {
            let recipe = catallax_recipe::Recipe::parse(&std::fs::read_to_string(file)?)?;
            let limit = limit.parse::<std::num::NonZeroU64>()?.get();
            let rows = catallax_compare::independent_steps(
                &recipe,
                &catallax_features::ImportFiles::new(),
                limit,
            )?;
            let kernel = |k: &catallax_compare::KernelBuild| {
                use catallax_compare::Construction;
                let mut v = json!({"work":k.work,"limit":k.limit});
                match &k.construction {
                    Construction::Reference => {
                        v["status"] = json!("reference");
                    }
                    Construction::Built { faces, empty } => {
                        v["status"] = json!("built");
                        v["faces"] = json!(faces);
                        v["empty"] = json!(empty);
                    }
                    Construction::Refused {
                        blocked_at,
                        refusal,
                    } => {
                        v["status"] = json!("refused");
                        v["blocked_at"] = json!(blocked_at);
                        v["refusal"] = json!({"code":refusal.code,"detail":refusal.detail});
                    }
                }
                v
            };
            println!(
                "{}",
                serde_json::to_string_pretty(
                    &json!({"format":"catallax-independent-steps/1", "scope":"independent-prefix-construction", "agreement":"not-requested", "steps":rows.iter().map(|r| json!({"step":r.step,"operation":r.operation,"brep":kernel(&r.brep),"field":kernel(&r.field)})).collect::<Vec<_>>()})
                )?
            );
        }
        [command] if command == "stream" => {
            let stdin = io::stdin();
            let mut stdout = io::stdout().lock();
            for line in stdin.lock().lines() {
                let line = line?;
                let response =
                    match serde_json::from_str::<catallax_cli::verification::Request>(&line) {
                        Ok(request) => catallax_cli::verification::answer(&request)?,
                        Err(error) => json!({"format":"catallax-response/1","manifest":{
                        "format":"catallax-manifest/3","status":"invalid-question",
                        "refusal":{"code":"cli.invalid-request","detail":error.to_string()}},
                        "work":{"used":0,"limit":0,"exhausted":false},"wall_ms":0}),
                    };
                writeln!(stdout, "{response}")?;
                stdout.flush()?;
            }
        }
        _ => {
            eprintln!(
                "usage: catallax-verify [--threads N] stream | about | replay QUESTION --out DIR | steps RECIPE --max-work N\nstream: one catallax-request/1 JSON object per line; one catallax-response/1 per line\nreplay: one portable catallax-question/1, preserving its work allowance\n--threads: kernel threads per request; every available core by default, or CATALLAX_THREADS\nUse control-catallax.py --help for parts, batch, cache, compare and timing commands."
            );
            if args != ["--help"] {
                std::process::exit(64);
            }
        }
    }
    Ok(())
}
