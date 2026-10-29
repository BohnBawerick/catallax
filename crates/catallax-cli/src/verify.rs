//! Persistent verification command. See docs/VERIFY-INTERFACE.md.
use serde_json::json;
use std::io::{self, BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    match args.as_slice() {
        [command] if command == "about" => {
            println!(
                "{}",
                json!({"name":"catallax-verify","request":"catallax-request/1",
                "response":"catallax-response/1","manifest":"catallax-manifest/3",
                "certified_recipe":"catallax-certified-recipe/1","certified_manifest":"catallax-manifest/4",
                "operation_log":"catallax-operation-log/1","question":"catallax-question/1",
                "field_print":"catallax-print/1","iges":"catallax-iges/1","stand_in":"catallax-stand-in/1",
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
                "usage: catallax-verify stream | about | replay QUESTION --out DIR\nstream: one catallax-request/1 JSON object per line; one catallax-response/1 per line\nreplay: one portable catallax-question/1, preserving its work allowance\nUse control-catallax.py --help for parts, batch, cache, compare and timing commands. Fresh proof evidence is required."
            );
            if args != ["--help"] {
                std::process::exit(64);
            }
        }
    }
    Ok(())
}
