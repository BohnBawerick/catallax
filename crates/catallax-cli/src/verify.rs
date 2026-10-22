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
                "sources":["question","recipe","compare"],"default_work":catallax_cli::verification::DEFAULT_WORK,
                "kernel":catallax_cli::about()})
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
                "usage: catallax-verify stream | about\nstream: one catallax-request/1 JSON object per line; one catallax-response/1 per line\nUse control-catallax.py --help for batch, cache, compare and timing commands."
            );
            if args != ["--help"] {
                std::process::exit(64);
            }
        }
    }
    Ok(())
}
