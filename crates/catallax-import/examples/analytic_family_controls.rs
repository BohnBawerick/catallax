#[path = "../tests/common/mod.rs"]
mod common;
use common::Soup;
use std::{error::Error, fs, path::PathBuf};

fn main() -> Result<(), Box<dyn Error>> {
    let out = PathBuf::from(std::env::args().nth(1).ok_or("output directory required")?);
    if let Some(parent) = out.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&out)?;
    let mut band = Soup::new();
    band.half_torus(["0", "0", "0"], 2, "5", "1");
    let text = band
        .finish()
        .lines()
        .map(|line| {
            if line.contains("=TOROIDAL_SURFACE(") {
                let (prefix, _) = line.rsplit_once(',').expect("torus radius");
                format!("{prefix},1.000000001);")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(out.join("torus-band.step"), &text)?;
    println!(
        "torus-band.step {}",
        catallax_import::sha256_hex(text.as_bytes())
    );
    let mut horn = Soup::new();
    horn.torus(["0", "0", "0"], 2, "1", "1");
    let text = horn.finish();
    fs::write(out.join("horn.step"), &text)?;
    println!("horn.step {}", catallax_import::sha256_hex(text.as_bytes()));
    Ok(())
}
