#!/usr/bin/env python3
"""Add release policy and generated evidence links to GitHub's PR notes."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("generated", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("--tag", required=True)
parser.add_argument("--alpha", choices=("true", "false"), required=True)
args = parser.parse_args()
label = " alpha" if args.alpha == "true" else ""
args.out.write_text(
    f"# Catallax {args.tag}{label}\n\n"
    + ("The 0.x series is alpha.\n\n" if label else "")
    + "File formats and APIs may change between releases. Keep your source recipes and original files; do not rely on round-tripping saved kernel artifacts across versions.\n\n"
    + "macOS and Windows downloads have no developer identity signature. macOS universal executables carry only the ad-hoc signature required to run on Apple Silicon. GitHub build-provenance attestations are separate from platform code signing.\n\n"
    + "## Evidence and verification\n\n"
    + "Attached `replay-summary.json` audits 477 earlier cases against the selected baseline. Only committed, reviewed refusal upgrades and exact-fingerprint invalid-to-valid STEP corrections are allowed; their separate counts and approval-file hashes are recorded, and an approved difference is not byte equality. Candidate manifests and strict output bytes must still agree at one, two and four workers without exceptions. `target-parts.json` records this commit's target scoreboard, including refusals and checks not yet run. Neither is a claim that all CAD targets are supported.\n\n"
    + "Check downloads with `sha256sum --check SHA256SUMS` or an equivalent SHA-256 tool. The CycloneDX SBOM is `catallax-cli.cdx.xml`. Verify build provenance with `gh attestation verify FILE --repo BohnBawerick/catallax`.\n\n"
    + "## Merged pull requests\n\n"
    + args.generated.read_text()
)
