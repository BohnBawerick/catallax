# Catallax + SolidDraft site

One static site for Catallax and SolidDraft at https://kernel.zigrawbohnbawerick.org/. The Pages workflow publishes this directory on site changes pushed to `main`; its custom domain is also recorded in `CNAME`. Domain DNS is managed separately.

## Preview

From the repository root:

```sh
nice -n 10 python3 -B -m http.server 8873 --bind 127.0.0.1 --directory site
```

Open http://127.0.0.1:8873/ in Chromium. No package install or build step is needed to serve the site.

## Publication rules

Keep Catallax's works-today list beside the shared sponsorship links. Name only capabilities on `main`, retain the not-yet line and label aspirations "in design". SolidDraft's section must distinguish its current SMLib backend from the planned Catallax-only cutover. [SolidDraft is public under MIT](https://github.com/BohnBawerick/soliddraft). It is not released for download yet; a release with the features is coming soon. For now it is research grade only.

The SolidDraft image is real application output, not interface artwork. `media/soliddraft-bracket-a3.svg` is an unmodified copy of the committed CLI drawing fixture at the revision recorded in `evidence.json`. Its source exit test invokes `soliddraft draw bracket.json --template iso-a3` with a fixed `SOURCE_DATE_EPOCH` and compares the SVG byte for byte. It was not regenerated for the site and is not evidence of Catallax integration. Keep the backend caption, source revision and checksum when replacing it.

Both projects use the same GitHub Sponsors and Ko-fi links. The manifesto from the repository README opens the page. Keep the special-thanks section near the bottom, with the Fornjot non-affiliation statement and links only to work the repository actually cites, the open-to-kernel-work contact, the plain link to https://zigrawbohnbawerick.org/ and the custom domain. Never claim the paid apps fund the work.

## Evidence and video capture

`evidence.json` records a pinned kernel revision, counting methods, the scope of historical robustness figures and the SolidDraft drawing's provenance. Counters and robustness figures are explicitly historical snapshots, not live measurements of `main`. Rebuild the kernel evidence when updating its figures:

```sh
CARGO_TARGET_DIR=/tmp/catallax-site-target nice -n 10 node site/render/collect-evidence.mjs
```

This runs a real workspace release build, numbers library tests and a compiled workspace test inventory. Rust lines include tests, comments and blank lines. The test count excludes ignored tests and doctests; it counts each compiled binary target. The transcript makes paths repository-relative. No source or test result is invented for the video.

Serve the site, then capture both films through a dedicated `chrome-devtools-axi` session and ffmpeg:

```sh
CHROME_DEVTOOLS_AXI_SESSION=catallax-render nice -n 10 node site/render/capture.mjs video-a.html site/media/kernel-assembly.mp4
CHROME_DEVTOOLS_AXI_SESSION=catallax-render nice -n 10 node site/render/capture.mjs video-b.html site/media/release-build.mp4
nice -n 10 ffmpeg -y -ss 5 -i site/media/kernel-assembly.mp4 -frames:v 1 site/media/kernel-assembly-poster.jpg
nice -n 10 ffmpeg -y -ss 11 -i site/media/release-build.mp4 -frames:v 1 site/media/release-build-poster.jpg
```

Set `SITE_URL` to use another preview port. `capture.mjs` renders each frame deterministically; the terminal film selects real transcript lines and compresses their timing. It is labelled as edited for time and is not a full-suite recording. The assembly film and interactive wireframes are art, not claims of generated models or shipped features.

Check desktop and phone layouts with reduced motion both enabled and disabled. Confirm the hero shape controls, orbit drag, all anatomy layers, video playback, external links and the Perth clock. Keep browser screenshots and scratch output outside the repository.

Third-party files retain their own terms; see `THIRD-PARTY.txt`.

## Multi-tool comparison

`comparison-results.json` is a public projection of a completed five-tool benchmark, not manually entered figures. Generate the page after importing quiet-host results:

```sh
python3 -B site/render/comparison.py --import-run /tmp/h2h/five-tools/head-to-head.json
python3 -B site/render/comparison.py
python3 -B site/render/comparison.py --check
```

The generator refuses incomplete target coverage, missing engines or a different machine class. The log-scale speed bars, per-part correctness table and per-tool target-range table all come from committed per-sample observations. It also updates the homepage's `comparison-summary` paragraph with Catallax's count of targets whose construction and edit completed in every measured sample, without changing the historical table. `--check` rejects a stale page or homepage summary. `comparison.css` uses the existing brand tokens and fonts. Preview `/comparison.html` at desktop and phone widths and exercise its target selector.

The page states that OpenCASCADE is both a competitor and the STEP importer. Properties and expected membership use independent analytic definitions, with unavailable property formulas explicitly unmeasured. Pinned releases, the Pali Ijo source-distribution exception, rebuild modes, failed phases and name-check limits remain visible. See the [rerun guide](../.claude/skills/verify-catallax/scripts/head-to-head/README.md) for preparation and commands.
