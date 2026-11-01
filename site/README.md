# Catallax + SolidDraft site

One static site for Catallax and SolidDraft at https://kernel.zigrawbohnbawerick.org/. The Pages workflow publishes this directory on site changes pushed to `main`; its custom domain is also recorded in `CNAME`. Domain DNS is managed separately.

The bare-bones design is published here, with `why.html` and `thanks.html` as its long-form pages. The review copy under `.lavish/` has been removed. `bare.css`, `bare.js`, `field.js` and `js/wire.js` drive the homepage; `read.css` and `read.js` drive the long-form pages. `css/site.css` remains a dependency of the comparison page.

## Preview

From the repository root:

```sh
nice -n 10 python3 -B -m http.server 8873 --bind 127.0.0.1 --directory site
```

Open http://127.0.0.1:8873/ in Chromium. No package install or build step is needed to serve the site.

## Publication rules

The site states the intended direction of Catallax and SolidDraft. Keep the owner's product claims as written, including SolidDraft running natively on Catallax; do not edit the copy per release or require current-backend, "Not yet" or "in design" distinctions. [SolidDraft is public under MIT](https://github.com/BohnBawerick/soliddraft).

The SolidDraft image is real application output, not interface artwork. `media/soliddraft-bracket-a3.svg` is an unmodified copy of the committed CLI drawing fixture at the revision recorded in `evidence.json`. Its source exit test invokes `soliddraft draw bracket.json --template iso-a3` with a fixed `SOURCE_DATE_EPOCH` and compares the SVG byte for byte. It was not regenerated for the site and is not evidence of Catallax integration. Keep the source revision and checksum when replacing it.

Both projects use the same GitHub Sponsors and Ko-fi links. Keep the special-thanks section near the bottom, with the Fornjot non-affiliation statement and links only to work the repository actually cites, the open-to-kernel-work contact, the plain link to https://zigrawbohnbawerick.org/ and the custom domain. Never claim the paid apps fund the work.

## Motion rules

The page must be complete without `bare.js`. Motion only decorates content that is already readable and reachable.

- Reveals animate opacity and transform, never `visibility` or `autoAlpha`. Hidden copy must stay in the tab order and the accessibility tree.
- The hero headline's entrance is CSS, so it runs from first paint. Do not hide text after it has painted.
- The kernel anatomy section pins only at `min-height: 640px`. Below that it stays a plain list, which also serves no-JS and high zoom.
- The fixed `.motion` button pauses all automatic motion and video playback. Reduced motion disables automatic playback and exposes native video controls.

Design tokens and rules are in [DESIGN.md](DESIGN.md), and product facts are in [PRODUCT.md](PRODUCT.md).

## Evidence and video capture

`evidence.json` records a pinned kernel revision, counting methods, the scope of historical robustness figures and the SolidDraft drawing's provenance. Counters and robustness figures are explicitly historical snapshots, not live measurements of `main`. Rebuild the kernel evidence when updating its figures:

```sh
CARGO_TARGET_DIR=/tmp/catallax-site-target nice -n 10 node site/render/collect-evidence.mjs
```

This runs a real workspace release build, numbers library tests and a compiled workspace test inventory. It writes the evidence and transcript and updates the three repository counts in the homepage's `data-metric` elements. Rust lines include tests, comments and blank lines. The test count excludes ignored tests and doctests; it counts each compiled binary target. The transcript makes paths repository-relative. No source or test result is invented for the video. Do not run this collector merely to publish a design change.

Serve the site, then capture both films through a dedicated `chrome-devtools-axi` session and ffmpeg:

```sh
CHROME_DEVTOOLS_AXI_SESSION=catallax-render nice -n 10 node site/render/capture.mjs video-a.html site/media/kernel-assembly.mp4
CHROME_DEVTOOLS_AXI_SESSION=catallax-render nice -n 10 node site/render/capture.mjs video-b.html site/media/release-build.mp4
nice -n 10 ffmpeg -y -ss 5 -i site/media/kernel-assembly.mp4 -frames:v 1 site/media/kernel-assembly-poster.jpg
nice -n 10 ffmpeg -y -ss 11 -i site/media/release-build.mp4 -frames:v 1 site/media/release-build-poster.jpg
```

Set `SITE_URL` to use another preview port. `capture.mjs` renders each frame deterministically; the terminal film selects real transcript lines and compresses their timing. It is labelled as edited for time and is not a full-suite recording. The assembly film and interactive wireframes are art, not claims of generated models or shipped features.

Check desktop and phone layouts with reduced motion both enabled and disabled. Confirm the hero shape controls, orbit drag, B-rep/field wipe, all anatomy layers, the Why and Thanks navigation, the pause control, video playback, external links and the Perth clock. Keep browser screenshots and scratch output outside the repository.

Third-party files retain their own terms; see `THIRD-PARTY.txt`.

## Multi-tool comparison

`comparison-results.json` is a public projection of a completed five-tool benchmark, not manually entered figures. Generate the page after importing quiet-host results:

```sh
python3 -B site/render/comparison.py --import-run /tmp/h2h/five-tools/head-to-head.json
python3 -B site/render/comparison.py
python3 -B site/render/comparison.py --check
```

The generator refuses incomplete target coverage, missing engines or a different machine class. The log-scale speed bars, per-part correctness table and per-tool target-range table all come from committed per-sample observations. It also generates the homepage's `comparison-summary` figure, including each tool's completed construction and edit builds, proportional bars, the scheduled sample count and the full-comparison link. Completion is not correctness; an edit skipped after a failed construction remains incomplete. Never hand-edit generated numbers. `--check` rejects a stale page or homepage summary. `comparison.css` uses the existing brand tokens and fonts. Preview `/comparison.html` at desktop and phone widths and exercise its target selector.

The page states that OpenCASCADE is both a competitor and the STEP importer. Properties and expected membership use independent analytic definitions, with unavailable property formulas explicitly unmeasured. Pinned releases, the Pali Ijo source-distribution exception, rebuild modes, failed phases and name-check limits remain visible. See the [rerun guide](../.claude/skills/verify-catallax/scripts/head-to-head/README.md) for preparation and commands.
