# Kernel site

Static site for https://kernel.zigrawbohnbawerick.org/. The Pages workflow publishes this directory on site changes pushed to `main`; its custom domain is also recorded in `CNAME`. Domain DNS is managed separately.

## Preview

From the repository root:

```sh
nice -n 10 python3 -B -m http.server 8873 --bind 127.0.0.1 --directory site
```

Open http://127.0.0.1:8873/ in Chromium. No package install or build step is needed to serve the site.

## Evidence and video capture

`evidence.json` records a pinned kernel revision, counting methods and the scope of historical robustness figures. Counters are a snapshot, not live measurements. Rebuild the evidence when updating them:

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
