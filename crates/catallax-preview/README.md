# Catallax GPU preview mesher

`Mesher::new` takes the renderer's wgpu 30 device and queue. A call to `mesh`
triangulates each face's trimmed UV domain, uploads one binary32 copy of each
shared edge sample, and submits one compute pass. Use `PreviewVertex::layout()`
and each successful face's `Uint32` index range to draw the returned buffers on
that device. Positions are relative to the binary64 body origin. Production
code never maps either buffer.

Planes and cylinders evaluate their carriers in WGSL. Other carriers use
sampled grids for now; `FaceDraw::sampled_grid` distinguishes these draws. The
grid density is chosen per face. It does not split shared trim chords or promise
a geometric bound. `preview_or_cpu` calls the existing CPU preview on device
failure or face refusal. In a mixed result, use the CPU mesh only for the refused
face IDs. The checked CPU mesh and its proof mark remain a separate result.

The P5 flanged-housing harness reads the same compute output on each platform.
It checks face coverage, bit-identical shared chords and distance from the
checked surface. On Linux the software Vulkan run is a correctness check, not
hardware acceptance. For native acceptance on Windows (DX12) and macOS (Metal):

```sh
cargo test --locked --release -p catallax-preview --test gpu_mesher p5_native_hardware -- --ignored --nocapture
```

The same test file exposes `p5_browser_webgpu` to `wasm-bindgen-test`. Browser
acceptance must verify the browser's hardware adapter separately; a native pass
or a software adapter does not count.
