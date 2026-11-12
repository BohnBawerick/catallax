//! Writes a fit's deviation map as evidence: coloured PLY point files for both
//! directions, a PNG of three orthographic views per direction, and a
//! plain-text summary.

use std::fmt::Write as _;
use std::path::Path;

use catallax_field::reconstruction::scan::{FittedPart, ScanFit};

/// Blue (zero) through green and yellow to red (`top` and beyond).
fn colour(value: f64, top: f64) -> [u8; 3] {
    const STOPS: [[f64; 3]; 5] = [
        [48.0, 18.0, 160.0],
        [33.0, 145.0, 140.0],
        [94.0, 201.0, 98.0],
        [253.0, 231.0, 37.0],
        [220.0, 40.0, 30.0],
    ];
    let x = (value / top).clamp(0.0, 1.0) * 4.0;
    let i = (x.floor() as usize).min(3);
    let f = x - i as f64;
    std::array::from_fn(|c| (STOPS[i][c] + f * (STOPS[i + 1][c] - STOPS[i][c])).round() as u8)
}

fn ply(points: &[([f64; 3], f64)], top: f64) -> String {
    let mut s = format!(
        "ply\nformat ascii 1.0\ncomment deviation in mm; colour scale 0 to {top} mm\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nproperty float deviation\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n",
        points.len()
    );
    for (p, d) in points {
        let c = colour(*d, top);
        let _ = writeln!(
            s,
            "{:.4} {:.4} {:.4} {:.5} {} {} {}",
            p[0], p[1], p[2], d, c[0], c[1], c[2]
        );
    }
    s
}

/// An RGB image, written as an uncompressed PNG.
struct Image {
    width: usize,
    height: usize,
    rgb: Vec<u8>,
}

impl Image {
    fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            rgb: vec![255; width * height * 3],
        }
    }

    fn dot(&mut self, x: f64, y: f64, c: [u8; 3]) {
        for dy in 0..2 {
            for dx in 0..2 {
                let (px, py) = (x as usize + dx, y as usize + dy);
                if px < self.width && py < self.height {
                    let at = (py * self.width + px) * 3;
                    self.rgb[at..at + 3].copy_from_slice(&c);
                }
            }
        }
    }

    fn png(&self) -> Vec<u8> {
        fn crc(bytes: &[u8]) -> u32 {
            let mut c = !0u32;
            for &b in bytes {
                c ^= u32::from(b);
                for _ in 0..8 {
                    c = if c & 1 == 1 {
                        0xEDB8_8320 ^ (c >> 1)
                    } else {
                        c >> 1
                    };
                }
            }
            !c
        }
        fn chunk(out: &mut Vec<u8>, kind: &[u8], data: &[u8]) {
            out.extend_from_slice(&(data.len() as u32).to_be_bytes());
            let mut body = kind.to_vec();
            body.extend_from_slice(data);
            out.extend_from_slice(&body);
            out.extend_from_slice(&crc(&body).to_be_bytes());
        }
        let mut raw = Vec::with_capacity(self.height * (self.width * 3 + 1));
        for row in self.rgb.chunks(self.width * 3) {
            raw.push(0);
            raw.extend_from_slice(row);
        }
        // A zlib stream of stored (uncompressed) deflate blocks.
        let mut z = vec![0x78, 0x01];
        let blocks: Vec<&[u8]> = raw.chunks(65_535).collect();
        for (i, block) in blocks.iter().enumerate() {
            z.push(u8::from(i + 1 == blocks.len()));
            let n = block.len() as u16;
            z.extend_from_slice(&n.to_le_bytes());
            z.extend_from_slice(&(!n).to_le_bytes());
            z.extend_from_slice(block);
        }
        let (mut a, mut b) = (1u32, 0u32);
        for &byte in &raw {
            a = (a + u32::from(byte)) % 65_521;
            b = (b + a) % 65_521;
        }
        z.extend_from_slice(&((b << 16) | a).to_be_bytes());
        let mut header = Vec::new();
        header.extend_from_slice(&(self.width as u32).to_be_bytes());
        header.extend_from_slice(&(self.height as u32).to_be_bytes());
        header.extend_from_slice(&[8, 2, 0, 0, 0]);
        let mut out = vec![0x89, b'P', b'N', b'G', 0x0D, 0x0A, 0x1A, 0x0A];
        chunk(&mut out, b"IHDR", &header);
        chunk(&mut out, b"IDAT", &z);
        chunk(&mut out, b"IEND", &[]);
        out
    }
}

/// One row of three views: from +z (x right, y up), from -y (x right, z up)
/// and from +x (y right, z up), nearest drawn last, then the colour bar.
fn views(image: &mut Image, points: &[([f64; 3], f64)], top: f64, y0: f64) {
    const SCALE: f64 = 12.0;
    let views: [(usize, usize, usize, f64); 3] = [(0, 1, 2, 1.0), (0, 2, 1, -1.0), (1, 2, 0, 1.0)];
    let lo: [f64; 3] =
        std::array::from_fn(|k| points.iter().map(|p| p.0[k]).fold(f64::INFINITY, f64::min));
    let hi: [f64; 3] = std::array::from_fn(|k| {
        points
            .iter()
            .map(|p| p.0[k])
            .fold(f64::NEG_INFINITY, f64::max)
    });
    let mut x0 = 10.0;
    for (u, v, depth, toward) in views {
        let mut order: Vec<usize> = (0..points.len()).collect();
        order.sort_by(|&a, &b| {
            (toward * points[a].0[depth]).total_cmp(&(toward * points[b].0[depth]))
        });
        let height = (hi[v] - lo[v]) * SCALE;
        for i in order {
            let (p, d) = points[i];
            image.dot(
                x0 + (p[u] - lo[u]) * SCALE,
                y0 + height - (p[v] - lo[v]) * SCALE,
                colour(d, top),
            );
        }
        x0 += (hi[u] - lo[u]) * SCALE + 30.0;
    }
    for k in 0..=200 {
        for h in 0..20 {
            image.dot(
                x0 + k as f64,
                y0 + h as f64,
                colour(top * k as f64 / 200.0, top),
            );
        }
    }
}

/// Writes the map files into `dir`.
pub fn write(dir: &Path, fit: &ScanFit, part: &FittedPart) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let d = &part.deviation;
    let scan: Vec<([f64; 3], f64)> = fit
        .samples
        .iter()
        .copied()
        .zip(d.scan_map.iter().copied())
        .collect();
    // Every k-th part sample, so that the map stays a few megabytes.
    let stride = d.part_map.len().div_ceil(200_000).max(1);
    let part_samples: Vec<([f64; 3], f64)> = d
        .part_map
        .iter()
        .step_by(stride)
        .map(|s| (s.0, s.1))
        .collect();
    let (scan_top, part_top) = (0.6, 1.0);
    std::fs::write(dir.join("scan-to-part.ply"), ply(&scan, scan_top))?;
    std::fs::write(dir.join("part-to-scan.ply"), ply(&part_samples, part_top))?;
    let view_stride = part_samples.len().div_ceil(60_000).max(1);
    let part_view: Vec<_> = part_samples.iter().copied().step_by(view_stride).collect();
    // Two rows of three views: scan to part on the scan's samples (colour
    // scale 0 to 0.6 mm), then part to scan on the part's samples (0 to 1 mm).
    let mut image = Image::new(1650, 760);
    views(&mut image, &scan, scan_top, 10.0);
    views(&mut image, &part_view, part_top, 390.0);
    std::fs::write(dir.join("deviation-map.png"), image.png())?;
    let mut summary = String::new();
    let _ = writeln!(summary, "scan sha256 {}", fit.sha256);
    let _ = writeln!(summary, "noise deviation used {:.6} mm", fit.sigma_mm);
    for (i, s) in fit.surfaces.iter().enumerate() {
        let _ = writeln!(
            summary,
            "surface {i} {:?} samples {} rms {:.5} max {:.5} {:?}",
            s.role, s.samples, s.rms_mm, s.max_mm, s.surface
        );
    }
    for r in &fit.refusals {
        let _ = writeln!(summary, "refusal {}: {}", r.code, r.detail);
    }
    for (side, s) in [
        ("scan to part", &d.scan_to_part),
        ("part to scan", &d.part_to_scan),
    ] {
        let _ = writeln!(
            summary,
            "{side}: samples {} max {:.6} rms {:.6} p95 {:.6} mm",
            s.samples, s.max_mm, s.rms_mm, s.p95_mm
        );
    }
    let _ = writeln!(summary, "display mesh chord bound {:.6} mm", d.chord_mm);
    let _ = writeln!(
        summary,
        "deviation-map.png: top row scan to part on every scan sample, colour 0 to {scan_top} mm; bottom row part to scan on every {}th part sample, colour 0 to {part_top} mm; views from +z, -y and +x; blue is zero, red the top of the scale",
        stride * view_stride
    );
    summary.push_str(&part.recipe.to_text());
    std::fs::write(dir.join("summary.txt"), summary)
}
