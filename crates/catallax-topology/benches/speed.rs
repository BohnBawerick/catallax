#![allow(missing_docs)]
#[path = "../tests/common/mod.rs"]
mod common;
#[path = "../tests/common/prism.rs"]
mod prism;
use common::*;
use std::hint::black_box;
use std::time::Instant;

fn measure<T>(operation: &str, faces: usize, samples: u32, mut f: impl FnMut() -> T) {
    let start = Instant::now();
    for _ in 0..samples {
        black_box(f());
    }
    println!(
        "operation={operation} faces={faces} samples={samples} ns_per_call={}",
        start.elapsed().as_nanos() / u128::from(samples)
    );
}
fn main() {
    if std::env::var_os("CATALLAX_TOPOLOGY_RUN_BENCH").as_deref() != Some(std::ffi::OsStr::new("1"))
    {
        println!(
            "Disabled against interface bodies. Set CATALLAX_TOPOLOGY_RUN_BENCH=1 to measure piece 3."
        );
        return;
    }
    let pose = Rng(0x3beac4).pose();
    for (i, shape) in primitives().into_iter().enumerate() {
        for (name, pose) in [("identity", Pose::identity()), ("posed", pose.clone())] {
            // Each sample owns a fresh model. Measure construction only, excluding
            // model allocation and destruction; report the aggregate sample time.
            let mut total = std::time::Duration::ZERO;
            for _ in 0..100 {
                let mut m = Model::new();
                let start = Instant::now();
                black_box(shape.build(&mut m, black_box(&pose)).unwrap());
                total += start.elapsed();
                black_box(&m);
            }
            println!(
                "operation=build.{}.variant{i}.{name} samples=100 ns_per_call={}",
                shape.name(),
                total.as_nanos() / 100
            );
        }
    }
    for n in [3, 10, 100, 1000, 10_000] {
        let descriptor = prism::prism(n);
        let mut m = Model::new();
        let b = assemble(&mut m, &descriptor);
        // Correctness can fail; elapsed time never does.
        sound(&m, b);
        let counts = m.counts(b);
        assert_eq!(
            (counts.faces, counts.edges, counts.vertices),
            (n + 2, 3 * n, 2 * n)
        );
        let samples = if n >= 1000 { 3 } else { 20 };
        measure("check", n + 2, samples, || m.check(black_box(b)));
        measure("counts", n + 2, samples, || m.counts(black_box(b)));
        measure("faces_of_body", n + 2, samples, || {
            m.faces_of_body(black_box(b))
        });
        measure("edges_of_body", n + 2, samples, || {
            m.edges_of_body(black_box(b))
        });
        measure("vertices_of_body", n + 2, samples, || {
            m.vertices_of_body(black_box(b))
        });
        if n == 10_000 {
            let edges = m.edges_of_body(b);
            let coedges: Vec<_> = edges.iter().map(|&e| m.coedges_of_edge(e)[0]).collect();
            let mut index = 0;
            measure("coedges_of_edge", n + 2, 30_000, || {
                let e = edges[index % edges.len()];
                index += 1;
                m.coedges_of_edge(black_box(e))
            });
            measure("owner_of_coedge", n + 2, 30_000, || {
                let c = coedges[index % coedges.len()];
                index += 1;
                m.owner_of_coedge(black_box(c))
            });
        }
    }
}
