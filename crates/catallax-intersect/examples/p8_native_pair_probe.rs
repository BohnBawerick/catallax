#[path = "../tests/fixtures/p8_native_curve_patch.rs"]
mod fixture;
use catallax_intersect::{CurveSurface, intersect_curve_surface};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::{Curve, GeomStore, Plane, Surface};

fn main() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (source, cut, target) = fixture::patches(&mut dec);
        let section = source.sections(&cut, &mut dec).unwrap();
        println!("source-branches={}", section.branches().len());
        let seam = Rational::parse("1/59").unwrap();
        let profile_end = Rational::parse("1/4").unwrap();
        let controls = source.carrier().controls();
        let origin = controls[3].clone();
        let a = controls[7].minus(&origin);
        let b = controls[11].minus(&origin);
        let supporting = Plane {
            origin,
            normal: a.cross(&b),
            ref_dir: a,
        };
        let mut selected = 0;
        for (index, branch) in section.branches().iter().enumerate() {
            let ends = branch.parameter_endpoints();
            let on_seam = ends.iter().any(|uv| uv[1].rational_ref() == Some(&seam));
            let on_profile = ends
                .iter()
                .any(|uv| uv[0].rational_ref() == Some(&profile_end));
            println!(
                "branch={index} closed={} seam-end={on_seam} profile-end={on_profile}",
                branch.is_closed()
            );
            if !on_seam || !on_profile {
                continue;
            }
            selected += 1;
            for uv in &ends {
                if uv[1].rational_ref() != Some(&seam) {
                    continue;
                }
                let left = source.point_at(uv, &mut dec).unwrap();
                let right = target.point_at(uv, &mut dec).unwrap();
                assert!(dec.sign(&left.minus(&right).norm_squared()).is_coincident());
                println!("branch={index} exact-shared-boundary-incidence=true");
            }
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::PatchSection(branch.clone()));
            let plane = geom.add_surface(Surface::Plane(supporting.clone()));
            let before = dec.stats().work;
            match intersect_curve_surface(&mut geom, curve, plane, &mut dec) {
                Ok(answer) => match answer.value {
                    CurveSurface::Meets(hits) => println!(
                        "branch={index} exact-support-plane-hits={} work={}",
                        hits.len(),
                        dec.stats().work - before
                    ),
                    other => println!("branch={index} support-plane={other:?}"),
                },
                Err(error) => println!(
                    "branch={index} support-plane-refusal={error:?} work={}",
                    dec.stats().work - before
                ),
            }
        }
        assert_eq!(
            selected, 1,
            "fixture must identify the actual captured component uniquely"
        );
    });
    println!("work={work:?}");
    assert!(!work.exhausted);
}
