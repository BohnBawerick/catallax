#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use common::*;
use std::time::Instant;

fn measured_roots(name: &str, polynomial: &UPoly) -> Vec<RealRoot> {
    let mut decider = Decider::default();
    let start = Instant::now();
    let outcome = decider.real_roots(polynomial);
    eprintln!(
        "{name}: {:?}, work={}, elapsed={:.6}s",
        if matches!(outcome, Outcome::Proven(_)) {
            "proven"
        } else {
            "refused"
        },
        decider.stats().work,
        start.elapsed().as_secs_f64()
    );
    stats_sum(decider.stats());
    match outcome {
        Outcome::Proven(roots) => roots,
        other => panic!("expected complete roots, got {other:?}"),
    }
}

fn check_four_simple_roots(polynomial: &UPoly, roots: &[RealRoot]) {
    assert_eq!(roots.len(), 4);
    let coefficients: Vec<_> = polynomial.coefficients().iter().map(oq).collect();
    for root in roots {
        assert_eq!(root.multiplicity, 1);
        assert!(root.value.as_rational().is_none());
        check_irrational_certificate(&root.value);
        let (lo, hi) = root.value.isolating_interval();
        assert!(oracle_eval(&coefficients, &oq(&lo)) * oracle_eval(&coefficients, &oq(&hi)) < 0);
    }
    for pair in roots.windows(2) {
        assert!(pair[0].value.isolating_interval().1 < pair[1].value.isolating_interval().0);
    }
}

#[test]
fn case40_scaled_large_leading_quartic() {
    let n = r("1e30") * Rational::from_integer(Integer::one().mul_pow2(607) - Integer::one());
    let n2 = &n * &n;
    let polynomial = UPoly::new(vec![r("6"), r("0"), -(&r("5") * &n2), r("0"), &n2 * &n2]);
    let roots = measured_roots("case40", &polynomial);
    check_four_simple_roots(&polynomial, &roots);
    for (index, root) in roots.iter().enumerate() {
        let (lo, hi) = root.value.isolating_interval();
        let (lo, hi) = (oq(&(&lo * &n)), oq(&(&hi * &n)));
        let radicand = if index == 0 || index == 3 { 3 } else { 2 };
        if index < 2 {
            assert!(hi < 0 && lo.clone().square() > radicand && hi.square() < radicand);
        } else {
            assert!(lo > 0 && lo.square() < radicand && hi.square() > radicand);
        }
    }
}

#[test]
fn case40_tiny_bite_quartic_zero() {
    let e = r("1e-30");
    let a = &r("4") * &e - &e * &e;
    let b = &r("8") * &e - &e * &e - r("12");
    let polynomial = UPoly::new(vec![&a * &b, r("0"), -(&a + &b), r("0"), r("1")]);
    let roots = measured_roots("case40 tiny bite", &polynomial);
    assert_eq!(roots.len(), 2);
    for (index, root) in roots.into_iter().enumerate() {
        assert_eq!(root.multiplicity, 1);
        check_irrational_certificate(&root.value);
        let (lo, hi) = root.value.isolating_interval();
        assert_eq!(
            lo.sign(),
            if index == 0 {
                Sign::Negative
            } else {
                Sign::Positive
            }
        );
        assert!((oq(&lo).square() - oq(&a)) * (oq(&hi).square() - oq(&a)) < 0);
        let degree = root.value.defining_polynomial().degree();
        let x = Real::from_algebraic(root.value);
        let inner = (real("4") - x.square()).sqrt() - real("3") + Real::from_rational(e.clone());
        let expression = real("1") - inner.square();
        let mut decider = Decider::default();
        let start = Instant::now();
        let outcome = decider.sign(&expression);
        eprintln!(
            "case40 tiny bite sign: {outcome:?}, work={}, elapsed={:?}",
            decider.stats().work,
            start.elapsed()
        );
        assert!(matches!(outcome, Outcome::ProvenCoincident(_)));
        assert_eq!(degree, Some(2));
        stats_sum(decider.stats());
    }
}

#[test]
fn non_biquadratic_quartic_quadratic_fields() {
    for (left, right) in [
        (poly(&[-1, 1, 1]), poly(&[-1, 2, 1])),
        (poly(&[-1, 1, 1]), poly(&[-3, 1, 1])),
        (poly(&[-1, 1, 1]), poly(&[-1, -1, 1])),
        (
            UPoly::new(vec![r("-1"), r("1/3"), r("1")]),
            UPoly::new(vec![r("-1"), r("2/3"), r("1")]),
        ),
    ] {
        let polynomial = &left * &right;
        let roots = measured_roots("non-biquadratic split", &polynomial);
        check_four_simple_roots(&polynomial, &roots);
        let left = UPoly::new(
            left.primitive_integer()
                .into_iter()
                .map(Rational::from_integer)
                .collect(),
        );
        let right = UPoly::new(
            right
                .primitive_integer()
                .into_iter()
                .map(Rational::from_integer)
                .collect(),
        );
        let mut counts = [0, 0];
        for root in roots {
            let defining = root.value.defining_polynomial();
            let index = if defining == left {
                0
            } else {
                assert_eq!(defining, right);
                1
            };
            counts[index] += 1;
            let x = Real::from_algebraic(root.value);
            proven_sign(&mut Decider::default(), &defining.eval_real(&x), 0);
        }
        assert_eq!(counts, [2, 2]);
    }
}

#[test]
fn irreducible_quartic_with_rational_resolvent_root_stays_quartic() {
    // Its real roots are ±sqrt(2) ±sqrt(3), not quadratic over Q.
    let polynomial = poly(&[1, 0, -10, 0, 1]);
    let roots = measured_roots("irreducible quartic", &polynomial);
    check_four_simple_roots(&polynomial, &roots);
    for root in roots {
        assert_eq!(root.value.defining_polynomial(), polynomial);
    }
}

#[test]
fn quartic_split_discards_nonreal_factor_and_preserves_multiplicity() {
    let real_factor = poly(&[-1, 2, 1]);
    let nonreal_factor = poly(&[1, 1, 1]);
    let quartic = &real_factor * &nonreal_factor;
    let roots = measured_roots("mixed real quartic", &(&quartic * &quartic));
    assert_eq!(roots.len(), 2);
    for root in roots {
        assert_eq!(root.multiplicity, 2);
        assert_eq!(root.value.defining_polynomial(), real_factor);
        check_irrational_certificate(&root.value);
    }
}

#[test]
fn case23_scaled_probe_plane_quartic() {
    // Increasing-degree coefficients from the scale 10^-30 probe-plane intersection.
    let polynomial = UPoly::new([
        "438465779342886310724550832791405980021223566338247198645426210144441769063860108512211564675503270571839796871238097594388303971902816208828433697392425013833746186412006149474416616133880200718854143417531948902525107665245758763583734405776143900239042704581240857785911081580178802089234054108161478574249/27646996131816418365991313272262769272359906074912746952991702567090043810554730210014742819375866606396450974737801670028744456421030725650144261652062214074618745928582020696675555512084045549585575088160263742934091646783978631991571317877548602406832709176386694599840151668816518615185584006010925744128000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        "-1094911715548042674306492485585595591491810164521228595980707408923830539204465712847496908613340912502562713027245435264712366752660707113957918727612432520725538557634973822785408904503368733590911835937735349502293132858996906932192340895235908105678906855688182407102464/55376024831528252148370813628103832758473042995237322141286914892068108198251181447000587646828339794668723003600699192831367065561522432748992094806034590758341774166155520958482631418587506045182439572600267219873080368883111107596169145564433321056748141083845866676566306912338570329399304955187010332887858727693802481963081252303027213201858103275299072265625",
        "-18171085501262581693049804937597187610473271579636478481379026476037961385007297493347961935928805383144016897754853019747946499969674564988345680379142145271261503576316364108503352362019123309873657539184588894725808031153338311160016887531463267150249454134651166638708945066913738941385641908296677613767/2281417161268054054615494112799027347182043030595827263210741080975692091788940139400630633044199929922363386098968985466239166569899508122497255966405524501274691436489434325067465274190529149355450288036662388943291742337154486721960719102200055569704456958403003607115715640639644358382013523933518774784000000000000000000000000000000000000000000000000000000000000",
        "62054782366352144811118993726372680086287112794140191598566720321162068886886251681700024467744454985047384756959023086424426350212561240064/12519733959184334248445786156668192021277548149710899608887788732455997867183464638077467062181819868157742847986035246515244911693052927112759848264506645500659942626953125",
        "1",
    ].map(r).to_vec());
    let roots = measured_roots("case23", &polynomial);
    check_four_simple_roots(&polynomial, &roots);
    // The quartic has degree four; one certified root in each disjoint interval
    // proves completeness independently of the kernel's root counter.
    for (index, root) in roots.iter().enumerate() {
        let (lo, hi) = root.value.isolating_interval();
        if index < 2 {
            assert!(lo > r("-3e-30") && hi < r("-1e-30"));
        } else {
            assert!(lo > r("1e-30") && hi < r("3e-30"));
        }
    }
}

#[test]
fn unenumerated_rational_root_remains_exact() {
    let denominator = Integer::one().mul_pow2(607) - Integer::one();
    let expected = Rational::new(Integer::one(), denominator.clone()).unwrap();
    let linear = UPoly::new(vec![r("-1"), Rational::from_integer(denominator)]);
    let polynomial = &linear * &poly(&[-2, 0, 1]);
    let roots = roots(&mut Decider::default(), &polynomial);
    assert_eq!(roots.len(), 3);
    assert_eq!(roots[1].value.as_rational(), Some(&expected));
    for root in [&roots[0], &roots[2]] {
        assert_eq!(root.multiplicity, 1);
        check_irrational_certificate(&root.value);
        assert_eq!(root.value.defining_polynomial(), poly(&[-2, 0, 1]));
    }
}

#[test]
fn cross_factor_separation_and_scale_relative_refinement() {
    // Both factors initially isolate their positive roots in overlapping
    // intervals. The second factor has multiplicity two in the input.
    let a = poly(&[-2, 0, 1]);
    let b = poly(&[-3, 0, 1]);
    let polynomial = &a * &(&b * &b);
    for policy in [
        BudgetPolicy::default(),
        BudgetPolicy {
            max_precision_bits: 0,
            ..BudgetPolicy::default()
        },
    ] {
        let roots = roots(&mut Decider::new(policy), &polynomial);
        assert_eq!(roots.len(), 4);
        for (index, root) in roots.iter().enumerate() {
            assert_eq!(
                root.multiplicity,
                if index == 0 || index == 3 { 2 } else { 1 }
            );
            check_irrational_certificate(&root.value);
        }
        for pair in roots.windows(2) {
            assert!(pair[0].value.isolating_interval().1 < pair[1].value.isolating_interval().0);
        }
    }
    for scale in [r("1e-30"), r("1"), r("1e30")] {
        let polynomial = UPoly::new(vec![-(&r("2") * &(&scale * &scale)), r("0"), r("1")]);
        let roots = roots(&mut Decider::default(), &polynomial);
        assert_eq!(roots.len(), 2);
        let target = &scale * &r("1e-80");
        for root in roots {
            let refined = root.value.refined(&target);
            check_irrational_certificate(&refined);
            let (lo, hi) = refined.isolating_interval();
            assert!(&hi - &lo <= target);
            let (lo, hi) = (oq(&(&lo / &scale)), oq(&(&hi / &scale)));
            assert!((lo.square() - 2) * (hi.square() - 2) < 0);
        }
    }
}

#[test]
fn rational_roots_fit_a_small_work_budget_without_precision_escalation() {
    let polynomial = poly(&[2, -9, 9]);
    let mut decider = Decider::new(BudgetPolicy {
        max_work: 10_000,
        ..BudgetPolicy::default()
    });
    let result = roots(&mut decider, &polynomial);
    assert_eq!(
        result
            .iter()
            .map(|root| (root.value.as_rational().cloned(), root.multiplicity))
            .collect::<Vec<_>>(),
        vec![(Some(r("1/3")), 1), (Some(r("2/3")), 1)]
    );
}

#[test]
fn rational_translation_of_irreducible_quartic_keeps_all_real_roots() {
    let denominator = Integer::one().mul_pow2(127) - Integer::one();
    let shift = r("22496") + Rational::new(Integer::one(), denominator).unwrap();
    let original = poly(&[1, 1, -10, 0, 1]);
    let polynomial = original.compose(&UPoly::new(vec![-&shift, r("1")]));
    let result = measured_roots("translated irreducible quartic", &polynomial);
    check_four_simple_roots(&polynomial, &result);
    for root in &result {
        assert_eq!(root.value.defining_polynomial().degree(), Some(4));
    }
}

#[test]
fn rotated_probe_plane_quartic_fits_default_budget() {
    let polynomial = UPoly::new([
        "744084836416383378744392617980469422776554284710885694386263293101584717734754069663073393987810109044527605505797867775036065994667233439809657946280670003069787163632158530300461016576232983087330863370055338801593835501217536304536222011483836028813089797234428275922961/2905497773183080566409804419226053137231986284713011403554479294236198037675029065467544636344402336202278975440815198860424841812221516193340622254569449768595205252795395008303071576176900905467686785586827103351448327726796353816232284751344435200000000",
        "-148619203228761990890098077104870402769054181018552855229035385984756424211226141166331858960739864462366368701933203934515118503006957295519024083199478934135675511682296313462681164235828481648832266011041134820877488427619764394316214580271/3263724614418660749721991502821041696201203618434426049652506705651668123248556673724728942127741410731962226907549716596820254897205937158141875974633615540377832286684177173485729665905269853093305720639107749326390789406720000",
        "445265105950935789564328753109520507517855831691251515754279718031788896110190814048270818442840025545248464552290545736920825832737116714347595170347412012725189412622012782983407040077656613229234743113085274381351/146644729272570545953553899689245679559850566569510625556609986147828696465462960792630151510814151964682082657823401739016751783219432407010018206060363934636420528620248287884648144622933562697113927680000",
        "-32912442371175564680622141741796550370886883304751833241629600964979576759196115753036037256413724967969856314300910268704737349455759699691440355568140364781550503726817/365762892784659978683698046687335340147162968552198136687841195042720380640054622367960362407222190808240595386012896755413743775164619031841046070607227453683170880",
        "1",
    ].map(r).to_vec());
    let result = measured_roots("rotated probe quartic", &polynomial);
    check_four_simple_roots(&polynomial, &result);
}

#[test]
fn translated_clustered_quadratic_fits_default_budget() {
    let denominator = Integer::one().mul_pow2(767) - Integer::one();
    let center = r("434139/7") + Rational::new(Integer::one(), denominator).unwrap();
    let scale = r("1e-20");
    let constant = &(&center * &center) - &(r("2") * (&scale * &scale));
    let polynomial = UPoly::new(vec![constant, -center.mul_pow2(1), r("1")]);
    let result = measured_roots("translated clustered quadratic", &polynomial);
    assert_eq!(result.len(), 2);
    for (index, root) in result.into_iter().enumerate() {
        assert_eq!(root.multiplicity, 1);
        assert_eq!(root.value.defining_polynomial().degree(), Some(2));
        check_irrational_certificate(&root.value);
        let (lo, hi) = root.value.isolating_interval();
        let a = oq(&(&(&lo - &center) / &scale));
        let b = oq(&(&(&hi - &center) / &scale));
        if index == 0 {
            assert!(b < 0 && a.square() > 2 && b.square() < 2);
        } else {
            assert!(a > 0 && a.square() < 2 && b.square() > 2);
        }
    }
}

#[test]
fn translated_clustered_cubic_fits_default_budget() {
    let denominator = Integer::one().mul_pow2(767) - Integer::one();
    let center = r("434139/7") + Rational::new(Integer::one(), denominator).unwrap();
    let scale = r("1e-20");
    let square = &center * &center;
    let polynomial = UPoly::new(vec![
        -(&square * &center) - r("2") * (&(&scale * &scale) * &scale),
        r("3") * square,
        -(&r("3") * &center),
        r("1"),
    ]);
    let result = measured_roots("translated clustered cubic", &polynomial);
    assert_eq!(result.len(), 1);
    let root = &result[0];
    assert_eq!(root.multiplicity, 1);
    assert_eq!(root.value.defining_polynomial().degree(), Some(3));
    check_irrational_certificate(&root.value);
    let (lo, hi) = root.value.isolating_interval();
    let a = oq(&(&(&lo - &center) / &scale));
    let b = oq(&(&(&hi - &center) / &scale));
    assert!(a.clone().square() * a < 2);
    assert!(b.clone().square() * b > 2);
}

#[test]
fn translated_clustered_rational_cubic_stays_rational() {
    let denominator = Integer::one().mul_pow2(767) - Integer::one();
    let center = r("434139/7") + Rational::new(Integer::one(), denominator).unwrap();
    let scale = r("1e-20");
    let square = &center * &center;
    let polynomial = UPoly::new(vec![
        -(&square * &center) - (&(&scale * &scale) * &scale),
        r("3") * square,
        -(&r("3") * &center),
        r("1"),
    ]);
    let result = measured_roots("translated clustered rational cubic", &polynomial);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].multiplicity, 1);
    assert_eq!(result[0].value.as_rational(), Some(&(center + scale)));
}

#[test]
fn repeated_root_queries_keep_distinct_certificates() {
    let polynomials: Vec<_> = (1..=257)
        .map(|t| UPoly::new(vec![-r(&(t * t + 1).to_string()), r("0"), r("1")]))
        .collect();
    let mut decider = Decider::default();
    // Distinct t²+1 lie strictly between consecutive integer squares. Query
    // more polynomials than memo slots, then revisit them in reverse order.
    for polynomial in polynomials.iter().chain(polynomials.iter().rev()) {
        let expected = oq(&(-&polynomial.coefficients()[0]));
        let result = roots(&mut decider, polynomial);
        assert_eq!(result.len(), 2);
        for (index, root) in result.into_iter().enumerate() {
            assert_eq!(root.multiplicity, 1);
            assert!(root.value.as_rational().is_none());
            let value = root.value.refined(&r("1/4"));
            check_irrational_certificate(&value);
            let (lo, hi) = value.isolating_interval();
            let (lo, hi) = (oq(&lo), oq(&hi));
            if index == 0 {
                assert!(hi < 0 && lo.square() > expected && hi.square() < expected);
            } else {
                assert!(lo > 0 && lo.square() < expected && hi.square() > expected);
            }
        }
    }
    stats_sum(decider.stats());
}

#[test]
fn repeated_root_query_preserves_work_and_obeys_expired_deadline() {
    let polynomial = UPoly::new(vec![r("-2"), r("0"), r("1")]);
    let mut decider = Decider::default();
    roots(&mut decider, &polynomial);
    let original_work = decider.stats().work;
    roots(&mut decider, &polynomial);
    assert_eq!(decider.stats().work, original_work * 2);
    decider.set_deadline(Some(Instant::now()));
    assert!(matches!(
        decider.real_roots(&polynomial),
        Outcome::Refused(error) if error.code == "numbers.deadline"
    ));
}

#[test]
fn certified_narrow_root_encloses_without_rechecking_its_polynomial() {
    let scale = Rational::from_integer(Integer::one().mul_pow2(607) - Integer::one());
    let polynomial = UPoly::new(vec![r("-2"), r("0"), &scale * &scale]);
    let root = roots(&mut Decider::default(), &polynomial).pop().unwrap();
    let value = Real::from_algebraic(root.value);
    let mut decider = Decider::new(BudgetPolicy {
        max_precision_bits: 512,
        max_work: 10_000,
    });
    let width = r("1e-220");
    let enclosure = decider
        .enclose(&value, &width)
        .expect("an already-isolated root needs only its certified endpoints");
    assert!(&enclosure.hi() - &enclosure.lo() <= width);
    let lo = oq(&(&enclosure.lo() * &scale));
    let hi = oq(&(&enclosure.hi() * &scale));
    assert!(lo > 0 && lo.square() < 2 && hi.square() > 2);
    stats_sum(decider.stats());
}

#[test]
fn shared_root_enclosures_preserve_limited_budgets_and_deadlines() {
    let polynomial = UPoly::new(vec![r("-2"), r("0"), r("0"), r("1")]);
    let root = roots(&mut Decider::default(), &polynomial).pop().unwrap();
    let value = Real::from_algebraic(root.value);
    let policy = BudgetPolicy {
        max_precision_bits: 512,
        max_work: 1_000,
    };
    let width = r("1e-30");
    let mut cold = Decider::new(policy);
    let cold_ball = cold.enclose(&value, &width).unwrap();
    let cold_stats = cold.take_stats();
    let certificate = Decider::default().enclose_at(&value, 64).unwrap();
    assert!(&certificate.hi() - &certificate.lo() < &cold_ball.hi() - &cold_ball.lo());
    let lo = oq(&certificate.lo());
    let hi = oq(&certificate.hi());
    assert!(lo > 0 && lo.clone().square() * lo < 2 && hi.clone().square() * hi > 2);
    let mut warm = Decider::new(policy);
    let warm_ball = warm.enclose(&value, &width).unwrap();
    assert_eq!(warm_ball.mid(), cold_ball.mid());
    assert_eq!(warm_ball.rad(), cold_ball.rad());
    assert_eq!(warm.take_stats(), cold_stats);
    warm.set_deadline(Some(Instant::now()));
    assert!(matches!(warm.enclose_at(&value, 64), Err(error) if error.code == "numbers.deadline"));
}

#[test]
fn bisection_prefixes_preserve_rounding_work_and_budget_boundaries() {
    let polynomial = UPoly::new(vec![r("-2"), r("0"), r("0"), r("1")]);
    let fresh = || {
        Real::from_algebraic(
            roots(&mut Decider::default(), &polynomial)
                .pop()
                .unwrap()
                .value,
        )
    };
    let warmed = fresh();
    let width = r("1e-90");
    for bits in [32, 64, 96, 128, 256, 512, 1024, 2048, 4096, 64] {
        Decider::default().enclose_at(&warmed, bits).unwrap();
        for max_work in [0, 100, 1_000, 10_000, 1_000_000] {
            let policy = BudgetPolicy {
                max_precision_bits: 512,
                max_work,
            };
            let mut cold = Decider::new(policy);
            let mut warm = Decider::new(policy);
            let before = cold.enclose(&fresh(), &width);
            let after = warm.enclose(&warmed, &width);
            match (before, after) {
                (Ok(a), Ok(b)) => {
                    assert_eq!(a.mid(), b.mid());
                    assert_eq!(a.rad(), b.rad());
                    assert_eq!(a.precision(), b.precision());
                    let lo = oq(&b.lo());
                    let hi = oq(&b.hi());
                    assert!(lo.clone().square() * lo < 2);
                    assert!(hi.clone().square() * hi > 2);
                }
                (Err(a), Err(b)) => assert_eq!(a.code, b.code),
                outcomes => panic!("cache history changed the outcome: {outcomes:?}"),
            }
            assert_eq!(cold.take_stats(), warm.take_stats());
        }
    }
}
