use super::*;

#[derive(Clone, Debug)]
pub struct Case {
    pub id: &'static str,
    pub number: usize,
    pub a: Primitive,
    pub b: Primitive,
}
pub const IDS: &[&str] = &[
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "09-overlap",
    "10",
    "11",
    "12",
    "13",
    "13-bite",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "19-gap",
    "19-bite",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "27-tube",
    "28",
    "29",
    "29-bite",
    "30",
    "30-irrational",
    "31",
    "32",
    "33",
    "34",
    "35",
    "36",
    "37",
    "37-gap",
    "37-bite",
    "38",
    "39",
    "40",
];
pub fn case(id: &'static str) -> Case {
    let number = id[..2].parse().unwrap();
    let box10 = || block_at([10; 3], [0; 3]);
    let (a, mut b) = match number {
        1 => (block_at([2; 3], [0; 3]), block_at([2; 3], [1; 3])),
        2 => (block_at([1; 3], [0; 3]), block_at([1; 3], [1, 0, 0])),
        3 => (block_at([2; 3], [0; 3]), block_at([2; 3], [2, 1, 0])),
        4 => (block_at([1; 3], [0; 3]), block_at([1; 3], [1, 1, 0])),
        5 => (block_at([1; 3], [0; 3]), block_at([1; 3], [1; 3])),
        6 => (block_at([4; 3], [0; 3]), block_at([2, 2, 6], [1, 1, -1])),
        7 => (block_at([2; 3], [0; 3]), block_at([2; 3], [1, 0, 0])),
        8 => (block_at([2; 3], [0; 3]), block_at([2; 3], [0; 3])),
        9 => (block_at([1; 3], [0; 3]), block_at([1; 3], [1, 0, 0])),
        10 => (box10(), cyl_at(3, 12, [5, 5, -1], false)),
        11 => (box10(), cyl_at(3, 5, [5, 5, 6], false)),
        12 | 13 => (box10(), cyl_at(2, 20, [-5, 5, 12], true)),
        14 => (box10(), cyl_at(5, 12, [5, 5, -1], false)),
        15 => (box10(), cyl_at(3, 12, [0, 5, -1], false)),
        16 => (box10(), cyl_at(2, 4, [5, 5, 3], false)),
        17 => (box10(), ball_at(3, [10; 3])),
        18 => (box10(), ball_at(3, [5, 5, 10])),
        19 => (box10(), ball_at(3, [5, 5, 13])),
        20 => (box10(), ball_at(3, [5, 5, 7])),
        21 => (box10(), ball_at(5, [5; 3])),
        22 => (box10(), ball_at(3, [10, 5, 5])),
        23 => (
            cyl_at(2, 10, [0, 0, -5], false),
            cyl_at(1, 10, [-5, 0, 0], true),
        ),
        24 => (
            cyl_at(1, 4, [0, 0, -2], false),
            cyl_at(1, 4, [-2, 0, 0], true),
        ),
        25 => (
            cyl_at(5, 10, [0; 3], false),
            cyl_at(5, 10, [8, 0, 0], false),
        ),
        26 => (
            cyl_at(5, 10, [0; 3], false),
            cyl_at(5, 10, [10, 0, 0], false),
        ),
        27 if id == "27-tube" => (
            cyl_at(3, 10, [0; 3], false),
            cyl_at(1, 12, [0, 0, -1], false),
        ),
        27 => (cyl_at(2, 5, [0; 3], false), cyl_at(2, 5, [0, 0, 5], false)),
        28 => (
            cyl_at(2, 10, [0, 0, -5], false),
            cyl_at(1, 10, [-5, 1, 0], true),
        ),
        29 => (
            cyl_at(2, 10, [0, 0, -5], false),
            cyl_at(1, 10, [-5, 3, 0], true),
        ),
        30 => (
            cyl_at(3, 20, [0, 0, -10], false),
            ball_at(if id == "30-irrational" { 6 } else { 5 }, [0; 3]),
        ),
        31 => (cyl_at(5, 20, [0, 0, -10], false), ball_at(5, [0; 3])),
        32 => (cyl_at(1, 6, [0, 0, -3], false), ball_at(2, [1, 0, 0])),
        33 => (cyl_at(1, 6, [0, 0, -3], false), ball_at(2, [2, 0, 0])),
        34 => (cyl_at(1, 6, [0, 0, -3], false), ball_at(2, [3, 0, 0])),
        35 => (cyl_at(1, 6, [0, 0, -3], false), ball_at(1, [1, 0, 0])),
        36 => (ball_at(5, [0; 3]), ball_at(5, [6, 0, 0])),
        37 => (ball_at(5, [0; 3]), ball_at(5, [10, 0, 0])),
        38 => (ball_at(5, [0; 3]), ball_at(2, [3, 0, 0])),
        39 => (ball_at(5, [0; 3]), ball_at(3, [0; 3])),
        40 => (ball_at(5, [0; 3]), ball_at(5, [0; 3])),
        _ => panic!("unknown case"),
    };
    let delta = catallax_numbers::Real::from(q("1e-30"));
    match id {
        "09" => b.pose.translation.x = &b.pose.translation.x + &delta,
        "09-overlap" | "37-bite" => b.pose.translation.x = &b.pose.translation.x - &delta,
        "13" | "19-gap" => b.pose.translation.z = &b.pose.translation.z + &delta,
        "13-bite" | "19-bite" => b.pose.translation.z = &b.pose.translation.z - &delta,
        "29-bite" => b.pose.translation.y = &b.pose.translation.y - &delta,
        "37-gap" => b.pose.translation.x = &b.pose.translation.x + &delta,
        _ => {}
    }
    Case { id, number, a, b }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Counts {
    pub lumps: usize,
    pub shells: usize,
    pub faces: usize,
    pub edges: usize,
    pub rings: usize,
    pub vertices: usize,
    pub surfaces: [usize; 3],
    pub curves: [usize; 4],
}
impl Counts {
    pub fn empty() -> Self {
        Self {
            lumps: 0,
            shells: 0,
            faces: 0,
            edges: 0,
            rings: 0,
            vertices: 0,
            surfaces: [0; 3],
            curves: [0; 4],
        }
    }
}
#[derive(Clone, Debug)]
pub enum Expected {
    Counts(Counts),
    Refused(&'static str, Stage),
}
pub fn operation(code: &str) -> (BoolOp, bool) {
    match code {
        "U" => (BoolOp::Unite, false),
        "I" => (BoolOp::Intersect, false),
        "D" => (BoolOp::Subtract, false),
        "R" => (BoolOp::Subtract, true),
        _ => panic!("operation"),
    }
}

// Columns: case operation lumps shells F E R V planes cylinders spheres
// lines circles ellipses branches. Values are the specification's counts.
// Unprinted kind counts follow the boundary decomposition: caps are planes,
// cylinder/sphere cuts are their retained side patches, and each closed-form
// section has the kind of its carrier, not CurveKind::Intersection.
// Zero-thickness intersections are regularised to Empty.
pub const TABLE: &str = "
01 U 1 1 12 30 0 20 12 0 0 30 0 0 0
01 I 1 1 6 12 0 8 6 0 0 12 0 0 0
01 D 1 1 9 21 0 14 9 0 0 21 0 0 0
01 R 1 1 9 21 0 14 9 0 0 21 0 0 0
02 U 1 1 6 12 0 8 6 0 0 12 0 0 0
02 I E
02 D 1 1 6 12 0 8 6 0 0 12 0 0 0
03 U 1 1 10 24 0 16 10 0 0 24 0 0 0
03 I E
03 D 1 1 6 12 0 8 6 0 0 12 0 0 0
04 U 2 2 12 24 0 16 12 0 0 24 0 0 0
04 I E
04 D 1 1 6 12 0 8 6 0 0 12 0 0 0
05 U 2 2 12 24 0 16 12 0 0 24 0 0 0
05 I E
06 U 1 1 16 36 0 24 16 0 0 36 0 0 0
06 I 1 1 6 12 0 8 6 0 0 12 0 0 0
06 D 1 1 10 24 0 16 10 0 0 24 0 0 0
07 U 1 1 6 12 0 8 6 0 0 12 0 0 0
07 I 1 1 6 12 0 8 6 0 0 12 0 0 0
07 D 1 1 6 12 0 8 6 0 0 12 0 0 0
08 U 1 1 6 12 0 8 6 0 0 12 0 0 0
08 I 1 1 6 12 0 8 6 0 0 12 0 0 0
08 D E
09 U 2 2 12 24 0 16 12 0 0 24 0 0 0
09 I E
09 D 1 1 6 12 0 8 6 0 0 12 0 0 0
09-overlap U 1 1 6 12 0 8 6 0 0 12 0 0 0
09-overlap I 1 1 6 12 0 8 6 0 0 12 0 0 0
09-overlap D 1 1 6 12 0 8 6 0 0 12 0 0 0
10 U 1 1 10 16 4 8 8 2 0 12 4 0 0
10 I 1 1 3 2 2 0 2 1 0 0 2 0 0
10 D 1 1 7 14 2 8 6 1 0 12 2 0 0
11 D 1 1 8 14 2 8 7 1 0 12 2 0 0
12 U 2 2 9 14 2 8 8 1 0 12 2 0 0
12 I E
12 D 1 1 6 12 0 8 6 0 0 12 0 0 0
13 U 2 2 9 14 2 8 8 1 0 12 2 0 0
13 I E
13 D 1 1 6 12 0 8 6 0 0 12 0 0 0
13-bite U 1 1 10 20 2 12 9 1 0 16 4 0 0
13-bite I 1 1 4 6 0 4 3 1 0 4 2 0 0
13-bite D 1 1 8 18 0 12 7 1 0 16 2 0 0
14 U 1 1 16 30 2 16 14 2 0 20 10 0 0
14 I 1 1 3 2 2 0 2 1 0 0 2 0 0
14 D 4 4 20 36 0 24 16 4 0 28 8 0 0
15 D 1 1 8 18 0 12 7 1 0 16 2 0 0
15 I 1 1 4 6 0 4 3 1 0 4 2 0 0
16 U 1 1 6 12 0 8 6 0 0 12 0 0 0
16 I 1 1 3 2 2 0 2 1 0 0 2 0 0
16 D 1 2 9 14 2 8 8 1 0 12 2 0 0
17 U 1 1 7 15 0 10 6 0 1 12 3 0 0
17 I 1 1 4 6 0 4 3 0 1 3 3 0 0
17 D 1 1 7 15 0 10 6 0 1 12 3 0 0
18 U 1 1 7 13 1 8 6 0 1 12 1 0 0
18 I 1 1 2 1 1 0 1 0 1 0 1 0 0
18 D 1 1 7 13 1 8 6 0 1 12 1 0 0
19 U 2 2 7 12 0 8 6 0 1 12 0 0 0
19 I E
19 D 1 1 6 12 0 8 6 0 0 12 0 0 0
19-gap U 2 2 7 12 0 8 6 0 1 12 0 0 0
19-gap I E
19-gap D 1 1 6 12 0 8 6 0 0 12 0 0 0
19-bite U 1 1 7 13 1 8 6 0 1 12 1 0 0
19-bite I 1 1 2 1 1 0 1 0 1 0 1 0 0
19-bite D 1 1 7 13 1 8 6 0 1 12 1 0 0
20 U 1 1 6 12 0 8 6 0 0 12 0 0 0
20 I 1 1 1 0 0 0 0 0 1 0 0 0 0
20 D 1 2 7 12 0 8 6 0 1 12 0 0 0
21 U 1 1 6 12 0 8 6 0 0 12 0 0 0
21 I 1 1 1 0 0 0 0 0 1 0 0 0 0
21 D 1 2 7 12 0 8 6 0 1 12 0 0 0
22 D 1 1 7 13 1 8 6 0 1 12 1 0 0
22 I 1 1 2 1 1 0 1 0 1 0 1 0 0
23 U 1 1 7 6 6 0 4 3 0 0 4 0 2
23 I 1 1 3 2 2 0 0 3 0 0 0 0 2
23 D 1 1 4 4 4 0 2 2 0 0 2 0 2
23 R 2 2 6 4 4 0 2 4 0 0 2 0 2
24 U 1 1 8 8 4 2 4 4 0 0 4 4 0
24 I 1 1 4 4 0 2 0 4 0 0 0 4 0
24 D 2 2 6 6 2 4 2 4 0 0 2 4 0
25 U 1 1 4 6 0 4 2 2 0 2 4 0 0
25 I 1 1 4 6 0 4 2 2 0 2 4 0 0
25 D 1 1 4 6 0 4 2 2 0 2 4 0 0
26 U 2 2 6 4 4 0 4 2 0 0 4 0 0
26 I E
27 U 1 1 3 2 2 0 2 1 0 0 2 0 0
27 I E
27-tube U 1 1 7 6 6 0 4 3 0 0 6 0 0
27-tube I 1 1 3 2 2 0 2 1 0 0 2 0 0
27-tube D 1 1 4 4 4 0 2 2 0 0 4 0 0
28 U 1 1 7 6 4 1 4 3 0 0 4 0 2
28 I 1 1 3 2 0 1 0 3 0 0 0 0 2
28 R 2 2 6 4 2 2 2 4 0 0 2 0 2
28 D REFUSE
29 U 2 2 6 4 4 0 4 2 0 0 4 0 0
29-bite U 1 1 6 5 5 0 4 2 0 0 4 0 1
30 U 1 1 5 4 4 0 2 2 1 0 4 0 0
30 I 1 1 3 2 2 0 0 1 2 0 2 0 0
30 D 2 2 6 4 4 0 2 2 2 0 4 0 0
30 R 1 1 2 2 2 0 0 1 1 0 2 0 0
30-irrational U 1 1 5 4 4 0 2 2 1 0 4 0 0
30-irrational I 1 1 3 2 2 0 0 1 2 0 2 0 0
30-irrational D 2 2 6 4 4 0 2 2 2 0 4 0 0
30-irrational R 1 1 2 2 2 0 0 1 1 0 2 0 0
31 U 1 1 3 2 2 0 2 1 0 0 2 0 0
31 I 1 1 1 0 0 0 0 0 1 0 0 0 0
31 D 2 2 6 4 4 0 2 2 2 0 4 0 0
31 R E
32 U 1 1 5 4 2 1 2 2 1 0 2 0 2
32 I 1 1 3 2 0 1 0 1 2 0 0 0 2
32 D 2 2 6 4 2 2 2 2 2 0 2 0 2
32 R REFUSE
33 U 1 1 4 3 3 0 2 1 1 0 2 0 1
33 I 1 1 2 1 1 0 0 1 1 0 0 0 1
33 D 1 1 4 3 3 0 2 1 1 0 2 0 1
33 R 1 1 2 1 1 0 0 1 1 0 0 0 1
34 U 2 2 4 2 2 0 2 1 1 0 2 0 0
35 U 1 1 4 3 3 0 2 1 1 0 2 0 1
35 I 1 1 2 1 1 0 0 1 1 0 0 0 1
35 D 1 1 4 3 3 0 2 1 1 0 2 0 1
35 R 1 1 2 1 1 0 0 1 1 0 0 0 1
36 U 1 1 2 1 1 0 0 0 2 0 1 0 0
36 I 1 1 2 1 1 0 0 0 2 0 1 0 0
36 D 1 1 2 1 1 0 0 0 2 0 1 0 0
37 U 2 2 2 0 0 0 0 0 2 0 0 0 0
37 I E
37-gap U 2 2 2 0 0 0 0 0 2 0 0 0 0
37-gap I E
37-bite U 1 1 2 1 1 0 0 0 2 0 1 0 0
38 U 1 1 1 0 0 0 0 0 1 0 0 0 0
38 I 1 1 1 0 0 0 0 0 1 0 0 0 0
38 D 1 2 2 0 0 0 0 0 2 0 0 0 0
38 R E
39 U 1 1 1 0 0 0 0 0 1 0 0 0 0
39 I 1 1 1 0 0 0 0 0 1 0 0 0 0
39 D 1 2 2 0 0 0 0 0 2 0 0 0 0
40 U 1 1 1 0 0 0 0 0 1 0 0 0 0
40 I 1 1 1 0 0 0 0 0 1 0 0 0 0
40 D E
";
pub fn expectations(id: &str) -> Vec<(BoolOp, bool, Expected)> {
    TABLE
        .lines()
        .filter_map(|line| {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.first().copied() != Some(id) {
                return None;
            }
            let (op, reverse) = operation(fields[1]);
            let expected = match fields[2] {
                "E" => Expected::Counts(Counts::empty()),
                "REFUSE" => Expected::Refused("boolean.self-touching-shell", Stage::Sew),
                _ => {
                    let n: Vec<usize> = fields[2..].iter().map(|s| s.parse().unwrap()).collect();
                    assert_eq!(n.len(), 13);
                    Expected::Counts(Counts {
                        lumps: n[0],
                        shells: n[1],
                        faces: n[2],
                        edges: n[3],
                        rings: n[4],
                        vertices: n[5],
                        surfaces: [n[6], n[7], n[8]],
                        curves: [n[9], n[10], n[11], n[12]],
                    })
                }
            };
            Some((op, reverse, expected))
        })
        .collect()
}

// Every explicitly listed containment query in cases 1..40, with additional
// contact witnesses where the prose gives the contact's coordinates.
// Boundary expectations use the regularised result, not a Boolean of closed sets.
pub fn witnesses(id: &str, op: BoolOp) -> Vec<(RPoint3, Location)> {
    let code = match op {
        BoolOp::Unite => "U",
        BoolOp::Intersect => "I",
        BoolOp::Subtract => "D",
        _ => unreachable!(),
    };
    let table = "
01 U i 1/2 1/2 1/2
01 U i 5/2 5/2 5/2
01 U i 3/2 3/2 3/2
01 U o 5/2 1/2 1/2
01 U b 2 1/2 1/2
01 I i 3/2 3/2 3/2
01 I o 1/2 1/2 1/2
01 D i 1/2 1/2 1/2
01 D o 3/2 3/2 3/2
01 D b 3/2 3/2 1
02 U b 1 1/2 0
03 U i 2 3/2 1
04 U b 1 1 1/2
05 U b 1 1 1
06 D o 2 2 2
06 D i 1/2 2 2
10 D o 5 5 5
10 D o 5 5 10
10 D i 1 1 1
10 D i 5 17/2 5
10 D b 5 8 5
10 D b 5 8 0
10 U i 5 5 21/2
10 U i 5 5 5
10 I b 5 5 10
11 D o 5 5 7
11 D i 5 5 5
12 U b 5 5 10
14 D i 1 1 5
14 D b 0 5 5
15 D o 1 5 5
15 D i 4 5 5
19 U b 5 5 10
20 D b 5 5 10
21 D b 5 5 10
23 D o 0 0 0
23 D i 0 3/2 0
24 I b 0 1 0
24 I b 0 -1 0
28 I b 0 2 0
29 U b 0 2 0
32 I b -1 0 0
36 I i 3 0 0
36 U b 3 0 4
37 U b 5 0 0
";
    let mut result: Vec<_> = table
        .lines()
        .filter_map(|line| {
            let f: Vec<_> = line.split_whitespace().collect();
            if f.len() != 6 || f[0] != id || f[1] != code {
                return None;
            }
            Some((
                p(f[3], f[4], f[5]),
                match f[2] {
                    "i" => Location::Inside,
                    "o" => Location::Outside,
                    "b" => Location::On,
                    _ => unreachable!(),
                },
            ))
        })
        .collect();
    if id == "09-overlap" && op == BoolOp::Intersect {
        result.push((
            RPoint3::new(r(1) - q("1e-30") / r(2), q("1/2"), q("1/2")),
            Location::Inside,
        ));
    }
    if id == "13-bite" && op == BoolOp::Intersect {
        result.push((
            RPoint3::new(r(5), r(5), r(10) - q("1e-31")),
            Location::Inside,
        ));
    }
    result
}
