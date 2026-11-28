type Glyph = &'static [&'static [[i8; 2]]];

const GLYPHS: &[(char, Glyph)] = &[
    ('A', &[&[[0, 0], [2, 6], [4, 0]], &[[1, 3], [3, 3]]]),
    (
        'B',
        &[
            &[[0, 0], [0, 6], [3, 6], [4, 5], [4, 4], [3, 3], [0, 3]],
            &[[3, 3], [4, 2], [4, 1], [3, 0], [0, 0]],
        ],
    ),
    (
        'C',
        &[&[
            [4, 1],
            [3, 0],
            [1, 0],
            [0, 1],
            [0, 5],
            [1, 6],
            [3, 6],
            [4, 5],
        ]],
    ),
    (
        'H',
        &[&[[0, 0], [0, 6]], &[[4, 0], [4, 6]], &[[0, 3], [4, 3]]],
    ),
    (
        'I',
        &[&[[1, 0], [3, 0]], &[[2, 0], [2, 6]], &[[1, 6], [3, 6]]],
    ),
    (
        'O',
        &[&[
            [1, 0],
            [0, 1],
            [0, 5],
            [1, 6],
            [3, 6],
            [4, 5],
            [4, 1],
            [3, 0],
            [1, 0],
        ]],
    ),
    (
        'P',
        &[&[[0, 0], [0, 6], [3, 6], [4, 5], [4, 4], [3, 3], [0, 3]]],
    ),
    (
        'R',
        &[
            &[[0, 0], [0, 6], [3, 6], [4, 5], [4, 4], [3, 3], [0, 3]],
            &[[2, 3], [4, 0]],
        ],
    ),
    (
        'S',
        &[&[
            [4, 5],
            [3, 6],
            [1, 6],
            [0, 5],
            [0, 4],
            [1, 3],
            [3, 3],
            [4, 2],
            [4, 1],
            [3, 0],
            [1, 0],
            [0, 1],
        ]],
    ),
    (
        'a',
        &[
            &[[4, 4], [4, 0]],
            &[
                [4, 3],
                [3, 4],
                [1, 4],
                [0, 3],
                [0, 1],
                [1, 0],
                [3, 0],
                [4, 1],
            ],
        ],
    ),
    (
        '0',
        &[
            &[
                [1, 0],
                [0, 1],
                [0, 5],
                [1, 6],
                [3, 6],
                [4, 5],
                [4, 1],
                [3, 0],
                [1, 0],
            ],
            &[[1, 1], [3, 5]],
        ],
    ),
    ('1', &[&[[1, 5], [2, 6], [2, 0]], &[[1, 0], [3, 0]]]),
    (
        '2',
        &[&[[0, 5], [1, 6], [3, 6], [4, 5], [4, 4], [0, 0], [4, 0]]],
    ),
    (
        '5',
        &[&[
            [4, 6],
            [0, 6],
            [0, 3],
            [3, 3],
            [4, 2],
            [4, 1],
            [3, 0],
            [1, 0],
            [0, 1],
        ]],
    ),
    (
        '6',
        &[&[
            [4, 5],
            [3, 6],
            [1, 6],
            [0, 5],
            [0, 1],
            [1, 0],
            [3, 0],
            [4, 1],
            [4, 2],
            [3, 3],
            [0, 3],
        ]],
    ),
    ('7', &[&[[0, 6], [4, 6], [1, 0]]]),
    (
        '8',
        &[&[
            [1, 3],
            [0, 4],
            [0, 5],
            [1, 6],
            [3, 6],
            [4, 5],
            [4, 4],
            [3, 3],
            [1, 3],
            [0, 2],
            [0, 1],
            [1, 0],
            [3, 0],
            [4, 1],
            [4, 2],
            [3, 3],
        ]],
    ),
    ('.', &[&[[1, 0], [1, 1]]]),
    ('|', &[&[[2, -1], [2, 6]]]),
    (
        '⌀',
        &[
            &[
                [1, 1],
                [0, 2],
                [0, 4],
                [1, 5],
                [3, 5],
                [4, 4],
                [4, 2],
                [3, 1],
                [1, 1],
            ],
            &[[0, 0], [4, 6]],
        ],
    ),
];

pub fn label(wire: &[Vec<[f64; 2]>], finish: bool) -> (String, f64) {
    let symbol = &wire[0];
    let strokes = &wire[if finish { 2 } else { 1 }..];
    let (step, origin) = if finish {
        let step = (symbol[2][1] - symbol[1][1]) / 18.0;
        let origin = [0, 1].map(|axis| {
            strokes
                .iter()
                .flatten()
                .map(|p| p[axis])
                .fold(f64::INFINITY, f64::min)
        });
        (step, origin)
    } else {
        let step = (symbol[2][1] - symbol[0][1]) / 10.0;
        (step, symbol[0].map(|coordinate| coordinate + 2.0 * step))
    };
    let mut pen = origin[0];
    let mut remaining = strokes;
    let mut text = String::new();
    while !remaining.is_empty() {
        let (gap, character, lines) = [0.0, 4.0]
            .into_iter()
            .find_map(|gap| {
                GLYPHS
                    .iter()
                    .filter(|(_, expected)| {
                        remaining.get(..expected.len()).is_some_and(|actual| {
                            actual.iter().zip(*expected).all(|(path, glyph)| {
                                path.len() == glyph.len()
                                    && path.iter().zip(*glyph).all(|([x, y], [u, v])| {
                                        ((x - pen) / step - gap - f64::from(*u)).abs() < 1e-6
                                            && ((y - origin[1]) / step - f64::from(*v)).abs() < 1e-6
                                    })
                            })
                        })
                    })
                    .max_by_key(|(_, glyph)| glyph.len())
                    .map(|(character, glyph)| (gap, *character, glyph.len()))
            })
            .unwrap_or_else(|| panic!("unreadable glyph after {text:?}: {:?}", remaining[0]));
        if gap != 0.0 {
            text.push(' ');
        }
        text.push(character);
        pen += (gap + if character == '.' { 3.0 } else { 6.0 }) * step;
        remaining = &remaining[lines..];
    }
    (text, step)
}
