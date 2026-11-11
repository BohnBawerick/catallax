#![allow(dead_code)]
//! Shared reader and comparator for the SMLib differential gates.
//!
//! A fixture line is `operation kind | inputs | outputs`. Inputs are 16-digit
//! hexadecimal binary64 bit patterns or decimal integers. Outputs are tokens
//! `d<hex>` for binary64 values and `i<decimal>` for integers, statuses,
//! booleans and enumerations. Integers compare exactly. A binary64 value
//! passes when its bits match, when both are NaN, or when it is within
//! the operation's unit-in-the-last-place budget. Opposite signed zeros
//! differ. There is no absolute-error fallback near zero.

use catallax_trace::Fault;
use catallax_trace::extent::Extent1;
use catallax_trace::spline::{BSplineCurve, ControlPoint, NO_WEIGHT};
use catallax_trace::vector::{Vector2, Vector3};
use std::collections::BTreeMap;

/// SMLib's success status, `SM_SUCCESS`.
pub const SUCCESS: i64 = 1000;

/// SMLib's status value for a port fault.
pub fn status_code(fault: Fault) -> i64 {
    match fault {
        Fault::Degenerate => 1001,
        Fault::InvalidInput => 1007,
    }
}

/// Reads one fixture line's inputs in order.
pub struct Input<'a> {
    tokens: std::str::SplitWhitespace<'a>,
}

impl Input<'_> {
    fn next(&mut self) -> &str {
        self.tokens.next().expect("fixture input exhausted")
    }

    pub fn d(&mut self) -> f64 {
        f64::from_bits(u64::from_str_radix(self.next(), 16).expect("hex binary64 input"))
    }

    pub fn i(&mut self) -> i64 {
        self.next().parse().expect("integer input")
    }

    pub fn v2(&mut self) -> Vector2 {
        Vector2::new(self.d(), self.d())
    }

    pub fn v3(&mut self) -> Vector3 {
        Vector3::new(self.d(), self.d(), self.d())
    }

    pub fn interval(&mut self) -> Extent1 {
        Extent1::new(self.d(), self.d())
    }

    /// `degree count rational knots... (x y z w)...`, Euclidean coordinates
    /// with weights, stored weighted as SMLib stores them.
    pub fn curve(&mut self) -> BSplineCurve {
        let degree = self.i() as usize;
        let count = self.i() as usize;
        let rational = self.i() != 0;
        let knots: Vec<f64> = (0..count + degree + 1).map(|_| self.d()).collect();
        let points = (0..count)
            .map(|_| {
                let (x, y, z, w) = (self.d(), self.d(), self.d(), self.d());
                ControlPoint {
                    coordinates: if rational {
                        [x * w, y * w, z * w, w]
                    } else {
                        [x, y, z, NO_WEIGHT]
                    },
                }
            })
            .collect();
        BSplineCurve::new(degree, knots, points).expect("valid test curve")
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Token {
    D(f64),
    I(i64),
}

/// The port's outputs for one case, in the oracle's order.
#[derive(Default)]
pub struct Output {
    tokens: Vec<Token>,
}

impl Output {
    pub fn d(&mut self, value: f64) {
        self.tokens.push(Token::D(value));
    }

    pub fn i(&mut self, value: i64) {
        self.tokens.push(Token::I(value));
    }

    pub fn flag(&mut self, value: bool) {
        self.i(i64::from(value));
    }

    pub fn v2(&mut self, v: Vector2) {
        self.d(v.x);
        self.d(v.y);
    }

    pub fn v3(&mut self, v: Vector3) {
        self.d(v.x);
        self.d(v.y);
        self.d(v.z);
    }

    /// Emits SMLib's status for `result` and passes its value on.
    pub fn status<T>(&mut self, result: Result<T, Fault>) -> Option<T> {
        match result {
            Ok(value) => {
                self.i(SUCCESS);
                Some(value)
            }
            Err(fault) => {
                self.i(status_code(fault));
                None
            }
        }
    }
}

fn parse_expected(text: &str) -> Vec<Token> {
    text.split_whitespace()
        .map(|token| match token.split_at(1) {
            ("d", hex) => Token::D(f64::from_bits(
                u64::from_str_radix(hex, 16).expect("hex output"),
            )),
            ("i", int) => Token::I(int.parse().expect("integer output")),
            _ => panic!("bad output token {token}"),
        })
        .collect()
}

/// Distance in representable binary64 values between two finite numbers.
pub fn ulps(a: f64, b: f64) -> u64 {
    let key = |x: f64| {
        let bits = x.to_bits() as i64;
        if bits < 0 { i64::MIN - bits } else { bits }
    };
    key(a).abs_diff(key(b))
}

fn close(port: f64, oracle: f64, budget: u64) -> bool {
    if port.to_bits() == oracle.to_bits() || (port.is_nan() && oracle.is_nan()) {
        return true;
    }
    if port.is_nan() || oracle.is_nan() || port.is_infinite() || oracle.is_infinite() {
        return false;
    }
    if port == 0.0 && oracle == 0.0 {
        return false; // signed zeros differ
    }
    ulps(port, oracle) <= budget
}

/// Outcome of one gate run.
pub struct Report {
    /// Ordered output tags and bits for exact worker-count comparisons.
    pub candidate_bits: Vec<u64>,
    pub cases: usize,
    pub kinds: BTreeMap<String, usize>,
    pub failures: Vec<String>,
    pub max_ulps: u64,
    pub exact_values: usize,
    pub values: usize,
    /// Every binary64 output that is not bit-identical, for review.
    pub inexact: Vec<String>,
}

impl Report {
    /// Fails the test on any mismatch or a missing required case kind.
    pub fn assert_passed(&self, required_kinds: &[&str]) {
        eprintln!(
            "{} cases {:?}; {} of {} binary64 outputs bit-identical; largest gap {} ulp",
            self.cases, self.kinds, self.exact_values, self.values, self.max_ulps
        );
        for value in &self.inexact {
            eprintln!("  not bit-identical: {value}");
        }
        for kind in required_kinds {
            assert!(self.kinds.contains_key(*kind), "gate lacks {kind} cases");
        }
        assert!(
            self.failures.is_empty(),
            "{} of {} cases differ from SMLib:\n{}",
            self.failures.len(),
            self.cases,
            self.failures.join("\n")
        );
    }
}

/// Runs every fixture line through `evaluate` and compares with SMLib.
pub fn run_gate(
    fixture: &str,
    budget: impl Fn(&str) -> u64,
    evaluate: impl Fn(&str, &mut Input, &mut Output),
) -> Report {
    let mut report = Report {
        candidate_bits: Vec::new(),
        cases: 0,
        kinds: BTreeMap::new(),
        failures: Vec::new(),
        max_ulps: 0,
        exact_values: 0,
        values: 0,
        inexact: Vec::new(),
    };
    for (number, line) in fixture.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut parts = line.splitn(3, '|');
        let (head, inputs, outputs) = (
            parts.next().expect("head"),
            parts.next().expect("inputs"),
            parts.next().expect("outputs"),
        );
        let mut head = head.split_whitespace();
        let (op, kind) = (head.next().expect("operation"), head.next().expect("kind"));
        report.cases += 1;
        *report.kinds.entry(kind.to_owned()).or_default() += 1;
        let mut input = Input {
            tokens: inputs.split_whitespace(),
        };
        let mut output = Output::default();
        evaluate(op, &mut input, &mut output);
        for token in &output.tokens {
            match token {
                Token::D(value) => report.candidate_bits.extend([0, value.to_bits()]),
                Token::I(value) => report.candidate_bits.extend([1, *value as u64]),
            }
        }
        let expected = parse_expected(outputs);
        let mut problem = None;
        if expected.len() != output.tokens.len() {
            problem = Some(format!(
                "{} outputs, SMLib {} ({:?} vs {:?})",
                output.tokens.len(),
                expected.len(),
                output.tokens,
                expected
            ));
        }
        for (index, (port, oracle)) in output.tokens.iter().zip(&expected).enumerate() {
            match (port, oracle) {
                (Token::I(a), Token::I(b)) if a == b => {}
                (Token::D(a), Token::D(b)) => {
                    report.values += 1;
                    if a.to_bits() == b.to_bits() {
                        report.exact_values += 1;
                    } else {
                        if a.is_finite() && b.is_finite() {
                            report.max_ulps = report.max_ulps.max(ulps(*a, *b));
                        }
                        report.inexact.push(format!(
                            "line {} {op} {kind} output {index}: port {a:e}, SMLib {b:e}, {} ulp",
                            number + 1,
                            ulps(*a, *b)
                        ));
                    }
                    if !close(*a, *b, budget(op)) {
                        problem.get_or_insert(format!("output {index}: port {a:e}, SMLib {b:e}"));
                    }
                }
                _ => {
                    problem
                        .get_or_insert(format!("output {index}: port {port:?}, SMLib {oracle:?}"));
                }
            }
        }
        if let Some(problem) = problem {
            report
                .failures
                .push(format!("line {} {op} {kind}: {problem}", number + 1));
        }
    }
    report
}
