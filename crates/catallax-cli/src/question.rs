//! The question file, format `exact-brep-question/2`.

use catallax_boolean::BoolOp;
use catallax_measure::Accuracy;
use catallax_numbers::Rational;
use catallax_shapes::Pose;
use catallax_step::StepOptions;

/// A shape and its exact parameters, in its local frame.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum ShapeSpec {
    /// The box `[0, sx] x [0, sy] x [0, sz]`.
    Block {
        /// `[sx, sy, sz]`.
        size: [Rational; 3],
    },
    /// A cylinder on the `z` axis from 0 to `height`.
    Cylinder {
        /// Radius.
        radius: Rational,
        /// Height.
        height: Rational,
    },
    /// A sphere at the origin.
    Sphere {
        /// Radius.
        radius: Rational,
    },
    /// A cone frustum on the `z` axis.
    Cone {
        /// Radius at `z = 0`.
        base_radius: Rational,
        /// Radius at `z = height`; 0 for an apex.
        top_radius: Rational,
        /// Height.
        height: Rational,
    },
    /// A torus at the origin about the `z` axis.
    Torus {
        /// Radius of the core circle.
        major: Rational,
        /// Radius of the tube.
        minor: Rational,
    },
}

/// One named operand.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Operand {
    /// The shape.
    pub shape: ShapeSpec,
    /// Where it is placed.
    pub pose: Pose,
}

/// One boolean step.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Step {
    /// Step id.
    pub id: String,
    /// The operation.
    pub op: BoolOp,
    /// Names of the two arguments: operands or earlier steps.
    pub args: [String; 2],
}

/// Which properties to measure, and how tightly.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MeasureRequest {
    /// Measure the volume.
    pub volume: bool,
    /// Measure the area.
    pub area: bool,
    /// Measure the centre of mass.
    pub centroid: bool,
    /// Requested width.
    pub accuracy: Accuracy,
}

/// Whether and how to write STEP.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ExportRequest {
    /// Write `answer.step`.
    pub step: bool,
    /// Export options; the product name is the question id.
    pub options: StepOptions,
}

/// Limits on the whole answer.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Wall-clock limit in milliseconds, counted from the program's start; default 50,000.
    pub deadline_ms: u64,
}

/// A parsed, valid question.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Question {
    /// Question id.
    pub id: String,
    /// Named operands, in file order.
    pub operands: Vec<(String, Operand)>,
    /// Steps in order.
    pub steps: Vec<Step>,
    /// Name of the operand or step that is the answer.
    pub result: String,
    /// Measurement request.
    pub measure: MeasureRequest,
    /// Export request.
    pub export: ExportRequest,
    /// Limits.
    pub limits: Limits,
}

/// Why a question is invalid. Becomes a `cli.invalid-question` refusal.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct QuestionError {
    /// The question id, if it could be read.
    pub question_id: Option<String>,
    /// Free text for people.
    pub detail: String,
}

/// Parses and validates a question file.
pub fn parse_question(text: &str) -> Result<Question, QuestionError> {
    parse(text).map_err(|detail| QuestionError {
        question_id: read_id(text),
        detail,
    })
}

fn read_id(text: &str) -> Option<String> {
    struct Visitor<'a>(&'a mut Option<String>);
    impl<'de> serde::de::Visitor<'de> for Visitor<'_> {
        type Value = ();

        fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("a question object")
        }

        fn visit_map<A: serde::de::MapAccess<'de>>(self, mut map: A) -> Result<(), A::Error> {
            while let Some(key) = map.next_key::<String>()? {
                if key == "id" {
                    *self.0 = Some(map.next_value::<String>()?);
                } else {
                    map.next_value::<serde::de::IgnoredAny>()?;
                }
            }
            Ok(())
        }
    }
    let mut id = None;
    let mut deserializer = serde_json::Deserializer::from_str(text);
    let _ = serde::Deserializer::deserialize_map(&mut deserializer, Visitor(&mut id));
    id
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct RawQuestion {
    format: String,
    id: String,
    operands: Operands,
    steps: Vec<RawStep>,
    result: String,
    #[serde(default, deserialize_with = "present")]
    measure: Option<RawMeasure>,
    #[serde(default, deserialize_with = "present")]
    export: Option<RawExport>,
    #[serde(default, deserialize_with = "present")]
    limits: Option<RawLimits>,
}

struct Number(Rational);

impl<'de> serde::Deserialize<'de> for Number {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let text = <String as serde::Deserialize>::deserialize(deserializer)?;
        Rational::parse(&text)
            .map(Self)
            .map_err(serde::de::Error::custom)
    }
}

struct Operands(Vec<(String, RawOperand)>);

impl<'de> serde::Deserialize<'de> for Operands {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Visitor;
        impl<'de> serde::de::Visitor<'de> for Visitor {
            type Value = Operands;

            fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str("an object of uniquely named operands")
            }

            fn visit_map<A: serde::de::MapAccess<'de>>(
                self,
                mut map: A,
            ) -> Result<Operands, A::Error> {
                let mut values = Vec::new();
                let mut names = std::collections::BTreeSet::new();
                while let Some((name, operand)) = map.next_entry::<String, RawOperand>()? {
                    if !names.insert(name.clone()) {
                        return Err(serde::de::Error::custom("duplicate operand name"));
                    }
                    values.push((name, operand));
                }
                Ok(Operands(values))
            }
        }
        deserializer.deserialize_map(Visitor)
    }
}

#[derive(serde::Deserialize)]
#[serde(tag = "shape", rename_all = "lowercase", deny_unknown_fields)]
enum RawOperand {
    Block {
        size: [Number; 3],
        #[serde(default, deserialize_with = "present")]
        pose: Option<RawPose>,
    },
    Cylinder {
        radius: Number,
        height: Number,
        #[serde(default, deserialize_with = "present")]
        pose: Option<RawPose>,
    },
    Sphere {
        radius: Number,
        #[serde(default, deserialize_with = "present")]
        pose: Option<RawPose>,
    },
    Cone {
        base_radius: Number,
        top_radius: Number,
        height: Number,
        #[serde(default, deserialize_with = "present")]
        pose: Option<RawPose>,
    },
    Torus {
        major: Number,
        minor: Number,
        #[serde(default, deserialize_with = "present")]
        pose: Option<RawPose>,
    },
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPose {
    rotation: [Number; 4],
    translation: [Number; 3],
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct RawStep {
    id: String,
    op: String,
    args: [String; 2],
}

#[derive(serde::Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct RawMeasure {
    #[serde(default, deserialize_with = "present")]
    properties: Option<Vec<String>>,
    #[serde(default, deserialize_with = "present")]
    relative: Option<Number>,
    #[serde(default, deserialize_with = "present")]
    absolute: Option<Number>,
}

#[derive(serde::Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct RawExport {
    #[serde(default, deserialize_with = "present")]
    step: Option<bool>,
    #[serde(default, deserialize_with = "present")]
    curve_tolerance: Option<Number>,
    #[serde(default, deserialize_with = "present")]
    significant_digits: Option<u32>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct RawLimits {
    deadline_ms: u64,
}

fn present<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: serde::Deserialize<'de>,
{
    T::deserialize(deserializer).map(Some)
}

fn parse(text: &str) -> Result<Question, String> {
    let raw: RawQuestion = serde_json::from_str(text).map_err(|error| error.to_string())?;
    if raw.format != crate::QUESTION_FORMAT {
        return Err("unsupported question format".into());
    }
    if raw.id.is_empty()
        || raw.id.len() > 200
        || !raw
            .id
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c))
    {
        return Err("invalid question id".into());
    }
    if raw.operands.0.is_empty() {
        return Err("at least one operand is required".into());
    }
    let mut names = std::collections::BTreeSet::new();
    let mut operands = Vec::with_capacity(raw.operands.0.len());
    for (name, raw_operand) in raw.operands.0 {
        validate_name(&name)?;
        names.insert(name.clone());
        let (shape, pose) = match raw_operand {
            RawOperand::Block { size, pose } => {
                let size = size.map(|n| n.0);
                for n in &size {
                    positive(n)?;
                }
                (ShapeSpec::Block { size }, pose)
            }
            RawOperand::Cylinder {
                radius,
                height,
                pose,
            } => {
                positive(&radius.0)?;
                positive(&height.0)?;
                (
                    ShapeSpec::Cylinder {
                        radius: radius.0,
                        height: height.0,
                    },
                    pose,
                )
            }
            RawOperand::Sphere { radius, pose } => {
                positive(&radius.0)?;
                (ShapeSpec::Sphere { radius: radius.0 }, pose)
            }
            RawOperand::Cone {
                base_radius,
                top_radius,
                height,
                pose,
            } => {
                positive(&base_radius.0)?;
                nonnegative(&top_radius.0)?;
                positive(&height.0)?;
                if base_radius.0 == top_radius.0 {
                    return Err("cone radii must differ".into());
                }
                (
                    ShapeSpec::Cone {
                        base_radius: base_radius.0,
                        top_radius: top_radius.0,
                        height: height.0,
                    },
                    pose,
                )
            }
            RawOperand::Torus { major, minor, pose } => {
                positive(&minor.0)?;
                if minor.0 >= major.0 {
                    return Err("torus requires minor < major".into());
                }
                (
                    ShapeSpec::Torus {
                        major: major.0,
                        minor: minor.0,
                    },
                    pose,
                )
            }
        };
        let pose = if let Some(pose) = pose {
            let [w, x, y, z] = pose.rotation.map(|n| n.0);
            let rotation = catallax_shapes::Rotation::from_quaternion(w, x, y, z)
                .map_err(|error| error.detail)?;
            let [x, y, z] = pose.translation.map(|n| n.0);
            Pose::new(rotation, catallax_shapes::RVector3::new(x, y, z))
        } else {
            Pose::identity()
        };
        operands.push((name, Operand { shape, pose }));
    }
    let mut steps = Vec::with_capacity(raw.steps.len());
    for step in raw.steps {
        validate_name(&step.id)?;
        if names.contains(&step.id) {
            return Err("operand and step names must be unique".into());
        }
        if step.args.iter().any(|name| !names.contains(name)) {
            return Err("step arguments must name operands or earlier steps".into());
        }
        let op = match step.op.as_str() {
            "unite" => BoolOp::Unite,
            "subtract" => BoolOp::Subtract,
            "intersect" => BoolOp::Intersect,
            _ => return Err("unknown boolean operation".into()),
        };
        names.insert(step.id.clone());
        steps.push(Step {
            id: step.id,
            op,
            args: step.args,
        });
    }
    if !names.contains(&raw.result) {
        return Err("result must name an operand or step".into());
    }
    let measure = raw.measure.unwrap_or_default();
    let mut accuracy = Accuracy::default();
    if let Some(n) = measure.relative {
        nonnegative(&n.0)?;
        accuracy.relative = n.0;
    }
    if let Some(n) = measure.absolute {
        nonnegative(&n.0)?;
        accuracy.absolute = n.0;
    }
    let mut request = MeasureRequest {
        volume: false,
        area: false,
        centroid: false,
        accuracy,
    };
    let properties = measure
        .properties
        .unwrap_or_else(|| vec!["volume".into(), "area".into(), "centroid".into()]);
    for property in properties {
        match property.as_str() {
            "volume" => request.volume = true,
            "area" => request.area = true,
            "centroid" => request.centroid = true,
            _ => return Err("unknown measured property".into()),
        }
    }
    let export = raw.export.unwrap_or_default();
    let mut options = StepOptions {
        product_name: raw.id.clone(),
        ..StepOptions::default()
    };
    if let Some(n) = export.curve_tolerance {
        positive(&n.0)?;
        options.curve_tolerance = n.0;
    }
    if let Some(n) = export.significant_digits {
        if !(6..=40).contains(&n) {
            return Err("significant_digits must be from 6 to 40".into());
        }
        options.significant_digits = n;
    }
    Ok(Question {
        id: raw.id,
        operands,
        steps,
        result: raw.result,
        measure: request,
        export: ExportRequest {
            step: export.step.unwrap_or(true),
            options,
        },
        limits: Limits {
            deadline_ms: raw.limits.map_or(50_000, |limits| limits.deadline_ms),
        },
    })
}

fn validate_name(name: &str) -> Result<(), String> {
    if name.is_empty()
        || name.len() > 64
        || !name
            .bytes()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || b"_-".contains(&c))
    {
        Err("invalid operand or step name".into())
    } else {
        Ok(())
    }
}

fn positive(n: &Rational) -> Result<(), String> {
    if n.sign() == catallax_numbers::Sign::Positive {
        Ok(())
    } else {
        Err("expected a positive geometric value".into())
    }
}

fn nonnegative(n: &Rational) -> Result<(), String> {
    if n.sign() == catallax_numbers::Sign::Negative {
        Err("expected a nonnegative value".into())
    } else {
        Ok(())
    }
}
