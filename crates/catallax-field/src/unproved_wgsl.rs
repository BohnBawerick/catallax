//! Display-only WGSL translation of the compiled field tape.
use crate::{Instruction, Plan};
use catallax_numbers::Refusal;
use catallax_recipe::BooleanKind;

pub(crate) fn number(x: f64) -> String {
    let rounded = x as f32;
    if !rounded.is_finite() || (x != 0.0 && rounded == 0.0) {
        return "invalid_range".into();
    }
    format!("{rounded:.9e}")
}
pub(crate) fn vector(x: &[f64]) -> String {
    format!(
        "vec{}<f32>({})",
        x.len(),
        x.iter().map(|x| number(*x)).collect::<Vec<_>>().join(",")
    )
}
pub(crate) fn frame(matrix: &[[f64; 3]; 3], translation: &[f64; 3]) -> String {
    format!(
        "vec3<f32>(dot({},p),dot({},p),dot({},p))+{}",
        vector(&matrix[0]),
        vector(&matrix[1]),
        vector(&matrix[2]),
        vector(translation)
    )
}
pub(crate) fn unsupported() -> Refusal {
    Refusal::new(
        "field.preview-wgsl-unsupported",
        "this leaf has no unproved WGSL translation",
    )
}
impl Plan {
    /// Translates supported tape leaves to unproved binary32 WGSL for display.
    /// This source carries no error bound and must never feed an export or proof.
    /// Unsupported leaves refuse the whole translation, without a partial shader.
    pub fn unproved_preview_wgsl(&self) -> Result<String, Refusal> {
        let mut source = String::new();
        for (i, leaf) in self.leaves.iter().enumerate() {
            source.push_str(&format!(
                "fn leaf_{i}(p: vec3<f32>) -> f32 {{ {} }}\n",
                leaf.wgsl().map_err(|mut refusal| {
                    refusal.detail = format!(
                        "leaf {i} of {}, {} result face records: {}",
                        self.leaves.len(),
                        self.faces
                            .iter()
                            .filter(|face| face.leaf == i && face.in_result)
                            .count(),
                        refusal.detail
                    );
                    refusal
                })?
            ));
        }
        source.push_str("fn field(p: vec3<f32>) -> f32 {\n");
        for (i, instruction) in self.instructions.iter().enumerate() {
            let expression = match instruction {
                Instruction::Empty => "1e30".into(),
                Instruction::Leaf(leaf) => format!("leaf_{leaf}(p)"),
                Instruction::Boolean { kind, a, b } => match kind {
                    BooleanKind::Union => format!("min(v{a},v{b})"),
                    BooleanKind::Intersect => format!("max(v{a},v{b})"),
                    BooleanKind::Subtract => format!("max(v{a},-v{b})"),
                },
            };
            source.push_str(&format!("let v{i} = {expression};\n"));
        }
        source.push_str(&format!(
            "return select(v{0},0.0,v{0}==0.0); }}\n",
            self.root
        ));
        if source.contains("invalid_range") {
            return Err(Refusal::new(
                "field.preview-wgsl-range",
                "field constants exceed finite nonzero binary32",
            ));
        }
        Ok(source)
    }
}
