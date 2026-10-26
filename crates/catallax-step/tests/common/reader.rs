use super::*;

// This reader recognizes the contract's STEP subset, not arbitrary STEP input.
// Quoted #, comma, apostrophe and parentheses must not affect graph traversal.
pub fn split(s: &str) -> Vec<&str> {
    let bytes = s.as_bytes();
    let (mut depth, mut quoted, mut start, mut i) = (0_i32, false, 0, 0);
    let mut out = Vec::new();
    while i < bytes.len() {
        match bytes[i] {
            b'\'' if quoted && bytes.get(i + 1) == Some(&b'\'') => {
                i += 2;
                continue;
            }
            b'\'' => quoted = !quoted,
            b'(' if !quoted => depth += 1,
            b')' if !quoted => {
                depth -= 1;
                assert!(depth >= 0);
            }
            b',' if !quoted && depth == 0 => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
        i += 1;
    }
    assert!(!quoted && depth == 0, "unbalanced STEP argument: {s}");
    if !s.is_empty() {
        out.push(&s[start..]);
    }
    out
}
pub fn list(s: &str) -> Vec<&str> {
    split(s.strip_prefix('(').unwrap().strip_suffix(')').unwrap())
}
pub fn reference(s: &str) -> usize {
    s.strip_prefix('#').unwrap().parse().unwrap()
}
pub fn refs(s: &str) -> Vec<usize> {
    let b = s.as_bytes();
    let (mut i, mut quoted) = (0, false);
    let mut out = Vec::new();
    while i < b.len() {
        if b[i] == b'\'' {
            if quoted && b.get(i + 1) == Some(&b'\'') {
                i += 2;
                continue;
            }
            quoted = !quoted;
        } else if !quoted && b[i] == b'#' {
            i += 1;
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            out.push(s[start..i].parse().unwrap());
            continue;
        }
        i += 1;
    }
    assert!(!quoted);
    out
}
pub fn string(s: &str) -> String {
    let mut out = String::from("'");
    for c in s.chars() {
        match c {
            '\'' => out.push_str("''"),
            '\\' => out.push_str("\\\\"),
            ' '..='~' => out.push(c),
            _ => {
                out.push_str("\\X2\\");
                for u in c.encode_utf16(&mut [0; 2]) {
                    out.push_str(&format!("{u:04X}"));
                }
                out.push_str("\\X0\\");
            }
        }
    }
    out.push('\'');
    out
}
pub fn number(s: &str) -> Rational {
    assert!(!s.contains('e') && !s.starts_with('+'));
    assert!(
        s.contains('.') || s.contains('E'),
        "real without STEP decimal: {s}"
    );
    let normalized = if let Some((m, e)) = s.split_once('E') {
        format!("{}e{e}", m.trim_end_matches('.'))
    } else {
        s.trim_end_matches('.').to_owned()
    };
    r(&normalized)
}
pub fn decimal(q: &Rational, digits: u32, upper: bool) -> String {
    // Piece 1 is finished; its documented directed rounding is the contract.
    let (lo, hi) = q.to_decimal_bounds(digits);
    let s = if upper { hi } else { lo };
    if let Some((m, e)) = s.split_once('e') {
        format!(
            "{}E{e}",
            if m.contains('.') {
                m.to_owned()
            } else {
                format!("{m}.")
            }
        )
    } else if s.contains('.') {
        s
    } else {
        format!("{s}.")
    }
}
pub fn rounded(s: &str, exact: &Rational, digits: u32) {
    assert!(
        s == decimal(exact, digits, false) || s == decimal(exact, digits, true),
        "{s} is not either directed {digits}-digit rounding of {exact}"
    );
}

#[derive(Debug)]
pub struct Entity {
    pub text: String,
    pub kind: String,
    pub args: Vec<String>,
}
#[derive(Debug)]
pub struct Step {
    pub entities: BTreeMap<usize, Entity>,
}
impl Step {
    pub fn parse(bytes: &[u8], name: &str) -> Self {
        assert!(bytes.is_ascii());
        let text = std::str::from_utf8(bytes).unwrap();
        assert!(text.ends_with("ENDSEC;\nEND-ISO-10303-21;\n"));
        let header = format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('catallax export'),'2;1');\nFILE_NAME({},'1970-01-01T00:00:00',(''),(''),'catallax','catallax','');\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF {{ 1 0 10303 442 1 1 4 }}'));\nENDSEC;\nDATA;\n",
            string(&format!("{name}.step"))
        );
        let body = text
            .strip_prefix(&header)
            .expect("fixed header")
            .strip_suffix("ENDSEC;\nEND-ISO-10303-21;\n")
            .unwrap();
        let mut entities = BTreeMap::new();
        for (i, line) in body.lines().enumerate() {
            let (label, value) = line.split_once('=').unwrap();
            assert_eq!(reference(label), i + 1);
            let text = value.strip_suffix(';').unwrap().to_owned();
            let mut quoted = false;
            for c in text.chars() {
                if c == '\'' {
                    quoted = !quoted;
                }
                if !quoted {
                    assert!(!c.is_ascii_whitespace(), "whitespace outside string");
                }
            }
            assert!(!quoted);
            let (kind, args) = if text.starts_with('(') {
                (String::new(), Vec::new())
            } else {
                let (k, a) = text.split_once('(').unwrap();
                (
                    k.to_owned(),
                    split(a.strip_suffix(')').unwrap())
                        .into_iter()
                        .map(str::to_owned)
                        .collect(),
                )
            };
            entities.insert(i + 1, Entity { text, kind, args });
        }
        let out = Self { entities };
        for e in out.entities.values() {
            for id in refs(&e.text) {
                assert!(out.entities.contains_key(&id), "dangling #{id}");
            }
        }
        let mut seen = (1..=20).collect::<BTreeSet<_>>();
        let mut order = Vec::new();
        for id in out.items() {
            out.walk(id, &mut seen, &mut order);
        }
        assert_eq!(
            order,
            (21..=out.entities.len()).collect::<Vec<_>>(),
            "depth-first numbering or orphan entity"
        );
        out
    }
    fn walk(&self, id: usize, seen: &mut BTreeSet<usize>, order: &mut Vec<usize>) {
        if seen.insert(id) {
            order.push(id);
            for child in refs(&self.entities[&id].text) {
                self.walk(child, seen, order);
            }
        }
    }
    pub fn e(&self, id: usize, kind: &str) -> &Entity {
        let e = &self.entities[&id];
        assert_eq!(e.kind, kind, "#{id}");
        e
    }
    pub fn ids(&self, kind: &str) -> Vec<usize> {
        self.entities
            .iter()
            .filter_map(|(id, e)| (e.kind == kind).then_some(*id))
            .collect()
    }
    pub fn count(&self, kind: &str) -> usize {
        self.ids(kind).len()
    }
    pub fn items(&self) -> Vec<usize> {
        list(&self.e(16, "ADVANCED_BREP_SHAPE_REPRESENTATION").args[1])
            .into_iter()
            .map(reference)
            .collect()
    }
    pub fn xyz_text(&self, id: usize) -> Vec<&str> {
        list(&self.e(id, "CARTESIAN_POINT").args[1])
    }
    pub fn xyz(&self, id: usize) -> [Rational; 3] {
        self.xyz_text(id)
            .into_iter()
            .map(number)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }
    pub fn direction(&self, id: usize) -> [Rational; 3] {
        list(&self.e(id, "DIRECTION").args[1])
            .into_iter()
            .map(number)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }
    pub fn vertex(&self, id: usize) -> usize {
        reference(&self.e(id, "VERTEX_POINT").args[1])
    }
    pub fn placement(&self, id: usize) -> ([Rational; 3], [Rational; 3], [Rational; 3]) {
        let a = &self.e(id, "AXIS2_PLACEMENT_3D").args;
        (
            self.xyz(reference(&a[1])),
            self.direction(reference(&a[2])),
            self.direction(reference(&a[3])),
        )
    }
    pub fn fixed(&self, name: &str, uncertainty: &str) {
        let n = string(name);
        let items = self
            .items()
            .iter()
            .map(|id| format!("#{id}"))
            .collect::<Vec<_>>()
            .join(",");
        let fixed = [
            "APPLICATION_CONTEXT('managed model based 3d engineering')".into(),
            "APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1)".into(),
            "PRODUCT_CONTEXT('',#1,'mechanical')".into(),
            format!("PRODUCT({n},{n},'',(#3))"),
            "PRODUCT_DEFINITION_FORMATION('','',#4)".into(),
            "PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design')".into(),
            "PRODUCT_DEFINITION('design','',#5,#6)".into(),
            "PRODUCT_DEFINITION_SHAPE('','',#7)".into(),
            "SHAPE_DEFINITION_REPRESENTATION(#8,#16)".into(),
            "(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))".into(),
            "(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))".into(),
            "(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())".into(),
            format!("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty}),#10,'distance_accuracy_value','maximum deviation from the exact model')"),
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''))".into(),
            "PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#4))".into(),
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION({n},({items}),#14)"),
            "AXIS2_PLACEMENT_3D('',#18,#19,#20)".into(),
            "CARTESIAN_POINT('',(0.,0.,0.))".into(),
            "DIRECTION('',(0.,0.,1.))".into(),
            "DIRECTION('',(1.,0.,0.))".into(),
        ];
        for (i, text) in fixed.iter().enumerate() {
            assert_eq!(&self.entities[&(i + 1)].text, text);
        }
    }
}
