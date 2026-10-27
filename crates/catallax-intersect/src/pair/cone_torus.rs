//! Parallel offset cone/torus pairs, clipped on the torus's minor-angle circle.

use super::*;
use catallax_shapes::{Cone, Torus};

#[derive(Clone, Debug)]
pub(crate) struct Coordinates {
    pub origin: Point3,
    pub e1: Vector3,
    pub e2: Vector3,
    pub e3: Vector3,
    pub major: Real,
    pub minor: Real,
    inverse_minor: Real,
    inverse_twice_radii: Real,
    squared_radii: Real,
}
impl Coordinates {
    pub fn minor_angle(&self, point: &Point3) -> [Real; 2] {
        let delta = point.minus(&self.origin);
        // On a ring torus, |p|² = R²+r²+2Rr cos(phi). No radial root is needed.
        [
            &(&delta.norm_squared() - &self.squared_radii) * &self.inverse_twice_radii,
            delta.dot(&self.e3) * self.inverse_minor.clone(),
        ]
    }
}

pub(crate) fn quarter_parameter(
    angle: &[Real; 2],
    quarter: usize,
    dec: &mut Decider,
) -> Result<Option<Real>, Failure> {
    let [cosine, sine] = angle;
    let (cosine, sine) = match quarter {
        0 => (cosine.clone(), sine.clone()),
        1 => (sine.clone(), -cosine),
        2 => (-cosine, -sine),
        3 => (-sine, cosine.clone()),
        _ => unreachable!("minor-angle quarter"),
    };
    if sign(&cosine, dec)? == Sign::Negative || sign(&sine, dec)? == Sign::Negative {
        return Ok(None);
    }
    Ok(Some(sine / (s(1) + cosine)))
}

#[derive(Clone, Debug)]
pub(crate) struct Polynomial {
    pub degree: usize,
    pub bernstein: [Real; 3],
}
#[derive(Clone, Debug)]
pub(crate) struct Sector {
    pub quarter: usize,
    pub lo: Real,
    pub span: Real,
    pub orders: [usize; 2],
    pub factors: [Polynomial; 4],
}
#[derive(Clone, Debug)]
pub(crate) struct Domain {
    pub start: usize,
    pub count: usize,
    /// Zero joins both sheets at the endpoints; +/-1 keeps a full-circle sheet.
    pub sheet: i64,
}
#[derive(Clone, Debug)]
pub(crate) struct ConeTorusFrame {
    pub coordinates: Arc<Coordinates>,
    pub distance: Real,
    pub radius_base: Real,
    pub radius_sine: Real,
    pub sectors: Vec<Sector>,
    pub domains: Vec<Domain>,
    pub sense: i64,
}
impl ConeTorusFrame {
    fn contains_angle(
        &self,
        domain: &Domain,
        angle: &[Real; 2],
        dec: &mut Decider,
    ) -> Result<bool, Failure> {
        for offset in 0..domain.count {
            let sector = &self.sectors[(domain.start + offset) % self.sectors.len()];
            let Some(value) = quarter_parameter(angle, sector.quarter, dec)? else {
                continue;
            };
            if compare(&value, &sector.lo, dec)? != Sign::Negative
                && compare(&value, &(&sector.lo + &sector.span), dec)? != Sign::Positive
            {
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub fn component(&self, witness: &Point3, dec: &mut Decider) -> Result<usize, Failure> {
        let angle = self.coordinates.minor_angle(witness);
        for (index, domain) in self.domains.iter().enumerate() {
            if domain.sheet != 0
                && sign(
                    &(witness
                        .minus(&self.coordinates.origin)
                        .dot(&self.coordinates.e2)
                        * s(domain.sheet)),
                    dec,
                )? != Sign::Positive
            {
                continue;
            }
            if self.contains_angle(domain, &angle, dec)? {
                return Ok(index);
            }
        }
        Err(Refusal::new(
            "intersect.invalid-branch-witness",
            "the cone/torus witness is outside every retained component",
        )
        .into())
    }
    pub fn contains(
        &self,
        component: usize,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Failure> {
        let domain = &self.domains[component];
        if domain.sheet != 0
            && sign(
                &(point
                    .minus(&self.coordinates.origin)
                    .dot(&self.coordinates.e2)
                    * s(domain.sheet)),
                dec,
            )? == Sign::Negative
        {
            return Ok(false);
        }
        self.contains_angle(domain, &self.coordinates.minor_angle(point), dec)
    }
    fn position(&self, quarter: usize, parameter: &Real, sheet: i64) -> Point3 {
        let [cosine, sine] = angle(quarter, parameter);
        let rho = &self.coordinates.major + &(&self.coordinates.minor * &cosine);
        let radius = &self.radius_base + &(&self.radius_sine * &sine);
        let x = (rho.square() + self.distance.square() - radius.square())
            / (s(2) * self.distance.clone());
        let y = if sheet == 0 {
            s(0)
        } else {
            root_real(&(rho.square() - x.square())) * s(sheet)
        };
        let z = &self.coordinates.minor * &sine;
        self.coordinates
            .origin
            .plus(&self.coordinates.e1.scaled(&x))
            .plus(&self.coordinates.e2.scaled(&y))
            .plus(&self.coordinates.e3.scaled(&z))
    }
}

fn angle(quarter: usize, parameter: &Real) -> [Real; 2] {
    let square = parameter.square();
    let denominator = s(1) + square.clone();
    let cosine = (s(1) - square) / denominator.clone();
    let sine = s(2) * parameter.clone() / denominator;
    match quarter {
        0 => [cosine, sine],
        1 => [-sine, cosine],
        2 => [-cosine, -sine],
        3 => [sine, -cosine],
        _ => unreachable!("minor-angle quarter"),
    }
}
fn coefficients(constant: &Real, cosine: &Real, sine: &Real, quarter: usize) -> [Real; 3] {
    let (cosine, sine) = match quarter {
        0 => (cosine.clone(), sine.clone()),
        1 => (sine.clone(), -cosine),
        2 => (-cosine, -sine),
        3 => (-sine, cosine.clone()),
        _ => unreachable!("minor-angle quarter"),
    };
    [constant + &cosine, s(2) * sine, constant - &cosine]
}
fn value(polynomial: &[Real; 3], at: &Real) -> Real {
    &polynomial[0] + &(at * &(&polynomial[1] + &(at * &polynomial[2])))
}
#[derive(Clone)]
struct Event {
    at: Real,
    roots: [usize; 4],
}
fn insert_event(events: &mut Vec<Event>, event: Event, dec: &mut Decider) -> Result<(), Refusal> {
    for (index, existing) in events.iter_mut().enumerate() {
        match compare(&event.at, &existing.at, dec)? {
            Sign::Negative => {
                events.insert(index, event);
                return Ok(());
            }
            Sign::Zero => {
                for (existing, incoming) in existing.roots.iter_mut().zip(event.roots) {
                    *existing = (*existing).max(incoming);
                }
                return Ok(());
            }
            Sign::Positive => {}
        }
    }
    events.push(event);
    Ok(())
}
fn add_root(
    events: &mut Vec<Event>,
    index: usize,
    at: Real,
    multiplicity: usize,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if sign(&at, dec)? != Sign::Negative && compare(&at, &s(1), dec)? != Sign::Positive {
        let mut roots = [0; 4];
        roots[index] = multiplicity;
        insert_event(events, Event { at, roots }, dec)?;
    }
    Ok(())
}
fn events(
    polynomials: &[[Real; 3]; 4],
    dec: &mut Decider,
) -> Result<(Vec<Event>, [usize; 4]), Refusal> {
    let mut events = Vec::with_capacity(12);
    events.push(Event {
        at: s(0),
        roots: [0; 4],
    });
    events.push(Event {
        at: s(1),
        roots: [0; 4],
    });
    let mut degrees = [0; 4];
    for (index, polynomial) in polynomials.iter().enumerate() {
        let [constant, linear, quadratic] = polynomial;
        let quadratic_sign = sign(quadratic, dec)?;
        if quadratic_sign == Sign::Zero {
            if sign(linear, dec)? != Sign::Zero {
                degrees[index] = 1;
                if index < 3 {
                    add_root(&mut events, index, &(-constant) / linear, 1, dec)?;
                }
            }
            continue;
        }
        degrees[index] = 2;
        // The fourth factor is rho+a+d = (rho+a-d)+2d, strictly positive
        // throughout the retained domain; its roots cannot bound a component.
        if index < 3 {
            let discriminant = linear.square() - s(4) * quadratic.clone() * constant.clone();
            match sign(&discriminant, dec)? {
                Sign::Negative => {}
                Sign::Zero => add_root(
                    &mut events,
                    index,
                    -linear / (s(2) * quadratic.clone()),
                    2,
                    dec,
                )?,
                Sign::Positive => {
                    let root = root_real(&discriminant);
                    for side in [1, -1] {
                        add_root(
                            &mut events,
                            index,
                            (-linear + s(side) * root.clone()) / (s(2) * quadratic.clone()),
                            1,
                            dec,
                        )?;
                    }
                }
            }
        }
        // A positive quadratic split at its vertex has positive Bernstein
        // coefficients on each retained subinterval. This keeps interval jets
        // finite even when direct power-basis evaluation would lose positivity.
        if quadratic_sign == Sign::Positive {
            add_root(
                &mut events,
                index,
                -linear / (s(2) * quadratic.clone()),
                0,
                dec,
            )?;
        }
    }
    Ok((events, degrees))
}
fn factor(
    coefficients: &[Real; 3],
    lo: &Event,
    hi: &Event,
    index: usize,
    mut degree: usize,
) -> Polynomial {
    let mut coefficients = coefficients.clone();
    for endpoint in [lo, hi] {
        for _ in 0..endpoint.roots[index] {
            match degree {
                2 => {
                    coefficients[0] = &coefficients[1] + &(&endpoint.at * &coefficients[2]);
                    coefficients[1] = coefficients[2].clone();
                    coefficients[2] = s(0);
                }
                1 => {
                    coefficients[0] = coefficients[1].clone();
                    coefficients[1] = s(0);
                }
                _ => unreachable!("endpoint multiplicity cannot exceed polynomial degree"),
            }
            degree -= 1;
        }
    }
    if hi.roots[index] % 2 == 1 {
        for coefficient in &mut coefficients[..=degree] {
            *coefficient = -&*coefficient;
        }
    }
    let first = value(&coefficients, &lo.at);
    let bernstein = match degree {
        0 => [first, s(0), s(0)],
        1 => [first, value(&coefficients, &hi.at), s(0)],
        2 => {
            let derivative = &coefficients[1] + &(&(&s(2) * &coefficients[2]) * &lo.at);
            let middle = &first + &((&hi.at - &lo.at) * derivative / s(2));
            [first, middle, value(&coefficients, &hi.at)]
        }
        _ => unreachable!("quadratic factor"),
    };
    Polynomial { degree, bernstein }
}

struct Boundary {
    quarter: usize,
    parameter: Real,
    order: usize,
    valid: bool,
    point: Option<usize>,
}
struct Span {
    start: usize,
    end: usize,
    sector: Option<usize>,
}

pub(super) fn parallel(
    cone: &Cone,
    torus: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !vector_zero(&cone.axis.cross(&torus.axis), dec)? {
        return Err(Refusal::new(
            "intersect.traced-pair-required",
            "skew cone/torus intersection requires certified branch tracing",
        ));
    }
    let e3 = unit(&torus.axis);
    let delta = cone.apex.minus(&torus.center);
    let axial = delta.dot(&e3);
    let offset = delta.minus(&e3.scaled(&axial));
    let distance = root_real(&offset.norm_squared());
    let e1 = offset.scaled(&(s(1) / distance.clone()));
    let e2 = e3.cross(&e1);
    let cone_sense = if sign(&cone.axis.dot(&e3), dec)? == Sign::Positive {
        s(1)
    } else {
        s(-1)
    };
    let slope = &cone.slope * &cone_sense;
    let radius_base = -(&slope * &axial);
    let radius_sine = &slope * &torus.minor;
    let coordinates = Arc::new(Coordinates {
        origin: torus.center.clone(),
        e1,
        e2,
        e3,
        major: torus.major.clone(),
        minor: torus.minor.clone(),
        inverse_minor: s(1) / torus.minor.clone(),
        inverse_twice_radii: s(1) / (s(2) * torus.major.clone() * torus.minor.clone()),
        squared_radii: torus.major.square() + torus.minor.square(),
    });
    // Four quarter seams, at most two roots per inequality, and one quadratic
    // vertex per factor per quarter. These are reserve sizes, not work limits.
    const BOUNDARIES: usize = 4 + 2 * 3 + 4 * 4;
    let mut frame = ConeTorusFrame {
        coordinates,
        distance,
        radius_base,
        radius_sine,
        sectors: Vec::with_capacity(BOUNDARIES),
        domains: Vec::with_capacity(6),
        sense: if sign(&sigma, dec)? == Sign::Positive {
            1
        } else {
            -1
        },
    };
    let mut boundaries: Vec<Boundary> = Vec::with_capacity(BOUNDARIES);
    let mut spans = Vec::with_capacity(BOUNDARIES);
    for quarter in 0..4 {
        let polynomials = [
            coefficients(
                &(&torus.major + &frame.radius_base - frame.distance.clone()),
                &torus.minor,
                &frame.radius_sine,
                quarter,
            ),
            coefficients(
                &(&torus.major - &frame.radius_base + frame.distance.clone()),
                &torus.minor,
                &(-&frame.radius_sine),
                quarter,
            ),
            coefficients(
                &(&frame.radius_base - &torus.major + frame.distance.clone()),
                &(-&torus.minor),
                &frame.radius_sine,
                quarter,
            ),
            coefficients(
                &(&torus.major + &frame.radius_base + frame.distance.clone()),
                &torus.minor,
                &frame.radius_sine,
                quarter,
            ),
        ];
        let (events, degrees) = events(&polynomials, dec)?;
        let mut indices = Vec::with_capacity(events.len());
        for (index, event) in events.iter().enumerate() {
            let boundary_index = if index == 0 && quarter != 0 {
                boundaries.len() - 1
            } else if quarter == 3 && index + 1 == events.len() {
                0
            } else {
                let mut valid = true;
                for (index, polynomial) in polynomials[..3].iter().enumerate() {
                    if event.roots[index] == 0
                        && sign(&value(polynomial, &event.at), dec)? == Sign::Negative
                    {
                        valid = false;
                        break;
                    }
                }
                boundaries.push(Boundary {
                    quarter,
                    parameter: event.at.clone(),
                    order: event.roots.iter().sum(),
                    valid,
                    point: None,
                });
                boundaries.len() - 1
            };
            boundaries[boundary_index].order = boundaries[boundary_index]
                .order
                .max(event.roots.iter().sum());
            indices.push(boundary_index);
        }
        for (index, endpoints) in events.windows(2).enumerate() {
            let [lo, hi] = endpoints else { unreachable!() };
            let middle = (&lo.at + &hi.at) / s(2);
            let mut valid = true;
            for polynomial in &polynomials[..3] {
                if sign(&value(polynomial, &middle), dec)? != Sign::Positive {
                    valid = false;
                    break;
                }
            }
            let sector = if valid {
                let factors = [
                    factor(&polynomials[0], lo, hi, 0, degrees[0]),
                    factor(&polynomials[1], lo, hi, 1, degrees[1]),
                    factor(&polynomials[2], lo, hi, 2, degrees[2]),
                    factor(&polynomials[3], lo, hi, 3, degrees[3]),
                ];
                frame.sectors.push(Sector {
                    quarter,
                    lo: lo.at.clone(),
                    span: &hi.at - &lo.at,
                    orders: [lo.roots.iter().sum(), hi.roots.iter().sum()],
                    factors,
                });
                Some(frame.sectors.len() - 1)
            } else {
                None
            };
            spans.push(Span {
                start: indices[index],
                end: indices[index + 1],
                sector,
            });
        }
    }
    let mut points = Vec::with_capacity(6);
    for (index, boundary) in boundaries.iter_mut().enumerate() {
        if boundary.valid && boundary.order != 0 {
            let adjacent = spans
                .iter()
                .any(|span| span.sector.is_some() && (span.start == index || span.end == index));
            if boundary.order > 1 || !adjacent {
                boundary.point = Some(points.len());
                points.push(frame.position(boundary.quarter, &boundary.parameter, 0));
            }
        }
    }
    let all_valid = spans.iter().all(|span| span.sector.is_some());
    let roots = boundaries
        .iter()
        .filter(|boundary| boundary.valid && boundary.order != 0)
        .count();
    if all_valid && roots <= 1 {
        for sheet in [1, -1] {
            frame.domains.push(Domain {
                start: 0,
                count: frame.sectors.len(),
                sheet,
            });
        }
    } else if !frame.sectors.is_empty() {
        let first = spans
            .iter()
            .position(|span| span.sector.is_some() && boundaries[span.start].order != 0)
            .expect("a proper retained interval starts at a boundary root");
        let mut domain: Option<Domain> = None;
        for offset in 0..spans.len() {
            let span = &spans[(first + offset) % spans.len()];
            if let Some(sector) = span.sector {
                let current = domain.get_or_insert(Domain {
                    start: sector,
                    count: 0,
                    sheet: 0,
                });
                current.count += 1;
                if boundaries[span.end].order != 0 {
                    frame.domains.push(domain.take().unwrap());
                }
            } else {
                assert!(domain.is_none(), "a retained interval ends at a root");
            }
        }
        assert!(domain.is_none(), "cyclic interval partition closes");
    }
    let mut curves = Vec::with_capacity(frame.domains.len());
    for domain in &frame.domains {
        let first = &frame.sectors[domain.start];
        let witness = frame.position(
            first.quarter,
            &(&first.lo + &(&first.span / &s(2))),
            if domain.sheet == 0 { 1 } else { domain.sheet },
        );
        let mut singular = Vec::new();
        for span in &spans {
            let Some(sector) = span.sector else { continue };
            if (sector + frame.sectors.len() - domain.start) % frame.sectors.len() >= domain.count {
                continue;
            }
            for endpoint in [span.start, span.end] {
                if let Some(point) = boundaries[endpoint].point
                    && !singular.contains(&point)
                {
                    singular.push(point);
                }
            }
        }
        curves.push(Spec {
            form: CurveForm::Branch,
            witness: Some(witness),
            contact: Contact::Crossing,
            singular,
        });
    }
    if curves.is_empty() && points.is_empty() {
        return Ok(Plan::disjoint());
    }
    Ok(Plan {
        terminal: None,
        points,
        frame: (!curves.is_empty()).then_some(Frame::ConeTorus(frame)),
        curves,
        divider: None,
    })
}
