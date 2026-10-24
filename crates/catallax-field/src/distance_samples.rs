//! Query-only candidates for distance features. No part-construction path calls this module.

use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_recipe::{BooleanKind, Origin};
use catallax_shapes::{Point3, RPoint3, RVector3, Surface, Vector3};

use crate::{Instruction, Plan, geometry::Leaf, samples::{Bounds, circle, fraction, placement_error}};

fn real_bound(value: &catallax_numbers::Real, dec: &mut Decider) -> Result<[Rational;2],Refusal> {
    if let Some(q) = value.rational_ref() { return Ok([q.clone(),q.clone()]); }
    let b = dec.enclose_budgeted_at(value,128)?;
    Ok([b.lo(),b.hi()])
}
fn points_bound(points: &[Point3], dec: &mut Decider) -> Result<Option<Bounds>,Refusal> {
    let mut bounds: Option<Bounds> = None;
    for p in points {
        let coordinate = [real_bound(&p.x,dec)?,real_bound(&p.y,dec)?,real_bound(&p.z,dec)?];
        if let Some(bounds) = &mut bounds {
            for i in 0..3 { bounds[i][0] = bounds[i][0].clone().min(coordinate[i][0].clone()); bounds[i][1] = bounds[i][1].clone().max(coordinate[i][1].clone()); }
        } else { bounds = Some(coordinate); }
    }
    Ok(bounds)
}
fn leaf_bound(leaf: &Leaf, dec: &mut Decider) -> Result<Option<Bounds>,Refusal> {
    match leaf {
        Leaf::Primitive(p) => return p.sample_bounds().map(Some),
        Leaf::Feature(f) => return f.sample_bounds(),
        Leaf::Cell(_) => (),
    }
    let constraints = leaf.constraints();
    if constraints.iter().all(|c|matches!(c.surface(),Surface::Plane(_))) {
        let faces = constraints.iter().map(|c| {
            let Surface::Plane(plane) = c.surface() else { unreachable!() };
            catallax_features::polyhedron::PlaneFace {plane:plane.clone(),origins:Vec::new()}
        }).collect();
        let cell = catallax_features::polyhedron::ConvexCell::new(faces,dec)?;
        return points_bound(cell.vertices(),dec);
    }
    let surfaces: Vec<_> = constraints.iter().map(|c|(c.surface().clone(),Vec::new())).collect();
    if let Some(skin) = catallax_features::cylinder_skin::CylinderSkin::new(&surfaces,dec)? {
        let axis = &skin.cylinder.axis;
        let norm = axis.norm_squared();
        let radius = skin.cylinder.radius();
        let components = [&axis.x,&axis.y,&axis.z];
        let mut result = std::array::from_fn(|_|[Rational::zero(),Rational::zero()]);
        for i in 0..3 {
            let radial = &radius * &(&catallax_numbers::Real::one() - &(&components[i].square()/&norm)).sqrt();
            let radial = real_bound(&radial,dec)?[1].clone();
            let at = |p: &Point3| [&p.x,&p.y,&p.z][i].clone();
            let a = real_bound(&at(&skin.ends[0]),dec)?;
            let b = real_bound(&at(&skin.ends[1]),dec)?;
            result[i] = [&a[0].clone().min(b[0].clone()) - &radial,&a[1].clone().max(b[1].clone()) + &radial];
        }
        return Ok(Some(result));
    }
    if let [constraint] = constraints {
        if let Surface::Sphere(s) = constraint.surface() {
            let radius = real_bound(&s.radius,dec)?[1].clone();
            let c = points_bound(std::slice::from_ref(&s.center),dec)?.unwrap();
            return Ok(Some(c.map(|v|[&v[0]-&radius,&v[1]+&radius])));
        }
    }
    let Leaf::Cell(cell) = leaf else { unreachable!() };
    cell.sample_bounds()
}
fn combine(kind: BooleanKind, a: &Option<Bounds>, b: &Option<Bounds>) -> Option<Bounds> {
    match kind {
        BooleanKind::Subtract => a.clone(),
        BooleanKind::Union => match (a,b) {
            (Some(a),Some(b)) => Some(std::array::from_fn(|i|[a[i][0].clone().min(b[i][0].clone()),a[i][1].clone().max(b[i][1].clone())])),
            _ => None,
        },
        BooleanKind::Intersect => match (a,b) {
            (Some(a),Some(b)) => Some(std::array::from_fn(|i|[a[i][0].clone().max(b[i][0].clone()),a[i][1].clone().min(b[i][1].clone())])),
            (a,b) => a.clone().or_else(||b.clone()),
        },
    }
}
fn rational_vector(v: &Vector3) -> Result<RVector3,Refusal> {
    let Some([x,y,z]) = v.rational_components() else { return Err(placement_error()); };
    Ok(RVector3::new(x.clone(),y.clone(),z.clone()))
}
fn frame(axis: &Vector3, reference: &Vector3) -> Result<[RVector3;3],Refusal> {
    let axis = rational_vector(axis)?;
    let reference = rational_vector(reference)?;
    let a = crate::profile::sqrt(&axis.norm_squared()).ok_or_else(placement_error)?;
    let u = crate::profile::sqrt(&reference.norm_squared()).ok_or_else(placement_error)?;
    let axis = axis.scaled(&a.recip().unwrap());
    let reference = reference.scaled(&u.recip().unwrap());
    Ok([reference.clone(),axis.cross(&reference),axis])
}
fn within(p: &RPoint3, bounds: &Bounds) -> bool {
    [&p.x,&p.y,&p.z].iter().enumerate().all(|(i,q)| *q >= &bounds[i][0] && *q <= &bounds[i][1])
}
fn roots(a: &Rational,b: &Rational,c: &Rational,cuts: &mut Vec<Rational>) {
    let mut add = |q: Rational| { if q > Rational::zero() && q < Rational::one() { cuts.push(q); } };
    if a.is_zero() {
        if !b.is_zero() { add(-c / b.clone()); }
        return;
    }
    let discriminant = b*b - Rational::from_i64(4)*a.clone()*c.clone();
    if discriminant < Rational::zero() { return; }
    // These enclosures only choose rational candidates. They never accept a
    // geometric distance or substitute a rounded boundary for a definition.
    let radical = Ball::from_rational(&discriminant,64).sqrt().unwrap();
    let b = Ball::from_rational(b,64);
    let inverse = Ball::from_rational(&(Rational::from_i64(2)*a.clone()),64).recip().unwrap();
    for q in [(&(-&b) + &radical)*inverse.clone(),(-b-radical)*inverse] { add(q.lo());add(q.hi()); }
}
fn cylinder_candidates(surface: &catallax_shapes::Cylinder,bounds: &Bounds,n: u32,visit: &mut impl FnMut(RPoint3)->bool) -> Result<bool,Refusal> {
    let [u,v,axis] = frame(&surface.axis,&surface.ref_dir)?;
    let origin = surface.origin.as_rational().ok_or_else(placement_error)?;
    let radius = surface.radius().rational_ref().cloned().ok_or_else(placement_error)?;
    let mut axial = [Rational::zero(),Rational::zero()];
    for corner in 0..8 {
        let [x,y,z] = std::array::from_fn(|i|bounds[i][(corner>>i)&1].clone());
        let value = RPoint3::new(x,y,z).minus(&origin).dot(&axis);
        if corner == 0 { axial = [value.clone(),value]; }
        else { axial[0] = axial[0].clone().min(value.clone());axial[1] = axial[1].clone().max(value); }
    }
    for i in 0..n {
        let height = &axial[0]+&((&axial[1]-&axial[0])*fraction(i,n));
        let centre = origin.plus(&axis.scaled(&height));
        for quarter in 0..4 {
            let (u,v) = match quarter {
                0 => (u.clone(),v.clone()),
                1 => (v.scaled(&(-Rational::one())),u.clone()),
                2 => (u.scaled(&(-Rational::one())),v.scaled(&(-Rational::one()))),
                _ => (v.clone(),u.scaled(&(-Rational::one()))),
            };
            let mut cuts = vec![Rational::zero(),Rational::one()];
            for j in 1..n { cuts.push(Rational::from_ratio(j.into(),n.into()).unwrap()); }
            let c = [&centre.x,&centre.y,&centre.z];
            let us = [&u.x,&u.y,&u.z];
            let vs = [&v.x,&v.y,&v.z];
            for k in 0..3 {
                for bound in &bounds[k] {
                    let h = c[k]-bound;
                    let ru = &radius*us[k];
                    roots(&(&h-&ru),&(Rational::from_i64(2)*radius.clone()*vs[k].clone()),&(&h+&ru),&mut cuts);
                }
            }
            cuts.sort();cuts.dedup();
            for pair in cuts.windows(2) {
                let t = (&pair[0]+&pair[1])/Rational::from_i64(2);
                let squared = &t*&t;
                let denominator = Rational::one()+squared.clone();
                let cosine = (Rational::one()-squared)/denominator.clone();
                let sine = Rational::from_i64(2)*t/denominator;
                let point = centre.plus(&u.scaled(&(&radius*&cosine))).plus(&v.scaled(&(&radius*&sine)));
                if within(&point,bounds) && !visit(point) { return Ok(false); }
            }
        }
    }
    Ok(true)
}

impl Plan {
    /// Visits additional rational candidates for distance-feature results. Plane
    /// grids split at field-native leaf bounds; cylinders split their exact
    /// rational charts at box crossings. The caller must certify every trim,
    /// origin and normal. Sampling is not a whole-boundary certificate.
    /// This query-only path is separate from compilation and point evaluation.
    pub fn visit_distance_candidates(&self, divisions: u32, mut visit: impl FnMut(RPoint3,&Origin)->bool) -> Result<bool,Refusal> {
        if !(1..=64).contains(&divisions) { return Err(Refusal::new("field.sample-options","distance sampling divisions must be in 1..=64")); }
        if self.is_proven_empty() { return Ok(true); }
        let mut dec = Decider::default();
        let mut bounds: Vec<Option<Bounds>> = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            let b = match *instruction {
                Instruction::Empty => Some(std::array::from_fn(|_|[Rational::zero(),Rational::zero()])),
                Instruction::Leaf(i) => leaf_bound(&self.leaves[i],&mut dec)?,
                Instruction::Boolean {kind,a,b} => combine(kind,&bounds[a],&bounds[b]),
            };
            bounds.push(b);
        }
        let result = bounds[self.root].as_ref().ok_or_else(||Refusal::new("field.sample-unbounded","distance result has no finite candidate bound"))?;
        let mut cuts: [Vec<Rational>;3] = std::array::from_fn(|i|result[i].to_vec());
        for register in self.leaf_registers.iter().flatten() {
            if let Some(b) = &bounds[*register] {
                for i in 0..3 { for q in &b[i] { if q > &result[i][0] && q < &result[i][1] { cuts[i].push(q.clone()); } } }
            }
        }
        for i in 0..3 {
            for j in 1..divisions { cuts[i].push(&result[i][0]+&((&result[i][1]-&result[i][0])*Rational::from_ratio(j.into(),divisions.into()).unwrap())); }
            cuts[i].sort();cuts[i].dedup();
        }
        for face in &self.faces {
            if !face.in_result || self.leaf_registers[face.leaf].is_none() { continue; }
            let mut emit = |p|visit(p,&face.origin);
            match &face.surface {
                Surface::Plane(plane) => {
                    let origin = plane.origin.as_rational().ok_or_else(placement_error)?;
                    let normal = rational_vector(&plane.normal)?;
                    let ns = [&normal.x,&normal.y,&normal.z];
                    let axis = (0..3).max_by_key(|i|ns[*i].abs()).unwrap();
                    if ns[axis].is_zero() { return Err(placement_error()); }
                    let other = [(axis+1)%3,(axis+2)%3];
                    let h = origin.minus(&RPoint3::origin()).dot(&normal);
                    for a in cuts[other[0]].windows(2) {
                        for b in cuts[other[1]].windows(2) {
                            let mut p = std::array::from_fn(|_|Rational::zero());
                            p[other[0]] = (&a[0]+&a[1])/Rational::from_i64(2);
                            p[other[1]] = (&b[0]+&b[1])/Rational::from_i64(2);
                            let residual = &h-&(&(&p[other[0]]*ns[other[0]])+&(&p[other[1]]*ns[other[1]]));
                            p[axis] = &residual/ns[axis];
                            let [x,y,z] = p;
                            let p = RPoint3::new(x,y,z);
                            if within(&p,result) && !emit(p) { return Ok(false); }
                        }
                    }
                }
                Surface::Cylinder(cylinder) => { if !cylinder_candidates(cylinder,result,divisions,&mut emit)? { return Ok(false); } }
                Surface::Sphere(sphere) => {
                    let [u,v,axis] = frame(&sphere.axis,&sphere.ref_dir)?;
                    let origin = sphere.center.as_rational().ok_or_else(placement_error)?;
                    let radius = sphere.radius.rational_ref().ok_or_else(placement_error)?;
                    let circle = circle(divisions);
                    for longitude in &circle { for latitude in &circle {
                        let p = origin.plus(&u.scaled(&(radius*&(&longitude[0]*&latitude[0]))))
                            .plus(&v.scaled(&(radius*&(&longitude[1]*&latitude[0])))).plus(&axis.scaled(&(radius*&latitude[1])));
                        if within(&p,result) && !emit(p) { return Ok(false); }
                    } }
                }
                Surface::Torus(torus) => {
                    let [u,v,axis] = frame(&torus.axis,&torus.ref_dir)?;
                    let origin = torus.center.as_rational().ok_or_else(placement_error)?;
                    let major = torus.major.rational_ref().ok_or_else(placement_error)?;
                    let minor = torus.minor.rational_ref().ok_or_else(placement_error)?;
                    let circle = circle(divisions);
                    for longitude in &circle { for tube in &circle {
                        let radial = major+&(minor*&tube[0]);
                        let p = origin.plus(&u.scaled(&(&radial*&longitude[0])))
                            .plus(&v.scaled(&(&radial*&longitude[1]))).plus(&axis.scaled(&(minor*&tube[1])));
                        if within(&p,result) && !emit(p) { return Ok(false); }
                    } }
                }
                _ => (),
            }
        }
        Ok(true)
    }
}
