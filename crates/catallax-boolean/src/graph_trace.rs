//! Temporary firstmate-authorized attribution of one unchanged tray run.

use std::{cell::RefCell, rc::Rc};
use catallax_numbers::{Real, current_operation_work};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::{FaceId, Model};

const CLASSES: [&str; 4] = ["plane-plane", "plane-cylinder", "cylinder-cylinder", "with-sphere"];

#[derive(Clone)]
pub(crate) struct Trace(Rc<RefCell<Ledger>>);
struct Ledger {
    start: u64,
    phase_start: u64,
    phase: &'static str,
    phases: Vec<(&'static str, u64)>,
    pairs: Vec<Pair>,
}
struct Pair {
    faces: [FaceId; 2],
    labels: [String; 2],
    class: usize,
    rejected: bool,
    late_rejected: bool,
    reached: bool,
    calls: u64,
    work: [u64; 4],
}
pub(crate) struct Charge {
    trace: Trace,
    pub(crate) index: usize,
    category: usize,
    start: u64,
}
fn number(x: &Real) -> String {
    x.rational_ref().map(ToString::to_string).unwrap_or_else(|| "non-rational-definition".into())
}
fn point(p: &Point3) -> [String; 3] { [&p.x, &p.y, &p.z].map(number) }
fn vector(p: &Vector3) -> [String; 3] { [&p.x, &p.y, &p.z].map(number) }
fn label(model: &Model, face: FaceId) -> String {
    let shape = match model.geom().surface(model.face(face).surface) {
        Surface::Plane(p) => format!("plane origin={:?} normal={:?}", point(&p.origin), vector(&p.normal)),
        Surface::Cylinder(c) => format!("cylinder origin={:?} axis={:?} radius_sq={}", point(&c.origin), vector(&c.axis), number(&c.radius_sq)),
        Surface::Sphere(s) => format!("sphere center={:?} radius={}", point(&s.center), number(&s.radius)),
        _ => unreachable!("the experiment has only planar, cylinder and sphere carriers"),
    };
    format!("{shape} origins={:?}", model.face_origins(face))
}
impl Trace {
    pub(crate) fn new() -> Self {
        let start = current_operation_work().used;
        Self(Rc::new(RefCell::new(Ledger {start, phase_start:start, phase:"setup", phases:Vec::new(), pairs:Vec::new()})))
    }
    pub(crate) fn checkpoint(&self, next: &'static str) {
        let mut ledger = self.0.borrow_mut();
        let now = current_operation_work().used;
        let phase = ledger.phase;
        let work = now - ledger.phase_start;
        ledger.phases.push((phase, work));
        ledger.phase = next;
        ledger.phase_start = now;
    }
    pub(crate) fn pair(&self, model: &Model, a: FaceId, b: FaceId) -> Charge {
        let sa = model.geom().surface(model.face(a).surface);
        let sb = model.geom().surface(model.face(b).surface);
        let class = match (sa, sb) {
            (Surface::Sphere(_), _) | (_, Surface::Sphere(_)) => 3,
            (Surface::Plane(_), Surface::Plane(_)) => 0,
            (Surface::Cylinder(_), Surface::Cylinder(_)) => 2,
            (Surface::Plane(_), Surface::Cylinder(_)) | (Surface::Cylinder(_), Surface::Plane(_)) => 1,
            _ => unreachable!("only the authorized planar experiment is traced"),
        };
        let mut ledger = self.0.borrow_mut();
        let index = ledger.pairs.len();
        ledger.pairs.push(Pair {faces:[a,b], labels:[label(model,a),label(model,b)],class,rejected:false,late_rejected:false,reached:false,calls:0,work:[0;4]});
        drop(ledger);
        self.charge(index, 0)
    }
    pub(crate) fn charge(&self, index: usize, category: usize) -> Charge {
        Charge {trace:self.clone(),index,category,start:current_operation_work().used}
    }
}
impl Charge {
    pub(crate) fn rejected(&self) { self.trace.0.borrow_mut().pairs[self.index].rejected=true; }
    pub(crate) fn late_rejected(&self) { self.trace.0.borrow_mut().pairs[self.index].late_rejected=true; }
    pub(crate) fn reached(&self) { self.trace.0.borrow_mut().pairs[self.index].reached=true; }
    pub(crate) fn called(&self) { self.trace.0.borrow_mut().pairs[self.index].calls+=1; }
}
impl Drop for Charge {
    fn drop(&mut self) {
        self.trace.0.borrow_mut().pairs[self.index].work[self.category] += current_operation_work().used-self.start;
    }
}
impl Drop for Ledger {
    fn drop(&mut self) {
        let now=current_operation_work().used;
        self.phases.push((self.phase,now-self.phase_start));
        for (phase,work) in &self.phases { eprintln!("PAIR_TRACE_PHASE phase={phase} work={work}"); }
        let mut owned=0u64;
        for (class,name) in CLASSES.iter().enumerate() {
            let rows:Vec<_>=self.pairs.iter().filter(|p|p.class==class).collect();
            let direct:u64=rows.iter().map(|p|p.work[0]).sum();
            let downstream:u64=rows.iter().map(|p|p.work[1..].iter().sum::<u64>()).sum();
            owned+=direct+downstream;
            eprintln!("PAIR_TRACE_CLASS class={name} considered={} rejected={} late_rejected={} reached={} calls={} direct={direct} downstream={downstream} total={}",rows.len(),rows.iter().filter(|p|p.rejected).count(),rows.iter().filter(|p|p.late_rejected).count(),rows.iter().filter(|p|p.reached).count(),rows.iter().map(|p|p.calls).sum::<u64>(),direct+downstream);
        }
        let total=now-self.start;
        eprintln!("PAIR_TRACE_TOTAL work={total} pair_owned={owned} shared={}",total-owned);
        let mut ranked:Vec<_>=self.pairs.iter().collect();
        ranked.sort_by_key(|p| std::cmp::Reverse(p.work.iter().sum::<u64>()));
        for (rank,p) in ranked.into_iter().take(5).enumerate() {
            eprintln!("PAIR_TRACE_HEAVY rank={} class={} faces={:?} direct={} hits={} curve_pairs={} trim={} total={} first={:?} second={:?}",rank+1,CLASSES[p.class],p.faces,p.work[0],p.work[1],p.work[2],p.work[3],p.work.iter().sum::<u64>(),p.labels[0],p.labels[1]);
        }
    }
}
