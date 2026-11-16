//! STEP geometry declarations, before unit conversion or admission.
use catallax_numbers::Rational;

/// Original Part 21 instance number, without the leading `#`.
pub type EntityId = u64;

/// A SELECT value used to trim a source curve.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum TrimSelect {
    /// A Cartesian point instance.
    Point(EntityId),
    /// A parameter in the basis curve's STEP parameterization.
    Parameter(Rational),
}

/// STEP's declared preference when both trim representations are supplied.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TrimmingPreference {
    /// Use Cartesian trim points first.
    Cartesian,
    /// Use numeric trim parameters first.
    Parameter,
    /// No representation is preferred.
    Unspecified,
}

/// Source knot encoding, retained without assuming a parameter remapping.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum KnotData {
    /// Explicit knot multiplicities, values and knot specification.
    Explicit(Vec<usize>, Vec<Rational>, String),
    /// Bezier end multiplicities.
    Bezier,
    /// Uniform knots, including the exterior knots.
    Uniform,
    /// Clamped uniform knots.
    QuasiUniform,
}

/// A source B-spline curve, including its declarations and control references.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SplineCurve {
    /// Written degree.
    pub degree: usize,
    /// Cartesian point references in source order; they may be 2D or 3D.
    pub controls: Vec<EntityId>,
    /// Written weights, or implicit unit weights.
    pub weights: Option<Vec<Rational>>,
    /// Knot encoding.
    pub knots: KnotData,
    /// Declared curve form, not a geometric certificate.
    pub form: String,
    /// Declared closure; `None` means `.U.`.
    pub closed: Option<bool>,
    /// Declared self-intersection; `None` means `.U.`.
    pub self_intersect: Option<bool>,
}

/// A source B-spline surface with its u-major control grid.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SplineSurface {
    /// Written u and v degrees.
    pub degree: [usize; 2],
    /// Cartesian control point references.
    pub controls: Vec<Vec<EntityId>>,
    /// Written weight grid, or implicit unit weights.
    pub weights: Option<Vec<Vec<Rational>>>,
    /// U and v knot encodings.
    pub knots: [KnotData; 2],
    /// Declared surface form.
    pub form: String,
    /// Declared u and v closure.
    pub closed: [Option<bool>; 2],
    /// Declared self-intersection.
    pub self_intersect: Option<bool>,
}

/// A Cartesian transformation as written, including omitted axes and scales.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Transformation {
    /// Declared Cartesian dimension; absent on the unspecialized base operator.
    pub dimension: Option<usize>,
    /// Local Cartesian origin.
    pub origin: EntityId,
    /// Optional source directions; the third is absent for 2D.
    pub axes: [Option<EntityId>; 3],
    /// Uniform scale followed by optional nonuniform second and third scales.
    pub scales: [Option<Rational>; 3],
}

/// A source curve. References preserve wrappers and their parameter conventions.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Curve {
    /// Origin point and VECTOR instance, including its magnitude.
    Line(EntityId, EntityId),
    /// Placement and radius.
    Circle(EntityId, Rational),
    /// Placement and the two written semi-axes, without reordering them.
    Ellipse(EntityId, Rational, Rational),
    /// Placement, real semi-axis and imaginary semi-axis.
    Hyperbola(EntityId, Rational, Rational),
    /// Placement and focal distance.
    Parabola(EntityId, Rational),
    /// Finite polynomial or rational spline source data.
    BSpline(Box<SplineCurve>),
    /// Ordered Cartesian points; source parameter knots are point indices.
    Polyline(Vec<EntityId>),
    /// Basis, trim SELECT lists, sense agreement and representation preference.
    Trimmed(EntityId, [Vec<TrimSelect>; 2], bool, TrimmingPreference),
    /// Parent curve and Cartesian transformation.
    Replica(EntityId, EntityId),
    /// Basis, signed distance, self-intersection declaration and reference direction.
    Offset3d(EntityId, Rational, Option<bool>, EntityId),
    /// 3D curve, associated geometry, representation preference and seam declaration.
    SurfaceCurve(EntityId, Vec<EntityId>, String, bool),
    /// Basis surface and definitional representation of the 2D curve.
    Pcurve(EntityId, EntityId),
    /// Ordered composite segments and self-intersection declaration.
    Composite(Vec<EntityId>, Option<bool>),
}

/// A source surface, preserving its STEP parameterization and wrappers.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Surface {
    /// Plane placement.
    Plane(EntityId),
    /// Placement and radius.
    Cylinder(EntityId, Rational),
    /// Placement, reference radius and signed semi-angle.
    Cone(EntityId, Rational, Rational),
    /// Placement and radius.
    Sphere(EntityId, Rational),
    /// Placement and the written major and minor radii.
    Torus(EntityId, Rational, Rational),
    /// Polynomial or rational spline source data.
    BSpline(Box<SplineSurface>),
    /// Swept curve and extrusion VECTOR instance.
    Extrusion(EntityId, EntityId),
    /// Swept curve and AXIS1_PLACEMENT instance.
    Revolution(EntityId, EntityId),
    /// Basis, [u1, u2, v1, v2] limits and u/v sense flags.
    RectangularTrim(EntityId, [Rational; 4], [bool; 2]),
    /// Parent surface and Cartesian transformation.
    Replica(EntityId, EntityId),
    /// Basis, signed distance and self-intersection declaration.
    Offset(EntityId, Rational, Option<bool>),
    /// Basis, ordered boundary curves and implicit-outer declaration.
    CurveBounded(EntityId, Vec<EntityId>, bool),
    /// U-major grid of SURFACE_PATCH references.
    RectangularComposite(Vec<Vec<EntityId>>),
}

/// One typed geometry instance. All scalar values retain their source units.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Geometry {
    /// One, two or three Cartesian coordinates.
    Point(Vec<Rational>),
    /// Two or three direction ratios.
    Direction(Vec<Rational>),
    /// Direction reference and written magnitude.
    Vector(EntityId, Rational),
    /// Origin and optional axis direction.
    Axis1Placement(EntityId, Option<EntityId>),
    /// Origin and optional 2D reference direction.
    Axis2Placement2d(EntityId, Option<EntityId>),
    /// Origin, optional axis and optional reference direction.
    Axis2Placement3d(EntityId, Option<EntityId>, Option<EntityId>),
    /// Cartesian transformation definition.
    Transformation(Transformation),
    /// Curve definition.
    Curve(Curve),
    /// Surface definition.
    Surface(Surface),
    /// Transition code, sense agreement and parent curve.
    CompositeSegment(String, bool, EntityId),
    /// Parent surface, u/v transition codes and u/v senses.
    SurfacePatch(EntityId, [String; 2], [bool; 2]),
    /// Items and representation context for an associated parameter curve.
    DefinitionalRepresentation(Vec<EntityId>, EntityId),
}

/// A representation's source items and context, without placement or admission.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Representation {
    /// Original item references, in source order.
    pub items: Vec<EntityId>,
    /// Representation context carrying dimensions and units.
    pub context: EntityId,
}

/// The dimensional and unit declarations of one representation context.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RepresentationContext {
    /// Declared Cartesian dimension, when present.
    pub dimension: Option<usize>,
    /// Unit instance references, without inventing a unit for missing data.
    pub units: Vec<EntityId>,
    /// Declared uncertainty references.
    pub uncertainties: Vec<EntityId>,
}

/// Topology declarations. Shell closure and orientation remain source claims.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Topology {
    /// Outer connected face set and oriented inner-shell references.
    Solid {
        /// Outer shell or connected face set.
        outer: EntityId,
        /// Oriented inner-shell instances in source order.
        voids: Vec<EntityId>,
    },
    /// Shell-based surface model, retaining open or closed shell references.
    SurfaceModel(Vec<EntityId>),
    /// Face-based surface model, retaining connected face sets.
    FaceModel(Vec<EntityId>),
    /// Face references and closure declaration; absent for a connected face set.
    Shell {
        /// Source face instances.
        faces: Vec<EntityId>,
        /// Declared closure, absent for an unspecialized face set.
        closed: Option<bool>,
    },
    /// Oriented shell reference and its declared orientation.
    OrientedShell(EntityId, bool),
    /// Bounds, optional support geometry and declared same-sense flag.
    Face {
        /// Face-bound instances in source order.
        bounds: Vec<EntityId>,
        /// Support surface, absent on a topological FACE.
        surface: Option<EntityId>,
        /// Declared sense against the support geometry.
        same_sense: Option<bool>,
    },
    /// Oriented face reference and its declared orientation.
    OrientedFace(EntityId, bool),
    /// Boundary loop, orientation and explicit outer-bound designation.
    FaceBound {
        /// Boundary LOOP instance.
        loop_id: EntityId,
        /// Declared traversal orientation.
        orientation: bool,
        /// Whether the entity explicitly names an outer bound.
        outer: bool,
    },
    /// Ordered oriented-edge references.
    EdgeLoop(Vec<EntityId>),
    /// Ordered Cartesian points of a polygonal loop.
    PolyLoop(Vec<EntityId>),
    /// A degenerate loop represented by one vertex.
    VertexLoop(EntityId),
    /// Endpoint vertices, optional curve and same-sense declaration.
    Edge {
        /// Source start and end vertices.
        vertices: [EntityId; 2],
        /// Support curve, absent on a topological EDGE.
        curve: Option<EntityId>,
        /// Declared sense against the support geometry.
        same_sense: Option<bool>,
    },
    /// Oriented edge reference and its declared orientation.
    OrientedEdge(EntityId, bool),
    /// A vertex's Cartesian point reference.
    Vertex(EntityId),
}

/// Representation placement records, before composing transformations.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PlacementLink {
    /// Mapping origin and mapped representation.
    RepresentationMap(EntityId, EntityId),
    /// Representation map and target transformation.
    MappedItem(EntityId, EntityId),
    /// Two representation references and their transformation operator.
    Relationship([EntityId; 2], EntityId),
    /// Source and target representation items defining a transformation.
    ItemDefinedTransformation([EntityId; 2]),
}
