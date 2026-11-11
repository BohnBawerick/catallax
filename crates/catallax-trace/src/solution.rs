// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: translated from SMLib SmSolutionArray.h and
// SmSolutionArray.cpp (SmSolutionEnd, SmSolution, SmSolutionArray) and checked
// against rust-brep crates/kernel/src/solution_array.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, which ports the same routines for
// a frozen fixture domain of five native objects. Object and tree-node
// pointers become opaque identities that are compared and never followed; the
// object kind and the parameter domain that SMLib reads through `IsKindOf` and
// virtual getters come from a caller-supplied `SolutionObjects`. SMLib's
// borrowed arrays, capacity bookkeeping, assertion macros, `Dump`, `Draw`,
// `AssertValid`, `SelfTest` and `GetPoint` (which evaluates geometry that is
// ported in later steps) are not carried over.
//! Search solutions: the values, objects and tree nodes a solver reports.
//!
//! A [`Solution`] is an unproved candidate. Its parameters steer later
//! refinement; no Catallax proof reads its [`SolutionEnd::value`] as a bound.

use crate::Fault;
use crate::extent::Extent1;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQRT, are_same};

/// `SM_SA_MAX_TREES`: the most objects or tree nodes one solution names.
pub const MAX_OBJECTS: usize = 5;

/// `SM_SA_MAX_TREES * 2`: the parameter slots of one solution end.
pub const MAX_PARAMETERS: usize = MAX_OBJECTS * 2;

/// Identity of an `SmObject` a solution names. It is compared, never followed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ObjectId(pub u32);

/// Identity of a tree node a solution names, as an index into the node arena.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct NodeId(pub u32);

/// `SmSolutionType`. The discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolutionType {
    /// `SM_ST_NODES_ONLY`: the solution holds only tree nodes.
    NodesOnly = 0,
    /// `SM_ST_SINGLE_VALUE`: one point or parameter set.
    SingleValue = 1,
    /// `SM_ST_RANGE_OF_VALUES`: a range between the start and the end.
    RangeOfValues = 2,
    /// `SM_ST_POLY_SINGLE_VALUE`: a poly-solver point solution.
    PolySingleValue = 3,
    /// `SM_ST_UNKNOWN`: a solution not yet filled in.
    Unknown = 4,
}

impl SolutionType {
    /// The solution type with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        [
            Self::NodesOnly,
            Self::SingleValue,
            Self::RangeOfValues,
            Self::PolySingleValue,
            Self::Unknown,
        ]
        .get(code)
        .copied()
    }
}

/// `SmSolutionPairType`. The discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolutionPair {
    /// `SM_SP_IDENTICAL`.
    Identical = 0,
    /// `SM_SP_DIFFCOUNTS`: different type, object count, variable count or value.
    DifferentCounts = 1,
    /// `SM_SP_DIFFOBJECTS`: different objects.
    DifferentObjects = 2,
    /// `SM_SP_DIFFPARAMETERS`: only the parameters differ.
    DifferentParameters = 3,
}

/// Which of two compared solutions an identical pair prefers to keep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preferred {
    /// No preference (SMLib's `NULL`).
    Neither,
    /// The first solution compared.
    First,
    /// The second solution compared.
    Second,
}

/// `SmSortKeyType`: how `add_sorted_solution` orders answers. The
/// discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortKey {
    /// `SM_SK_BY_SOLUTION_VALUE`.
    ByValue = 0,
    /// `SM_SK_BY_FIRST_PARAMETER`.
    ByFirstParameter = 1,
    /// `SM_SK_BY_FIRST_TWO_PARAMETERS`.
    ByFirstTwoParameters = 2,
    /// `SM_SK_BY_VALUE_WITH_PARAMETER`.
    ByValueWithParameter = 3,
    /// `SM_SK_BY_ALL_PARAMETERS`.
    ByAllParameters = 4,
}

impl SortKey {
    /// The sort key with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        [
            Self::ByValue,
            Self::ByFirstParameter,
            Self::ByFirstTwoParameters,
            Self::ByValueWithParameter,
            Self::ByAllParameters,
        ]
        .get(code)
        .copied()
    }
}

/// The classes of `SmObject` a solution can name, as far as the solution
/// code distinguishes them (`IsKindOf`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObjectKind {
    /// `SmCurve`: one parameter.
    Curve,
    /// `SmEdge`: one parameter.
    Edge,
    /// `SmFace`: two parameters.
    Face,
    /// `SmSurface`: two parameters.
    Surface,
    /// `SmVolume`: three parameters.
    Volume,
    /// `SmVertex`, `SmPoly*` and `SmCPolyFace`: no parameters.
    Parameterless,
    /// Any other object: not a solution domain.
    Other,
}

impl ObjectKind {
    /// How many solution parameters objects of this kind take in
    /// `GetIndex`, which counts curves, edges, faces and surfaces only.
    fn index_width(self) -> usize {
        match self {
            Self::Curve | Self::Edge => 1,
            Self::Face | Self::Surface => 2,
            _ => 0,
        }
    }
}

/// What a solution needs to know about the objects it names.
///
/// SMLib reads these through virtual methods on the objects themselves; the
/// caller of the port supplies them here. An unknown identity is
/// [`ObjectKind::Other`].
pub trait SolutionObjects {
    /// The class of the object.
    fn kind(&self, id: ObjectId) -> ObjectKind;

    /// The parameter interval in one domain direction: `GetInterval` for an
    /// edge, `GetNaturalInterval` for a curve, and the `u`, `v` or `w`
    /// interval of the face, surface or volume domain for directions 0, 1, 2.
    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1;

    /// `SmTol::GetZoneTol3d`: the zone tolerance, zero for an object that
    /// has none.
    fn zone_tolerance(&self, id: ObjectId) -> f64;
}

/// A [`SolutionObjects`] that knows no object: every identity is
/// [`ObjectKind::Other`].
#[derive(Clone, Copy, Debug, Default)]
pub struct NoObjects;

impl SolutionObjects for NoObjects {
    fn kind(&self, _id: ObjectId) -> ObjectKind {
        ObjectKind::Other
    }

    fn parameter_interval(&self, _id: ObjectId, _direction: usize) -> Extent1 {
        Extent1::default()
    }

    fn zone_tolerance(&self, _id: ObjectId) -> f64 {
        0.0
    }
}

/// `SmSolutionEnd`: the value and parameters of one end of a solution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolutionEnd {
    /// `m_dSolutionValue`: typically a distance; the accuracy of the answer.
    pub value: f64,
    /// `m_adParameters`: the concatenated parameters of the solved objects.
    pub parameters: [f64; MAX_PARAMETERS],
}

impl Default for SolutionEnd {
    fn default() -> Self {
        Self {
            value: 0.0,
            parameters: [0.0; MAX_PARAMETERS],
        }
    }
}

impl SolutionEnd {
    /// `operator==`: the value and every parameter, used or not, agree
    /// within SMLib's scaled zero.
    pub fn is_same_as(&self, other: &Self) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        }
        are_same(self.value, other.value)
            && self
                .parameters
                .iter()
                .zip(&other.parameters)
                .all(|(a, b)| are_same(*a, *b))
    }
}

/// `SmSolution`: the answer to one search problem.
///
/// `variables` and `object_count` are capped at [`MAX_PARAMETERS`] and
/// [`MAX_OBJECTS`] wherever they index an array; SMLib leaves a larger
/// count undefined.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Solution {
    /// `m_eSolutionType`.
    pub kind: SolutionType,
    /// `m_lNumVariables`: how many parameters the ends carry.
    pub variables: usize,
    /// `m_vStart`: the point, or the start of a range.
    pub start: SolutionEnd,
    /// `m_vEnd`: the end of a range.
    pub end: SolutionEnd,
    /// `m_lNumObjects`: how many of `objects` are meaningful.
    pub object_count: usize,
    /// `m_apObjects`: the objects the solution came from.
    pub objects: [Option<ObjectId>; MAX_OBJECTS],
    /// `m_apNodes`: the tree nodes the solution came from.
    pub nodes: [Option<NodeId>; MAX_OBJECTS],
}

impl Default for Solution {
    /// The default constructor: type unknown, nothing set.
    fn default() -> Self {
        Self {
            kind: SolutionType::Unknown,
            variables: 0,
            start: SolutionEnd::default(),
            end: SolutionEnd::default(),
            object_count: 0,
            objects: [None; MAX_OBJECTS],
            nodes: [None; MAX_OBJECTS],
        }
    }
}

fn max3(a: f64, b: f64, c: f64) -> f64 {
    if a > b {
        if a > c { a } else { c }
    } else if b > c {
        b
    } else {
        c
    }
}

fn are_sequential(a: f64, b: f64, c: f64) -> bool {
    (a < b && b < c) || (a > b && b > c)
}

impl Solution {
    /// A solution of all-zero bytes, which is `SM_ST_NODES_ONLY`, as SMLib
    /// leaves the slots that `SetSize` and `InsertAt` open up.
    pub fn zeroed() -> Self {
        Self {
            kind: SolutionType::NodesOnly,
            ..Self::default()
        }
    }

    fn variable_count(&self) -> usize {
        self.variables.min(MAX_PARAMETERS)
    }

    fn named_objects(&self) -> usize {
        self.object_count.min(MAX_OBJECTS)
    }

    /// `operator==`: the same type, counts, ends and object identities.
    pub fn is_same_as(&self, other: &Self) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        }
        self.kind == other.kind
            && self.variables == other.variables
            && self.start.is_same_as(&other.start)
            && self.end.is_same_as(&other.end)
            && self.object_count == other.object_count
            && (0..self.named_objects()).all(|i| self.objects[i] == other.objects[i])
    }

    /// `ClassifySolutionEnd`: where one end of `other` falls against this
    /// solution, `compare_start` choosing the start or the end of `other`.
    ///
    /// Returns 0 for a match with a single value or a point inside a range,
    /// -1 and +1 for a match with the minimum or maximum end of a range,
    /// -2 and +2 for below the minimum and above the maximum, and +3 when
    /// the two are inconsistent and cannot be classified.
    pub fn classify_solution_end(&self, other: &Self, compare_start: bool) -> i64 {
        if self.variables != other.variables || self.object_count != other.object_count {
            return 3;
        }
        for k in 0..self.named_objects() {
            if self.objects[k] != other.objects[k] {
                return 3;
            }
        }
        let to_classify = if compare_start || other.kind == SolutionType::SingleValue {
            &other.start
        } else {
            &other.end
        };

        if self.kind == SolutionType::SingleValue {
            let mut ret = 3;
            if EFF_ZERO_SQRT < (self.start.value - to_classify.value).abs() {
                return ret;
            }
            for i in 0..self.variable_count() {
                let scaled_zero = EFF_ZERO * (1.0 + self.start.parameters[i].abs());
                if (self.start.parameters[i] - to_classify.parameters[i]).abs() < scaled_zero {
                    if ret != 3 && ret != 0 {
                        continue;
                    }
                    ret = 0;
                    continue;
                }
                if ret == 0 {
                    return 3;
                }
                if to_classify.parameters[0] < self.start.parameters[0] {
                    ret = -2;
                } else {
                    ret = 2;
                }
                break;
            }
            return ret;
        }

        let mut ret = 3;
        for i in 0..self.variable_count() {
            let at = to_classify.parameters[i];
            let scaled_zero = (1.0 + at.abs()).abs() * EFF_ZERO;
            let (start, end) = (self.start.parameters[i], self.end.parameters[i]);

            if (start - at).abs() < scaled_zero && (end - at).abs() < scaled_zero {
                continue;
            }
            if (start - at).abs() < scaled_zero {
                if ret == -2 || ret == 2 || ret == 1 {
                    return 3;
                }
                if ret != 3 && ret != -1 {
                    continue;
                }
                if (self.start.value - to_classify.value).abs() > EFF_ZERO_SQRT {
                    return 3;
                }
                ret = -1;
                continue;
            }
            if (end - at).abs() < scaled_zero {
                if ret == -2 || ret == 2 || ret == -1 {
                    return 3;
                }
                if ret != 3 && ret != 1 {
                    continue;
                }
                if (self.end.value - to_classify.value).abs() > EFF_ZERO_SQRT {
                    return 3;
                }
                ret = 1;
                continue;
            }
            if are_sequential(start, at, end) {
                if ret != 3 && ret != 0 {
                    return 3;
                }
                ret = 0;
                continue;
            }
            if are_sequential(at, start, end) {
                if ret != 3 && ret != -2 {
                    return 3;
                }
                ret = -2;
                continue;
            }
            if are_sequential(start, end, at) {
                if ret != 3 && ret != 2 {
                    return 3;
                }
                ret = 2;
                continue;
            }
            return 3;
        }
        ret
    }

    /// `GetIndex`: the index of the object that is, or with `want_object`
    /// false is not, `target`, and the index of its first solution parameter.
    ///
    /// On failure SMLib leaves the object index at the object count and the
    /// parameter index at the count so far, and this returns them in `Err`.
    pub fn index_of_object(
        &self,
        objects: &dyn SolutionObjects,
        target: Option<ObjectId>,
        want_object: bool,
    ) -> Result<(usize, usize), (usize, usize)> {
        let mut parameter_index = 0;
        for ii in 0..self.named_objects() {
            let object = self.objects[ii];
            if (want_object && target == object) || (!want_object && target != object) {
                return Ok((ii, parameter_index));
            }
            if let Some(id) = object {
                parameter_index += objects.kind(id).index_width();
            }
        }
        Err((self.object_count, parameter_index))
    }

    /// `GetObjectIndexForParameterIndex`: the object and the domain direction
    /// (0 for `s` or `u`, 1 for `v`, 2 for `w`) that solution parameter
    /// `parameter_index` belongs to. Both are zero when it returns `Err`.
    ///
    /// SMLib looks up the kind of `objects[0]` on every pass of its loop and
    /// not the kind of the object being counted, because it only updates the
    /// object index when it finds the answer. A mixed pair of curve and
    /// surface therefore answers as SMLib does, including where that differs
    /// from what the loop intends.
    pub fn object_for_parameter_index(
        &self,
        objects: &dyn SolutionObjects,
        parameter_index: usize,
    ) -> Result<(usize, usize), Fault> {
        if parameter_index >= self.variables || self.object_count == 0 {
            return Err(Fault::Degenerate);
        }
        let mut walked = 0;
        for ii in 0..self.named_objects() {
            let Some(id) = self.objects[0] else {
                return Err(Fault::Degenerate);
            };
            match objects.kind(id) {
                ObjectKind::Edge | ObjectKind::Curve => {
                    if walked + 1 > parameter_index {
                        return Ok((ii, 0));
                    }
                    walked += 1;
                }
                ObjectKind::Face | ObjectKind::Surface => {
                    if walked + 2 > parameter_index {
                        return Ok((ii, parameter_index - walked));
                    }
                    walked += 2;
                }
                ObjectKind::Volume => {
                    if walked + 3 > parameter_index {
                        return Ok((ii, parameter_index - walked));
                    }
                    walked += 3;
                }
                ObjectKind::Parameterless => {}
                ObjectKind::Other => return Err(Fault::Degenerate),
            }
        }
        Ok((0, 0))
    }

    /// `GetParameterInterval`: the interval of the domain that solution
    /// parameter `parameter_index` lies in. Fails when SMLib has no objects
    /// to search or the object has no parameters.
    pub fn parameter_interval(
        &self,
        objects: &dyn SolutionObjects,
        parameter_index: usize,
    ) -> Result<Extent1, Fault> {
        let (object_index, direction) =
            self.object_for_parameter_index(objects, parameter_index)?;
        let id = self.objects[object_index].ok_or(Fault::Degenerate)?;
        match objects.kind(id) {
            ObjectKind::Edge | ObjectKind::Curve => Ok(objects.parameter_interval(id, 0)),
            ObjectKind::Face | ObjectKind::Surface | ObjectKind::Volume => {
                Ok(objects.parameter_interval(id, direction))
            }
            ObjectKind::Parameterless | ObjectKind::Other => Err(Fault::Degenerate),
        }
    }

    /// `GetObjectTolerance`: the zone tolerance of the object at
    /// `object_index`, zero for a bad index or a missing object.
    pub fn object_tolerance(&self, objects: &dyn SolutionObjects, object_index: usize) -> f64 {
        if object_index >= self.object_count || object_index >= MAX_OBJECTS {
            return 0.0;
        }
        self.objects[object_index].map_or(0.0, |id| objects.zone_tolerance(id))
    }

    /// `IsEqual`: the same counts, type, value, objects and parameters, the
    /// parameters within the scaled zero or a ten-thousandth of their domain.
    pub fn is_equal(&self, objects: &dyn SolutionObjects, other: &Self) -> bool {
        if self.object_count != other.object_count
            || self.variables != other.variables
            || self.kind != other.kind
            || EFF_ZERO_SQRT < (self.start.value - other.start.value).abs()
        {
            return false;
        }
        if (0..self.named_objects()).any(|j| self.objects[j] != other.objects[j]) {
            return false;
        }
        for i in 0..self.variable_count() {
            let length1 = self
                .parameter_interval(objects, i)
                .map_or(0.0, |interval| 0.0001 * interval.length());
            let length2 = other
                .parameter_interval(objects, i)
                .map_or(0.0, |interval| 0.0001 * interval.length());
            let scaled_zero = max3(
                (1.0 + self.start.parameters[i].abs()).abs() * EFF_ZERO,
                length1,
                length2,
            );
            if (self.start.parameters[i] - other.start.parameters[i]).abs() > scaled_zero
                || (self.kind == SolutionType::RangeOfValues
                    && (self.end.parameters[i] - other.end.parameters[i]).abs() > scaled_zero)
            {
                return false;
            }
        }
        true
    }

    /// `SmGlobalSolver::IdenticalSolutions`, the default definition of an
    /// identical pair: the same counts, type, value, objects and parameters.
    /// Parameters match within `100 * SM_EFF_ZERO` scaled, or a
    /// ten-thousandth of their domain. An identical pair prefers the
    /// solution with the smaller value.
    pub fn compare_for_identity(
        &self,
        objects: &dyn SolutionObjects,
        other: &Self,
    ) -> (SolutionPair, Preferred) {
        if self.object_count != other.object_count
            || self.variables != other.variables
            || self.kind != other.kind
            || EFF_ZERO_SQRT < (self.start.value - other.start.value).abs()
        {
            return (SolutionPair::DifferentCounts, Preferred::Neither);
        }
        if (0..self.named_objects()).any(|j| self.objects[j] != other.objects[j]) {
            return (SolutionPair::DifferentObjects, Preferred::Neither);
        }
        for i in 0..self.variable_count() {
            let length1 = self
                .parameter_interval(objects, i)
                .map_or(0.0, |interval| 0.0001 * interval.length());
            let length2 = other
                .parameter_interval(objects, i)
                .map_or(0.0, |interval| 0.0001 * interval.length());
            let scaled_zero = max3(
                (1.0 + self.start.parameters[i].abs()).abs() * 100.0 * EFF_ZERO,
                length1,
                length2,
            );
            if (self.start.parameters[i] - other.start.parameters[i]).abs() > scaled_zero
                || (self.kind == SolutionType::RangeOfValues
                    && (self.end.parameters[i] - other.end.parameters[i]).abs() > scaled_zero)
            {
                return (SolutionPair::DifferentParameters, Preferred::Neither);
            }
        }
        let preferred = if self.start.value < other.start.value {
            Preferred::First
        } else {
            Preferred::Second
        };
        (SolutionPair::Identical, preferred)
    }
}

/// `SmSolutionArray`: the solutions of one search, in order.
///
/// The size of the array is its length. Slots that a grow opens up hold
/// [`Solution::zeroed`], as SMLib's zero-filled allocation does.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SolutionArray {
    items: Vec<Solution>,
}

impl std::ops::Index<usize> for SolutionArray {
    type Output = Solution;

    fn index(&self, index: usize) -> &Solution {
        &self.items[index]
    }
}

impl std::ops::IndexMut<usize> for SolutionArray {
    fn index_mut(&mut self, index: usize) -> &mut Solution {
        &mut self.items[index]
    }
}

impl SolutionArray {
    /// An empty array.
    pub fn new() -> Self {
        Self::default()
    }

    /// `GetSize`.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the array holds no solution.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// The solutions in order.
    pub fn as_slice(&self) -> &[Solution] {
        &self.items
    }

    /// The solutions in order, for editing in place.
    pub fn as_mut_slice(&mut self) -> &mut [Solution] {
        &mut self.items
    }

    /// `SetSize`: truncates, or grows with zeroed solutions.
    pub fn set_size(&mut self, size: usize) {
        self.items.resize(size, Solution::zeroed());
    }

    /// `ReSet` and `RemoveAll`: empties the array.
    pub fn clear(&mut self) {
        self.items.clear();
    }

    /// `GetAt`.
    pub fn get_at(&self, index: usize) -> Result<Solution, Fault> {
        self.items.get(index).copied().ok_or(Fault::InvalidInput)
    }

    /// `SetAt`: replaces an existing solution.
    pub fn set_at(&mut self, index: usize, solution: Solution) -> Result<(), Fault> {
        *self.items.get_mut(index).ok_or(Fault::InvalidInput)? = solution;
        Ok(())
    }

    /// `Add`: appends and returns the index of the new solution.
    pub fn add(&mut self, solution: Solution) -> usize {
        self.items.push(solution);
        self.items.len() - 1
    }

    /// `Append`: adds every solution of `source` at the end and returns the
    /// size the array had before.
    pub fn append(&mut self, source: &Self) -> usize {
        let old = self.items.len();
        self.items.extend_from_slice(&source.items);
        old
    }

    /// `Copy`: becomes a copy of `source`.
    pub fn copy_from(&mut self, source: &Self) {
        self.items.clone_from(&source.items);
    }

    /// `SetAtGrow`: sets one solution, growing the array to hold it.
    pub fn set_at_grow(&mut self, index: usize, solution: Solution) {
        if index >= self.items.len() {
            self.set_size(index + 1);
        }
        self.items[index] = solution;
    }

    /// `InsertAt`: inserts `count` copies of `solution` at `index`. An index
    /// at or past the end grows the array to `index + count` and fills the
    /// gap with zeroed solutions.
    pub fn insert_at(&mut self, index: usize, solution: Solution, count: usize) {
        if index >= self.items.len() {
            self.set_size(index + count);
        } else {
            let tail = self.items.split_off(index);
            self.items
                .extend(std::iter::repeat_n(Solution::zeroed(), count));
            self.items.extend(tail);
        }
        for slot in &mut self.items[index..index + count] {
            *slot = solution;
        }
    }

    /// `InsertAt(start, array)`: inserts every solution of `source` at `start`.
    pub fn insert_array_at(&mut self, start: usize, source: &Self) {
        if let Some(first) = source.items.first() {
            self.insert_at(start, *first, source.items.len());
            for (offset, solution) in source.items.iter().enumerate() {
                self.items[start + offset] = *solution;
            }
        }
    }

    /// `RemoveAt`: removes `count` solutions from `index` on.
    pub fn remove_at(&mut self, index: usize, count: usize) -> Result<(), Fault> {
        if index
            .checked_add(count)
            .is_none_or(|end| end > self.items.len())
        {
            return Err(Fault::InvalidInput);
        }
        self.items.drain(index..index + count);
        Ok(())
    }

    /// `Swap`: exchanges two solutions.
    pub fn swap(&mut self, a: usize, b: usize) -> Result<(), Fault> {
        if a >= self.items.len() || b >= self.items.len() {
            return Err(Fault::InvalidInput);
        }
        self.items.swap(a, b);
        Ok(())
    }

    /// `RemoveDuplicateSolutions`: drops later solutions whose first, then
    /// second, start parameter is exactly equal to an earlier one's.
    ///
    /// SMLib does not step back after a removal, so the solution that slides
    /// into the removed place is not compared with the current one. This
    /// keeps that behaviour.
    pub fn remove_duplicate_solutions(&mut self) {
        for parameter in 0..2 {
            let mut ii = 0;
            while ii < self.items.len() {
                let value = self.items[ii].start.parameters[parameter];
                let mut jj = ii + 1;
                while jj < self.items.len() {
                    if value == self.items[jj].start.parameters[parameter] {
                        self.items.remove(jj);
                    }
                    jj += 1;
                }
                ii += 1;
            }
        }
    }

    /// `SmSolutionArray::AddSortedSolution`: inserts `new` at its place in
    /// the order `key` gives, without combining solutions. Returns the index
    /// of the added solution, or `None` for a rejected one: a value of
    /// `SM_BIG_DOUBLE`, or one equal to a solution already held (which the
    /// held one replaces when the new value is smaller in magnitude).
    pub fn add_sorted_solution(
        &mut self,
        objects: &dyn SolutionObjects,
        new: Solution,
        key: SortKey,
    ) -> Option<usize> {
        if new.start.value == BIG_DOUBLE {
            return None;
        }
        let mut found_place = false;
        let mut insertion_index = self.items.len();
        let mut ii = 0;
        while ii < self.items.len() && !found_place {
            let old = self.items[ii];
            if old.is_equal(objects, &new) {
                if new.start.value.abs() < old.start.value.abs() {
                    self.items[ii] = new;
                }
                return None;
            }
            match key {
                SortKey::ByValue => {
                    if old.start.value > new.start.value {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByFirstParameter => {
                    if old.start.parameters[0] > new.start.parameters[0] {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByFirstTwoParameters => {
                    let scale = (1.0 + old.start.parameters[0]) * EFF_ZERO_SQRT;
                    if old.start.parameters[0] > new.start.parameters[0] + scale
                        || ((old.start.parameters[0] - new.start.parameters[0]).abs() < scale
                            && old.start.parameters[1] > new.start.parameters[1])
                    {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByValueWithParameter | SortKey::ByAllParameters => {
                    let mut by_parameters = true;
                    if key == SortKey::ByValueWithParameter {
                        let scale = (1.0 + old.start.value) * EFF_ZERO_SQRT;
                        if (old.start.value - new.start.value).abs() > scale {
                            if old.start.value > new.start.value {
                                insertion_index = ii;
                                found_place = true;
                            }
                            by_parameters = false;
                        }
                    }
                    if by_parameters {
                        for jj in 0..old.variable_count() {
                            let scale = (1.0 + old.start.parameters[jj]) * EFF_ZERO_SQRT;
                            if old.start.parameters[jj] > new.start.parameters[jj] + scale {
                                insertion_index = ii;
                                found_place = true;
                                break;
                            } else if old.start.parameters[jj] < new.start.parameters[jj] - scale {
                                break;
                            }
                        }
                    }
                }
            }
            ii += 1;
        }
        if found_place {
            self.insert_at(insertion_index, new, 1);
            Some(insertion_index)
        } else {
            Some(self.add(new))
        }
    }
}
