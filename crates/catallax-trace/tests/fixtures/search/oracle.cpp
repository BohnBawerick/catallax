// Test-only oracle driver for the catallax-trace search gate (step 2:
// SmGlobalSolver and SmSolutionArray). Reads one request line on stdin: an
// operation name and its inputs, doubles as 16-digit hexadecimal binary64 bit
// patterns and integers in decimal. Prints one response line: tokens `d<hex>`
// for binary64 values and `i<decimal>` for integers, statuses, booleans and
// enumerations. Catallax never compiles or links this file.
//
// Objects. Solutions name objects by index into a fixed table, -1 for none:
//   0 SmLine on [0, 1]          4 SmVertex (no parameters)
//   1 SmLine on [-2, 3.5]       5 SmEdge on [0.25, 0.75]
//   2 SmPlane on [0,2]x[-1,1]   6 SmFace on [1,3]x[2,5]
//   3 SmPlane on [0.5,4]x[2,3]  7 SmObject (no kind)
//
// Solution tokens: type nvars nobjs obj0..obj4 value start0..start3
//                  value end0..end3 (parameters 4..9 are zero)
// Tree tokens:     N, then N nodes `bbox(6 doubles) aux child1 child2`
//                  (children are node indices, -1 for none).
#include "StdAfx.h"
#include <SmContext.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmGlobalSolver.h>
#include <SmLine.h>
#include <SmPlane.h>
#include <SmSolutionArray.h>
#include <SmTree.h>
#include <SmVertex.h>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

struct VertexProbe : SmVertex { VertexProbe() : SmVertex() {} };
struct EdgeProbe : SmEdge { EdgeProbe() : SmEdge() {} };
struct FaceProbe : SmFace { FaceProbe() : SmFace() {} };

// SmGlobalSolver with its protected methods exposed to the driver.
struct Probe : SmGlobalSolver {
    SmSolutionPairType identical(const SmSolution &a, const SmSolution &b, const SmSolution *&p) const {
        return IdenticalSolutions(a, b, p);
    }
    SmBoolean branch(SmTreeNode **b) { return BranchMayContainAnswers(b); }
    SmBoolean ready(SmTreeNode **b) const { return AreReadyForLocalSolve(b); }
    SmStatus findBest(SmTreeNode **b, ULONG &i, SmBoolean &swap) const { return FindBestSubdivisionNode(b, i, swap); }
    SmStatus subdivide(SmTreeNode **o, SmTreeNode **c1, SmTreeNode **c2) { return Subdivide(o, c1, c2); }
    void setVectors(const SmVector3d *v) { m_cpOptVectors = v; }
    void setTree(ULONG i, SmTree *t) { m_apTrees[i] = t; }
    double best() const { return m_dBestAnswerSoFarSq; }
};

struct ArrayProbe : SmSolutionArray { using SmSolutionArray::SetAtGrow; };

struct Harness : SmTree {
    Harness() : SmTree((SmObject *)NULL) {}
    void top(SmTreeNode *n) { m_pTopNode = n; }
};

static std::istringstream in;

static double rd() {
    std::string token;
    if (!(in >> token)) { std::fprintf(stderr, "missing input\n"); std::exit(3); }
    uint64_t bits = std::strtoull(token.c_str(), nullptr, 16);
    double value;
    std::memcpy(&value, &bits, sizeof value);
    return value;
}
static long ri() {
    long value;
    if (!(in >> value)) { std::fprintf(stderr, "missing integer\n"); std::exit(3); }
    return value;
}
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }

static SmContext *context;
static std::vector<SmObject *> table;

static long objectId(const SmObject *p) {
    if (!p) return -1;
    for (size_t i = 0; i < table.size(); ++i) if (table[i] == p) return (long)i;
    std::fprintf(stderr, "unknown object\n");
    std::exit(3);
}
static SmObject *objectAt(long id) { return id < 0 ? NULL : table.at((size_t)id); }

static void buildTable() {
    SmLine *a = new (*context) SmLine(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmExtent1d(0.0, 1.0), 1.0, 3, context);
    SmLine *b = new (*context) SmLine(SmPoint3d(0, 1, 0), SmVector3d(0, 1, 0), SmExtent1d(-2.0, 3.5), 2.0, 3, context);
    SmPlane *s1 = new (*context) SmPlane(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), SmVector2d(1, 1), SmExtent2d(0, -1, 2, 1), context);
    SmPlane *s2 = new (*context) SmPlane(SmPoint3d(0, 0, 1), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), SmVector2d(1, 1), SmExtent2d(0.5, 2, 4, 3), context);
    VertexProbe *v = new (a) VertexProbe();
    EdgeProbe *e = new (a) EdgeProbe();
    e->SetInterval(SmExtent1d(0.25, 0.75));
    FaceProbe *f = new (a) FaceProbe();
    f->SetUVDomain(SmExtent2d(SmPoint2d(1, 2), SmPoint2d(3, 5)));
    SmObject *base = ::new SmObject(FALSE);
    table = {a, b, s1, s2, v, e, f, base};
}

static std::vector<SmTreeNode *> nodes;

static long nodeId(const SmTreeNode *p) {
    if (!p) return -1;
    for (size_t i = 0; i < nodes.size(); ++i) if (nodes[i] == p) return (long)i;
    std::fprintf(stderr, "unknown node\n");
    std::exit(3);
}

static SmSolution rsol() {
    SmSolution s;
    s.m_eSolutionType = (SmSolutionType)ri();
    s.m_lNumVariables = (ULONG)ri();
    s.m_lNumObjects = (ULONG)ri();
    for (int k = 0; k < SM_SA_MAX_TREES; ++k) s.m_apObjects[k] = objectAt(ri());
    s.m_vStart.m_dSolutionValue = rd();
    for (int k = 0; k < 4; ++k) s.m_vStart.m_adParameters[k] = rd();
    s.m_vEnd.m_dSolutionValue = rd();
    for (int k = 0; k < 4; ++k) s.m_vEnd.m_adParameters[k] = rd();
    return s;
}

static void osol(const SmSolution &s) {
    oi((long)s.m_eSolutionType);
    oi((long)s.m_lNumVariables);
    oi((long)s.m_lNumObjects);
    for (int k = 0; k < SM_SA_MAX_TREES; ++k) oi(objectId(s.m_apObjects[k]));
    for (int k = 0; k < SM_SA_MAX_TREES; ++k) oi(nodeId(s.m_apNodes[k]));
    od(s.m_vStart.m_dSolutionValue);
    for (int k = 0; k < SM_SA_MAX_TREES * 2; ++k) {
        if (k >= 4 && s.m_vStart.m_adParameters[k] != 0.0) { std::fprintf(stderr, "stray parameter\n"); std::exit(3); }
        if (k < 4) od(s.m_vStart.m_adParameters[k]);
    }
    od(s.m_vEnd.m_dSolutionValue);
    for (int k = 0; k < SM_SA_MAX_TREES * 2; ++k) {
        if (k >= 4 && s.m_vEnd.m_adParameters[k] != 0.0) { std::fprintf(stderr, "stray parameter\n"); std::exit(3); }
        if (k < 4) od(s.m_vEnd.m_adParameters[k]);
    }
}

// The solver gate compares node pairs, so it prints each solution briefly.
static void osolBrief(const SmSolution &s) {
    oi((long)s.m_eSolutionType);
    oi((long)s.m_lNumVariables);
    oi((long)s.m_lNumObjects);
    for (int k = 0; k < SM_SA_MAX_TREES; ++k) oi(nodeId(s.m_apNodes[k]));
    od(s.m_vStart.m_dSolutionValue);
}

static void oarray(SmSolutionArray &a) {
    oi((long)a.GetSize());
    for (ULONG i = 0; i < a.GetSize(); ++i) osol(a[i]);
}

static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }

static SmExtent3d rbox() {
    double a[6];
    for (double &v : a) v = rd();
    return SmExtent3d(a[0], a[1], a[2], a[3], a[4], a[5]);
}

// Reads N nodes and links them.
static void rtree() {
    long n = ri();
    std::vector<long> c1, c2;
    for (long i = 0; i < n; ++i) {
        SmTreeNode *node = new SmTreeNode();
        node->m_sBBox = rbox();
        node->m_eAuxDataType = (SmAuxDataType)ri();
        c1.push_back(ri());
        c2.push_back(ri());
        nodes.push_back(node);
    }
    for (long i = 0; i < n; ++i) {
        if (c1[i] >= 0) nodes[i]->m_pChild1 = nodes[c1[i]];
        if (c2[i] >= 0) nodes[i]->m_pChild2 = nodes[c2[i]];
    }
}

static SmTreeNode *nodeAt(long id) { return id < 0 ? NULL : nodes.at((size_t)id); }

int main() {
    SmContext ctx;
    context = &ctx;
    buildTable();
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    Probe probe;
    SmVector3d vectors[3];
    if (op == "solcmp") {
        SmSolution a = rsol(), b = rsol();
        oi(a == b);
        oi(a.ClassifySolutionEnd(b, TRUE));
        oi(a.ClassifySolutionEnd(b, FALSE));
        oi(a.IsEqual(b));
        const SmSolution *preferred = NULL;
        oi((long)probe.identical(a, b, preferred));
        oi(preferred == NULL ? 0 : preferred == &a ? 1 : 2);
    } else if (op == "solidx") {
        SmSolution s = rsol();
        long p = ri(), target = ri(), flag = ri();
        ULONG objectIndex, paramIndex;
        oi((long)s.GetIndex(objectAt(target), flag ? TRUE : FALSE, objectIndex, paramIndex));
        oi((long)objectIndex);
        oi((long)paramIndex);
        ULONG object, direction;
        SmStatus st = s.GetObjectIndexForParameterIndex((ULONG)p, object, direction);
        oi((long)st);
        oi((long)object);
        oi((long)direction);
        SmExtent1d interval;
        st = s.GetParameterInterval((ULONG)p, interval);
        oi((long)st);
        if (st == SM_SUCCESS) { od(interval.GetMin()); od(interval.GetMax()); }
        od(s.GetObjectTolerance((ULONG)p));
    } else if (op == "addsorted") {
        long operation = ri(), requested = ri(), category = ri();
        double tolerance = rd(), best = rd();
        long variables = ri(), key = ri(), combine = ri();
        probe.SetSolverOperation((SmSolverOperationType)operation);
        probe.SetSolutionRequested((SmSolutionRequestedType)requested);
        probe.SetOperationCategory((SmOperationCategory)category);
        probe.Set3DTolerance(tolerance);
        probe.SetBestAnswerSoFarSq(best);
        probe.SetNumVariables((ULONG)variables);
        SmSolutionArray array;
        probe.SetSolutions(&array);
        long initial = ri();
        for (long i = 0; i < initial; ++i) array.Add(rsol());
        long adds = ri();
        for (long i = 0; i < adds; ++i) {
            SmSolution s = rsol();
            probe.AddSortedSolution(s, (SmSortKeyType)key, combine ? TRUE : FALSE);
            oi((long)array.GetSize());
            od(probe.best());
        }
        oarray(array);
    } else if (op == "array") {
        ArrayProbe array;
        long ops = ri();
        for (long i = 0; i < ops; ++i) {
            static const char *const names[] = {"add", "ins", "rem", "size", "set", "grow", "swap", "dedup",
                                                "sorted", "append", "copy", "insarr", "reset"};
            std::string name = names[ri()];
            if (name == "add") {
                oi((long)array.Add(rsol()));
            } else if (name == "ins") {
                long index = ri(), count = ri();
                SmSolution s = rsol();
                array.InsertAt((ULONG)index, s, (ULONG)count);
            } else if (name == "rem") {
                long index = ri(), count = ri();
                array.RemoveAt((ULONG)index, (ULONG)count);
            } else if (name == "size") {
                array.SetSize((ULONG)ri());
            } else if (name == "set") {
                long index = ri();
                SmSolution s = rsol();
                array.SetAt((ULONG)index, s);
            } else if (name == "grow") {
                long index = ri();
                SmSolution s = rsol();
                array.SetAtGrow((ULONG)index, s);
            } else if (name == "swap") {
                long a = ri(), b = ri();
                array.Swap((ULONG)a, (ULONG)b);
            } else if (name == "dedup") {
                array.RemoveDuplicateSolutions();
            } else if (name == "sorted") {
                long key = ri();
                SmSolution s = rsol();
                oi((long)array.AddSortedSolution(s, (SmSortKeyType)key));
            } else if (name == "append" || name == "copy" || name == "insarr") {
                long index = name == "insarr" ? ri() : 0;
                long count = ri();
                SmSolutionArray other;
                for (long k = 0; k < count; ++k) other.Add(rsol());
                if (name == "append") array.Append(other);
                else if (name == "copy") array.Copy(other);
                else array.InsertAt((ULONG)index, other);
            } else if (name == "reset") {
                array.ReSet();
            } else {
                std::fprintf(stderr, "unknown array operation %s\n", name.c_str());
                return 2;
            }
            oi((long)array.GetSize());
        }
        oarray(array);
    } else if (op == "branch") {
        long operation = ri(), trees = ri();
        double tolerance = rd(), best = rd(), atDistance = rd();
        for (auto &v : vectors) v = r3();
        SmTreeNode *branch[SM_GS_MAX_TREES];
        for (int i = 0; i < SM_GS_MAX_TREES; ++i) {
            SmTreeNode *node = new SmTreeNode();
            node->m_sBBox = rbox();
            nodes.push_back(node);
            branch[i] = node;
        }
        probe.SetSolverOperation((SmSolverOperationType)operation);
        probe.SetNumTrees((ULONG)trees);
        probe.Set3DTolerance(tolerance);
        probe.SetBestAnswerSoFarSq(best);
        probe.SetAtDistance(atDistance);
        probe.setVectors(vectors);
        oi(probe.branch(branch));
    } else if (op == "tree") {
        long operation = ri(), trees = ri();
        double tolerance = rd(), best = rd();
        for (auto &v : vectors) v = r3();
        rtree();
        SmTreeNode *branch[SM_GS_MAX_TREES];
        for (int i = 0; i < SM_GS_MAX_TREES; ++i) branch[i] = nodeAt(ri());
        probe.SetSolverOperation((SmSolverOperationType)operation);
        probe.SetNumTrees((ULONG)trees);
        probe.Set3DTolerance(tolerance);
        probe.SetBestAnswerSoFarSq(best);
        probe.setVectors(vectors);
        oi(probe.ready(branch));
        ULONG index = 0;
        SmBoolean swap = FALSE;
        SmStatus st = probe.findBest(branch, index, swap);
        oi((long)st);
        if (st == SM_SUCCESS) { oi((long)index); oi(swap); }
        SmTreeNode *c1[SM_GS_MAX_TREES] = {NULL, NULL, NULL};
        SmTreeNode *c2[SM_GS_MAX_TREES] = {NULL, NULL, NULL};
        st = probe.subdivide(branch, c1, c2);
        oi((long)st);
        if (st == SM_SUCCESS) for (long i = 0; i < trees; ++i) { oi(nodeId(c1[i])); oi(nodeId(c2[i])); }
    } else if (op == "solve") {
        long operation = ri(), requested = ri(), trees = ri(), variables = ri();
        double tolerance = rd(), best = rd();
        for (auto &v : vectors) v = r3();
        SmSolutionArray array;
        long prior = ri();
        for (long i = 0; i < prior; ++i) array.Add(rsol());
        rtree();
        std::vector<Harness *> harness;
        for (int i = 0; i < SM_GS_MAX_TREES; ++i) {
            Harness *h = ::new Harness();
            SmTreeNode *top = nodeAt(ri());
            h->top(top);
            harness.push_back(h);
        }
        probe.SetNumTrees((ULONG)trees);
        probe.SetNumVariables((ULONG)variables);
        for (long i = 0; i < trees; ++i) probe.setTree((ULONG)i, harness[i]);
        SmStatus st = probe.SolveIt((SmSolverOperationType)operation, (SmSolutionRequestedType)requested,
                                    tolerance, best, vectors, array);
        oi((long)st);
        oi((long)array.GetSize());
        for (ULONG i = 0; i < array.GetSize(); ++i) osolBrief(array[i]);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    std::fflush(stdout);
    _exit(0);
}
