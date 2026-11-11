// Test-only oracle driver for the catallax-trace local solver gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. The test functions are defined
// identically in ../../local_solve.rs. Catallax never compiles or links this.
#include "StdAfx.h"
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtentNd.h>
#include <SmLocalSolve1d.h>
#include <SmLocalSolveNd.h>
#include <SmMatrix.h>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

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

// One-variable families. p[0..3] are parameters.
class Family1 : public SmEvalFunctionObject {
public:
    long kind;
    double p[4];
    SmStatus Evaluate(double t, double &f, double &fp, SmBoolean &found, SmBoolean) override {
        found = FALSE;
        switch (kind) {
        case 1: // cubic with roots p0, p1, p2
            f = (t - p[0]) * (t - p[1]) * (t - p[2]);
            fp = (t - p[1]) * (t - p[2]) + (t - p[0]) * (t - p[2]) + (t - p[0]) * (t - p[1]);
            break;
        case 2: // parabola (t - p0)^2 + p1
            f = (t - p[0]) * (t - p[0]) + p[1];
            fp = 2.0 * (t - p[0]);
            break;
        case 3: { // signed square root about p0, derivative softened by p1
            double d = t - p[0];
            f = d >= 0.0 ? std::sqrt(d) : -std::sqrt(-d);
            fp = 0.5 / std::sqrt(std::fabs(d) + p[1]);
            break;
        }
        case 4: // line t - p0 that cannot be evaluated beyond p1
            if (t > p[1]) return SM_ERR;
            f = t - p[0];
            fp = 1.0;
            break;
        case 5: // quartic (t - p0)^2 (t - p1)^2 - p2, a double well
            f = (t - p[0]) * (t - p[0]) * (t - p[1]) * (t - p[1]) - p[2];
            fp = 2.0 * (t - p[0]) * (t - p[1]) * (t - p[1]) + 2.0 * (t - p[0]) * (t - p[0]) * (t - p[1]);
            break;
        default:
            return SM_ERR;
        }
        if (f == 0.0) found = TRUE;
        return SM_SUCCESS;
    }
};

// N-variable families.
class FamilyN : public SmEvalNFunctionsObject {
public:
    long kind;
    double p[4];
    SmStatus Evaluate(const SmTArray<double> &x, SmTArray<double> &f, SmMatrix *jacobian, SmBoolean &found) override {
        found = FALSE;
        switch (kind) {
        case 1: // circle x^2 + y^2 = p0 and line y = p1 x + p2
        case 2: // the same with a finite-difference Jacobian
            f[0] = x[0] * x[0] + x[1] * x[1] - p[0];
            f[1] = x[1] - p[1] * x[0] - p[2];
            if (jacobian && kind == 1) {
                (*jacobian)[0][0] = 2.0 * x[0];
                (*jacobian)[0][1] = 2.0 * x[1];
                (*jacobian)[1][0] = -p[1];
                (*jacobian)[1][1] = 1.0;
            }
            if (jacobian && kind == 2) {
                SmTArray<double> xx = x;
                if (ComputeJacobian(xx, f, *jacobian) != SM_SUCCESS) return SM_ERR;
            }
            break;
        case 3: // two unit circles centred at (p0, 0) and (p1, 0)
            f[0] = (x[0] - p[0]) * (x[0] - p[0]) + x[1] * x[1] - 1.0;
            f[1] = (x[0] - p[1]) * (x[0] - p[1]) + x[1] * x[1] - 1.0;
            if (jacobian) {
                (*jacobian)[0][0] = 2.0 * (x[0] - p[0]);
                (*jacobian)[0][1] = 2.0 * x[1];
                (*jacobian)[1][0] = 2.0 * (x[0] - p[1]);
                (*jacobian)[1][1] = 2.0 * x[1];
            }
            break;
        case 4: // a rank-deficient linear system in three variables
            f[0] = x[0] + x[1] + x[2] - p[0];
            f[1] = 2.0 * x[0] + 2.0 * x[1] + 2.0 * x[2] - p[1];
            f[2] = x[0] - x[2] - p[2];
            if (jacobian) {
                double rows[3][3] = {{1, 1, 1}, {2, 2, 2}, {1, 0, -1}};
                for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) (*jacobian)[i][j] = rows[i][j];
            }
            break;
        default:
            return SM_ERR;
        }
        return SM_SUCCESS;
    }
};

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "solve1") {
        Family1 family;
        family.kind = ri();
        for (double &v : family.p) v = rd();
        double guess = rd();
        long has_interval = ri();
        double lo = rd(), hi = rd();
        long periodic = ri(), algorithm = ri();
        double acceptable = rd(), desired = rd();
        SmLocalSolve1d solver = has_interval ? SmLocalSolve1d(family, SmExtent1d(lo, hi), periodic ? TRUE : FALSE)
                                             : SmLocalSolve1d(family);
        solver.SetSolverAlgorithm((SmLocalSolverAlgorithmType)algorithm);
        if (desired > 0.0) solver.SetDesiredAccuracy(desired);
        SmBoolean found = FALSE;
        double t = 7.0;
        SmStatus status = solver.SolveIt(guess, acceptable, found, t);
        oi((long)status);
        if (status == SM_SUCCESS) {
            oi(found); od(t); od(solver.GetFoundAccuracy()); oi((long)solver.GetTerminationReason());
        }
    } else if (op == "solven") {
        FamilyN family;
        family.kind = ri();
        for (double &v : family.p) v = rd();
        long n = ri();
        SmTArray<double> guess;
        for (long i = 0; i < n; ++i) guess.Add(rd());
        long has_intervals = ri();
        SmExtentNd intervals(n);
        SmTArray<SmBoolean> periodic;
        for (long i = 0; i < n; ++i) {
            double lo = rd(), hi = rd();
            intervals[i] = SmExtent1d(lo, hi);
            periodic.Add(ri() ? TRUE : FALSE);
        }
        long handler = ri(), max_hits = ri(), cutting = ri();
        double acceptable = rd();
        SmLocalSolveNd solver(family, has_intervals ? &intervals : NULL, has_intervals ? &periodic : NULL);
        solver.SetBoundaryHandler((SmBoundaryHandlerType)handler, (ULONG)max_hits);
        solver.SetStepCuttingEnabled(cutting ? TRUE : FALSE);
        SmBoolean found = FALSE;
        SmTArray<double> solution;
        SmStatus status = solver.SolveIt(guess, acceptable, found, solution);
        oi((long)status);
        if (status == SM_SUCCESS) {
            oi(found);
            for (ULONG i = 0; i < solution.GetSize(); ++i) od(solution[i]);
            oi((long)solver.GetTerminationReason());
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
