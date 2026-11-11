#!/usr/bin/env python3
"""Regenerate spline object fixtures with the pinned strict SMLib driver."""
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
RNG = random.Random(2026111103)
CASES = []
NO_Z = sys.float_info.max
NO_W = -sys.float_info.max


def h(value):
    return struct.pack('>d', float(value)).hex()


def full_knots(degree, count):
    inner = count - degree - 1
    return [0.] * (degree + 1) + [float(i + 1) for i in range(inner)] + [float(inner + 1)] * (degree + 1)


def packed(points, rational, dim=3):
    values = []
    for point in points:
        weight = RNG.uniform(.5, 2) if rational else 1.
        values += [h(point[0] * weight), h(point[1] * weight),
                   h(point[2] * weight if dim == 3 else NO_Z), h(weight if rational else NO_W)]
    return values


def curve(degree, count, rational, dim, kind, t, derivatives, left, nonzero, points=None, knots=None):
    knots = full_knots(degree, count) if knots is None else knots
    if points is None:
        points = [[RNG.uniform(-4, 4) for _ in range(3)] for _ in range(count)]
    if kind == 'singular':
        points = [[2., 1., 3.]] * count
    elif kind == 'tangent':
        points[1] = points[0][:]
    elif kind == 'touching':
        points[-1] = points[0][:]
    elif kind == 'near-singular' and knots[-1] > 1e-12:
        points[1] = [x + 1e-13 for x in points[0]]
    tokens = [str(degree), str(count), str(dim)] + list(map(h, knots)) + packed(points, rational, dim)
    tokens += [h(1e-10), h(t), str(derivatives), str(int(left)), str(int(nonzero))]
    CASES.append(('curve', kind, ' '.join(tokens)))


for degree in (1, 2, 3, 4):
    for rational in (False, True):
        for dim in (2, 3):
            for count in (degree + 1, degree + 3):
                end = count - degree
                for derivatives in range(5):
                    curve(degree, count, rational, dim, 'regular', RNG.uniform(.05, end - .05), derivatives, True, False)
                for t in (0., float(end), -1., end + 1.):
                    for derivatives in (0, 1, 2):
                        curve(degree, count, rational, dim, 'boundary', t, derivatives, False, True)
                for kind in ('singular', 'tangent', 'touching', 'near-singular'):
                    for derivatives in (1, 2):
                        curve(degree, count, rational, dim, kind, 0., derivatives, True, True)

for degree in (7,15,31):
    for rational in (False,True):
        for dim in (2,3):
            for count in (degree+1,degree+3):
                for derivatives in (0,1,4):
                    curve(degree,count,rational,dim,'regular',.37,derivatives,True,False)

# Different one-sided derivatives at an interior knot of full multiplicity.
for rational in (False, True):
    for left in (False, True):
        for derivatives in range(5):
            curve(2, 5, rational, 3, 'touching', 1., derivatives, left, False,
                  knots=[0, 0, 0, 1, 1, 2, 2, 2])


for rational in (False,True):
    for t in (-2.,-1.,0.,.5,1.,2.,3.,4.,5.):
        for derivatives in (0,1,2):
            curve(2,5,rational,3,'boundary',t,derivatives,True,False,
                  knots=[-2,-1,0,1,2,3,4,5])


for width in (1e-24,2e-24,1e-13):
    for derivatives in (0,1,2):
        curve(2,3,False,3,'near-singular',width/2,derivatives,True,False,
              points=[[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]],knots=[0.,0.,0.,width,width,width])


def surface(p, q, n, m, rational, kind):
    u, v = full_knots(p, n), full_knots(q, m)
    points = [[float(i), float(j), RNG.uniform(-2, 2)] for i in range(n) for j in range(m)]
    if kind == 'singular':
        points = [[1., 2., 3.]] * (n * m)
    elif kind == 'tangent':
        points[:m] = points[m:2*m]
    elif kind == 'touching':
        points[-m:] = points[:m]
    elif kind == 'near-singular':
        points[:m] = [[x + 1e-13 for x in point] for point in points[m:2*m]]
    tokens = list(map(str, (p, q, n, m))) + list(map(h, u + v)) + packed(points, rational)
    tokens += list(map(h, (1e-10, 1., 1., 1.))) + [str(0 if kind=='boundary' else 10 if kind=='tangent' else 4)]
    if kind == 'boundary':
        tokens += list(map(h, (.25, .25, u[-1] - .25, v[-1] - .25)))
    else:
        tokens += list(map(h, (u[0], v[0], u[-1], v[-1])))
    CASES.append(('surface', kind, ' '.join(tokens)))


for p, q in ((1, 1), (2, 1), (2, 2), (3, 2), (3, 3)):
    for rational in (False, True):
        for extra in (0, 2):
            for kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary'):
                surface(p, q, p + 1 + extra, q + 1 + extra, rational, kind)


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2])
    raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request = f'{op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        (raw / f'{index:04d}-{op}.txt').write_text(
            f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode:
            sys.exit(f'case {index} {op}: oracle exit {run.returncode}: {run.stderr}')
        response = run.stdout.split()
        assert response[0] == op
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    (HERE / 'cases.txt').write_text('# Generated by generate.py from the pinned strict SMLib; see README.md.\n' + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
