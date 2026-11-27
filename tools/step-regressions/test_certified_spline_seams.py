"""Raw two-chart spline prism STEP files against independent definition and file oracles.

The scenarios, formulas and limits are the designer's frozen contract in
fixtures/certified-spline-seams.json. Neither oracle reads a kernel property or
an OCCT number: the first integrates the retained ellipse law, the second the
written Bezier polynomials. OCCT is then tested against those answers.
"""
import copy
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
from importlib.metadata import version
import json
from math import comb, dist, isfinite
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_FORWARD, TopAbs_IN, TopAbs_OUT,
                        TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX)
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_loopless_faces import BINARY, DRIVER, ROOT, unique
from test_p5_flanged_housing import ordinary_shape
from raw_step import raw_shapes

CONTRACT = json.loads((Path(__file__).parent / 'fixtures/certified-spline-seams.json').read_text())
EXAMPLE = ROOT / '.claude/skills/verify-catallax/references/examples/certified-recipe.json'
ORACLE = CONTRACT['property_oracle']
EPSILONS = (1e-6, 1e-8, 1e-10, 1e-12)
ROUNDOFF = Decimal('1e-35')
# Accepted base geometry and certificate digests, after exact chart reassembly.
FINGERPRINTS = {
    'offset': ('30fd8b0c12cfeb0f123abf5affe1491c80cc3fbad43145f96186ff3502abd477',
               '77259346ecb7cb4508c17d33842114850148a646e7e4fa03c8e5f45419f611c9'),
    'second': ('a69cfe0f941f6ed11502316403d84c83d0c458097f047bb31040db1765883ea9',
               '5aaf9c6ec621c994ebbb4377bec71b0b701ad320ee5292afdf433857701f6f3c'),
    'after_trim': ('ed1ec54082e4776bd689ebef58772649b79f22b18c5a5d2389af0d837a704f4f',
                   '91c5409dce8e40add112b56ebb65968cc80286818e8b22a5f3cdbab97f050395'),
}
# Orthonormal frame of the sweep: n, then the ellipse's a=2 and b=8/5 axes.
AXIS = (Fraction(3, 5), Fraction(0), Fraction(4, 5))
MAJOR = (Fraction(0), Fraction(1), Fraction(0))
MINOR = (Fraction(-4, 5), Fraction(0), Fraction(3, 5))


def problems():
    offset = json.loads(EXAMPLE.read_text())
    second = copy.deepcopy(offset)
    second['recipe']['steps'].append(
        {'id': 'second', 'op': 'offset', 'body': 'offset', 'distance': '-1/4'})
    second['recipe']['result'] = 'second'
    trimmed = copy.deepcopy(offset)
    trimmed['recipe']['steps'].extend([
        {'id': 'trim_window', 'op': 'block', 'size': ['8', '8', '1'],
         'pose': {'rotation': ['3', '0', '1', '0'], 'translation': ['-7/5', '-4', '24/5']}},
        {'id': 'trimmed', 'op': 'boolean', 'kind': 'intersect',
         'target': 'offset', 'tool': 'trim_window'},
        {'id': 'after_trim', 'op': 'offset', 'body': 'trimmed', 'distance': '-1/4'},
    ])
    trimmed['recipe']['result'] = 'after_trim'
    trimmed['field_points'][0] = ['21/10', '0', '14/5']
    return {'offset': offset, 'second': second, 'after_trim': trimmed}


def step_only(problem):
    problem = copy.deepcopy(problem)
    problem['outputs'] = {'step': problem['outputs']['step']}
    problem['properties'] = None
    problem['field_points'] = []
    problem['edit'] = None
    return problem


def run(problem, out, allowance):
    source = out.with_suffix('.json')
    source.write_text(json.dumps(problem))
    subprocess.run(
        [str(DRIVER), 'run', str(source), '--binary', str(BINARY), '--max-work', str(allowance),
         '--out', str(out), '--fresh', '--workers', '1', '--json'],
        check=True, capture_output=True, text=True, timeout=CONTRACT['hang_guard_seconds'])
    return json.loads((out / 'run.json').read_text())['cases'][0]


# Definition oracle: exact intervals from the retained ellipse law alone.

def atan_interval(q, terms):
    total = sum(((-1) ** k * q ** (2 * k + 1) / (2 * k + 1) for k in range(terms)), Fraction(0))
    next_total = total + (-1) ** terms * q ** (2 * terms + 1) / (2 * terms + 1)
    return min(total, next_total), max(total, next_total)


def pi_interval():
    five = atan_interval(Fraction(1, 5), 100)
    other = atan_interval(Fraction(1, 239), 30)
    return 16 * five[0] - 4 * other[1], 16 * five[1] - 4 * other[0]


def target(distance, height):
    pi = pi_interval()
    m = Fraction(9, 25)
    series = sum((Fraction(comb(2 * k, k), 4 ** k) ** 2 / (2 * k - 1) * m ** k
                  for k in range(1, 65)), Fraction(0))
    tail = m ** 65 / (1 - m)
    perimeter = 4 * pi[0] * (1 - series - tail), 4 * pi[1] * (1 - series)
    areas = [(Fraction(16, 5) + distance ** 2) * p - distance * length
             for p in pi for length in perimeter]
    lengths = [length - 2 * p * distance for p in pi for length in perimeter]
    area = min(areas), max(areas)
    length = min(lengths), max(lengths)
    return {'pi': pi, 'area': area, 'length': length, 'height': height,
            'volume': (area[0] * height, area[1] * height),
            'boundary': (2 * area[0] + height * length[0], 2 * area[1] + height * length[1])}


def volume_envelope(goal, epsilon):
    """Volumes of the eroded and dilated target, outward at every interval endpoint."""
    pi = goal['pi']
    height = goal['height']
    inner = ((goal['area'][0] - goal['length'][1] * epsilon + pi[0] * epsilon ** 2)
             * (height - 2 * epsilon))
    outer = (goal['volume'][1] + goal['boundary'][1] * epsilon
             + (pi[1] * height + pi[1] * goal['length'][1] / 2) * epsilon ** 2
             + 4 * pi[1] * epsilon ** 3 / 3)
    return inner, outer


# File oracle: exact arithmetic on the numerals actually written.

def fields(value):
    out, depth, start, quoted = [], 0, 0, False
    for index, character in enumerate(value):
        if character == "'":
            quoted = not quoted
        elif not quoted and character == '(':
            depth += 1
        elif not quoted and character == ')':
            depth -= 1
        elif not quoted and character == ',' and depth == 0:
            out.append(value[start:index])
            start = index + 1
    return out + [value[start:]]


def items(value):
    return fields(value[1:-1])


def references(value):
    return [int(item) for item in re.findall(r'#(\d+)', value)]


def entities(text):
    return {int(number): (name, fields(body))
            for number, name, body in re.findall(r'#(\d+)=(\w+)\((.*?)\);', text)}


def named(es, name):
    return {number: body for number, (kind, body) in es.items() if kind == name}


def coordinates(es, number):
    return tuple(Fraction(item) for item in items(es[number][1][1]))


def chart(es, surface):
    rows = [references(row) for row in items(surface[3])]
    degree_flags = surface[1:3] + surface[5:8]
    if degree_flags != ['3', '1', '.F.', '.F.', '.F.']:
        raise AssertionError(f'not an open cubic prism chart: {degree_flags}')
    multiplicities = [int(item) for item in items(surface[8])]
    knots = [Fraction(item) for item in items(surface[10])]
    # Interior multiplicity 3 makes every group of four controls one Bezier span.
    if (multiplicities != [4] + [3] * (len(knots) - 2) + [4]
            or len(rows) != 3 * (len(knots) - 1) + 1):
        raise AssertionError(f'unsupported knot cover: {multiplicities}')
    if ([int(item) for item in items(surface[9])] != [2, 2]
            or [Fraction(item) for item in items(surface[11])] != [0, 1]):
        raise AssertionError('v direction is not one exact linear span')
    return {'ids': rows, 'knots': knots,
            'lower': [coordinates(es, row[0]) for row in rows],
            'upper': [coordinates(es, row[1]) for row in rows]}


def charts(es):
    return sorted((chart(es, surface) for surface in named(es, 'B_SPLINE_SURFACE_WITH_KNOTS').values()),
                  key=lambda value: value['knots'][0])


def polynomial(points):
    return [[p[0], 3 * (p[1] - p[0]), 3 * (p[2] - 2 * p[1] + p[0]),
             p[3] - 3 * p[2] + 3 * p[1] - p[0]] for p in zip(*points)]


def product(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def vector_area(rows):
    """1/2 integral(r cross r') over Bezier spans, exactly."""
    total = [Fraction(0)] * 3
    for start in range(0, len(rows) - 1, 3):
        p = polynomial(rows[start:start + 4])
        derivative = [[i * value for i, value in enumerate(row)][1:] for row in p]
        for axis in range(3):
            first = product(p[(axis + 1) % 3], derivative[(axis + 2) % 3])
            second = product(p[(axis + 2) % 3], derivative[(axis + 1) % 3])
            total[axis] += sum(((a - b) / (2 * (i + 1))
                                for i, (a, b) in enumerate(zip(first, second))), Fraction(0))
    return total


def to_decimal(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def norm(vector):
    return sum((value * value for value in vector), Decimal(0)).sqrt()


def length_interval(rows):
    """Chord and control-polygon sums of repeated half splits bracket the length."""
    segments = [tuple(tuple(to_decimal(value) for value in row) for row in rows[start:start + 4])
                for start in range(0, len(rows) - 1, 3)]
    width = Decimal(ORACLE['file_length_width_mm'])
    half = Decimal(1) / 2

    def middle(a, b):
        return tuple((x + y) * half for x, y in zip(a, b))

    def distance(a, b):
        return norm([x - y for x, y in zip(a, b)])

    for _ in range(18):
        lower = sum((distance(s[0], s[3]) for s in segments), Decimal(0))
        upper = sum((distance(s[0], s[1]) + distance(s[1], s[2]) + distance(s[2], s[3])
                     for s in segments), Decimal(0))
        if upper - lower <= width:
            return lower - ROUNDOFF, upper + ROUNDOFF
        finer = []
        for s in segments:
            a, b, c = middle(s[0], s[1]), middle(s[1], s[2]), middle(s[2], s[3])
            d, e = middle(a, b), middle(b, c)
            f = middle(d, e)
            finer.extend(((s[0], a, d, f), (f, e, c, s[3])))
        segments = finer
    raise AssertionError('profile length interval did not reach its width')


def file_oracle(sides):
    profile = sides[0]['lower'] + sides[1]['lower'][1:]
    sweep = tuple(b - a for a, b in zip(sides[0]['lower'][0], sides[0]['upper'][0]))
    for side in sides:
        for lower, upper in zip(side['lower'], side['upper']):
            if tuple(b - a for a, b in zip(lower, upper)) != sweep:
                raise AssertionError('a side row left the exact common sweep')
    area = vector_area(profile)
    with localcontext() as context:
        context.prec = 60
        height = norm([to_decimal(value) for value in sweep])
        lengths = [length_interval(side['lower']) for side in sides]
        return {
            'volume': sum((a * b for a, b in zip(area, sweep)), Fraction(0)),
            'cap_area': norm([to_decimal(value) for value in area]),
            'sides': [(height * lower, height * upper) for lower, upper in lengths],
        }


def reassembled(sides):
    """The one global spline the two charts claim to be, for the base fingerprint."""
    rows = [[lower, upper] for side in sides[:1]
            for lower, upper in zip(side['lower'], side['upper'])]
    rows += [[lower, upper] for lower, upper in zip(sides[1]['lower'][1:], sides[1]['upper'][1:])]
    knots = sides[0]['knots'] + sides[1]['knots'][1:]
    value = {'degree': [3, 1],
             'controls': [[[str(value) for value in point] for point in row] for row in rows],
             'multiplicities': [[4] + [3] * (len(knots) - 2) + [4], [2, 2]],
             'knots': [[str(knot) for knot in knots], ['0', '1']]}
    return canonical(value)


def canonical(value):
    text = json.dumps(value, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(text.encode()).hexdigest()


def volumes(shape):
    """Adaptive integration at each fixed epsilon, with spans enabled."""
    series = []
    for epsilon in EPSILONS:
        measured = GProp_GProps()
        error = BRepGProp.VolumePropertiesGK_s(shape, measured, epsilon, True, True, True,
                                               False, False)
        centre = measured.CentreOfMass()
        series.append((measured.Mass(), error, (centre.X(), centre.Y(), centre.Z())))
    return series


def areas(face):
    series = []
    for epsilon in EPSILONS:
        measured = GProp_GProps()
        error = BRepGProp.SurfaceProperties_s(face, measured, epsilon, False)
        series.append((measured.Mass(), error))
    return series


class CertifiedSplineSeams(unittest.TestCase):
    def assert_file_entities(self, es, manifest, sides):
        counts = CONTRACT['file']
        faces = named(es, 'ADVANCED_FACE')
        solid = named(es, 'MANIFOLD_SOLID_BREP')
        shell = named(es, 'CLOSED_SHELL')
        self.assertEqual((len(solid), len(shell)), (counts['solids'], counts['shells']))
        self.assertEqual(references(next(iter(solid.values()))[1]), list(shell))
        self.assertEqual(sorted(references(next(iter(shell.values()))[1])), sorted(faces))
        self.assertEqual(len(named(es, 'PLANE')), counts['planes'])
        self.assertEqual(len(sides), counts['splines'])
        self.assertEqual(len(named(es, 'EDGE_CURVE')), counts['edges'])
        self.assertEqual(len(named(es, 'VERTEX_POINT')), counts['vertices'])
        step = manifest['step']
        self.assertEqual(step['added_edges'], counts['added_edges'])
        self.assertEqual(step['added_vertices'], counts['added_vertices'])
        # Native topology stays three faces and two vertex-free rings.
        native = CONTRACT['native']
        result = manifest['result']
        self.assertEqual(result['faces_by_kind'], {'bspline': 1, 'plane': 2})
        self.assertEqual(result['faces'], native['faces'])
        self.assertEqual((result['edges'], result['vertices']),
                         (native['physical_rings'], native['vertices']))

    def assert_incidence(self, es):
        uses = {}
        faces = named(es, 'ADVANCED_FACE')
        for number, face in faces.items():
            for bound in references(face[1]):
                _, bound_fields = es[bound]
                for oriented in references(es[references(bound_fields[1])[0]][1][1]):
                    entry = es[oriented][1]
                    same = (entry[4] == '.T.') == (bound_fields[2] == '.T.')
                    uses.setdefault(references(entry[3])[0], []).append((number, same))
        self.assertEqual(sorted(uses), sorted(named(es, 'EDGE_CURVE')))
        for edge in uses.values():
            self.assertEqual(len(edge), 2)
            self.assertNotEqual(edge[0][0], edge[1][0], 'an edge used twice by one face')
            self.assertNotEqual(edge[0][1], edge[1][1], 'two same-direction uses of one edge')

    def assert_parameter_cover(self, es, sides):
        """Curves share the surfaces' own control rows, and pcurves trace each chart's rectangle."""
        domains = [[[str(Fraction(side['knots'][0])), str(Fraction(side['knots'][-1]))],
                    ['0', '1']] for side in sides]
        self.assertEqual(domains, CONTRACT['file_side_domains'])
        surfaces = {number: chart(es, body)
                    for number, body in named(es, 'B_SPLINE_SURFACE_WITH_KNOTS').items()}
        traced = {number: set() for number in surfaces}
        for wrapper in named(es, 'SURFACE_CURVE').values():
            curve = es[references(wrapper[1])[0]][1]
            ids = references(curve[2])
            self.assertEqual(wrapper[3], '.CURVE_3D.')
            self.assertEqual(curve[4:6], ['.F.', '.F.'])
            for pcurve in references(wrapper[2]):
                surface, representation = (references(es[pcurve][1][1])[0],
                                           references(es[pcurve][1][2])[0])
                line = es[references(es[representation][1][1])[0]][1]
                origin = coordinates(es, references(line[1])[0])
                vector = es[references(line[2])[0]][1]
                self.assertEqual(vector[2], '1.')
                direction = coordinates(es, references(vector[1])[0])
                traced[surface].add((origin, direction))
                rows = surfaces[surface]['ids']
                if direction == (1, 0):
                    self.assertEqual(ids, [row[int(origin[1])] for row in rows])
                    self.assertEqual(origin[0], 0)
                else:
                    self.assertEqual(direction, (0, 1))
                    self.assertEqual(origin[1], 0)
                    index = {surfaces[surface]['knots'][0]: 0,
                             surfaces[surface]['knots'][-1]: -1}[origin[0]]
                    self.assertEqual(ids, rows[index])
                    self.assertEqual([int(item) for item in items(curve[6])], [2, 2])
                    self.assertEqual(curve[1], '1')
        for number, seen in traced.items():
            low, high = surfaces[number]['knots'][0], surfaces[number]['knots'][-1]
            self.assertEqual(seen, {((Fraction(0), Fraction(0)), (1, 0)),
                                    ((Fraction(0), Fraction(1)), (1, 0)),
                                    ((low, Fraction(0)), (0, 1)),
                                    ((high, Fraction(0)), (0, 1))})
        for edge in named(es, 'EDGE_CURVE').values():
            for vertex in references(edge[1]) + references(edge[2]):
                self.assertEqual(es[vertex][0], 'VERTEX_POINT')
            self.assertEqual(edge[4], '.T.')

    def assert_closed_c1(self, sides):
        first, second = sides
        self.assertEqual(first['lower'][-1], second['lower'][0])
        self.assertEqual(second['lower'][-1], first['lower'][0])
        for before, after in ((first, second), (second, first)):
            left = before['knots'][-1] - before['knots'][-2]
            right = after['knots'][1] - after['knots'][0]
            incoming = tuple((a - b) / left
                             for a, b in zip(before['lower'][-1], before['lower'][-2]))
            outgoing = tuple((a - b) / right
                             for a, b in zip(after['lower'][1], after['lower'][0]))
            self.assertEqual(incoming, outgoing)

    def assert_regions(self, es, manifest):
        step = manifest['step']
        regions = step['regions']
        faces = named(es, 'ADVANCED_FACE')
        names = sorted(face[0].strip("'") for face in faces.values())
        self.assertEqual(sorted(region['file_region'] for region in regions), names)
        self.assertEqual(len(regions), 4)
        by_model = {}
        for region in regions:
            by_model.setdefault((region['face'], region['model_region']), []).append(region)
            self.assertEqual(region['file_to_model_bound_mm'], step['file_to_model_bound_mm'])
            self.assertEqual(region['file_to_definition_bound_mm'],
                             step['file_to_definition_bound_mm'])
        self.assertEqual(sorted(len(rows) for rows in by_model.values()), [1, 1, 2])
        for (face, _), rows in by_model.items():
            if len(rows) == 1:
                self.assertEqual(rows[0]['file_region'], f'face_{face}')
                self.assertNotIn('file_domain', rows[0])
            else:
                self.assertEqual([row['file_region'] for row in rows],
                                 [f'face_{face}_patch_{index}' for index in range(2)])
                self.assertEqual([row['file_domain'] for row in rows],
                                 CONTRACT['file_side_domains'])

    def assert_raw_topology(self, shape):
        counts = CONTRACT['file']
        self.assertEqual(len(unique(shape, TopAbs_SOLID)), counts['solids'])
        self.assertEqual(len(unique(shape, TopAbs_SHELL)), counts['shells'])
        faces = unique(shape, TopAbs_FACE)
        edges = unique(shape, TopAbs_EDGE)
        vertices = unique(shape, TopAbs_VERTEX)
        self.assertEqual((len(faces), len(edges), len(vertices)),
                         (counts['planes'] + counts['splines'], counts['edges'],
                          counts['vertices']))
        kinds = [BRepAdaptor_Surface(TopoDS.Face(face)).GetType() for face in faces]
        self.assertEqual(kinds.count(GeomAbs_Plane), counts['planes'])
        self.assertEqual(kinds.count(GeomAbs_BSplineSurface), counts['splines'])
        for face in faces:
            self.assertTrue(BRepCheck_Analyzer(face).IsValid())
        for edge in edges:
            signs = []
            for face in faces:
                explorer = TopExp_Explorer(face, TopAbs_EDGE)
                while explorer.More():
                    if explorer.Current().IsSame(edge):
                        signs.append(explorer.Current().Orientation() == TopAbs_FORWARD)
                    explorer.Next()
                self.assertFalse(BRep_Tool.IsClosed_s(TopoDS.Edge(edge), TopoDS.Face(face)),
                                 'a two-chart prism needs no self-seam')
            self.assertEqual(sorted(signs), [False, True])
            self.assertLessEqual(BRep_Tool.Tolerance_s(TopoDS.Edge(edge)),
                                 float(CONTRACT['raw']['max_tolerance']))
        for vertex in vertices:
            incident = [edge for edge in edges
                        if any(vertex.IsSame(end) for end in unique(edge, TopAbs_VERTEX))]
            adjacent = [{index for index, face in enumerate(faces)
                         if any(edge.IsSame(item) for item in unique(face, TopAbs_EDGE))}
                        for edge in incident]
            around = set().union(*adjacent)
            # Three edges and three faces, each face between exactly two edges.
            self.assertEqual((len(incident), len(around)), (3, 3))
            self.assertTrue(all(sum(index in pair for pair in adjacent) == 2 for index in around))

    def assert_converged(self, shape, file, midpoint):
        series = volumes(shape)
        for value, error, _ in series:
            self.assertTrue(isfinite(value) and value > 0 and error >= 0, series)
        self.assertLessEqual(abs(series[-1][0] - series[-2][0]),
                             float(ORACLE['last_two_volume_agreement_mm3']))
        self.assertLessEqual(abs(Fraction(series[-1][0]) - file['volume']),
                             Fraction(ORACLE['volume_numeric_allowance_mm3']))
        self.assertLessEqual(dist(series[-1][2], midpoint), 1e-8)
        allowance = Decimal(ORACLE['area_numeric_allowance_mm2'])
        for face in unique(shape, TopAbs_FACE):
            surface = BRepAdaptor_Surface(TopoDS.Face(face))
            if surface.GetType() == GeomAbs_Plane:
                bounds = (file['cap_area'] - ROUNDOFF, file['cap_area'] + ROUNDOFF)
            else:
                index = 0 if surface.FirstUParameter() < 1 else 1
                bounds = file['sides'][index]
            series = areas(face)
            for value, error in series:
                self.assertTrue(isfinite(value) and value > 0 and error >= 0, series)
            self.assertLessEqual(abs(series[-1][0] - series[-2][0]),
                                 float(ORACLE['last_two_area_agreement_mm2']))
            self.assertGreaterEqual(Decimal(series[-1][0]), bounds[0] - allowance)
            self.assertLessEqual(Decimal(series[-1][0]), bounds[1] + allowance)

    def assert_probes(self, shape, scenario):
        centre = tuple(Fraction(value) for value in scenario['centroid'])
        distance = Fraction(scenario['distance'])
        lower, upper = (Fraction(value) for value in scenario['axial'])
        margin = Fraction(1, 1000)

        def along(base, direction, length):
            return tuple(x + length * y for x, y in zip(base, direction))

        samples = [(centre, True)]
        for direction, semi_axis in ((MAJOR, 2 - distance), (MINOR, Fraction(8, 5) - distance)):
            for sign in (-1, 1):
                samples.append((along(centre, direction, sign * (semi_axis - margin)), True))
                samples.append((along(centre, direction, sign * (semi_axis + margin)), False))
        for level, inside in ((lower - margin, False), (lower + margin, True),
                              (upper - margin, True), (upper + margin, False)):
            samples.append((tuple(level * value for value in AXIS), inside))
        self.assertEqual(len(samples), 13)
        for point, inside in samples:
            state = BRepClass3d_SolidClassifier(
                shape, gp_Pnt(*(float(value) for value in point)), 1e-7).State()
            self.assertEqual(state, TopAbs_IN if inside else TopAbs_OUT, point)

    def assert_kernel_centroid(self, manifest, scenario):
        enclosure = manifest['properties']['values']['centroid']
        for value, interval in zip(scenario['centroid'], enclosure):
            self.assertLessEqual(Fraction(interval['lo']), Fraction(value))
            self.assertGreaterEqual(Fraction(interval['hi']), Fraction(value))

    def test_supported_inputs_write_one_raw_valid_covered_solid(self):
        inputs = problems()
        self.assertEqual(f"cadquery-ocp=={version('cadquery-ocp')}", CONTRACT['raw']['reader'])
        for scenario in CONTRACT['positive']:
            with self.subTest(scenario=scenario['id']), \
                    tempfile.TemporaryDirectory(prefix='catallax-spline-seam-') as temporary:
                root = Path(temporary)
                problem = inputs[scenario['id']]
                caps = CONTRACT['work_caps']
                writer = run(step_only(problem), root / 'writer', caps['writer_only_application'])
                full = run(problem, root / 'full', caps['replay'])
                for case in (writer, full):
                    self.assertFalse(case['response']['work']['exhausted'])
                    self.assertEqual(case['response']['manifest']['step']['status'], 'answered')
                path = Path(writer['path']) / 'answer.step'
                self.assertEqual(path.read_bytes(),
                                 (Path(full['path']) / 'answer.step').read_bytes())
                manifest = writer['response']['manifest']
                es = entities(path.read_text())
                sides = charts(es)
                self.assert_file_entities(es, manifest, sides)
                self.assert_incidence(es)
                self.assert_parameter_cover(es, sides)
                self.assert_closed_c1(sides)
                self.assert_regions(es, manifest)
                self.assertEqual(full['response']['manifest']['step'], manifest['step'])

                geometry, certificate = FINGERPRINTS[scenario['id']]
                self.assertEqual(reassembled(sides), geometry)
                self.assertEqual(canonical(manifest['certificate']), certificate)
                self.assertEqual(canonical(full['response']['manifest']['certificate']), certificate)

                lower, upper = (Fraction(value) for value in scenario['axial'])
                goal = target(Fraction(scenario['distance']), upper - lower)
                file = file_oracle(sides)
                inner, outer = volume_envelope(
                    goal, Fraction(manifest['step']['file_to_definition_bound_mm']))
                self.assertLessEqual(inner, file['volume'])
                self.assertLessEqual(file['volume'], outer)

                origins = [coordinates(es, references(es[references(plane[1])[0]][1][1])[0])
                           for plane in named(es, 'PLANE').values()]
                midpoint = tuple(float((a + b) / 2) for a, b in zip(*origins))

                shapes = raw_shapes(path)
                self.assertEqual(len(shapes), 1)
                raw = shapes[0]
                self.assertTrue(BRepCheck_Analyzer(raw).IsValid())
                self.assert_raw_topology(raw)
                self.assert_probes(raw, scenario)
                self.assert_kernel_centroid(full['response']['manifest'], scenario)
                ordinary = ordinary_shape(path)
                self.assertTrue(BRepCheck_Analyzer(ordinary).IsValid())
                for shape in (raw, ordinary):
                    self.assert_converged(shape, file, midpoint)


if __name__ == '__main__':
    unittest.main()
