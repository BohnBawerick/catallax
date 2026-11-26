"""Compare Catallax STEP files whose only change is affine chart boundary encoding."""

from collections import defaultdict
import re

ENTITY = re.compile(rb'^#(\d+)=(.*);$')
TOKEN = re.compile(rb"'(?:[^']|'')*'|#(\d+)")
PCURVE = re.compile(rb"PCURVE\('',#(\d+),#(\d+)\)")
REPRESENTATION = re.compile(rb"DEFINITIONAL_REPRESENTATION\('',\(#(\d+)\),#(\d+)\)")
SPLINE = re.compile(rb"B_SPLINE_CURVE_WITH_KNOTS\('',(\d+),\((#[\d,#]+)\),"
                    rb"\.UNSPECIFIED\.,(\.F\.|\.T\.),\.F\.,\(([^()]+)\),\(([^()]+)\),"
                    rb"\.UNSPECIFIED\.\)")
POINT = re.compile(rb"CARTESIAN_POINT\('',\(([^,]+),([^,]+)\)\)")


def projection(path):
    outside = []
    entities = {}
    incoming = defaultdict(set)
    for line in path.read_bytes().splitlines():
        match = ENTITY.fullmatch(line)
        if not match:
            outside.append(line)
            continue
        key, text = match.groups()
        if key in entities:
            raise ValueError('duplicate STEP entity')
        entities[key] = text
        for ref in TOKEN.finditer(text):
            if ref[1] is not None:
                incoming[ref[1]].add(key)
    curves, curve_ids, controls = {}, set(), set()
    for key, text in entities.items():
        pcurve = PCURVE.fullmatch(text)
        if not pcurve:
            continue
        rep_id = pcurve[2]
        rep = REPRESENTATION.fullmatch(entities[rep_id])
        if not rep:
            continue
        curve_id, context = rep.groups()
        spline = SPLINE.fullmatch(entities[curve_id])
        if not spline or incoming[curve_id] != {rep_id} or incoming[rep_id] != {key}:
            continue
        degree, point_list, closed, multiplicities, knots = spline.groups()
        point_ids = [point[1:] for point in point_list.split(b',')]
        coords = [POINT.fullmatch(entities[point]) for point in point_ids]
        if not all(coords):
            continue
        curves[rep_id] = (int(degree), [point.groups() for point in coords], closed,
                          multiplicities.split(b','), knots.split(b','))
        curve_ids.add(curve_id)
        controls.update(point_ids)
        entities[rep_id] = b"DEFINITIONAL_REPRESENTATION('',(UV_CURVE),#" + context + b')'
    hidden = curve_ids | {point for point in controls if incoming[point] <= curve_ids}
    visible = [key for key in entities if key not in hidden]
    ids = {key: str(i).encode() for i, key in enumerate(visible)}
    graph = [TOKEN.sub(lambda m: b'#' + ids[m[1]] if m[1] else m[0], entities[key])
             for key in visible]
    return outside, graph, {ids[key]: value for key, value in curves.items()}


def affine_pcurves_only(before, after):
    """Require identical written non-UV entities, endpoints and parameter limits."""
    try:
        h1, g1, curve1 = projection(before)
        h2, g2, curve2 = projection(after)
    except (OSError, KeyError, ValueError):
        return False
    if (h1, g1, curve1.keys()) != (h2, g2, curve2.keys()):
        return False
    changed = False
    for key, old in curve1.items():
        new = curve2[key]
        if old == new:
            continue
        degree, points, closed, mults, knots = new
        if old[0] != 7 or degree != 1 or len(points) != 2 or mults != [b'2', b'2']:
            return False
        if (len(knots) != 2 or closed != old[2] or points != [old[1][0], old[1][-1]]
                or knots != [old[4][0], old[4][-1]]):
            return False
        changed = True
    return changed
