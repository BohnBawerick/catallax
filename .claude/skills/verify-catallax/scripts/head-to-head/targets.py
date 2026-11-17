"""Authored target recipes in millimetres, evaluated by each external tool adapter.

These define whole targets, not successful prefixes. Picks are points on the
specified edges. A missing operation must raise NotImplementedError.
"""
from pathlib import Path


def rounded_rectangle(hw, hh, r):
    """Lines and quarter arcs counterclockwise from the right side's lower end, as the
    P9 recipe's first section: ('line', a, b) or ('arc', a, midpoint, b) in (x, y)."""
    m = r * (1 - 2 ** -0.5)
    return [('line', (hw, -hh + r), (hw, hh - r)), ('arc', (hw, hh - r), (hw - m, hh - m), (hw - r, hh)),
            ('line', (hw - r, hh), (-hw + r, hh)), ('arc', (-hw + r, hh), (-hw + m, hh - m), (-hw, hh - r)),
            ('line', (-hw, hh - r), (-hw, -hh + r)), ('arc', (-hw, -hh + r), (-hw + m, -hh + m), (-hw + r, -hh)),
            ('line', (-hw + r, -hh), (hw - r, -hh)), ('arc', (hw - r, -hh), (hw - m, -hh + m), (hw, -hh + r))]


def build(api, name, edit=None, case=None):
    box, cyl, sphere = api.box, api.cylinder, api.sphere
    union, cut, common = api.union, api.cut, api.common
    blend = api.blend
    if name.startswith('census/'):
        # Real part: the tool's own STEP import, then one drilled hole whose radius is the edit.
        cutter = case['census']['cutter']
        part = api.import_step(case['census']['source'])
        return cut(part, cyl(edit or float(cutter['radius']), float(cutter['height']),
                             tuple(map(float, cutter['base'])), tuple(cutter['axis'])))
    if name.startswith('phase2/'):
        # Phase-2 checkpoint targets; tools/checkpoint/phase2/README.md gives each definition.
        family = Path(name).stem
        if family == 'P8-spring-threaded-boss':
            part = union(box((50, 50, 4), (-25, -25, 0)), cyl(8, 20, (0, 0, 4)))
            part = union(part, api.coil(8, 7, 2880, 12, 0.6))
            return union(part, api.coil(17, 2.5, 1800, edit or 15, 1))
        if family == 'P9-lofted-duct':
            # Wall = outer loft minus an inner loft inset by 2, 1 mm past each end. The
            # edit opens the inner outlet radius, thinning the round end's wall.
            outer = api.loft([('segments', rounded_rectangle(20, 10, 5), 0), ('circle', 12, 40)])
            inner = api.loft([('segments', rounded_rectangle(18, 8, 3), -1), ('circle', edit or 10, 41)])
            return cut(outer, inner)
        if family == 'curved-housing':
            housing = api.loft([('circle', 10, 0), ('circle', edit or 14, 15), ('circle', 10, 30)])
            tool = api.bezier_prism([(-20, 34), (-10, 20), (10, 20), (20, 34)], 40, 20, -40)
            return api.blend_above(cut(housing, tool), 1, 2)
        raise ValueError(name)
    if name.startswith('P0/'):
        family = Path(name).stem
        height = edit if family not in ('block-sphere', 'general-torus') else None
        if family == 'block-sphere':
            return union(box((edit or 10, 10, 10)), sphere(6, (10, 5, 5)))
        if family == 'crossing-pipes':
            a = cut(cyl(5, height or 20, (0, 0, -10)), cyl(3, 22, (0, 0, -11)))
            b = cut(cyl(5, 20, (-10, 0, 0), (1, 0, 0)), cyl(3, 22, (-11, 0, 0), (1, 0, 0)))
            return union(a, b)
        if family == 'equal-cylinder-sphere':
            return union(cyl(5, height or 20, (0, 0, -10)), sphere(5))
        if family in ('line-tangent', 'near-miss', 'seam'):
            x = {'line-tangent': 10, 'near-miss': 10.0000001, 'seam': 4}[family]
            return union(cyl(5, height or 10), cyl(5, 10, (x, 0, 0)))
        if family in ('general-cone', 'general-torus'):
            if family == 'general-cone':
                a, b = api.cone(5, 2, height or 10), api.cone(5, 2, 10)
                shift = (2, 1, 1)
            else:
                a, b = api.torus(edit or 8, 2), api.torus(8, 2)
                shift = (1, 2, 1)
            return union(a, api.rotate_translate(b, (2, 1, 1, 0), shift))
        raise ValueError(name)
    if name.startswith('P1-'):
        part = union(box((60, 40, 5)), box((60, 5, 25), (0, 0, 5)))
        part = union(part, cyl(edit or 8, 10, (30, 25, 5)))
        part = cut(part, cyl(4, 17, (30, 25, -1)))
        for x in (10, 50):
            for y in (12, 32):
                part = cut(part, cyl(2, 7, (x, y, -1)))
        # Exact capsule: the union of the original overlapping slot cutters.
        slot = union(cyl(3, 7, (23, 6, 18), (0, -1, 0)), cyl(3, 7, (37, 6, 18), (0, -1, 0)))
        slot = union(slot, box((14, 7, 6), (23, -1, 15)))
        return cut(part, slot)
    if name.startswith('P2-'):
        part = cut(cyl(12, 12, axis=(0, 1, 0)), cyl(6, 12, axis=(0, 1, 0)))
        return blend(part, [(12, 0, 0), (12, 12, 0)], edit or 1, chamfer=True, curve='circle')
    if name.startswith('P3-'):
        part = union(box((40, 30, 10)), cyl(8, 12, (20, 15, 10)))
        part = blend(part, [(28, 15, 10)], edit or 2, curve='circle')
        return blend(part, [(20, 0, 10), (40, 15, 10), (20, 30, 10), (0, 15, 10)], 1)
    if name.startswith('P4-'):
        part = api.shell_top(box((40, 30, 12), (-20, -15, 0)), edit or 2)
        return common(part, cyl(23, 14, (0, 0, -1)))
    if name.startswith('P5-'):
        part = union(cyl(22, 5, axis=(0, 1, 0)), cyl(12, 30, axis=(0, 1, 0)))
        part = cut(part, cyl(8, 30, axis=(0, 1, 0)))
        for x, z in ((17, 0), (0, -17), (-17, 0), (0, 17)):
            tool = union(cyl(2, 7, (x, 6, z), (0, -1, 0)), api.cone(4, 2, 2, (x, 5, z), (0, -1, 0)))
            part = cut(part, tool)
        return blend(part, [(12, 5, 0)], edit or 2, curve='circle')
    raise ValueError(f'unknown target {name}')
