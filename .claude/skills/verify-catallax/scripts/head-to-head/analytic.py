"""Independent scalar definitions. No CAD kernel imports or exported-body references.

Property formulas cover eight target families. Four intersecting curved P0
families have analytic membership but no property oracle here. P5's junction
fillet reaches the countersink openings; a full torus is not an exact reference
for that trimmed blend. Its property and membership checks remain unmeasured.
Simpson quadrature integrates smooth elementary fillet sections, not meshes,
with two resolutions checked against each other.
"""
import math
from pathlib import Path

PI = math.pi


def integral(fn, lo=0, hi=PI/2):
    def simpson(n):
        h = (hi-lo)/n
        return h/3 * (fn(lo)+fn(hi)+sum((4 if i % 2 else 2)*fn(lo+i*h) for i in range(1, n)))
    a, b = simpson(1024), simpson(2048)
    if abs(a-b) > max(1e-10, abs(b)*1e-12):
        raise ValueError('analytic section quadrature did not converge')
    return b


def fillet_add(radius, blend, start):
    def rho(t):
        return radius+blend-blend*math.sin(t)
    def section(t):
        return PI*(rho(t)**2-radius**2)*blend*math.sin(t)
    volume = integral(section)
    moment = integral(lambda t: section(t)*(start+blend-blend*math.cos(t)))
    area = 2*PI*blend*integral(rho)-PI*((radius+blend)**2-radius**2)-2*PI*radius*blend
    return volume, moment, area


def properties(name, edit):
    family = Path(name).stem
    if family in ('near-miss', 'line-tangent', 'seam'):
        h, r = edit or 10, 5
        distance = 4 if family == 'seam' else 10
        overlap = 2*r*r*math.acos(distance/(2*r))-distance*math.sqrt(4*r*r-distance*distance)/2
        union_area = 2*PI*r*r-overlap
        union_perimeter = 4*r*(PI-math.acos(distance/(2*r)))
        extra = h-10
        v = union_area*10+PI*r*r*extra
        moment = union_area*50+PI*r*r*extra*(10+h)/2
        # The two end planes together expose exactly twice the union section.
        a = 2*union_area+union_perimeter*10+2*PI*r*extra
        x = (distance*union_area*5)/v
        return {'volume': v, 'area': a, 'centroid': [x, 0, moment/v]}
    if family == 'equal-cylinder-sphere':
        h = edit or 20
        upper = -10+h
        cap_v = PI*(2*5**3/3 - (25*upper-upper**3/3)) if upper < 5 else 0
        cap_m = PI*(625/4-(25*upper**2/2-upper**4/4)) if upper < 5 else 0
        v = 25*PI*h+cap_v
        area = 10*PI*h+50*PI
        if upper < 5:
            area += 10*PI*(5-upper)-PI*(25-upper**2)
        return {'volume': v, 'area': area, 'centroid': [0, 0, (25*PI*h*(-10+h/2)+cap_m)/v]}
    if name.startswith('P1-'):
        r = edit or 8
        slot = 5*(84+9*PI)
        pieces = [(12000, (30, 20, 2.5)), (7500, (30, 2.5, 17.5)),
                  (10*PI*r*r, (30, 25, 10)), (-240*PI, (30, 25, 7.5)),
                  (-80*PI, (30, 22, 2.5)), (-slot, (30, 2.5, 18))]
        v = sum(p[0] for p in pieces)
        return {'volume': v, 'area': 9022+PI*(20*r+148),
                'centroid': [sum(m*c[i] for m, c in pieces)/v for i in range(3)]}
    if name.startswith('P2-'):
        c = edit or 1
        v = PI*(1296-24*c*c+2*c**3/3)
        a = 24*PI*(12-2*c)+144*PI+2*PI*((12-c)**2-36)+2*PI*(24-c)*math.sqrt(2)*c
        return {'volume': v, 'area': a, 'centroid': [0, 6, 0]}
    if name.startswith('P3-'):
        def section(t):
            d = 1-math.cos(t)
            return (40-2*d)*(30-2*d)
        v = 1200*9+integral(lambda t: section(t)*math.cos(t))
        m = 1200*81/2+integral(lambda t: (9+math.sin(t))*section(t)*math.cos(t))
        a = 1200+38*28+140*9+integral(lambda t: 140-8*(1-math.cos(t)))
        dv, dm, da = fillet_add(8, edit or 2, 10)
        v += PI*64*12+dv
        m += PI*64*12*16+dm
        a += 2*PI*8*12+da
        return {'volume': v, 'area': a, 'centroid': [20, 15, m/v]}
    if name.startswith('P4-'):
        t, r = edit or 2, 23
        x0, y0 = math.sqrt(r*r-15*15), math.sqrt(r*r-20*20)
        def primitive(x):
            return (x*math.sqrt(r*r-x*x)+r*r*math.asin(x/r))/2
        outer = 1200-4*(15*(20-x0)-primitive(20)+primitive(x0))
        perimeter = 4*(x0+y0+r*(math.asin(20/r)-math.asin(x0/r)))
        inner = (40-2*t)*(30-2*t)
        v = 12*outer-(12-t)*inner
        m = 72*outer-(144-t*t)*inner/2
        return {'volume': v, 'area': 2*outer+12*perimeter+2*(70-4*t)*(12-t), 'centroid': [0, 0, m/v]}
    return None


def definition(case, label):
    name = case['recipe']
    family = Path(name).stem
    edit = float(case['edit']['value']) if label == 'edit' else None
    def box(p, lo, hi):
        return max(*(lo[i]-p[i] for i in range(3)), *(p[i]-hi[i] for i in range(3)))
    def cylinder(p, r, h, base=(0, 0, 0), axis=(0, 0, 1)):
        q = [p[i]-base[i] for i in range(3)]
        z = sum(q[i]*axis[i] for i in range(3))
        radial = math.sqrt(sum((q[i]-z*axis[i])**2 for i in range(3)))
        return max(radial-r, -z, z-h)
    def sphere(p, r, center=(0, 0, 0)):
        return math.dist(p, center)-r
    def cone(p, h):
        z = p[2]
        return max(math.hypot(p[0], p[1])-(5-3*z/h), -z, z-h)
    def torus(p, major):
        return math.hypot(math.hypot(p[0], p[1])-major, p[2])-2
    def junction(p, radius, blend, start, axis=2, center=(0, 0)):
        radial_axes = [i for i in range(3) if i != axis]
        rho = math.hypot(p[radial_axes[0]]-center[0], p[radial_axes[1]]-center[1])
        z = p[axis]-start
        if z < 0 or z > blend:
            return max(-z, z-blend, rho-radius)
        outer = radius+blend-math.sqrt(max(0, blend*blend-(z-blend)**2))
        return max(rho-outer, -z, z-blend)
    if name.startswith('P0/'):
        h = edit or (20 if family in ('crossing-pipes', 'equal-cylinder-sphere') else 10)
        if family == 'block-sphere':
            fn = lambda p: min(box(p, (0, 0, 0), (edit or 10, 10, 10)), sphere(p, 6, (10, 5, 5)))
            bounds = (0, -1, -1, 16, 11, 11)
        elif family == 'crossing-pipes':
            fn = lambda p: min(max(cylinder(p, 5, h, (0, 0, -10)), -cylinder(p, 3, 22, (0, 0, -11))), max(cylinder(p, 5, 20, (-10, 0, 0), (1, 0, 0)), -cylinder(p, 3, 22, (-11, 0, 0), (1, 0, 0))))
            bounds = (-10, -5, -10, 10, 5, max(5, h-10))
        elif family == 'equal-cylinder-sphere':
            fn = lambda p: min(cylinder(p, 5, h, (0, 0, -10)), sphere(p, 5))
            bounds = (-5, -5, -10, 5, 5, max(5, h-10))
        elif family in ('line-tangent', 'near-miss', 'seam'):
            x = {'line-tangent': 10, 'near-miss': 10.0000001, 'seam': 4}[family]
            fn = lambda p: min(cylinder(p, 5, h), cylinder(p, 5, 10, (x, 0, 0)))
            bounds = (-5, -5, 0, x+5, 5, max(h, 10))
        else:
            shift = (2, 1, 1) if family == 'general-cone' else (1, 2, 1)
            # Inverse of the exact normalized homogeneous quaternion (2,1,1,0).
            matrix = ((2/3, 1/3, 2/3), (1/3, 2/3, -2/3), (-2/3, 2/3, 1/3))
            def fn(p):
                q = [p[i]-shift[i] for i in range(3)]
                inv = [sum(matrix[j][i]*q[j] for j in range(3)) for i in range(3)]
                return min(cone(p, h), cone(inv, 10)) if family == 'general-cone' else min(torus(p, edit or 8), torus(inv, 8))
            bounds = (-15, -15, -15, 15, 15, 15)
    elif name.startswith('P1-'):
        def fn(p):
            value = min(box(p, (0, 0, 0), (60, 40, 5)), box(p, (0, 0, 5), (60, 5, 30)), cylinder(p, edit or 8, 10, (30, 25, 5)))
            value = max(value, -cylinder(p, 4, 17, (30, 25, -1)))
            for x in (10, 50):
                for y in (12, 32):
                    value = max(value, -cylinder(p, 2, 7, (x, y, -1)))
            slot = min(cylinder(p, 3, 7, (23, 6, 18), (0, -1, 0)), cylinder(p, 3, 7, (37, 6, 18), (0, -1, 0)), box(p, (23, -1, 15), (37, 6, 21)))
            return max(value, -slot)
        bounds = (0, 0, 0, 60, 40, 30)
    elif name.startswith('P2-'):
        c = edit or 1
        fn = lambda p: max(-p[1], p[1]-12, 6-math.hypot(p[0], p[2]), math.hypot(p[0], p[2])-12+max(0, c-p[1], c-(12-p[1])))
        bounds = (-12, 0, -12, 12, 12, 12)
    elif name.startswith('P3-'):
        def fn(p):
            d = 1-math.sqrt(max(0, 1-(p[2]-9)**2)) if 9 <= p[2] <= 10 else 0
            block = box(p, (d, d, 0), (40-d, 30-d, 10))
            return min(block, cylinder(p, 8, 12, (20, 15, 10)), junction(p, 8, edit or 2, 10, center=(20, 15)))
        bounds = (0, 0, 0, 40, 30, 22)
    elif name.startswith('P4-'):
        t = edit or 2
        fn = lambda p: max(box(p, (-20, -15, 0), (20, 15, 12)), -box(p, (-20+t, -15+t, t), (20-t, 15-t, 25)), cylinder(p, 23, 14, (0, 0, -1)))
        bounds = (-20, -15, 0, 20, 15, 12)
    elif name.startswith('P5-'):
        # The rolling blend meets the countersink openings at its flange contact.
        # Do not substitute an untrimmed torus or a successful unfilleted prefix.
        fn = None
        bounds = (-22, 0, -22, 22, 30, 22)
    else:
        raise ValueError(name)
    return properties(name, edit), bounds, fn
