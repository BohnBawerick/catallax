"""Independent standard-library oracle. Does not load any kernel implementation."""
from fractions import Fraction as Q

# STEP case 1: geometry children are private; fixed #1..#20 are included.
counts = dict(MANIFOLD_SOLID_BREP=1, CLOSED_SHELL=1, ADVANCED_FACE=6,
              FACE_BOUND=6, EDGE_LOOP=6, ORIENTED_EDGE=24, EDGE_CURVE=12,
              VERTEX_POINT=8, CARTESIAN_POINT=1+8+12+6, LINE=12, VECTOR=12,
              DIRECTION=2+12+2*6, AXIS2_PLACEMENT_3D=1+6, PLANE=6)
# Fixed #17..#20 have counted types; sixteen other fixed entities remain.
assert sum(counts.values()) + 16 == 170
assert 3 * (Q(1, 3) - Q(33333, 100000)) ** 2 == Q(1, 30000000000)
assert Q(10**30, 3) - 333333333333333330000000000000 == Q(10**13, 3)

# Alternating arctan series: two adjacent partial sums enclose atan(3/4).
x = Q(3, 4)
s = sum((-1)**k * x**(2*k+1) / (2*k+1) for k in range(200))
t = s + x**401 / 401
assert Q('0.643501108793284386802809228717') < s < t < Q('0.643501108793284386802809228718')
assert Q('0.64350110879328438') < s < t < Q('0.64350110879328439')
assert Q('1.4142135623730950')**2 < 2 < Q('1.4142135623730951')**2

# Pose q=(2,1,0,0), using the homogeneous quaternion formula.
w, x, y, z = map(Q, [2, 1, 0, 0])
d = w*w+x*x+y*y+z*z
R = [[(w*w+x*x-y*y-z*z)/d, 2*(x*y-w*z)/d, 2*(x*z+w*y)/d],
     [2*(x*y+w*z)/d, (w*w-x*x+y*y-z*z)/d, 2*(y*z-w*x)/d],
     [2*(x*z-w*y)/d, 2*(y*z+w*x)/d, (w*w-x*x-y*y+z*z)/d]]
assert R == [[1, 0, 0], [0, Q(3,5), Q(-4,5)], [0, Q(4,5), Q(3,5)]]
for i in range(3):
    for j in range(3):
        assert sum(R[k][i]*R[k][j] for k in range(3)) == int(i == j)

# Intersections: cylinder/sphere centres and near-tangent lines.
assert 6**2 - 3**2 == 27
e = Q(1, 10**30)
assert 2**2 - (2-e)**2 == 4*e-e**2
# Branch 17 lies on x²+y²=4 and y²+z²=1.
assert 2**2+0**2 == 4 and 0**2+1**2 == 1
# Branch 25 at x=1,y=0 has z²=3 on the sphere centred at (2,0,0).
assert (1-2)**2 + 3 == 2**2
print('verified: entity count, rounding bounds, arctan enclosure, sqrt(2), quaternion, intersection equations')
