from pathlib import Path
import tempfile
import unittest

from step_pcurve import affine_pcurves_only


def step(degree, *, v='0.', end='2.', radius='1.'):
    points = '#7,#8' if degree == 1 else '#7,#9,#9,#9,#9,#9,#9,#8'
    mult = '2,2' if degree == 1 else '8,8'
    interior = '' if degree == 1 else "#9=CARTESIAN_POINT('',(1.,0.));\n"
    return f"""ISO-10303-21;
DATA;
#1=CYLINDRICAL_SURFACE('',#2,{radius});
#2=AXIS2_PLACEMENT_3D('',#3,$,$);
#3=CARTESIAN_POINT('',(0.,0.,0.));
#4=PCURVE('',#1,#5);
#5=DEFINITIONAL_REPRESENTATION('',(#6),#10);
#6=B_SPLINE_CURVE_WITH_KNOTS('',{degree},({points}),.UNSPECIFIED.,.F.,.F.,({mult}),(0.,{end}),.UNSPECIFIED.);
#7=CARTESIAN_POINT('',(0.,0.));
#8=CARTESIAN_POINT('',(2.,{v}));
{interior}#10=REPRESENTATION_CONTEXT('','');
ENDSEC;
END-ISO-10303-21;
"""


class AffinePcurveTests(unittest.TestCase):
    def test_only_one_affine_replacement_is_allowed(self):
        with tempfile.TemporaryDirectory() as temp:
            old, new = [Path(temp) / file for file in ('old.step', 'new.step')]
            old.write_text(step(7))
            for text, expected in (
                (step(1), True),
                (step(7), False),
                (step(1, v='0.1'), False),
                (step(1, end='3.'), False),
                (step(1, radius='2.'), False),
                (step(1).replace('(2,2)', '(1,3)'), False),
                (step(1).replace('(#6)', '(#99)'), False),
            ):
                with self.subTest(candidate=text):
                    new.write_text(text)
                    self.assertEqual(affine_pcurves_only(old, new), expected)
