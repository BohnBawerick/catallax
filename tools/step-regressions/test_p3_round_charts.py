"""Public P3 block rounds and junction torus; raw translation bypasses caller-side processing."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Torus
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_OUT, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_loopless_faces import BINARY, DRIVER, ROOT, raw_solid, unique
from test_partial_periodic_trims import converged_volume

# Independent OpenCASCADE construction from the dimensions, recorded in
# docs/milestones/2026-10-27-p3-material.md for both junction radii.
VOLUMES = {'2': 14428.640126047918, '5/2': 14455.207561827563}


class P3RoundCharts(unittest.TestCase):
    def export(self, root, radius):
        recipe = json.loads((ROOT / 'parts/P3-filleted-boss.json').read_text())
        for step in recipe['steps']:
            if step['id'] == 'junction':
                step['radius'] = radius
        source = root / 'recipe.json'
        source.write_text(json.dumps(recipe))
        out = root / 'run'
        completed = subprocess.run(
            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary', str(BINARY),
             '--out', str(out), '--fresh', '--workers', '1', '--json'],
            capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        manifest = case['response']['manifest']
        self.assertTrue(manifest['exact'], manifest)
        self.assertEqual(manifest['marks'], [])
        # One STEP face per source face: the two complete bands gain one seam
        # each and the three rings one origin vertex each; nothing is split.
        self.assertEqual(manifest['result']['faces'], 13)
        self.assertEqual(manifest['step']['added_edges'], 2)
        self.assertEqual(manifest['step']['added_vertices'], 3)
        return Path(case['path']) / 'answer.step'

    def test_original_and_edited_junction_import_without_repair(self):
        # The edit moves the junction's outer radius at z=11 from 10-sqrt(3)
        # to 17/2, so this point changes from outside to inside.
        probes = [((20, 15, 21), TopAbs_IN, TopAbs_IN), ((20, 15, 23), TopAbs_OUT, TopAbs_OUT),
                  ((28.4, 15, 11), TopAbs_OUT, TopAbs_IN)]
        for edited, (radius, expected_volume) in enumerate(VOLUMES.items()):
            with self.subTest(radius=radius), \
                    tempfile.TemporaryDirectory(prefix='catallax-p3-charts-') as temporary:
                path = self.export(Path(temporary), radius)
                raw = raw_solid(path)
                reader = STEPControl_Reader()
                self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
                self.assertGreater(reader.TransferRoots(), 0)
                # Only raw_solid bypasses caller-side processing; internal
                # repair inside the reader remains unproven.
                for mode, shape in [('raw', raw), ('standard', reader.OneShape())]:
                    with self.subTest(reader=mode):
                        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                        self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                        self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                                               (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)),
                                         (13, 25, 15))
                        kinds = []
                        for face in unique(shape, TopAbs_FACE):
                            face = TopoDS.Face(face)
                            kinds.append(int(GeomAdaptor_Surface(BRep_Tool.Surface_s(face)).GetType()))
                            self.assertTrue(BRepCheck_Analyzer(face).IsValid())
                            # Finite trims, not the legacy unbounded cylinder V
                            # domain or full periodic torus chart.
                            u0, u1, v0, v1 = BRepTools.UVBounds_s(face)
                            self.assertLess(max(abs(u0), abs(u1), abs(v0), abs(v1)), 100)
                        self.assertEqual(sorted(kinds), sorted(
                            [int(GeomAbs_Plane)] * 7 + [int(GeomAbs_Cylinder)] * 5
                            + [int(GeomAbs_Torus)]))
                        volume, spread = converged_volume(shape)
                        area = GProp_GProps()
                        BRepGProp.SurfaceProperties_s(shape, area, 1e-12, False)
                        # Diagnostic tube comparison, not a certification.
                        self.assertLessEqual(abs(volume - expected_volume),
                                             area.Mass() * 1e-7 + spread)
                        for point, *states in probes:
                            classifier = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7)
                            self.assertEqual(classifier.State(), states[edited], point)


if __name__ == '__main__':
    unittest.main()
