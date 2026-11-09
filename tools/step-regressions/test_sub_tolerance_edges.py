"""Collapse a redundant short edge; refuse a slab that would lose its volume."""
import os
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID, TopAbs_VERTEX

from test_loopless_faces import ROOT, raw_solid, unique

BEFORE = Path(__file__).parent / 'fixtures/sub-tolerance'


class SubToleranceEdges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='catallax-sub-tolerance-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        subprocess.run(
            ['cargo', 'run', '--locked', '--release', '-p', 'catallax-step',
             '--example', 'sub_tolerance', '--', str(cls.directory)],
            cwd=ROOT, env={**os.environ, 'CARGO_BUILD_JOBS': '2'}, check=True,
            capture_output=True, text=True)

    def test_redundant_split_edge_collapses_and_reads_as_a_cube(self):
        self.assertFalse(BRepCheck_Analyzer(raw_solid(BEFORE / 'split-before.step')).IsValid())
        path = self.directory / 'split.step'
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        for mode, shape in [('raw', raw_solid(path)), ('standard', reader.OneShape())]:
            with self.subTest(reader=mode):
                self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                                       (TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)),
                                 (1, 6, 12, 8))
                properties = GProp_GProps()
                BRepGProp.VolumeProperties_s(shape, properties)
                self.assertAlmostEqual(properties.Mass(), 1000, delta=1e-9)

    def test_a_slab_that_would_collapse_to_a_sheet_writes_nothing(self):
        self.assertFalse(BRepCheck_Analyzer(raw_solid(BEFORE / 'thin-before.step')).IsValid())
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(BEFORE / 'thin-before.step')), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        self.assertFalse(BRepCheck_Analyzer(reader.OneShape()).IsValid())
        self.assertFalse((self.directory / 'thin.step').exists())
        refusal = (self.directory / 'thin.refusal').read_text()
        self.assertTrue(refusal.startswith('step.sub-tolerance-edge:'), refusal)

    def test_a_readable_subdivision_remains_in_the_file(self):
        shape = raw_solid(self.directory / 'readable-split.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        self.assertEqual(len(unique(shape, TopAbs_EDGE)), 13)
        self.assertEqual(len(unique(shape, TopAbs_VERTEX)), 9)

    def test_a_thin_cylinder_refuses_its_generated_line_seam(self):
        path = BEFORE / 'thin-cylinder-before.step'
        self.assertFalse(BRepCheck_Analyzer(raw_solid(path)).IsValid())
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        self.assertFalse(BRepCheck_Analyzer(reader.OneShape()).IsValid())
        self.assertFalse((self.directory / 'thin-cylinder.step').exists())
        self.assertTrue((self.directory / 'thin-cylinder.refusal').read_text()
                        .startswith('step.sub-tolerance-edge:'))

    def test_a_readable_cylinder_seam_stays_in_the_file(self):
        path = self.directory / 'readable-cylinder.step'
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        for shape in (raw_solid(path), reader.OneShape()):
            self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
            self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
            self.assertEqual(len(unique(shape, TopAbs_EDGE)), 3)
            properties = GProp_GProps()
            BRepGProp.VolumeProperties_s(shape, properties)
            self.assertAlmostEqual(properties.Mass(), math.pi * 2e-7, delta=1e-14)


if __name__ == '__main__':
    unittest.main()
