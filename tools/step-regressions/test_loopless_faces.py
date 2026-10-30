"""Import complete periodic faces without ShapeFix, sewing or SameParameter."""
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import OCP
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.StepShape import StepShape_BrepWithVoids, StepShape_ManifoldSolidBrep
from OCP.StepToTopoDS import StepToTopoDS_Builder
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.Transfer import Transfer_TransientProcess

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'


def unique(shape, kind):
    shapes = []
    explorer = TopExp_Explorer(shape, kind)
    while explorer.More():
        item = explorer.Current()
        if not any(item.IsSame(previous) for previous in shapes):
            shapes.append(item)
        explorer.Next()
    return shapes


def raw_solid(path):
    solids = raw_solids(path)
    if len(solids) != 1:
        raise AssertionError(f'expected one STEP solid, got {len(solids)}')
    return solids[0]


def raw_solids(path):
    reader = STEPControl_Reader()
    # Translate each STEP solid directly, without TransferRoots or ShapeFix.
    settings = {'read.stdsameparameter.mode': 0, 'read.maxprecision.mode': 0}
    previous = {name: Interface_Static.IVal_s(name) for name in settings}
    try:
        for name, value in settings.items():
            if not Interface_Static.SetIVal_s(name, value):
                raise AssertionError(f'cannot disable {name}')
        if reader.ReadFile(str(path)) != IFSelect_RetDone:
            raise AssertionError(f'cannot read {path}')
        model = reader.StepModel()
        solids = [model.Value(i) for i in range(1, model.NbEntities() + 1)
                  if isinstance(model.Value(i), StepShape_ManifoldSolidBrep)]
        if not solids:
            raise AssertionError(f'no STEP solids in {path}')
        shapes = []
        for solid in solids:
            process = Transfer_TransientProcess()
            process.SetModel(model)
            builder = StepToTopoDS_Builder()
            builder.SetPrecision(1e-7)
            builder.SetMaxTol(1e-7)
            if isinstance(solid, StepShape_BrepWithVoids):
                builder.Init(theBRepWithVoids=solid, theTP=process)
            else:
                builder.Init(solid, process)
            if not builder.IsDone():
                raise AssertionError(f'cannot translate {path}')
            shapes.append(builder.Value())
        return shapes
    finally:
        for name, value in previous.items():
            Interface_Static.SetIVal_s(name, value)


class LooplessFaceImports(unittest.TestCase):
    def test_complete_sphere_and_torus(self):
        poses = [
            {'rotation': ['1', '0', '0', '0'], 'translation': ['0', '0', '0']},
            {'rotation': ['1', '2', '3', '4'], 'translation': ['7/3', '-11/7', '13/5']},
        ]
        fixtures = [
            ('sphere', {'radius': '5'}, (1, 3, 2, 2), 4 * math.pi * 5**3 / 3),
            ('torus', {'major': '3', 'minor': '1'}, (1, 2, 1, 0), 2 * math.pi**2 * 3),
        ]
        with tempfile.TemporaryDirectory(prefix='catallax-loopless-step-') as temporary:
            directory = Path(temporary)
            for name, dimensions, counts, volume in fixtures:
                for index, pose in enumerate(poses):
                    with self.subTest(solid=name, pose=index, ocp=OCP.__version__):
                        recipe = {'format': 'catallax-recipe/3', 'name': name,
                                  'steps': [{'id': 'body', 'op': name, **dimensions, 'pose': pose}],
                                  'result': 'body'}
                        source = directory / f'{name}-{index}.json'
                        source.write_text(json.dumps(recipe))
                        out = directory / f'{name}-{index}'
                        completed = subprocess.run(
                            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary',
                             str(BINARY), '--out', str(out), '--fresh', '--workers', '1', '--json'],
                            check=True, capture_output=True, text=True)
                        self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                        case = json.loads((out / 'run.json').read_text())['cases'][0]
                        manifest = case['response']['manifest']
                        self.assertTrue(manifest['exact'], manifest)
                        self.assertEqual(manifest['marks'], [])
                        shape = raw_solid(Path(case['path']) / 'answer.step')
                        self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), case['path'])
                        edges = unique(shape, TopAbs_EDGE)
                        actual = (len(unique(shape, TopAbs_FACE)), len(edges),
                                  len(unique(shape, TopAbs_VERTEX)),
                                  sum(BRep_Tool.Degenerated_s(TopoDS.Edge(edge)) for edge in edges))
                        self.assertEqual(actual, counts)
                        self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                        properties = GProp_GProps()
                        BRepGProp.VolumeProperties_s(shape, properties)
                        self.assertAlmostEqual(properties.Mass(), volume, delta=volume * 1e-10)


if __name__ == '__main__':
    unittest.main()
