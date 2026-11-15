"""N03 raw OCCT check on actual numerical STEP bytes, without healing.

Uses the repository's supported raw StepToTopoDS builder, never TransferRoots,
ShapeFix, sewing or SameParameter. No file is rewritten before validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import OCP
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_ON,
                        TopAbs_OUT, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE)
from test_loopless_faces import unique
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.StepData import StepData_Factors
from OCP.StepShape import StepShape_ManifoldSolidBrep
from OCP.StepToTopoDS import StepToTopoDS_Builder
from OCP.Transfer import Transfer_TransientProcess


def raw_solid(path):
    """Same raw builder as test_loopless_faces, using OCP 7.8's required factors."""
    reader = STEPControl_Reader()
    settings = {'read.stdsameparameter.mode': 0, 'read.maxprecision.mode': 0}
    previous = {key: Interface_Static.IVal_s(key) for key in settings}
    try:
        for key, value in settings.items():
            assert Interface_Static.SetIVal_s(key, value), key
        assert reader.ReadFile(str(path)) == IFSelect_RetDone, 'STEP read failed'
        model = reader.StepModel()
        solids = [model.Value(i) for i in range(1, model.NbEntities() + 1)
                  if isinstance(model.Value(i), StepShape_ManifoldSolidBrep)]
        assert len(solids) == 1, 'expected one explicit STEP solid'
        process = Transfer_TransientProcess()
        process.SetModel(model)
        # These N03 files explicitly use millimetres and radians. No scaling or
        # precision/healing policy is inferred from numerical source tolerance.
        factors = StepData_Factors()
        factors.InitializeFactors(1., 1., 1.)
        builder = StepToTopoDS_Builder()
        builder.SetPrecision(1e-7)
        builder.SetMaxTol(1e-7)
        builder.Init(solids[0], process, factors)
        assert builder.IsDone(), 'raw solid translation failed'
        return builder.Value()
    finally:
        for key, value in previous.items():
            Interface_Static.SetIVal_s(key, value)



def volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()


def topology(shape):
    faces = unique(shape, TopAbs_FACE)
    wires = [len(unique(face, TopAbs_WIRE)) for face in faces]
    vertices = len(unique(shape, TopAbs_VERTEX))
    edges = len(unique(shape, TopAbs_EDGE))
    return {'solids': len(unique(shape, TopAbs_SOLID)), 'faces': len(faces),
            'edges': edges, 'vertices': vertices, 'two_loop_faces': wires.count(2),
            'euler_with_face_holes': vertices - edges + sum(2 - n for n in wires)}


def state(shape, point):
    result = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State()
    return {TopAbs_IN: 'inside', TopAbs_OUT: 'outside', TopAbs_ON: 'boundary'}.get(result, 'unknown')


def check(path):
    before = path.read_bytes()
    text = before.decode('utf-8')
    assert 'catallax numerical approximate' in text
    assert 'complete_body_bound=unestablished' in text
    assert 'UNCERTAINTY_MEASURE_WITH_UNIT' not in text
    assert 'SI_UNIT(.MILLI.,.METRE.)' in text and 'SI_UNIT($,.RADIAN.)' in text
    raw = raw_solid(path)
    assert BRepCheck_Analyzer(raw).IsValid(), 'raw import invalid'
    raw_topology = topology(raw)
    assert raw_topology == {'solids': 1, 'faces': 10, 'edges': 24, 'vertices': 16,
                            'two_loop_faces': 2, 'euler_with_face_holes': 0}, raw_topology
    raw_volume = volume(raw)
    assert abs(raw_volume - 12.) < 1e-9, raw_volume
    probes = {'kept_wall': (.75, 2., .5), 'removed_wall': (.25, 2., .5),
              'hole': (2., 2., .5), 'outside': (-.25, 2., .5)}
    raw_states = {key: state(raw, p) for key, p in probes.items()}
    assert raw_states == {'kept_wall': 'inside', 'removed_wall': 'inside',
                          'hole': 'outside', 'outside': 'outside'}, raw_states
    # Remove x <= .5 through the entire tube height; expected material is 10 mm^3.
    tool = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), 1.5, 6., 3.).Shape()
    operation = BRepAlgoAPI_Cut(raw, tool)
    operation.Build()
    assert operation.IsDone(), 'subsequent cut failed'
    cut = operation.Shape()
    assert not cut.IsNull() and BRepCheck_Analyzer(cut).IsValid(), 'cut invalid'
    cut_topology = topology(cut)
    assert cut_topology == raw_topology, cut_topology
    cut_volume = volume(cut)
    assert abs(cut_volume - 10.) < 1e-9, cut_volume
    cut_states = {key: state(cut, p) for key, p in probes.items()}
    assert cut_states == {'kept_wall': 'inside', 'removed_wall': 'outside',
                          'hole': 'outside', 'outside': 'outside'}, cut_states
    assert state(cut, (.5, 2., .5)) == 'boundary'
    assert path.read_bytes() == before, 'reference check changed input file'
    return {'file': path.name, 'sha256': hashlib.sha256(before).hexdigest(),
            'category': 'numerical approximate', 'healing': False,
            'raw_builder': 'StepToTopoDS_Builder with explicit identity StepData_Factors',
            'same_parameter': False, 'unit_factors': [1., 1., 1.], 'read_maxprecision_mode': 0,
            'builder_precision_mm': '1e-7', 'builder_max_tolerance_mm': '1e-7',
            'raw_import': {'valid': True, 'topology': raw_topology,
                           'volume_mm3': raw_volume, 'material_probes': raw_states},
            'subsequent_cut': {'valid': True, 'topology': cut_topology,
                               'volume_mm3': cut_volume, 'material_probes': cut_states,
                               'new_boundary_probe': 'boundary'},
            'input_bytes_unchanged': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    report = {'format': 'catallax-n03-raw-reference/1', 'ocp_version': OCP.__version__,
              'python_version': platform.python_version(), 'cases': [check(p) for p in args.files]}
    text = json.dumps(report, indent=2) + '\n'
    if args.out:
        args.out.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
