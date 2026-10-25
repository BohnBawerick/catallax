#!/usr/bin/env python3
"""Independent OCCT read-back; no ShapeFix, sewing, or shape processing.

Generate fixtures with the ignored Rust writer test documented in README.md.
Run this script with the fixture directory. Missing bindings fail this check.
"""
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import OCP
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.IGESControl import IGESControl_Reader
from OCP.IGESToBRep import (IGESToBRep, IGESToBRep_AlgoContainer,
                           IGESToBRep_BRepEntity, IGESToBRep_ToolContainer)
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.Transfer import Transfer_TransientProcess

CASES = {
    'box': ({'Plane': 6}, 30.0),
    'bored-plate': ({'Plane': 6, 'Cylinder': 1}, 160 - 8 * math.pi),
    'P2-collar': ({'Plane': 2, 'Cylinder': 2, 'Cone': 2}, 3818 * math.pi / 3),
    'cylinder': ({'Plane': 2, 'Cylinder': 1}, 20 * math.pi),
    'cone': ({'Plane': 2, 'Cone': 1}, 26 * math.pi / 3),
    'rounded-box': ({'Plane': 6}, 10 / 231),
    'posed-cylinder': ({'Plane': 2, 'Cylinder': 1}, 20 * math.pi / 63),
}


def unique(shape, kind):
    out = []
    explorer = TopExp_Explorer(shape, kind)
    while explorer.More():
        current = explorer.Current()
        if not any(current.IsSame(previous) for previous in out):
            out.append(current)
        explorer.Next()
    return out


def check(directory, name, expected_kinds, expected_volume):
    path = directory / f'{name}-mm.igs'
    report = json.loads(path.with_suffix('.json').read_text())
    scale = Fraction(report['scale'])
    bound = (abs(scale - 1) * Fraction(report['coordinate_bound'])
             + scale * Fraction(report['trim_deviation']))
    assert bound == Fraction(report['max_deviation'])
    assert bound <= Fraction(1, 10_000_000)
    reader = IGESControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    model = reader.IGESModel()
    assert model.GlobalSection().UnitFlag() == 2
    assert model.GlobalSection().UnitName().ToCString() == 'MM'
    # IGESControl installs a boundary tool that invokes ShapeFix internally.
    # Replace it with the base IGESToBRep tool and bypass the control actor and
    # XSAlgo::ProcessShape entirely. This is stricter than an empty sequence.
    algo = IGESToBRep_AlgoContainer()
    algo.SetToolContainer(IGESToBRep_ToolContainer())
    IGESToBRep.SetAlgoContainer_s(algo)
    assert algo.ToolContainer().IGESBoundary().DynamicType().Name() == 'IGESToBRep_IGESBoundary'
    transfer = IGESToBRep_BRepEntity()
    transfer.SetModel(model)
    transfer.SetModeApprox(False)
    transfer.SetTransferProcess(Transfer_TransientProcess(model.NbEntities()))
    roots = [model.Entity(i) for i in range(1, model.NbEntities() + 1)
             if model.Entity(i).TypeNumber() == 186]
    assert len(roots) == 1
    shape = transfer.TransferManifoldSolid(roots[0])
    assert not shape.IsNull()
    assert len(unique(shape, TopAbs_SOLID)) == 1
    assert BRepCheck_Analyzer(shape).IsValid(), name
    faces = unique(shape, TopAbs_FACE)
    assert len(faces) == report['faces'] == sum(expected_kinds.values())
    kinds = Counter(BRepAdaptor_Surface((getattr(TopoDS, "Face_s", None) or TopoDS.Face)(f)).GetType().name.removeprefix('GeomAbs_')
                    for f in faces)
    assert kinds == expected_kinds, (name, kinds)
    edges = unique(shape, TopAbs_EDGE)
    assert len(edges) == report['edges'] == report['source_edges'] + report['added_edges']
    assert report['vertices'] == report['source_vertices'] + report['added_vertices']
    seams = 0
    for edge in edges:
        uses = 0
        owners = 0
        for face in faces:
            ex = TopExp_Explorer(face, TopAbs_EDGE)
            count = 0
            while ex.More():
                count += int(ex.Current().IsSame(edge))
                ex.Next()
            uses += count
            owners += int(count > 0)
        assert uses == 2, (name, uses)
        seams += int(owners == 1)
    assert seams == report['added_edges']
    # Signed volume checks orientation/units as a diagnostic. It is deliberately
    # NOT compared to the displacement certificate as if that bounded volume.
    properties = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, properties)
    expected = expected_volume * float(scale) ** 3
    assert math.isclose(properties.Mass(), expected, rel_tol=2e-10, abs_tol=1e-10), (name, properties.Mass(), expected)
    return {'case': name, 'valid': True, 'faces': len(faces), 'surface_kinds': dict(kinds),
            'shared_edges': len(edges), 'seams': seams, 'units': 'MM',
            'volume_diagnostic_mm3': properties.Mass(), 'max_deviation_mm': str(bound)}


if __name__ == '__main__':
    directory = Path(sys.argv[1])
    results = [check(directory, name, *expected) for name, expected in CASES.items()]
    print(json.dumps({'binding': 'cadquery-ocp', 'version': OCP.__version__,
                      'transfer': 'base IGESToBRep; no healing', 'cases': results}, sort_keys=True))
