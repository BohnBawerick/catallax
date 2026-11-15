"""Emit mixed-fillet mass references with OCCT 7.9.3 and both moment flags enabled."""

import json
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Trsf


def main():
    box = BRepPrimAPI_MakeBox(20.0, 20.0, 5.0).Shape()
    sanity = GProp_GProps()
    BRepGProp.VolumeProperties_s(box, sanity)
    assert abs(sanity.MatrixOfInertia().Value(1, 1) - 2000 * (400 + 25) / 12) < 1e-7
    fuse = BRepAlgoAPI_Fuse(box, BRepPrimAPI_MakeBox(20.0, 5.0, 10.0).Shape())
    fuse.SimplifyResult()
    stock = fuse.Shape()
    edges = []
    explorer = TopExp_Explorer(stock, TopAbs_EDGE)
    while explorer.More():
        edge = TopoDS.Edge_s(explorer.Current())
        if not any(edge.IsSame(other) for other in edges):
            edges.append(edge)
        explorer.Next()
    assert len(edges) == 18
    rows = []
    for radius, posed in [(1.0, False), (0.75, False), (0.75, True)]:
        fillet = BRepFilletAPI_MakeFillet(stock)
        for edge in edges:
            fillet.Add(radius, edge)
        fillet.Build()
        assert fillet.IsDone()
        shape = fillet.Shape()
        if posed:
            transform = gp_Trsf()
            transform.SetValues(-2/3, 2/15, 11/15, 7/3,
                                2/3, -1/3, 2/3, -11/7,
                                1/3, 14/15, 2/15, 13/5)
            shape = BRepBuilderAPI_Transform(shape, transform, True).Shape()
        assert BRepCheck_Analyzer(shape).IsValid()
        mass = GProp_GProps()
        BRepGProp.VolumePropertiesGK_s(
            shape, mass, Eps=1e-12, OnlyClosed=True, CGFlag=True, IFlag=True
        )
        area = GProp_GProps()
        BRepGProp.SurfaceProperties_s(shape, area)
        tensor = mass.MatrixOfInertia()
        assert min(tensor.Value(i, i) for i in range(1, 4)) > 1000
        if not posed:
            assert abs(mass.CentreOfMass().X() - 10) < 1e-9
        rows.append({'radius': radius, 'posed': posed, 'volume': mass.Mass(),
                     'area': area.Mass(), 'centroid': list(mass.CentreOfMass().Coord()),
                     'inertia': [[tensor.Value(i, j) for j in range(1, 4)] for i in range(1, 4)]})
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
