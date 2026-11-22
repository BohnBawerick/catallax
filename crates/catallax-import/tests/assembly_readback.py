#!/usr/bin/env python3
"""Read both files produced by assembly_step.rs using OCCT XDE."""
import math
from pathlib import Path
import sys

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone
from OCP.Quantity import Quantity_Color, Quantity_TOC_sRGB
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_ExtendedString
from OCP.TDataStd import TDataStd_Name
from OCP.TDF import TDF_Label
from OCP.TDocStd import TDocStd_Document
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCP.XCAFDoc import (
    XCAFDoc_ColorSurf,
    XCAFDoc_ColorTool,
    XCAFDoc_DimTolTool,
    XCAFDoc_DocumentTool,
    XCAFDoc_ShapeTool,
)
from OCP.gp import gp_Pnt
from OCP.collections import Sequence_TDF_Label


def name(label):
    attribute = TDataStd_Name()
    assert label.FindAttribute(TDataStd_Name.GetID_s(), attribute)
    return attribute.Get().ToExtString()


def components(label):
    seq = Sequence_TDF_Label()
    assert XCAFDoc_ShapeTool.GetComponents_s(label, seq)
    return {name(seq.Value(i)): seq.Value(i) for i in range(1, seq.Length() + 1)}


def definition(occurrence):
    label = TDF_Label()
    assert XCAFDoc_ShapeTool.GetReferredShape_s(occurrence, label)
    return label


def check_close(actual, expected):
    assert all(math.isclose(a, b, rel_tol=0.0, abs_tol=1e-7)
               for a, b in zip(actual, expected, strict=True)), (actual, expected)


def coordinates(point):
    return (point.X(), point.Y(), point.Z())


def read(path):
    reader = STEPCAFControl_Reader()
    reader.SetNameMode(True)
    reader.SetColorMode(True)
    reader.SetGDTMode(True)
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    document = TDocStd_Document(TCollection_ExtendedString("BinXCAF"))
    assert reader.Transfer(document)
    shapes = XCAFDoc_DocumentTool.ShapeTool_s(document.Main())
    roots = Sequence_TDF_Label()
    shapes.GetFreeShapes(roots)
    assert roots.Length() == 1
    root = roots.Value(1)
    assert name(root) == "Drive assembly"
    assert BRepCheck_Analyzer(XCAFDoc_ShapeTool.GetShape_s(root)).IsValid()
    instances = components(root)
    assert set(instances) == {"Fixed base", "Rotated carrier"}
    base_instance = instances["Fixed base"]
    carrier_instance = instances["Rotated carrier"]
    base = definition(base_instance)
    carrier = definition(carrier_instance)
    assert name(base) == "Base's plate"
    assert name(carrier) == "Carrier"
    nested = components(carrier)
    assert set(nested) == {"Left shaft", "Right shaft"}
    shaft = definition(nested["Left shaft"])
    assert name(shaft) == "Shaft #17"
    assert shaft.IsEqual(definition(nested["Right shaft"]))
    assert XCAFDoc_ShapeTool.GetShape_s(shaft).ShapeType() == TopAbs_SOLID
    base_motion = XCAFDoc_ShapeTool.GetLocation_s(base_instance).Transformation()
    check_close(coordinates(gp_Pnt().Transformed(base_motion)), (4, 6, 8))
    assembly_motion = XCAFDoc_ShapeTool.GetLocation_s(carrier_instance).Transformation()
    for label, expected in [("Left shaft", (10, 22, 30)), ("Right shaft", (7, 20, 30))]:
        local = XCAFDoc_ShapeTool.GetLocation_s(nested[label]).Transformation()
        point = gp_Pnt().Transformed(local).Transformed(assembly_motion)
        check_close(coordinates(point), expected)
        end = gp_Pnt(1, 0, 0).Transformed(local).Transformed(assembly_motion)
        check_close(coordinates(end), (expected[0], expected[1] + 1, expected[2]))
    for label, channels in [(base, (0.125, 0.25, 0.5)), (shaft, (0.75, 0.25, 0))]:
        colour = Quantity_Color()
        assert XCAFDoc_ColorTool.GetColor_s(label, XCAFDoc_ColorSurf, colour)
        check_close(colour.Values(Quantity_TOC_sRGB), channels)
    dim_tol = XCAFDoc_DocumentTool.DimTolTool_s(document.Main())
    datums = Sequence_TDF_Label()
    dim_tol.GetDatumLabels(datums)
    assert datums.Length() == 2
    owners = set()
    for i in (1, 2):
        first, second = Sequence_TDF_Label(), Sequence_TDF_Label()
        assert XCAFDoc_DimTolTool.GetRefShapeLabel_s(datums.Value(i), first, second)
        assert first.Length() == 1 and second.Length() == 0
        face = first.Value(1)
        assert XCAFDoc_ShapeTool.GetShape_s(face).ShapeType() == TopAbs_FACE
        owners.add(name(face.Father()))
    assert owners == {"Base's plate", "Shaft #17"}, owners
    for get_labels, expected_owners in [
        (dim_tol.GetDimensionLabels, {"Shaft #17"}),
        (dim_tol.GetGeomToleranceLabels, {"Base's plate", "Shaft #17"}),
    ]:
        labels = Sequence_TDF_Label()
        get_labels(labels)
        assert labels.Length() == len(expected_owners)
        owners = set()
        for i in range(1, labels.Length() + 1):
            first, second = Sequence_TDF_Label(), Sequence_TDF_Label()
            assert XCAFDoc_DimTolTool.GetRefShapeLabel_s(labels.Value(i), first, second)
            assert first.Length() == 1 and second.Length() == 0
            face = first.Value(1)
            assert XCAFDoc_ShapeTool.GetShape_s(face).ShapeType() == TopAbs_FACE
            owners.add(name(face.Father()))
        assert owners == expected_owners, owners
    print(f"{path.name}: nested poses, shared shaft, names, colours and five PMI attachments")


if __name__ == "__main__":
    for filename in ("first.step", "roundtrip.step"):
        read(Path(sys.argv[1]) / filename)
