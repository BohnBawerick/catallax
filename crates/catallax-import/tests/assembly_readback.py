#!/usr/bin/env python3
"""Read both files produced by assembly_step.rs using OCCT XDE."""
import math
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Quantity import Quantity_Color, Quantity_TOC_sRGB
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_ExtendedString
from OCP.TDataStd import TDataStd_Name
from OCP.TDF import TDF_AttributeIterator, TDF_Label
from OCP.TDocStd import TDocStd_Document
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCP.TopoDS import TopoDS
from OCP.XCAFDimTolObjects import (
    XCAFDimTolObjects_DimensionType_Size_Diameter,
    XCAFDimTolObjects_GeomToleranceType_Parallelism,
)
from OCP.XCAFDoc import (
    XCAFDoc_ColorSurf,
    XCAFDoc_ColorTool,
    XCAFDoc_Datum,
    XCAFDoc_DimTolTool,
    XCAFDoc_Dimension,
    XCAFDoc_DocumentTool,
    XCAFDoc_GeomTolerance,
    XCAFDoc_ShapeTool,
)
from OCP.gp import gp_Pnt
from OCP.collections import Sequence_TDF_Label


class WrongPmiFace(AssertionError):
    """The annotation resolves to a face other than the specified one."""


def reject_swapped_pmi(path):
    source = path.read_text()
    plate = re.search(r"(#\d+)=ADVANCED_BREP_SHAPE_REPRESENTATION\('Base''s plate',", source)
    assert plate is not None
    uses = re.compile(
        rf"GEOMETRIC_ITEM_SPECIFIC_USAGE\('','',#\d+,{plate[1]},(#\d+)\);"
    )
    first, second = list(uses.finditer(source))
    left, right = first.span(1), second.span(1)
    assert first[1] != second[1]
    swapped = (source[:left[0]] + second[1] + source[left[1]:right[0]]
               + first[1] + source[right[1]:])
    with TemporaryDirectory(dir=path.parent) as directory:
        edited = Path(directory) / path.name
        edited.write_text(swapped)
        try:
            read(edited)
        except WrongPmiFace:
            print(f"{path.name}: detected swapped PMI face references")
            return
    raise AssertionError(f"{path.name}: accepted the plate's swapped PMI face references")


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


def check_attachment(label, expected):
    first, second = Sequence_TDF_Label(), Sequence_TDF_Label()
    assert XCAFDoc_DimTolTool.GetRefShapeLabel_s(label, first, second)
    assert first.Length() == 1 and second.Length() == 0
    face_label = first.Value(1)
    shape = XCAFDoc_ShapeTool.GetShape_s(face_label)
    assert shape.ShapeType() == TopAbs_FACE
    owner = name(face_label.Father())
    kind, centre = expected[owner]
    face = BRepAdaptor_Surface(TopoDS.Face(shape))
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props)
    actual = coordinates(props.CentreOfMass())
    try:
        assert face.GetType() == kind
        check_close(actual, centre)
    except AssertionError as error:
        raise WrongPmiFace(
            f"{owner}: expected {kind} at {centre}, got {face.GetType()} at {actual}"
        ) from error
    return owner, face


def annotation_object(label, attribute_type):
    attributes = TDF_AttributeIterator(label)
    while attributes.More():
        attribute = attributes.Value()
        if isinstance(attribute, attribute_type):
            return attribute.GetObject()
        attributes.Next()
    raise AssertionError(f"missing {attribute_type.__name__} on {name(label)}")


def check_pmi(dim_tol):
    datums = Sequence_TDF_Label()
    dim_tol.GetDatumLabels(datums)
    assert datums.Length() == 2
    datum_labels = {}
    for i in (1, 2):
        label = datums.Value(i)
        owner, _ = check_attachment(label, {
            "Base's plate": (GeomAbs_Plane, (0, 4, 1.5)),
            "Shaft #17": (GeomAbs_Plane, (0, 0, 0)),
        })
        assert annotation_object(label, XCAFDoc_Datum).GetName().ToCString() == "A"
        datum_labels[owner] = label
    assert set(datum_labels) == {"Base's plate", "Shaft #17"}

    dimensions = Sequence_TDF_Label()
    dim_tol.GetDimensionLabels(dimensions)
    assert dimensions.Length() == 1
    label = dimensions.Value(1)
    _, shaft = check_attachment(label, {"Shaft #17": (GeomAbs_Cylinder, (0, 0, 2.5))})
    dimension = annotation_object(label, XCAFDoc_Dimension)
    assert dimension.GetType() == XCAFDimTolObjects_DimensionType_Size_Diameter
    check_close((shaft.Cylinder().Radius(), dimension.GetValue()), (1, 2))

    tolerances = Sequence_TDF_Label()
    dim_tol.GetGeomToleranceLabels(tolerances)
    assert tolerances.Length() == 2
    owners = set()
    for i in (1, 2):
        label = tolerances.Value(i)
        owner, _ = check_attachment(label, {
            "Base's plate": (GeomAbs_Plane, (6, 4, 1.5)),
            "Shaft #17": (GeomAbs_Plane, (0, 0, 5)),
        })
        owners.add(owner)
        tolerance = annotation_object(label, XCAFDoc_GeomTolerance)
        assert tolerance.GetType() == XCAFDimTolObjects_GeomToleranceType_Parallelism
        check_close((tolerance.GetValue(),), (0.01,))
        references = Sequence_TDF_Label()
        assert XCAFDoc_DimTolTool.GetDatumOfTolerLabels_s(label, references)
        assert references.Length() == 1 and references.Value(1).IsEqual(datum_labels[owner])
    assert owners == set(datum_labels)


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
    check_pmi(XCAFDoc_DocumentTool.DimTolTool_s(document.Main()))
    print(f"{path.name}: nested poses, shared shaft, names, colours and five PMI attachments")


if __name__ == "__main__":
    for filename in ("first.step", "roundtrip.step"):
        path = Path(sys.argv[1]) / filename
        read(path)
        reject_swapped_pmi(path)
