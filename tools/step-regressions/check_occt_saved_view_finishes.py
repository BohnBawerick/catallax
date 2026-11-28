"""Check P1 bore callouts in Open CASCADE's transferred saved views.

Generate the input files with CATALLAX_PMI_OUT set when running the
saved_views_write_their_callouts_and_read_back Rust test. Pass that directory
to this check with cadquery-ocp available.
"""

from pathlib import Path
import sys

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.collections import Sequence_TDF_Label
from OCP.GeomAbs import GeomAbs_Cylinder
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_ExtendedString
from OCP.TDF import TDF_AttributeIterator
from OCP.TDocStd import TDocStd_Document
from OCP.TopoDS import TopoDS
from OCP.XCAFDoc import (
    XCAFDoc_Dimension, XCAFDoc_DimTolTool, XCAFDoc_DocumentTool,
    XCAFDoc_GeomTolerance, XCAFDoc_ShapeTool, XCAFDoc_View,
)


def attribute(label, kind):
    attributes = TDF_AttributeIterator(label)
    while attributes.More():
        value = attributes.Value()
        if isinstance(value, kind):
            return value.GetObject()
        attributes.Next()
    return None


def check_callout(tool, view, name, normal, annotation):
    callouts = Sequence_TDF_Label()
    assert tool.GetRefGDTLabel(view, callouts), name
    found = False
    for callout in callouts:
        item = attribute(callout, XCAFDoc_Dimension) or attribute(callout, XCAFDoc_GeomTolerance)
        label = item.GetPresentationName() if item is not None else None
        if label is None or label.ToCString() != annotation:
            continue
        found = True
        assert item.HasPlane() and item.GetPlane().Direction().Coord() == normal, (name, annotation)
        assert not item.GetPresentation().IsNull(), (name, annotation)
        primary, secondary = Sequence_TDF_Label(), Sequence_TDF_Label()
        XCAFDoc_DimTolTool.GetRefShapeLabel_s(callout, primary, secondary)
        faces = [*primary, *secondary]
        assert len(faces) == 1, f"{name}: expected one {annotation} face, got {len(faces)}"
        face = TopoDS.Face(XCAFDoc_ShapeTool.GetShape_s(faces[0]))
        surface = BRepAdaptor_Surface(face)
        assert surface.GetType() == GeomAbs_Cylinder, (name, annotation)
        bore = surface.Cylinder()
        # P1 has one radius-4 bore at (30, 25) and four radius-2 holes.
        assert (*bore.Location().Coord()[:2], bore.Radius()) == (30, 25, 4), (name, annotation)
    assert found, f"{name}: {annotation} is missing"


def check(path, normals, annotation):
    reader = STEPCAFControl_Reader()
    reader.SetViewMode(True)
    reader.SetGDTMode(True)
    document = TDocStd_Document(TCollection_ExtendedString("BinXCAF"))
    assert reader.ReadFile(str(path)) == IFSelect_RetDone, path
    assert reader.Transfer(document), path
    tool = XCAFDoc_DocumentTool.ViewTool_s(document.Main())
    views = Sequence_TDF_Label()
    tool.GetViewLabels(views)
    remaining = dict(normals)
    for view in views:
        name = attribute(view, XCAFDoc_View).Name().ToCString()
        check_callout(tool, view, name, remaining.pop(name), annotation)
    assert not remaining, f"missing saved views: {remaining}"


if __name__ == "__main__":
    directory = Path(sys.argv[1])
    bore = {"Bore": (0, -1, 0)}
    top = {"Top": (0, 0, 1)}
    check(directory / "c4c-bore-view.step", bore, "Bore finish")
    check(directory / "c4c-top-view.step", top, "Bore finish")
    check(directory / "c4c-saved-views.step", bore | top, "Bore position")
