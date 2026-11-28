"""Check the Bore finish attachment in both standalone P1 STEP saved views.

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
    XCAFDoc_ShapeTool, XCAFDoc_View,
)


def attribute(label, kind):
    attributes = TDF_AttributeIterator(label)
    while attributes.More():
        value = attributes.Value()
        if isinstance(value, kind):
            return value.GetObject()
        attributes.Next()
    return None


def check(path, name, normal):
    reader = STEPCAFControl_Reader()
    reader.SetViewMode(True)
    reader.SetGDTMode(True)
    document = TDocStd_Document(TCollection_ExtendedString("BinXCAF"))
    assert reader.ReadFile(str(path)) == IFSelect_RetDone, path
    assert reader.Transfer(document), path
    tool = XCAFDoc_DocumentTool.ViewTool_s(document.Main())
    views = Sequence_TDF_Label()
    tool.GetViewLabels(views)
    view, = views
    assert attribute(view, XCAFDoc_View).Name().ToCString() == name
    callouts = Sequence_TDF_Label()
    assert tool.GetRefGDTLabel(view, callouts), name
    found = False
    for callout in callouts:
        dimension = attribute(callout, XCAFDoc_Dimension)
        label = dimension.GetPresentationName() if dimension is not None else None
        if label is None or label.ToCString() != "Bore finish":
            continue
        found = True
        assert dimension.HasPlane() and dimension.GetPlane().Direction().Coord() == normal, name
        primary, secondary = Sequence_TDF_Label(), Sequence_TDF_Label()
        XCAFDoc_DimTolTool.GetRefShapeLabel_s(callout, primary, secondary)
        faces = [*primary, *secondary]
        assert len(faces) == 1, f"{name}: expected one Bore finish face, got {len(faces)}"
        face = TopoDS.Face(XCAFDoc_ShapeTool.GetShape_s(faces[0]))
        surface = BRepAdaptor_Surface(face)
        assert surface.GetType() == GeomAbs_Cylinder, name
        bore = surface.Cylinder()
        # P1 has one radius-4 bore at (30, 25) and four radius-2 holes.
        assert (*bore.Location().Coord()[:2], bore.Radius()) == (30, 25, 4), name
    assert found, f"{name}: Bore finish is missing"


if __name__ == "__main__":
    directory = Path(sys.argv[1])
    check(directory / "c4c-bore-view.step", "Bore", (0, -1, 0))
    check(directory / "c4c-top-view.step", "Top", (0, 0, 1))
