"""Pin one drilled-hole edit for a census part from its source bounding box.

Run with the checker's OCP Python. The hole runs along the axis of the smallest
bounding-box extent, 1 mm past each side. Its radius is a tenth of the smaller
remaining extent; the edit makes it 1.5 times larger. The hole goes through the
centre of the other two extents when that removes material. Otherwise it moves
to the first of a fixed list of offsets, in eighths of those extents, where it
does: OCCT's common of the source and the cutter must exceed a thousandth of the
source volume at the drilled radius, and grow by as much again at the edited
radius. That is a thousand times the scoreboard's removal threshold, so the
choice does not hinge on how a kernel measures. All values are rounded to
0.001 mm and written as exact decimal strings, so every kernel receives the same
numbers. The source volume, also from OCCT, is the reference for judging whether
the drill removed material.
"""
from decimal import ROUND_HALF_EVEN, Decimal
from fractions import Fraction
import json
import sys

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader

MM = Decimal('0.001')
# A margin of 1000 times the scoreboard's REMOVAL fraction.
MARGIN = 1e-3
# Offsets from the centre in eighths of each remaining extent, nearest first. The
# largest, 3/8, keeps a radius of 1/10 of the smaller extent inside the box.
EIGHTHS = (0, 1, -1, 2, -2, 3, -3)
OFFSETS = sorted(((a, b) for a in EIGHTHS for b in EIGHTHS),
                 key=lambda ab: (max(map(abs, ab)), abs(ab[0]) + abs(ab[1]), -ab[0], -ab[1]))


def mm(value):
    return Decimal(repr(value)).quantize(MM, rounding=ROUND_HALF_EVEN)


def volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return abs(props.Mass())


def removed(source, base, axis, radius, height):
    frame = gp_Ax2(gp_Pnt(*map(float, base)), gp_Dir(*map(float, axis)))
    common = BRepAlgoAPI_Common(source, BRepPrimAPI_MakeCylinder(frame, float(radius), float(height)).Shape())
    if not common.IsDone():
        return None
    return volume(common.Shape())


def cutter(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone or reader.TransferRoots() == 0:
        raise ValueError('source did not transfer')
    source = reader.OneShape()
    source_volume = volume(source)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(source, box, False, False)
    lo, hi = box.CornerMin(), box.CornerMax()
    low = [mm(v) for v in (lo.X(), lo.Y(), lo.Z())]
    high = [mm(v) for v in (hi.X(), hi.Y(), hi.Z())]
    extent = [h - l for l, h in zip(low, high)]
    axis = min(range(3), key=lambda i: (extent[i], i))
    others = [i for i in range(3) if i != axis]
    radius = (min(extent[i] for i in others) / 10).quantize(MM, rounding=ROUND_HALF_EVEN)
    edited = (radius * Decimal('1.5')).quantize(MM, rounding=ROUND_HALF_EVEN)
    height = extent[axis] + 2
    direction = [int(i == axis) for i in range(3)]
    tried = []
    for offset in OFFSETS:
        base = [None] * 3
        for i, eighths in zip(others, offset):
            base[i] = ((low[i] + high[i]) / 2 + extent[i] * eighths / 8).quantize(MM, rounding=ROUND_HALF_EVEN)
        base[axis] = low[axis] - 1
        drilled = removed(source, base, direction, radius, height)
        widened = removed(source, base, direction, edited, height)
        tried.append({'offset_eighths': list(offset), 'removed_mm3': drilled, 'edited_removed_mm3': widened})
        if (drilled is not None and widened is not None and drilled > MARGIN * source_volume
                and widened - drilled > MARGIN * source_volume):
            break
    else:
        raise ValueError('no offset removes material at both radii')
    return {'axis': direction, 'base': [str(v) for v in base],
            'height': str(height), 'radius': str(radius), 'edited_radius': str(edited),
            'offset_eighths': list(offset), 'tried': tried,
            'rule': 'hole along the smallest bounding-box extent, through the centre of the others or the first '
                    f'listed offset in eighths of them where OCCT removes more than {Fraction(MARGIN)} of the source '
                    'volume at both radii; radius one tenth of the smaller remaining extent; edit 1.5 times; '
                    '0.001 mm rounding',
            'bounding_box': 'OCP BRepBndLib::AddOptimal without triangulation, source as transferred',
            'source_volume_mm3': source_volume}


if __name__ == '__main__':
    print(json.dumps(cutter(sys.argv[1]), sort_keys=True))
