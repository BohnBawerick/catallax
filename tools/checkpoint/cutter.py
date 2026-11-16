"""Pin one drilled-hole edit for a census part from its source bounding box.

Run with the checker's OCP Python. The hole runs along the axis of the smallest
bounding-box extent, through the centre of the other two, 1 mm past each side.
Its radius is a tenth of the smaller remaining extent; the edit makes it 1.5 times
larger. All values are rounded to 0.001 mm and written as exact decimal strings,
so every kernel receives the same numbers.
"""
from decimal import ROUND_HALF_EVEN, Decimal
import json
import sys

from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader

MM = Decimal('0.001')


def mm(value):
    return Decimal(repr(value)).quantize(MM, rounding=ROUND_HALF_EVEN)


def cutter(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone or reader.TransferRoots() == 0:
        raise ValueError('source did not transfer')
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(reader.OneShape(), box, False, False)
    lo, hi = box.CornerMin(), box.CornerMax()
    low = [mm(v) for v in (lo.X(), lo.Y(), lo.Z())]
    high = [mm(v) for v in (hi.X(), hi.Y(), hi.Z())]
    extent = [h - l for l, h in zip(low, high)]
    axis = min(range(3), key=lambda i: (extent[i], i))
    others = [i for i in range(3) if i != axis]
    radius = (min(extent[i] for i in others) / 10).quantize(MM, rounding=ROUND_HALF_EVEN)
    edited = (radius * Decimal('1.5')).quantize(MM, rounding=ROUND_HALF_EVEN)
    base = [None] * 3
    for i in others:
        base[i] = ((low[i] + high[i]) / 2).quantize(MM, rounding=ROUND_HALF_EVEN)
    base[axis] = low[axis] - 1
    return {'axis': [int(i == axis) for i in range(3)], 'base': [str(v) for v in base],
            'height': str(extent[axis] + 2), 'radius': str(radius), 'edited_radius': str(edited),
            'rule': 'hole along the smallest bounding-box extent, through the centre of the others; '
                    'radius one tenth of the smaller remaining extent; edit 1.5 times; 0.001 mm rounding',
            'bounding_box': 'OCP BRepBndLib::AddOptimal without triangulation, source as transferred'}


if __name__ == '__main__':
    print(json.dumps(cutter(sys.argv[1]), sort_keys=True))
