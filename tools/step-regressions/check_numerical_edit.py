"""N04 raw, healing-off readback of the resized and cut N02 tube."""
import argparse
import hashlib
import json
from pathlib import Path

import OCP
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from check_numerical_step import raw_solid, state, topology, volume


def check(path):
    before = path.read_bytes()
    assert b'catallax numerical approximate' in before
    raw = raw_solid(path)
    assert BRepCheck_Analyzer(raw).IsValid(), 'raw import invalid'
    counts = topology(raw)
    assert counts == {'solids': 1, 'faces': 10, 'edges': 24, 'vertices': 16,
                      'two_loop_faces': 2, 'euler_with_face_holes': 0}, counts
    assert abs(volume(raw) - 18.) < 1e-9, volume(raw)
    probes = {'kept_wall': (.75, 2., 1.25), 'later_removed_wall': (.25, 2., .75),
              'cut_away': (.75, 2., 1.75), 'hole': (2., 2., .75),
              'outside': (-.25, 2., .75), 'new_cap': (.75, 2., 1.5)}
    expected = {'kept_wall': 'inside', 'later_removed_wall': 'inside',
                'cut_away': 'outside', 'hole': 'outside', 'outside': 'outside',
                'new_cap': 'boundary'}
    observed = {key: state(raw, point) for key, point in probes.items()}
    assert observed == expected, observed
    tool = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), 1.5, 6., 4.).Shape()
    operation = BRepAlgoAPI_Cut(raw, tool)
    operation.Build()
    assert operation.IsDone()
    cut = operation.Shape()
    assert BRepCheck_Analyzer(cut).IsValid()
    assert abs(volume(cut) - 15.) < 1e-9, volume(cut)
    assert state(cut, probes['later_removed_wall']) == 'outside'
    assert state(cut, probes['kept_wall']) == 'inside'
    assert state(cut, probes['hole']) == 'outside'
    assert path.read_bytes() == before
    return {'format': 'catallax-n04-raw-reference/1', 'ocp_version': OCP.__version__,
            'sha256': hashlib.sha256(before).hexdigest(), 'healing': False,
            'same_parameter': False, 'raw_valid': True, 'topology': counts,
            'volume_mm3': volume(raw), 'material_probes': observed,
            'subsequent_cut': {'valid': True, 'volume_mm3': volume(cut)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = json.dumps(check(args.file), indent=2) + '\n'
    if args.out:
        args.out.write_text(result)
    print(result, end='')


if __name__ == '__main__':
    main()
