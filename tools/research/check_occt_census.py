"""Execute authored solid/surface/edge/failure controls for the OCCT protocol."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=Path(__file__).with_name('occt_census.py'))
    args = parser.parse_args()
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pln, gp_Pnt, gp_Dir
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCP.TopoDS import TopoDS_Shell

    with tempfile.TemporaryDirectory(prefix='catallax-occt-census-controls-') as folder:
        root = Path(folder)
        face = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0., 1., 0., 1.).Face()
        shell = TopoDS_Shell()
        builder = BRep_Builder()
        builder.MakeShell(shell)
        builder.Add(shell, face)
        shapes = [('solid', BRepPrimAPI_MakeBox(1., 2., 3.).Shape()), ('sheet', face),
            ('open-shell', shell), ('edge', BRepBuilderAPI_MakeEdge(gp_Pnt(0., 0., 0.), gp_Pnt(1., 0., 0.)).Shape())]
        for name, shape in shapes:
            writer = STEPControl_Writer()
            assert writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone
            assert writer.Write(str(root / (name + '.step'))) == IFSelect_RetDone
        (root / 'unreadable.step').write_text('unreadable control\n')
        rows = {}
        for name in ['solid', 'sheet', 'open-shell', 'edge', 'unreadable']:
            result = root / (name + '.json')
            with (root / (name + '.stdout')).open('wb') as stdout, (root / (name + '.stderr')).open('wb') as stderr:
                subprocess.run([sys.executable, '-B', str(args.protocol.resolve()), 'probe',
                    str(root / (name + '.step')), str(result)], check=True, timeout=60, stdout=stdout, stderr=stderr)
            rows[name] = json.loads(result.read_text())
        assert rows['solid']['outcome'] == 'healed_valid_solid'
        assert rows['solid']['solid_count'] == rows['solid']['valid_solid_count'] == 1
        # Standard STEP transfer wraps the free face in an open shell.
        assert rows['sheet']['outcome'] == 'open_shell' and rows['sheet']['solid_count'] == 0
        assert rows['edge']['outcome'] == 'opened_other' and rows['edge']['solid_count'] == 0
        assert rows['open-shell']['outcome'] == 'open_shell'
        assert rows['open-shell']['solid_count'] == 0 and rows['open-shell']['open_shell_count'] == 1
        assert rows['unreadable']['outcome'] == 'failed' and rows['unreadable']['opened'] is False
        assert all(row['source_unchanged'] for row in rows.values())
        print(json.dumps({'controls': {name: row['outcome'] for name, row in rows.items()}}), flush=True)


if __name__ == '__main__':
    main()
