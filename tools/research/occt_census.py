"""Same-source standard OCCT STEP transfer and default ShapeFix measurement."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def probe(source, output):
    import OCP
    from OCP.BRep import BRep_Tool
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.ShapeFix import ShapeFix_Shape
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer

    def unique(shape, kind, avoid=None):
        result = []
        explorer = TopExp_Explorer(shape, kind) if avoid is None else TopExp_Explorer(shape, kind, avoid)
        while explorer.More():
            item = explorer.Current()
            if not any(item.IsSame(previous) for previous in result):
                result.append(item)
            explorer.Next()
        return result

    row = {'ocp_version': OCP.__version__, 'source_sha256': sha(source),
           'opened': False, 'outcome': 'failed', 'stage': 'read'}
    started = time.monotonic()
    write(output, row)
    try:
        reader = STEPControl_Reader()
        row['reader_integer_defaults'] = {name: Interface_Static.IVal_s(name) for name in
            ['read.stdsameparameter.mode', 'read.maxprecision.mode', 'read.precision.mode', 'read.step.nonmanifold']}
        row['reader_real_defaults'] = {name: Interface_Static.RVal_s(name) for name in
            ['read.precision.val', 'read.maxprecision.val']}
        row['read_status'] = str(reader.ReadFile(str(source)))
        if row['read_status'] != str(IFSelect_RetDone):
            return
        row['stage'] = 'transfer'
        row['roots_available'] = reader.NbRootsForTransfer()
        write(output, row)
        row['roots_transferred'] = reader.TransferRoots()
        row['all_roots_transferred'] = row['roots_transferred'] == row['roots_available']
        row['result_shapes'] = reader.NbShapes()
        shape = reader.OneShape()
        row['opened'] = not shape.IsNull()
        if not row['opened']:
            return
        row['stage'] = 'check_transferred'
        write(output, row)
        row['transferred_whole_shape_valid'] = BRepCheck_Analyzer(shape).IsValid()
        row['stage'] = 'shape_fix'
        fixer = ShapeFix_Shape(shape)
        row['shape_fix_defaults'] = {'precision': fixer.Precision(),
            'min_tolerance': fixer.MinTolerance(), 'max_tolerance': fixer.MaxTolerance()}
        write(output, row)
        row['shape_fix_reported_change'] = fixer.Perform()
        healed = fixer.Shape()
        row['stage'] = 'check_healed'
        write(output, row)
        if healed.IsNull():
            row['failure'] = 'healing_returned_null'
            return
        row['healed_whole_shape_valid'] = BRepCheck_Analyzer(healed).IsValid()
        solids = unique(healed, TopAbs_SOLID)
        row['solid_count'] = len(solids)
        row['valid_solid_count'] = sum(BRepCheck_Analyzer(solid).IsValid() for solid in solids)
        shells = unique(healed, TopAbs_SHELL)
        row['shell_count'] = len(shells)
        row['open_shell_count'] = sum(not BRep_Tool.IsClosed_s(shell) for shell in shells)
        row['free_shell_count'] = len(unique(healed, TopAbs_SHELL, TopAbs_SOLID))
        for name, kind, owner in [('face', TopAbs_FACE, TopAbs_SHELL),
                ('wire', TopAbs_WIRE, TopAbs_FACE), ('edge', TopAbs_EDGE, TopAbs_WIRE),
                ('vertex', TopAbs_VERTEX, TopAbs_EDGE)]:
            row['free_' + name + '_count'] = len(unique(healed, kind, owner))
        # Occurrences avoid a quadratic face-deduplication pass on large assemblies.
        explorer = TopExp_Explorer(healed, TopAbs_FACE)
        row['face_occurrences'] = 0
        while explorer.More():
            row['face_occurrences'] += 1
            explorer.Next()
        row['outcome'] = ('healed_valid_solid' if row['healed_whole_shape_valid'] and solids
            and row['valid_solid_count'] == len(solids)
            and all(row['free_' + name + '_count'] == 0 for name in ['shell', 'face', 'wire', 'edge', 'vertex'])
            else 'open_shell' if row['open_shell_count'] and not solids else 'opened_other')
        row['stage'] = 'complete'
    except Exception as error:
        row['failure'] = type(error).__name__ + ': ' + str(error)
    finally:
        row['elapsed_seconds'] = time.monotonic() - started
        row['source_unchanged'] = sha(source) == row['source_sha256']
        write(output, row)


def collect(args):
    selection = json.loads(args.selection.read_text())
    expected = json.loads(args.hashes.read_text())
    assert len(selection) == len(expected) == 3000
    assert set(key for key, _ in selection) == set(expected)
    assert sha(args.selection) == '1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5'
    assert sha(args.hashes) == '3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7'
    args.out.mkdir()
    raw = args.out / 'raw'
    raw.mkdir()
    identity = {'format': 'catallax-occt-census/1', 'complete': False, 'denominator': 3000,
        'seed': 20261111, 'selection_sha256': sha(args.selection),
        'source_root': str(args.sources.resolve()), 'selection_path': str(args.selection.resolve()),
        'source_hash_manifest_sha256': sha(args.hashes), 'protocol_sha256': sha(Path(__file__)),
        'workers': args.workers, 'external_timeout_seconds': args.timeout,
        'reader_policy': 'fresh-process STEPControl_Reader default transfer; no setting overrides',
        'healing_policy': 'ShapeFix_Shape default Perform; no tolerance overrides',
        'catallax_work_limit_comparable': False}
    write(args.out / 'identity.json', identity)

    def one(item):
        key, member = item
        source = args.sources / member
        assert sha(source) == expected[key], (key, 'source hash before')
        output = raw / (key + '.json')
        started = time.monotonic()
        with (raw / (key + '.stdout')).open('wb') as stdout, (raw / (key + '.stderr')).open('wb') as stderr:
            try:
                process = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()),
                    'probe', str(source), str(output)], stdout=stdout, stderr=stderr, timeout=args.timeout)
                result = json.loads(output.read_text()) if output.exists() else {}
                if process.returncode:
                    result.update(outcome='failed', failure='process_exit_' + str(process.returncode))
            except subprocess.TimeoutExpired:
                result = json.loads(output.read_text()) if output.exists() else {}
                result.update(outcome='timeout', failure='external_deadline')
        result.update(model_id=key, source_sha256=expected[key], wall_seconds=time.monotonic()-started)
        result['source_unchanged'] = sha(source) == expected[key]
        assert result['source_unchanged'], (key, 'source hash after')
        write(output, result)
        return result

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(one, item) for item in selection]):
            rows.append(future.result())
            if len(rows) % 100 == 0:
                print(json.dumps({'completed': len(rows), 'outcomes': dict(Counter(r['outcome'] for r in rows))}), flush=True)
    rows.sort(key=lambda row: row['model_id'])
    assert len(rows) == 3000 and len({r['model_id'] for r in rows}) == 3000
    write(args.out / 'rows.json', rows)
    identity.update(complete=True, rows_sha256=sha(args.out / 'rows.json'),
        outcomes=dict(Counter(r['outcome'] for r in rows)),
        opened=sum(bool(r.get('opened')) for r in rows),
        transferred_whole_shape_valid=sum(bool(r.get('transferred_whole_shape_valid')) for r in rows),
        healed_whole_shape_valid=sum(bool(r.get('healed_whole_shape_valid')) for r in rows),
        ocp_versions=sorted({r['ocp_version'] for r in rows if 'ocp_version' in r}))
    write(args.out / 'identity.json', identity)
    print(json.dumps(identity, sort_keys=True), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    child = modes.add_parser('probe')
    child.add_argument('source', type=Path)
    child.add_argument('output', type=Path)
    parent = modes.add_parser('collect')
    for name in ['sources', 'selection', 'hashes', 'out']:
        parent.add_argument('--' + name, type=Path, required=True)
    parent.add_argument('--workers', type=int, default=4)
    parent.add_argument('--timeout', type=int, default=60)
    options = parser.parse_args()
    if options.mode == 'probe':
        probe(options.source, options.output)
    else:
        collect(options)
