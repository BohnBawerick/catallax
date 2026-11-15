"""Qualify complete STEP replacement on the ten fixed P06 capability sources.

Source and changed STEP bytes stay in the requested private output directory.
Only entity addresses are shifted; quoted text and comments are preserved.
Both source and replacement hashes are pinned by the existing correspondence
manifest before any replacement request runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def renumber(source):
    out = bytearray()
    i = 0
    quoted = comment = False
    while i < len(source):
        if comment:
            if source[i:i + 2] == b'*/':
                out.extend(b'*/')
                i += 2
                comment = False
            else:
                out.append(source[i])
                i += 1
        elif quoted:
            out.append(source[i])
            if source[i] == 39:
                if source[i:i + 2] == b"''":
                    out.append(39)
                    i += 1
                else:
                    quoted = False
            i += 1
        elif source[i:i + 2] == b'/*':
            comment = True
            out.extend(b'/*')
            i += 2
        elif source[i] == 39:
            quoted = True
            out.append(39)
            i += 1
        elif source[i] == 35:
            end = i + 1
            while end < len(source) and 48 <= source[end] <= 57:
                end += 1
            assert end > i + 1
            out.extend(b'#' + str(int(source[i + 1:end]) + 1_000_000_000).encode())
            i = end
        else:
            out.append(source[i])
            i += 1
    assert not quoted and not comment
    return bytes(out)


def digest(value):
    return hashlib.sha256(value).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--sources', required=True, type=Path,
                        help='private MODEL_ID.step sources')
    parser.add_argument('--out', required=True, type=Path,
                        help='new private output directory')
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[2]
                        / 'docs/phase1/import-correspondence/real-controls.json')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    assert len(manifest['rows']) == 10
    args.out.mkdir()
    rows = []
    for expected in manifest['rows']:
        key = expected['model_id']
        old = args.sources / (key + '.step')
        source = old.read_bytes()
        assert digest(source) == expected['source_sha256'], (key, 'source hash')
        changed = renumber(source)
        assert digest(changed) == expected['replacement_sha256'], (key, 'replacement hash')
        new = args.out / (key + '-renumbered.step')
        new.write_bytes(changed)
        with (args.out / (key + '.json')).open('wb') as out, \
                (args.out / (key + '.stderr')).open('wb') as err:
            result = subprocess.run([str(args.binary), str(old), str(new)],
                                    stdout=out, stderr=err, timeout=60, check=False)
        assert result.returncode == 0, (key, 'process failed')
        observed = json.loads((args.out / (key + '.json')).read_text())
        assert observed['status'] == 'replaced', (key, observed.get('code'))
        assert observed['source_sha256'] == expected['source_sha256']
        assert observed['replacement_sha256'] == expected['replacement_sha256']
        assert observed['unresolved_references'] == 0
        assert observed['old_revision_preserved']
        assert observed['save_reopen_checked'] and observed['old_revision_replay_checked']
        assert observed['work_used'] <= observed['work_limit'] == 20_000_000
        assert old.read_bytes() == source and new.read_bytes() == changed
        row = {'model_id': key, **observed}
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    report = {'schema': 'catallax-file-replacement-controls/1',
              'complete': True, 'count': 10,
              'binary_sha256': digest(args.binary.read_bytes()),
              'entity_id_offset': 1_000_000_000, 'geometry_changed': False,
              'source_hash_manifest_sha256': manifest['source_hash_manifest_sha256'],
              'rows': rows}
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
