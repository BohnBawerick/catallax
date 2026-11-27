"""Project certified face identities for answer records, not for wire manifests."""
from copy import deepcopy
import hashlib
import re


def face_key(revision, step, handle):
    """The legacy `catallax_features::certified::key` for label `face`."""
    digest = hashlib.sha256(b'catallax.elliptic-prism-certificate/1\0')
    digest.update(bytes.fromhex(revision))
    for text in (step.encode(), b'face'):
        digest.update(len(text).to_bytes(8, 'little'))
        digest.update(text)
    digest.update(handle.to_bytes(4, 'little'))
    return digest.hexdigest()


def _project_certificate(certificate, step):
    regions = certificate['regions']
    # Junction regions follow faces_of_body; prism regions are lower, upper, side.
    ranks = {row['face']: rank for rank, row in enumerate(regions)}
    if len(ranks) != len(regions):
        raise ValueError('certificate has duplicate result faces')
    keys = {}
    revision = certificate['identity']['revision']
    for row in regions:
        handle = row['face']
        if row['id'] != face_key(revision, step, handle):
            raise ValueError('certificate region key does not bind its face')
        keys[row['id']] = {'face_rank': ranks[handle]}
        row['face'] = ranks[handle]

    def key(value):
        return keys.get(value, value)

    def identity(row):
        row['representation'] = key(row['representation'])
        scope = row['scope']
        if 'entities' in scope:
            scope['entities'] = [key(entity) for entity in scope['entities']]

    for row in regions:
        row['id'] = key(row['id'])
        identity(row['identity'])
    for row in certificate['evidence_records']:
        identity(row['identity'])
    spline = certificate['spline']
    if isinstance(spline, list):
        for row in spline:
            row['face'] = ranks[row['face']]
    target = certificate['retained_target']
    for field in ('envelope_faces', 'saddle_faces', 'sphere_band_faces'):
        if field in target:
            target[field] = [ranks[face] for face in target[field]]
    return ranks, keys


def certificate_by_face_rank(certificate, step):
    result = deepcopy(certificate)
    _project_certificate(result, step)
    return result


def _generated_name(name, handle, rank):
    match = re.fullmatch(rf'face_{handle}(?:_patch_([0-9]+))?', name)
    if match:
        # A literal name, even one resembling a rank label, is still distinct.
        return {'face_rank': rank, 'patch': match[1]}
    return name


def _project_named_result(result, step):
    ranks, keys = _project_certificate(result['certificate'], step)
    for row in result['names']['faces']:
        row['face'] = ranks[row['face']]
    result['names']['faces'].sort(key=lambda row: row['face'])
    return ranks, keys


def _project_outputs(result, ranks, keys):
    for row in (result.get('display') or {}).get('faces', ()):
        handle = row['face']
        row['face'] = ranks[handle]
        row['group'] = _generated_name(row['group'], handle, row['face'])
    for row in (result.get('step') or {}).get('regions', ()):
        handle = row['face']
        row['face'] = ranks[handle]
        row['file_region'] = _generated_name(row['file_region'], handle, row['face'])
        row['model_region'] = keys.get(row['model_region'], row['model_region'])


def manifest_by_face_rank(manifest, request):
    if manifest.get('format') != 'catallax-manifest/4' or 'certificate' not in manifest:
        return manifest
    problem = request['problem']
    recipe = problem.get('recipe', problem)
    result = deepcopy(manifest)
    ranks, keys = _project_named_result(result, recipe['result'])
    _project_outputs(result, ranks, keys)
    edit = result.get('edit') or {}
    if 'certificate' in edit:
        edited_ranks, _ = _project_named_result(edit, recipe['result'])
        for row in edit['resolved']:
            if row['kind'] == 'face':
                row['old_handle'] = ranks[row['old_handle']]
                for match in row['matches']:
                    match['face'] = edited_ranks[match['face']]
                row['matches'].sort(key=lambda match: match['face'])
        edit['resolved'].sort(key=lambda row: (row['kind'] != 'face', row['old_handle']))
    return result
