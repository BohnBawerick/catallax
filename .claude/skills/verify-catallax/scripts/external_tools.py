"""Additional released-tool engines for the existing head-to-head command."""
import importlib.util
import json
from pathlib import Path

from head_to_head import DATA, process, save, sha


class Graph:
    """Serialize authored target operations for the native Dart runner."""
    def __init__(self):
        self.nodes = []

    def add(self, op, *args):
        self.nodes.append({'op': op, 'args': args})
        return len(self.nodes) - 1

    def box(self, size, corner=(0, 0, 0)):
        return self.add('box', size, corner)

    def cylinder(self, radius, height, base=(0, 0, 0), axis=(0, 0, 1)):
        return self.add('cylinder', radius, height, base, axis)

    def sphere(self, radius, center=(0, 0, 0)):
        return self.add('sphere', radius, center)

    def cone(self, radius, top, height, base=(0, 0, 0), axis=(0, 0, 1)):
        return self.add('cone', radius, top, height, base, axis)

    def torus(self, major, minor):
        return self.add('torus', major, minor)

    def rotate_translate(self, body, quaternion, offset):
        return self.add('rotate_translate', body, quaternion, offset)

    def union(self, a, b):
        return self.add('union', a, b)

    def cut(self, a, b):
        return self.add('cut', a, b)

    def common(self, a, b):
        return self.add('common', a, b)

    def blend(self, part, points, radius, chamfer=False, curve=None):
        return self.add('blend', part, points, radius, chamfer, curve)

    def shell_top(self, part, thickness):
        return self.add('shell_top', part, thickness)


def module(name):
    spec = importlib.util.spec_from_file_location(name, DATA / f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def runtime(path):
    value = json.loads(Path(path).read_text())
    for key, filename in (('pins_sha256', 'tools.json'), ('requirements_sha256', 'requirements.lock'),
                          ('runner_sha256', 'pali_parts.dart')):
        if value[key] != sha(DATA / filename):
            raise ValueError(f'tool installation has stale {filename}; rerun prepare.py')
    if value['pali_binary_sha256'] != sha(value['pali_binary']):
        raise ValueError('Pali benchmark executable checksum mismatch')
    prepare = module('prepare')
    for key, root in (('python', Path(value['python']).parents[1] / 'lib'),
                      ('checker', Path(value['checker']).parents[1] / 'lib')):
        if value[key + '_tree_sha256'] != prepare.tree_hash(root):
            raise ValueError(f'{key} installed files changed since checksum-verified preparation')
    return value


def external(engine, case_path, case, directory, args):
    receipt = args.external_runtime
    if engine == 'pali_ijo':
        targets = module('targets')
        request = {'part': case['recipe']}
        for label in ('construction', 'edit'):
            graph = Graph()
            targets.build(graph, case['recipe'], float(case['edit']['value']) if label == 'edit' else None)
            request[label] = graph.nodes
        save(directory / 'request.json', request)
        command = [receipt['pali_binary'], directory / 'request.json', directory]
    else:
        command = [receipt['python'], '-B', DATA / 'python_tools.py', engine, case_path.resolve(), directory]
    execution = process(command, directory, args.core, args.timeout)
    partial = directory / 'phases.json'
    result = json.loads(partial.read_text()) if partial.exists() else {'variants': {}}
    result['process'] = execution
    result['status'] = 'completed' if execution['exit'] == 0 else 'incomplete' if execution['status'] == 'timeout' else 'failed-rebuild'
    # A crash before a checkpoint is a failed rebuild, not an unsupported part.
    for label in ('construction', 'edit'):
        variant = result['variants'].setdefault(label, {})
        if 'build' not in variant:
            variant['build'] = {'status': 'incomplete' if execution['status'] == 'timeout' else 'command-failed',
                                'reason': 'process ended before complete build checkpoint'}
    return result
