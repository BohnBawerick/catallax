"""Record installed versions and hash the native code actually loaded."""
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import sys


def receipt():
    packages = {}
    for name in ('smlibpy', 'cadquery', 'cadquery-ocp'):
        distribution = metadata.distribution(name)
        native = {}
        for member in distribution.files or ():
            path = Path(distribution.locate_file(member))
            if '.so' in path.name and path.is_file():
                with path.open('rb') as stream:
                    native[str(member)] = hashlib.file_digest(stream, 'sha256').hexdigest()
        packages[name] = {'version': distribution.version, 'native_sha256': native}
    return {'python': sys.version, 'packages': packages}


if __name__ == '__main__':
    print(json.dumps(receipt(), sort_keys=True))
