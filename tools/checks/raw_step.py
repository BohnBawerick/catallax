#!/usr/bin/env python3
"""Audit replay STEP solids directly, without TransferRoots or shape repair."""
import json
import os
from pathlib import Path
import sys

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.StepShape import StepShape_BrepWithVoids, StepShape_ManifoldSolidBrep
from OCP.StepToTopoDS import StepToTopoDS_Builder
from OCP.Transfer import Transfer_TransientProcess

from checks import RAW_POLICY, ROOT, raw_step_errors, read, sha


def valid(path):
    reader = STEPControl_Reader()
    settings = {'read.stdsameparameter.mode': 0, 'read.maxprecision.mode': 0}
    previous = {name: Interface_Static.IVal_s(name) for name in settings}
    try:
        for name, value in settings.items():
            if not Interface_Static.SetIVal_s(name, value):
                return False
        if reader.ReadFile(str(path)) != IFSelect_RetDone:
            return False
        model = reader.StepModel()
        solids = [model.Value(i) for i in range(1, model.NbEntities() + 1)
                  if isinstance(model.Value(i), StepShape_ManifoldSolidBrep)]
        if not solids:
            return False
        for solid in solids:
            process = Transfer_TransientProcess()
            process.SetModel(model)
            builder = StepToTopoDS_Builder()
            builder.SetPrecision(1e-7)
            builder.SetMaxTol(1e-7)
            if isinstance(solid, StepShape_BrepWithVoids):
                builder.Init(theBRepWithVoids=solid, theTP=process)
            else:
                builder.Init(solid, process)
            if not builder.IsDone() or not BRepCheck_Analyzer(builder.Value()).IsValid():
                return False
        return True
    finally:
        for name, value in previous.items():
            Interface_Static.SetIVal_s(name, value)


def main(report):
    observations, checked = [], 0
    # OCCT can print entity diagnostics through native stdout/stderr. Keep them
    # off public logs, including paths and geometry from the private corpus.
    saved = [os.dup(fd) for fd in (1, 2)]
    try:
        with open(os.devnull, 'w') as sink:
            for fd in (1, 2):
                os.dup2(sink.fileno(), fd)
            for index, case in enumerate(json.loads(Path(report).read_text())['cases']):
                path = Path(case['path']) / 'answer.step'
                if path.exists():
                    checked += 1
                    try:
                        passed = valid(path)
                    except Exception:
                        passed = False
                    observations.append({'index': index, 'step_sha256': sha(path), 'valid': passed})
                else:
                    observations.append({'index': index, 'step_sha256': None, 'valid': None})
    finally:
        for fd, original in zip((1, 2), saved):
            os.dup2(original, fd)
            os.close(original)
    policy = read(ROOT / RAW_POLICY)
    errors = raw_step_errors(observations, policy)
    print(json.dumps({'raw_step_checked': checked,
                      'raw_step_invalid': [row['index'] for row in observations if row['valid'] is False],
                      'known_raw_invalid': policy['known_raw_invalid'],
                      'raw_step_failures': errors}))
    return bool(errors)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
