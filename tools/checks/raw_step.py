#!/usr/bin/env python3
"""Audit replay STEP solids directly, without TransferRoots or shape repair."""
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import sys

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.StepShape import StepShape_BrepWithVoids, StepShape_ManifoldSolidBrep
from OCP.StepToTopoDS import StepToTopoDS_Builder
from OCP.Transfer import Transfer_TransientProcess

from checks import CORRECTIONS, RAW_POLICY, ROOT, raw_step_errors, read, sha


def raw_shapes(path):
    reader = STEPControl_Reader()
    settings = {'read.stdsameparameter.mode': 0, 'read.maxprecision.mode': 0}
    previous = {name: Interface_Static.IVal_s(name) for name in settings}
    try:
        for name, value in settings.items():
            if not Interface_Static.SetIVal_s(name, value):
                return []
        if reader.ReadFile(str(path)) != IFSelect_RetDone:
            return []
        model = reader.StepModel()
        solids = [model.Value(i) for i in range(1, model.NbEntities() + 1)
                  if isinstance(model.Value(i), StepShape_ManifoldSolidBrep)]
        shapes = []
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
            if not builder.IsDone():
                return []
            shapes.append(builder.Value())
        return shapes
    finally:
        for name, value in previous.items():
            Interface_Static.SetIVal_s(name, value)


def valid(path):
    shapes = raw_shapes(path)
    return bool(shapes) and all(BRepCheck_Analyzer(shape).IsValid() for shape in shapes)


def corrected_properties_valid(path, manifest):
    if Fraction(manifest['step']['max_deviation']) > Fraction(1, 10_000_000):
        return False
    shapes = raw_shapes(path)
    if len(shapes) != 1 or not BRepCheck_Analyzer(shapes[0]).IsValid():
        return False
    shape = shapes[0]
    values = []
    for epsilon in (1e-6, 1e-7, 1e-8):
        properties = GProp_GProps()
        error = BRepGProp.VolumePropertiesGK_s(
            shape, properties, epsilon, True, False, False, False, False)
        if error < 0 or not math.isfinite(properties.Mass()):
            return False
        values.append(properties.Mass())
    if abs(values[-1] - values[-2]) > max(1, abs(values[-1])) * 1e-10:
        return False
    area = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, area, 1e-12, False)
    spread = max(values) - min(values)
    margin = area.Mass() * 1e-7 + spread
    interval = manifest['properties']['volume']
    return (math.isfinite(margin) and margin >= 0
            and float(interval['lo']) - margin <= values[-1] <= float(interval['hi']) + margin)


def main(report):
    observations, checked = [], 0
    corrections = {row['index']: row for row in read(ROOT / CORRECTIONS)['corrections']}
    property_checks = []
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
                    if index in corrections:
                        try:
                            properties_valid = (
                                sha(path) == corrections[index]['after_step_sha256']
                                and corrected_properties_valid(path, case['response']['manifest']))
                        except Exception:
                            properties_valid = False
                        property_checks.append({'index': index, 'passed': properties_valid})
                else:
                    observations.append({'index': index, 'step_sha256': None, 'valid': None})
                    if index in corrections:
                        property_checks.append({'index': index, 'passed': False})
    finally:
        for fd, original in zip((1, 2), saved):
            os.dup2(original, fd)
            os.close(original)
    policy = read(ROOT / RAW_POLICY)
    errors = raw_step_errors(observations, policy)
    errors += [f"index {row['index']}: corrected STEP property check failed"
               for row in property_checks if not row['passed']]
    print(json.dumps({'raw_step_checked': checked,
                      'raw_step_invalid': [row['index'] for row in observations if row['valid'] is False],
                      'known_raw_invalid': policy['known_raw_invalid'],
                      'correction_properties': property_checks,
                      'raw_step_failures': errors}))
    return bool(errors)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
