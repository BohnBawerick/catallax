#!/usr/bin/env python3
"""Prepare private M5 inputs from the pinned R8 corpus and unchanged P0 recipes."""

import argparse
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from r8_step_census import ARCHIVE_BYTES, ARCHIVE_MD5, SEED, choose, digest, restrict


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_SHA256 = "52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe"
# First ten complete admissions below 1M import work in the published R8 census.
MODEL_IDS = ("00000131", "00000173", "00000190", "00000405", "00000511",
             "00000549", "00000605", "00000632", "00000633", "00000695")
LENGTHS = {"block": ("size",), "cylinder": ("radius", "height"),
           "sphere": ("radius",), "cone": ("base_radius", "top_radius", "height"),
           "torus": ("major", "minor")}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def quaternion_product(a, b):
    w, x, y, z = map(Fraction, a)
    v, i, j, k = map(Fraction, b)
    return [str(n) for n in (w*v-x*i-y*j-z*k, w*i+x*v+y*k-z*j,
                             w*j-x*k+y*v+z*i, w*k+x*j-y*i+z*v)]


def p0_copy(recipe, variant):
    """Change geometry exactly, or move only the periodic charts of round stock."""
    result = copy.deepcopy(recipe)
    result["name"] += f" / {variant}"
    for step in result["steps"]:
        if step["op"] == "boolean":
            continue
        dimensions = LENGTHS[step["op"]]
        pose = step["pose"]
        if variant == "translated":
            pose["translation"] = [str(Fraction(x) + delta) for x, delta in
                                   zip(pose["translation"], (101, -37, 19))]
        elif variant == "scaled":
            for key in dimensions:
                value = step[key]
                step[key] = ([str(2 * Fraction(x)) for x in value] if isinstance(value, list)
                             else str(2 * Fraction(value)))
            pose["translation"] = [str(2 * Fraction(x)) for x in pose["translation"]]
        elif variant == "seam-shifted":
            if step["op"] != "block":
                pose["rotation"] = quaternion_product(pose["rotation"], (3, 0, 0, 4))
        else:
            raise ValueError(f"unknown P0 copy: {variant}")
    return result


def prepare_p0(scratch):
    inputs = scratch / "p0-inputs"
    inputs.mkdir()
    records = []
    for path in sorted((ROOT / "parts/P0").glob("*.json")):
        source = json.loads(path.read_text())
        for variant in ("base", "translated", "scaled", "seam-shifted"):
            name = f"{path.stem}-{variant}.json"
            if variant == "base":
                shutil.copyfile(path, inputs / name)
            else:
                write_json(inputs / name, p0_copy(source, variant))
            records.append({"input": name, "source": path.name,
                            "variant": variant, "sha256": digest(inputs / name)})
    write_json(scratch / "p0-selection.json", records)


def import_wrapper(data, model_id):
    sha = hashlib.sha256(data).hexdigest()
    return {"format": "catallax-import-recipe/1", "files": {sha: data.decode("utf-8")},
            "recipe": {"format": "catallax-recipe/3", "name": f"ABC {model_id}",
                       "steps": [{"id": "imported", "op": "import_step",
                                  "file": {"name": "part.step", "sha256": sha},
                                  "pose": {"rotation": ["1", "0", "0", "0"],
                                           "translation": ["0", "0", "0"]}}],
                       "result": "imported"}}


def prepare_real(archive, scratch):
    if (archive.stat().st_size != ARCHIVE_BYTES or digest(archive, "md5") != ARCHIVE_MD5
            or digest(archive) != ARCHIVE_SHA256):
        raise ValueError("pinned R8 archive integrity mismatch")
    entries = subprocess.check_output(["bsdtar", "-tf", str(archive)], text=True).splitlines()
    selected, total = choose(entries, 3000, SEED)
    selected = restrict(selected, MODEL_IDS)
    extracted = scratch / "extracted"
    inputs = scratch / "real-inputs"
    extracted.mkdir()
    inputs.mkdir()
    members = scratch / "extract-list-private.txt"
    members.write_text("".join(entry + "\n" for _, entry in selected))
    records = []
    try:
        subprocess.run(["bsdtar", "-xf", str(archive), "-C", str(extracted),
                        "-T", str(members)], check=True)
        for model_id, member in selected:
            path = extracted / member
            if not path.is_file() or path.is_symlink():
                raise ValueError("missing or non-regular corpus model")
            data = path.read_bytes()
            wrapper = import_wrapper(data, model_id)
            write_json(inputs / f"{model_id}.json", wrapper)
            records.append({"model_id": model_id, "source_sha256": hashlib.sha256(data).hexdigest()})
    finally:
        shutil.rmtree(extracted)
        members.unlink()
    write_json(scratch / "real-selection.json", {"archive_sha256": ARCHIVE_SHA256,
               "sample_seed": SEED, "chunk_models": total, "models": records})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("p0", "real"))
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    scratch = args.scratch.resolve()
    if scratch.is_relative_to(ROOT):
        parser.error("private inputs must stay outside the checkout")
    scratch.mkdir(parents=True, exist_ok=True)
    if args.kind == "p0":
        prepare_p0(scratch)
    elif args.archive is None:
        parser.error("real parts require --archive")
    else:
        prepare_real(args.archive.resolve(), scratch)


if __name__ == "__main__":
    main()
