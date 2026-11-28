"""Generate saved-view STEP fixtures and evaluate the full AP242 Edition 1 rules."""

import os
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory


def validate(path, directory, classpath):
    log = directory / f"{path.stem}.log"
    subprocess.run(
        ["java", "-cp", classpath, "jsdai.util.Validate", str(path), "-stdout", str(log)],
        cwd=directory,
        check=True,
    )
    text = log.read_text()
    counts = []
    for label in ("erroneous instances", "violated global rules", "violated uniqueness rules"):
        results = re.findall(rf"^count of {label}:\s*(\d+)\s*$", text, re.MULTILINE)
        if not results:
            raise AssertionError(f"{path.name}: missing validator result for {label}\n{text}")
        counts.extend(map(int, results))
    global_rules = re.findall(r'^count of global rules in schema instance .*?:\s*(\d+)\s*$', text, re.MULTILINE)
    if not global_rules or not all(map(int, global_rules)):
        raise AssertionError(f"{path.name}: no EXPRESS global rules evaluated\n{text}")
    return counts, text


def main():
    classpath = os.environ.get("JSDAI_AP242_E1_CLASSPATH")
    if not classpath:
        raise SystemExit("JSDAI_AP242_E1_CLASSPATH must provide JSDAI Validate and the compiled AP242 Edition 1 schema")
    classpath = os.pathsep.join(str(Path(entry).absolute()) for entry in classpath.split(os.pathsep))
    root = Path(__file__).resolve().parents[2]
    target = root / "target"
    target.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="ap242-views-", dir=target) as temp:
        directory = Path(temp)
        repos = directory / "repositories"
        repos.mkdir()
        (directory / "jsdai.properties").write_text(f"repositories={repos.as_posix()}\n")
        for package, test, case in (
            ("catallax-step", "pmi", "saved_views_write_their_callouts_and_read_back"),
            ("catallax-import", "assembly_step", "a_leaf_datum_target_round_trips_through_assembly_discovery"),
            ("catallax-import", "pmi_patch", "a_tolerance_value_changes_its_own_records_only"),
        ):
            subprocess.run(
                ["cargo", "test", "--locked", "-p", package, "--test", test, case, "--", "--exact"],
                cwd=root, env=dict(os.environ, CATALLAX_PMI_OUT=str(directory)), check=True,
            )
        graphical = root / "crates/catallax-step/tests/fixtures/graphical-only-pmi.step"
        invalid = directory / "unlinked-camera.step"
        invalid.write_text("\n".join(line for line in graphical.read_text().splitlines()
                                     if "=ITEM_DEFINED_TRANSFORMATION(" not in line) + "\n")
        counts, log = validate(invalid, directory, classpath)
        if not any(counts) or not re.search(r"camera_model[^\n]*violated", log, re.IGNORECASE):
            raise AssertionError(f"validator missed the deliberately unlinked camera\n{log}")
        for artifact in (
            graphical,
            directory / "c4c-bore-view.step",
            directory / "c4c-top-view.step",
            directory / "c4c-saved-views.step",
            directory / "c4c-assembly-views.step",
            directory / "c4d-source.step",
            directory / "c4d-flatness-0.02.step",
        ):
            if "{ 1 0 10303 442 1 1 4 }" not in artifact.read_text().split("\nDATA;\n", 1)[0]:
                raise AssertionError(f"{artifact.name}: expected AP242 Edition 1")
            counts, log = validate(artifact, directory, classpath)
            if any(counts):
                raise AssertionError(f"{artifact.name}: EXPRESS validation failed\n{log}")


if __name__ == "__main__":
    main()
