#!/usr/bin/env python3
"""Replace CycloneDX's machine-local workspace IDs with commit-qualified IDs."""
import argparse
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
REPOSITORY = "https://github.com/BohnBawerick/catallax"


def normalize(source, destination, commit, workspace=ROOT):
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("invalid SBOM source commit")
    tree = ET.parse(source)
    root = tree.getroot()
    namespace = root.tag.split("}")[0].removeprefix("{")
    if not namespace.startswith("http://cyclonedx.org/schema/bom/"):
        raise ValueError("expected a CycloneDX XML SBOM")
    ET.register_namespace("", namespace)
    prefix = "path+" + workspace.resolve().as_uri() + "/"
    replacement = "git+" + REPOSITORY + "@" + commit + "/"
    components = root.findall(f".//{{{namespace}}}component")
    if not components or root.findtext(f"{{{namespace}}}metadata/{{{namespace}}}component/{{{namespace}}}name") != "catallax-cli":
        raise ValueError("expected the catallax-cli dependency inventory")
    for element in root.iter():
        for key, value in element.attrib.items():
            if value.startswith("path+file:"):
                if not value.startswith(prefix):
                    raise ValueError("SBOM includes a path dependency outside the workspace")
                element.set(key, replacement + value.removeprefix(prefix))
    # Every dependency edge must still resolve after IDs and references change.
    ids = {element.get("bom-ref") for element in root.iter() if element.get("bom-ref")}
    for dependency in root.iter(f"{{{namespace}}}dependency"):
        if dependency.get("ref") not in ids:
            raise ValueError("SBOM has a dangling dependency reference")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tree.write(destination, encoding="utf-8", xml_declaration=True)
    print(f"Public CycloneDX SBOM: {len(components)} components, source {commit}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    normalize(args.source, args.destination, args.commit)
