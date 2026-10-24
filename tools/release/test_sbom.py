"""Publishing must not leak workspace paths or break SBOM dependency edges."""
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

from sbom import normalize


class SbomTests(unittest.TestCase):
    def fixture(self, workspace, dependency=None):
        prefix = "path+" + workspace.as_uri()
        dependency = dependency or prefix + "/crates/library#0.1.0"
        return f'''<bom xmlns="http://cyclonedx.org/schema/bom/1.3" version="1">
          <metadata><component type="application" bom-ref="{prefix}/crates/cli#0.1.0"><name>catallax-cli</name></component></metadata>
          <components><component type="library" bom-ref="{prefix}/crates/library#0.1.0"><name>library</name></component></components>
          <dependencies><dependency ref="{prefix}/crates/cli#0.1.0"><dependency ref="{dependency}"/></dependency></dependencies>
        </bom>'''

    def test_rewrites_component_ids_and_dependency_edges_together(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            source, out = workspace / "input.xml", workspace / "output.xml"
            source.write_text(self.fixture(workspace))
            normalize(source, out, "a" * 40, workspace)
            text = out.read_text()
            self.assertNotIn(workspace.as_uri(), text)
            tree = ET.fromstring(text)
            ids = {node.get("bom-ref") for node in tree.iter() if node.get("bom-ref")}
            refs = {node.get("ref") for node in tree.iter() if node.get("ref")}
            self.assertEqual(ids, refs)
            self.assertTrue(all("@" + "a" * 40 in reference for reference in ids))

    def test_outside_workspace_or_dangling_dependency_blocks_publish(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            for dependency in ("path+" + workspace.parent.as_uri() + "/external#0.1.0", "missing"):
                with self.subTest(dependency=dependency):
                    source = workspace / "input.xml"
                    source.write_text(self.fixture(workspace, dependency))
                    with self.assertRaises(ValueError):
                        normalize(source, workspace / "output.xml", "a" * 40, workspace)


if __name__ == "__main__":
    unittest.main()
