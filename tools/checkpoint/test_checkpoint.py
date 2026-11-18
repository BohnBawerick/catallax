import json
from pathlib import Path
import unittest

from checkpoint import census_markdown


class CensusMarkdownTests(unittest.TestCase):
    def test_collected_rows_do_not_imply_baseline_acceptance(self):
        report = json.loads(
            (Path(__file__).parent / 'fixtures/collected-unmeasured.json').read_text())
        markdown = census_markdown(report)
        self.assertIn('Collected: True. Evaluated: 3000/3000.', markdown)
        self.assertIn('Acceptance passed: False.', markdown)
        self.assertIn('| Tier losses | 0 |', markdown)
        self.assertIn('| Frozen valid-export losses | 0 |', markdown)
        self.assertIn('| Baseline not measured | 1 |', markdown)
        self.assertNotIn('Complete: True', markdown)


if __name__ == '__main__':
    unittest.main()
