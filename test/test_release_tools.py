from pathlib import Path
import sys
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from check_release_version import read_project_version, validate_release_tag  # noqa: E402
from update_bioconda_recipe import update_recipe_text  # noqa: E402


SAMPLE_RECIPE = """{% set name = "replidec" %}
{% set version = "0.3.5" %}

source:
  url: https://pypi.org/packages/source/{{ name[0] }}/{{ name }}/replidec-{{ version }}.tar.gz
  sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

build:
  number: 4

test:
  commands:
    - Replidec --help

about:
  doc_url: "https://github.com/pengSherryYel/Replidec/blob/v.{{ version }}/README.md"
"""


class ReleaseVersionTests(unittest.TestCase):
    def test_project_version_is_a_valid_release_version(self):
        version = read_project_version(PROJECT_ROOT)
        self.assertEqual(validate_release_tag(f"v{version}", version), version)

    def test_release_tag_validation(self):
        self.assertEqual(validate_release_tag("v0.3.6", "0.3.6"), "0.3.6")

    def test_release_tag_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not match"):
            validate_release_tag("v0.3.5", "0.3.6")


class BiocondaRecipeTests(unittest.TestCase):
    def test_recipe_update(self):
        sha256 = "b" * 64
        updated = update_recipe_text(SAMPLE_RECIPE, "0.3.6", sha256)

        self.assertIn('{% set version = "0.3.6" %}', updated)
        self.assertIn(f"sha256: {sha256}", updated)
        self.assertIn("number: 0", updated)
        self.assertIn("    - Replidec --version", updated)
        self.assertIn("/blob/v{{ version }}/README.md", updated)

    def test_recipe_update_preserves_command_indentation(self):
        sha256 = "b" * 64

        for indent in (" ", "    ", "\t"):
            with self.subTest(indent=repr(indent)):
                recipe = SAMPLE_RECIPE.replace(
                    "    - Replidec --help",
                    f"{indent}- Replidec --help",
                )
                updated = update_recipe_text(recipe, "0.3.6", sha256)

                self.assertIn(f"{indent}- Replidec --help", updated)
                self.assertIn(f"{indent}- Replidec --version", updated)

    def test_invalid_sha_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "SHA256"):
            update_recipe_text(SAMPLE_RECIPE, "0.3.6", "not-a-sha")


if __name__ == "__main__":
    unittest.main()
