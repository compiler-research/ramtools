# Sphinx configuration for the RAMTools documentation.
# Build locally with:  python3 -m sphinx -b html docs/source build/docs/html
# The API reference is generated from the headers by Doxygen (run below) and
# rendered by Breathe, so it needs doxygen on PATH.

import os
import subprocess

project = "RAMTools"
copyright = "2025, compiler-research"
author = "compiler-research"

extensions = ["sphinx.ext.todo", "breathe"]
templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"
html_title = "RAMTools"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/compiler-research/ramtools/",
    "source_branch": "develop",
    "source_directory": "docs/source/",
}

highlight_language = "bash"

# Keep punctuation as written: no automatic dashes or curly quotes.
smartquotes = False

# Doxygen XML for Breathe: the repository Doxyfile with XML output turned on,
# limited to the public headers, written under build/.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DOXYGEN_OUT = os.path.join(ROOT_DIR, "build", "doxygen-docs")
with open(os.path.join(ROOT_DIR, "Doxyfile")) as f:
    doxyfile = f.read()
doxyfile += "\n".join(
    [
        "",
        "INPUT = inc",
        "OUTPUT_DIRECTORY = " + DOXYGEN_OUT,
        "GENERATE_HTML = NO",
        "GENERATE_XML = YES",
        "QUIET = YES",
        "WARN_IF_UNDOCUMENTED = NO",
        "",
    ]
)
os.makedirs(DOXYGEN_OUT, exist_ok=True)
subprocess.run(["doxygen", "-"], input=doxyfile.encode(), cwd=ROOT_DIR, check=True)

breathe_projects = {"ramtools": os.path.join(DOXYGEN_OUT, "xml")}
breathe_default_project = "ramtools"
