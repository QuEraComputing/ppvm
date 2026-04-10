"""Auto-generate one API reference page per source module."""

from pathlib import Path

import mkdocs_gen_files

src = Path(__file__).parent.parent.parent / "src"
pkg = src / "ppvm"
nav = mkdocs_gen_files.Nav()

for path in sorted(pkg.rglob("*.py")):
    module_path = path.relative_to(pkg)
    parts = tuple(module_path.with_suffix("").parts)

    # Skip __init__, __main__, and private modules
    if parts[-1] in ("__init__", "__main__") or parts[-1].startswith("_"):
        continue

    module_name = "ppvm." + ".".join(parts)
    doc_path = module_path.with_suffix(".md")
    full_doc_path = Path("api") / "ppvm" / doc_path

    nav[parts] = str(doc_path)

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        fd.write(f"::: {module_name}\n")

# Write SUMMARY.md alongside the generated pages so literate-nav finds it
# when the nav entry is "- api/ppvm/"
with mkdocs_gen_files.open("api/ppvm/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
