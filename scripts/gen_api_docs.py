"""
Regenerate docs/API.md from the source tree.

Walks src/syndrumnet with ast, collects every public class and function with the
first line of its docstring, and emits the markdown reference. Parsing rather
than importing keeps this runnable without the pipeline's heavier dependencies
installed.

Usage
-----
    python scripts/gen_api_docs.py > docs/API.md
"""

import ast
import sys
from pathlib import Path

# Subpackages in the order the pipeline uses them.
PACKAGE_ORDER = [
    "io",
    "data",
    "metrics",
    "propagation",
    "scoring",
    "eval",
    "viz",
    "utils",
]

HEADER = """# API Reference

Every public class and function in `syndrumnet`, grouped by subpackage in the
order the pipeline uses them: load and parse data, build the network and
modules, compute metrics, propagate, score, evaluate, visualise.

This page is a map. Full parameter and return documentation lives in the
docstrings themselves, in NumPy style, and is the authoritative source:

```python
from syndrumnet.propagation.prince import PRINCE
help(PRINCE)
```

Regenerate this file after changing the public surface:

```bash
python scripts/gen_api_docs.py > docs/API.md
```

---
"""

FOOTER = """---

## Entry points

The four scripts in `scripts/` wire the above together. Each takes
`--config <path>` and accepts dotted overrides such as
`--propagation.alpha 0.7`.

| Script | Stage |
|---|---|
| `build_all_data.py` | Download sources, map IDs, build the network and modules |
| `run_pipeline.py` | Compute TQAB, PQAB and CQAB for every drug pair |
| `evaluate.py` | Score predictions against known synergies (AUC-ROC, AUC-PR) |
| `make_figures.py` | Render figures into `reports/figures/` |
"""


def first_docstring_line(node: ast.AST) -> str:
    """Return the first line of a node's docstring, or an empty string."""
    doc = ast.get_docstring(node)
    return doc.strip().split("\n")[0] if doc else ""


def collect(module_path: Path) -> list[str]:
    """Render the public surface of one module as markdown bullet lines."""
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    lines: list[str] = []

    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            continue
        if node.name.startswith("_"):
            continue

        summary = first_docstring_line(node)

        if isinstance(node, ast.ClassDef):
            lines.append(f"- **`{node.name}`** - {summary}")
            methods = [
                m.name
                for m in node.body
                if isinstance(m, ast.FunctionDef) and not m.name.startswith("_")
            ]
            if methods:
                rendered = ", ".join(f"`{m}()`" for m in methods)
                lines.append(f"  - methods: {rendered}")
        else:
            args = ", ".join(a.arg for a in node.args.args)
            lines.append(f"- `{node.name}({args})` - {summary}")

    return lines


def main() -> int:
    root = Path(__file__).resolve().parent.parent / "src" / "syndrumnet"
    if not root.is_dir():
        print(f"package not found at {root}", file=sys.stderr)
        return 1

    by_package: dict[str, list[tuple[str, list[str]]]] = {}

    for path in sorted(root.rglob("*.py")):
        if path.name == "__init__.py":
            continue
        parts = path.relative_to(root).parts
        if len(parts) < 2:
            continue  # top-level module, not part of a subpackage
        entries = collect(path)
        if entries:
            dotted = "syndrumnet." + ".".join(parts).removesuffix(".py")
            by_package.setdefault(parts[0], []).append((dotted, entries))

    out = [HEADER]
    for package in PACKAGE_ORDER:
        if package not in by_package:
            continue
        out.append(f"## `syndrumnet.{package}`\n")
        for dotted, entries in by_package[package]:
            out.append(f"### `{dotted}`\n")
            out.extend(entries)
            out.append("")
    out.append(FOOTER)

    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
