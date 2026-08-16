#!/usr/bin/env python3
"""Corrected launcher for scCRS T/B-state association scoring.

Compared with v2, this version also classifies plasmablast/plasma-cell and
antibody-secreting-cell labels as B lineage, so B-cell differentiation states
can be scored for labels such as 'IgA-secreting plasma cells'.
"""
from pathlib import Path

source = Path(__file__).with_name("scCRS_state_cytokine_association.py")
code = source.read_text(encoding="utf-8")
code = code.replace(
    'if "b cell" in x or x in {"b", "b cells"}: return "B"',
    'if "b cell" in x or x in {"b", "b cells"} or any(token in x for token in ("plasma", "plasmablast", "antibody-secreting", "asc")): return "B"'
)
code = code.replace(
    'if "t" in x or "mait" in x: return "T"',
    'if any(token in x for token in ("cd4", "cd8", "t cell", "mait", "gdt", "gamma delta")): return "T"'
)
exec(compile(code, str(source), "exec"))
