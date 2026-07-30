#!/usr/bin/env python3
r"""Consistency checks for response-to-reviewers.tex. Run before submitting.

Checks three things that have drifted before:
  1. Every item's \status{} is one of exactly COMPLETE / PARTIAL / NOT STARTED.
  2. The invariant: an item is COMPLETE if and only if it carries no \action and
     no \authors marker.
  3. The per-item status in the summary table agrees with the subsection status.

Exit code 0 if clean, 1 if anything is wrong.
Usage: python3 check-response-letter.py [path]
"""
import re, sys

VALID = {'COMPLETE', 'PARTIAL', 'NOT STARTED'}
TABLE_WORD = {'COMPLETE': 'complete', 'PARTIAL': 'partial', 'NOT STARTED': 'not started'}

path = sys.argv[1] if len(sys.argv) > 1 else 'response-to-reviewers.tex'
s = open(path).read()
fail = []

# ---- items, their status, and their markers ----
parts = re.split(r'(?m)^(\\subsection\{(?:JR-\d|R\d+\.\d+)\.[^}]*\}|\\section\{Additional corrections\})', s)
items = {}
for i in range(1, len(parts), 2):
    head, body = parts[i], parts[i + 1]
    name = re.search(r'\{((?:JR-\d|R\d+\.\d+)|Additional corrections)', head).group(1)
    m = re.search(r'\\status\{([^}]*)\}', body[:200])
    st = m.group(1) if m else None
    items[name] = (st, body.count(r'\action{'), body.count(r'\authors{'))

for name, (st, a, u) in items.items():
    if st is None:
        fail.append(f"{name}: no \\status{{}}")
        continue
    if st not in VALID:
        fail.append(f"{name}: status {st!r} is not one of {sorted(VALID)}")
        continue
    if st == 'COMPLETE' and a + u:
        fail.append(f"{name}: COMPLETE but carries {a} ACTION and {u} AUTHORS marker(s)")
    if st != 'COMPLETE' and a + u == 0:
        fail.append(f"{name}: {st} but carries no marker saying what remains")

# ---- summary table agreement ----
for line in s.split('\n'):
    m = re.match(r'^\\?t?e?x?t?b?f?\{?((?:JR-\d|R\d+\.\d+))\}?\s*&(.*?)&\s*(.*?)\s*\\\\\s*$', line)
    if not m:
        continue
    key, cell = m.group(1), m.group(3).replace(r'\textbf{', '').replace('}', '').lower()
    if key in items and items[key][0]:
        want = TABLE_WORD[items[key][0]]
        if not cell.startswith(want):
            fail.append(f"{key}: table says {m.group(3)!r}, subsection says {items[key][0]!r}")

print(f"checked {len(items)} items in {path}")
if fail:
    print(f"\n{len(fail)} PROBLEM(S):")
    for f in fail:
        print("  -", f)
    sys.exit(1)
print("all clean: closed vocabulary, invariant holds, table agrees with subsections")
