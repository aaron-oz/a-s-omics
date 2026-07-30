#!/usr/bin/env python3
r"""Consistency checks for response-to-reviewers.tex. Run before submitting.

Checks three things that have drifted before:
  1. Every item's \status{} is one of exactly DONE / INCOMPLETE / TODO.
  2. The invariant: a DONE item carries no marker; an INCOMPLETE item carries only
     \incomplete markers; a TODO item carries only \todo markers. The marker word
     always matches the item's status, so the two cannot disagree.
  3. The per-item status in the summary table agrees with the subsection status.

Exit code 0 if clean, 1 if anything is wrong.
Usage: python3 check-response-letter.py [path]
"""
import re, sys

VALID = {'DONE', 'INCOMPLETE', 'TODO'}
TABLE_WORD = {'DONE': 'done', 'INCOMPLETE': 'incomplete', 'TODO': 'to do'}

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
    items[name] = (st, body.count(r'\todo{'), body.count(r'\incomplete{'))

for name, (st, n_todo, n_inc) in items.items():
    if st is None:
        fail.append(f"{name}: no \\status{{}}")
        continue
    if st not in VALID:
        fail.append(f"{name}: status {st!r} is not one of {sorted(VALID)}")
        continue
    if st == 'DONE' and n_todo + n_inc:
        fail.append(f"{name}: DONE but carries {n_todo + n_inc} marker(s)")
    if st != 'DONE' and n_todo + n_inc == 0:
        fail.append(f"{name}: {st} but carries no marker saying what remains")
    if st == 'INCOMPLETE' and n_todo:
        fail.append(f"{name}: INCOMPLETE but carries {n_todo} \\todo marker(s); use \\incomplete")
    if st == 'TODO' and n_inc:
        fail.append(f"{name}: TODO but carries {n_inc} \\incomplete marker(s); use \\todo")

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

n_todo = sum(a for _, (st, a, b) in items.items())
n_inc = sum(b for _, (st, a, b) in items.items())
print(f"checked {len(items)} items in {path}")
print(f"outstanding markers: {n_todo} TODO (not started), {n_inc} INCOMPLETE (started, unfinished)")
if n_todo + n_inc:
    print("  " + ", ".join(sorted(n for n, (st, a, b) in items.items() if a + b)))
    print("  -> these must all be cleared before submission")
if fail:
    print(f"\n{len(fail)} PROBLEM(S):")
    for f in fail:
        print("  -", f)
    sys.exit(1)
print("all clean: closed vocabulary, invariant holds, table agrees with subsections")
