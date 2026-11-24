"""Build a streamlined change copy: paragraph-level alignment of the original submission (old.tex)
against the current manuscript (new.tex). Unchanged paragraphs are collapsed, lightly edited ones
get word-level latexdiff markup, and rewritten or new ones are shown whole (new in blue, the
original it replaced in red), so that no paragraph is a word-by-word interleaving of two texts.
Writes oldstar.tex and newstar.tex for latexdiff, then (second call, 'post') post-processes."""
import re, sys, difflib, json
import os
T = os.environ.get("CHANGE_TMP", os.path.expanduser("~/.cache/a-s-omics-change-copy"))   # scratch dir; see README-build.md section 2c
EDIT_MIN = 0.62      # share of words in common above which a paragraph counts as an edit of the old one
REPL_MIN = 0.25      # above this (below EDIT_MIN) the new paragraph is shown as replacing the old one

def body(lines):
    b = next(i for i, l in enumerate(lines) if "\\begin{document}" in l)
    e = next(i for i, l in enumerate(lines) if "\\end{document}" in l)
    return b, e
def is_prose(l):
    t = l.strip()
    if len(t) < 80 or t.startswith("%") or t.startswith("&"): return False
    return (not t.startswith("\\")) or t.startswith("\\caption") or t.startswith("\\item")
def words(l):
    return re.sub(r"\\[a-zA-Z]+\*?|[{}$~\[\]]", " ", l).split()
def esc(s):
    for a, b in (("%", "\\%"), ("&", "\\&"), ("_", "\\_"), ("#", "\\#"), ("^", "")): s = s.replace(a, b)
    return s

if sys.argv[1] == "pre":
    old = open(T + "/old2.tex").read().split("\n"); new = open(T + "/new2.tex").read().split("\n")
    ob, oe = body(old); nb, ne = body(new)
    def prose_idx(lines, b, e):
        out = []; depth = 0
        for i in range(b, e):
            l = lines[i]
            opens = len(re.findall(r"\\begin\{(equation|align\*?|alignat\*?|tabular|enumerate_never)\}", l))
            closes = len(re.findall(r"\\end\{(equation|align\*?|alignat\*?|tabular|enumerate_never)\}", l))
            if depth == 0 and opens == 0 and is_prose(l) and not l.rstrip().endswith("\\\\"): out.append(i)
            depth += opens - closes
        return out
    oprose = prose_idx(old, ob, oe); nprose = prose_idx(new, nb, ne)
    ow = {i: words(old[i]) for i in oprose}; nw = {i: words(new[i]) for i in nprose}
    # caption pairs by label
    def lab(lines, i):
        for k in range(i, min(i + 3, len(lines))):
            m = re.search(r"\\label\{([^}]*)\}", lines[k])
            if m: return m.group(1)
    olab = {lab(old, i): i for i in oprose if old[i].strip().startswith("\\caption")}
    cand = []
    for j in nprose:
        for i in oprose:
            sm = difflib.SequenceMatcher(None, ow[i], nw[j], autojunk=False)
            if sm.quick_ratio() < REPL_MIN: continue
            r = sm.ratio()
            if r >= REPL_MIN: cand.append((r, j, i))
    cand.sort(reverse=True)
    mo, mn = {}, {}
    for j in nprose:                      # captions first, by label
        if new[j].strip().startswith("\\caption"):
            i = olab.get(lab(new, j))
            if i is not None:
                r = difflib.SequenceMatcher(None, ow[i], nw[j], autojunk=False).ratio()
                mn[j] = (i, r); mo[i] = j
    for r, j, i in cand:
        if j in mn or i in mo: continue
        mn[j] = (i, r); mo[i] = j
    cls = {}
    for j in nprose:
        if j not in mn: cls[j] = "new"
        else:
            r = mn[j][1]
            cls[j] = "same" if old[mn[j][0]].strip() == new[j].strip() else ("edit" if r >= EDIT_MIN else "repl")
    removed = [i for i in oprose if i not in mo]
    oldheads = set(re.findall(r"\\(?:sub)*section\*\{(.*?)\}\s*(?:\\label|$)", "\n".join(old), re.M))
    # special non-prose lines whose change should show (title, affiliation, email)
    def findold(key): return next((l for l in old if key in l), None)
    special = {"Assessing differences on maps": findold("Morphogenic Field Analysis}"),
               "OZ Labs, State College": findold("Bucknell University"), "aaron@ozlabs.tech": findold("bucknell.edu")}
    def wrap(l, tag):
        t = l.strip()
        if t.startswith("\\caption{"):
            m = re.match(r"(\\caption\{\{\\bf [^}]*\}\s*)(.*)\}$", t)
            return ("    " + m.group(1) + "{" + tag + " " + m.group(2) + "}}") if m else l
        if t.startswith("\\item"): return "    \\item {" + tag + " " + t[5:].strip() + "}"
        return "{" + tag + " " + t + "\\par}"
    def oldtext(i):
        t = old[i].strip()
        t = re.sub(r"^\\caption\{\{\\bf [^}]*\}\s*", "", t); t = re.sub(r"\}$", "", t) if old[i].strip().startswith("\\caption") else t
        t = re.sub(r"^\\item\s*", "", t)
        return re.sub(r"\\label\{[^}]*\}", "", t)
    ostar, nstar = old[:ob], new[:nb]
    ostar = list(nstar)                    # same preamble on both sides: only the body is compared
    emitted = set(); stats = []; sec = "front matter"; meta = {}
    removed_sorted = sorted(removed)
    for j in range(nb, ne):
        l = new[j]
        hm = re.match(r"\s*\\((?:sub)*section)\*\{(.*?)\}\s*(\\label.*)?$", l)
        if hm:
            sec = hm.group(2)
            if hm.group(2) not in oldheads and not hm.group(1).startswith("subsub"):
                l = l.replace("{" + hm.group(2) + "}", "{" + hm.group(2) + " \\NEWHEAD}", 1)
            ostar.append(l); nstar.append(l); continue
        if j in cls:
            c = cls[j]
            if j in mn:                    # removed old paragraphs that came before this one's partner
                for i in removed_sorted:
                    if i < mn[j][0] and i not in emitted and c in ("same", "edit"):
                        emitted.add(i); blk = "{\\OLDPAR{Removed} " + oldtext(i) + "\\par}"
                        ostar.append(blk); nstar.append(blk); ostar.append(""); nstar.append("")
                        stats.append((sec, "removed", len(ow[i])))
            stats.append((sec, c, len(nw[j])))
            if c == "same":
                stub = "{\\SAMEPAR{" + str(len(new[j].split())) + "}{" + esc(" ".join(nw[j][:10])) + "}}"
                if l.strip().startswith("\\caption") or l.strip().startswith("\\item"): ostar.append(l); nstar.append(l)
                else: ostar.append(stub); nstar.append(stub)
            elif c == "edit":
                ostar.append(old[mn[j][0]]); nstar.append(l)
            elif c == "repl":
                i = mn[j][0]
                if l.strip().startswith("\\caption") or l.strip().startswith("\\item"):
                    w = wrap(l, "\\NEWPAR{Rewritten}").rstrip()
                    w = w[:-2] + " \\OLDINLINE{" + oldtext(i) + "}}}" if l.strip().startswith("\\caption") else w
                    ostar.append(w); nstar.append(w)
                else:
                    blk = "{\\OLDPAR{Original, replaced by the next paragraph} " + oldtext(i) + "\\par}"
                    w = wrap(l, "\\NEWPAR{Rewritten}")
                    for x in (blk, "", w): ostar.append(x); nstar.append(x)
            else:
                w = wrap(l, "\\NEWPAR{New}"); ostar.append(w); nstar.append(w)
            continue
        key = next((k for k in special if k in l), None)
        if key and special[key]: ostar.append(special[key]); nstar.append(l); continue
        ostar.append(l); nstar.append(l)
    rest = [i for i in removed_sorted if i not in emitted]
    tail = ["", "\\section*{Paragraphs removed by Claude with no counterpart in the current text}", ""]
    for i in rest:
        tail += ["{\\OLDPAR{Removed} " + oldtext(i) + "\\par}", ""]; stats.append(("(removed, listed at the end)", "removed", len(ow[i])))
    end = new[ne:]
    # put the removed list before the bibliography
    def place(lst):
        k = next((x for x in range(len(lst) - 1, -1, -1) if "\\bibliography{" in lst[x] or "\\nolinenumbers" in lst[x]), len(lst))
        return lst[:k] + tail + lst[k:]
    open(T + "/oldstar2.tex", "w").write("\n".join(place(ostar) + end))
    open(T + "/newstar2.tex", "w").write("\n".join(place(nstar) + end))
    json.dump(stats, open(T + "/stats2.json", "w"))
    from collections import Counter
    print(Counter(c for _, c, _ in stats))

else:   # post: fix up latexdiff output
    s = open(T + "/diffstar2.tex").read()
    s = "\n".join(l for l in s.split("\n") if not re.search(r"\\includegraphics|\\figimg\[", l))
    stats = json.load(open(T + "/stats2.json"))
    order = []; agg = {}
    for sec, c, n in stats:
        if sec not in agg: agg[sec] = {"same": 0, "edit": 0, "repl": 0, "new": 0, "removed": 0}; order.append(sec)
        agg[sec][c] += n
    rows = ["\\begin{center}\\footnotesize\\begin{tabular}{p{0.44\\linewidth}rrrrr}\\hline Section (words of prose) & unchanged & edited & rewritten & new & removed\\\\\\hline"]
    for sec in order:
        a = agg[sec]
        if sum(a.values()) < 20: continue
        rows.append(esc(re.sub(r"\\[a-zA-Z]+", "", sec))[:62] + " & " + " & ".join(f"{a[k]:,}" for k in ("same", "edit", "repl", "new", "removed")) + "\\\\")
    tot = {k: sum(agg[s][k] for s in agg) for k in ("same", "edit", "repl", "new", "removed")}
    rows.append("\\hline\\textbf{Total} & " + " & ".join(f"{tot[k]:,}" for k in ("same", "edit", "repl", "new", "removed")) + "\\\\\\hline\\end{tabular}\\end{center}")
    macros = r"""
\newcommand{\NEWPAR}[1]{\color{blue}\textbf{[#1]} }
\newcommand{\OLDPAR}[1]{\color{red}\small\textbf{[#1]} }
\newcommand{\OLDINLINE}[1]{{\color{red}\textbf{[Original caption]} #1}}
\newcommand{\NEWHEAD}{{\normalsize\normalfont\color{blue} [new section]}}
\newcommand{\SAMEPAR}[2]{\small\color{gray}[Unchanged paragraph, #1 words: ``#2 \dots'']\par}
"""
    note = r"""
\begin{center}\fbox{\parbox{0.94\linewidth}{\small\textbf{How to read this copy.} It shows only what Claude changed in \texttt{paper-plos-voice.tex} between Aaron's last saved edit on 2026-11-23 and the morning of 2026-11-24 (his bracketed comments addressed, the editing pass, and the follow-up fixes), paragraph by paragraph.
\textcolor{gray}{Gray one-liners} are paragraphs that did not change.
A paragraph in black with \textcolor{blue}{blue underlined} and \textcolor{red}{red struck-through} words is an edit of an original paragraph: black plus red is the sentence as Aaron left it, black plus blue is the current one.
A paragraph entirely in \textcolor{blue}{blue} is new, or was rewritten so heavily that marking single words would be unreadable; where it replaced an original paragraph, the version Aaron left is printed whole in \textcolor{red}{red} just before it.
Paragraphs that were removed are printed whole in red. Figure images, tables and equations are not compared. This is a reading aid, not the marked-up copy for the journal.}}\end{center}

\noindent{\small\textbf{Where the changes are.} Words of running prose in each section, by what happened to the paragraph they are in.}
"""
    s = s.replace("\\begin{document}", macros + "\\begin{document}", 1)
    s = s.replace("\\section*{Abstract}", note + "\n".join(rows) + "\n\n\\section*{Abstract}", 1)
    open("paper-plos-voice-claude-edits.tex", "w").write(s)
    print("written; add/del marks:", len(re.findall(r"\\DIF(add|del)\{", s)))
