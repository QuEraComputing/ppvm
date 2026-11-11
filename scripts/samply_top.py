#!/usr/bin/env python3
"""Summarize a samply profile as the top functions by self and inclusive sample
time — a headless substitute for the flame-graph UI.

Record with presymbolication so symbols are available without the UI:
    samply record --save-only --unstable-presymbolicate -o profile.json.gz -- <cmd>

Then:
    scripts/samply_top.py profile.json.gz [top_n]

If a `<profile>.syms.json` sidecar sits next to the profile, frames are
symbolicated from it (otherwise raw addresses are shown).
"""
import bisect
import gzip
import json
import os
import sys
from collections import Counter


def load(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return json.load(f)


def sidecar_path(profile_path):
    base = profile_path[:-3] if profile_path.endswith(".gz") else profile_path
    cand = base + ".syms.json"
    return cand if os.path.exists(cand) else None


def build_symbolizer(sidecar):
    """Return addr -> name. Primary: per-lib known_addresses (exact sampled
    addresses samply already resolved). Fallback: symbol_table rva ranges."""
    st = sidecar["string_table"]
    exact = {}
    ranges = []  # (start, end, name)
    for lib in sidecar["data"]:
        symtab = lib["symbol_table"]
        for e in symtab:
            ranges.append((e["rva"], e["rva"] + e["size"], st[e["symbol"]]))
        for addr, idx in lib.get("known_addresses", []):
            exact[addr] = st[symtab[idx]["symbol"]]
    ranges.sort()
    starts = [r[0] for r in ranges]

    def resolve(addr):
        if addr in exact:
            return exact[addr]
        i = bisect.bisect_right(starts, addr) - 1
        if i >= 0 and addr < ranges[i][1]:
            return ranges[i][2]
        return f"0x{addr:x}"

    return resolve


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: samply_top.py <profile.json.gz> [top_n]")
    path = sys.argv[1]
    top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    prof = load(path)

    sc = sidecar_path(path)
    resolve = build_symbolizer(load(sc)) if sc else (lambda a: f"0x{a:x}")

    self_c, incl_c = Counter(), Counter()
    total = 0.0
    for t in prof["threads"]:
        s = t["samples"]
        prefix = t["stackTable"]["prefix"]
        st_frame = t["stackTable"]["frame"]
        fr_addr = t["frameTable"]["address"]
        stacks = s["stack"]
        weights = s.get("weight") or [1] * s["length"]

        def name_of(stack_idx):
            return resolve(fr_addr[st_frame[stack_idx]])

        for i in range(s["length"]):
            node = stacks[i]
            if node is None:
                continue
            w = weights[i] or 1
            total += w
            self_c[name_of(node)] += w
            seen = set()
            while node is not None:
                nm = name_of(node)
                if nm not in seen:  # count each function once per sample
                    seen.add(nm)
                    incl_c[nm] += w
                node = prefix[node]

    if total == 0:
        sys.exit("no samples in profile")

    def short(nm):
        # Drop balanced <...> generic args, then keep the final path segment
        # (the function/method name) so distinct functions stay distinct.
        if nm.startswith("0x"):
            return nm
        out, depth = [], 0
        for ch in nm:
            if ch == "<":
                depth += 1
            elif ch == ">":
                depth = max(0, depth - 1)
            elif depth == 0:
                out.append(ch)
        segs = [s for s in "".join(out).split("::") if s and s != " as "]
        return segs[-1] if segs else nm

    def dump(title, counter):
        print(f"\n=== top {top_n} by {title} ===")
        for nm, c in counter.most_common(top_n):
            print(f"{100 * c / total:6.2f}%  {c:8.0f}  {short(nm)}")

    print(f"{os.path.basename(path)}: {total:.0f} weighted samples"
          f"{'  (symbolicated)' if sc else '  (RAW addrs — no sidecar)'}")
    dump("SELF time", self_c)
    dump("INCLUSIVE time", incl_c)


if __name__ == "__main__":
    main()
