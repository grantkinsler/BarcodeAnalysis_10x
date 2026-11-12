"""Step 2: extract lineage barcodes from step-1 read tables, cluster with starcode, assign to cells.

Adapted from SydShafferLab/BarcodeAnalysis_v2_10X scripts/ExtractBarcodes_10X.py for tagged lineage-barcode
constructs (barcode followed by a short tag that identifies the library). Expected read layout:

    [upstream anchor] [bclen-nt barcode] [tag] [downstream anchor]
    e.g. AAAGGACGAAACACCG  NNNNNNNNNNNNNNNNNNNN  GCCGAATC  GTTTTAGAGCTAGAAA

With an empty "tag_seqs" the construct is treated as untagged: [upstream] [barcode] [downstream].

Changes vs the upstream script (see CHANGES.md):
  1. Downstream-anchor-only reads: upstream takes the bclen nt preceding the downstream anchor, which for a
     tagged construct is part barcode + the tag. Here the barcode is the bclen nt preceding [tag][downstream].
  2. One record per read (upstream emits two records, one of them wrong, when both anchors are present).
  3. Tag check: the tag next to the barcode must be within Hamming <= tag_max_mm of a tag in tag_seqs; reads
     matching exclude_tag_seqs (other libraries' tags) or no tag are dropped; reads whose tag runs off the end
     of the read are kept only if the visible prefix matches (keep_partial_tag). When both anchors are present
     the insert must be exactly bclen + tag length.
  4. Starcode is run per GROUP (samples that should share one barcode namespace are pooled), and per-cell
     counts are keyed by (sample, CB): cell barcodes are not unique across 10x samples, so the upstream
     (CB, barcode) key merges different cells from different samples.
  5. Reverse-complemented reads are also checked (rare for unmapped reads, needed with step 1 --all-reads).

usage: python extract_barcodes_10x.py <config.json> [--group NAME]
Relative paths in the config are resolved against the config file's directory.

outputs, per group, in <out_dir>/<group>/:
  step1_output.txt         sample, CB, UMI, raw barcode, anchor (both/up/down), tag status
  sc_input.txt             raw barcodes, one per line (starcode input)
  sc_output.txt            starcode clusters (--seq-id)
  EB_output_modified.txt   step1_output with barcode replaced by its starcode centroid
  combo_counts.txt         sample, CB, barcode, n_UMI, n_reads
  extraction_stats.tsv     per-sample counts of every keep/drop reason
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")
DEFAULTS = {"starcode": "starcode", "threads": 4, "upstream_anchor": "AAAGGACGAAACACCG",
            "downstream_anchor": "GTTTTAGAGCTAGAAA", "bclen": 20, "sc_mm": 4, "tag_seqs": {},
            "exclude_tag_seqs": {}, "tag_max_mm": 1, "keep_partial_tag": True}


def rc(s):
    return s.translate(COMP)[::-1]


def ham(a, b):
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def classify_tag(obs, tags, excl, max_mm, tag_len):
    """-> 'match' | 'partial' | 'excluded' | 'mismatch' | 'absent' | 'untagged'"""
    if tag_len == 0:
        return "untagged"
    if len(obs) == 0:
        return "absent"
    if len(obs) < tag_len:  # tag runs off the end of the read: compare the visible prefix only
        mm = max_mm if len(obs) >= tag_len // 2 else 0
        return "partial" if any(ham(obs, t[:len(obs)]) <= mm for t in tags) else "mismatch"
    if any(ham(obs, t) <= max_mm for t in tags):
        return "match"
    if any(ham(obs, t) <= max_mm for t in excl):
        return "excluded"
    return "mismatch"


def extract(seq, cfg, tags, excl, tag_len):
    """-> (barcode or None, anchor, tag status / drop reason)"""
    U, D, L, T = cfg["upstream_anchor"], cfg["downstream_anchor"], cfg["bclen"], tag_len
    for s in (seq, rc(seq)):
        u, d = s.find(U), s.find(D)
        if u < 0 and d < 0:
            continue
        if u >= 0 and d > u:
            insert = s[u + len(U):d]
            if len(insert) != L + T:
                return None, "both", f"insert_len_{len(insert)}"
            bc, tag, anchor = insert[:L], insert[L:], "both"
        elif u >= 0:
            if u + len(U) + L > len(s):
                return None, "up", "bc_truncated"
            bc, tag, anchor = s[u + len(U):u + len(U) + L], s[u + len(U) + L:u + len(U) + L + T], "up"
        else:
            if d < L + T:
                return None, "down", "bc_truncated"
            bc, tag, anchor = s[d - T - L:d - T], s[d - T:d], "down"
        if "N" in bc:
            return None, anchor, "N_in_bc"
        status = classify_tag(tag, tags, excl, cfg["tag_max_mm"], T)
        keep = status in ("match", "untagged") or (status == "partial" and cfg["keep_partial_tag"])
        return (bc if keep else None), anchor, status
    return None, "none", "no_anchor"


def load_config(path):
    cfg = dict(DEFAULTS, **json.load(open(path)))
    base = os.path.dirname(os.path.abspath(path))
    resolve = lambda p: p if os.path.isabs(p) else os.path.normpath(os.path.join(base, p))
    cfg["out_dir"] = resolve(cfg["out_dir"])
    cfg["samples"] = {s: resolve(p) for s, p in cfg["samples"].items()}
    if os.sep in cfg["starcode"]:
        cfg["starcode"] = resolve(cfg["starcode"])
    lens = {len(t) for t in list(cfg["tag_seqs"].values()) + list(cfg["exclude_tag_seqs"].values())}
    if len(lens) > 1:
        sys.exit(f"all tags must be the same length, got lengths {sorted(lens)}")
    if cfg["exclude_tag_seqs"] and not cfg["tag_seqs"]:
        sys.exit("exclude_tag_seqs given but tag_seqs is empty")
    cfg["tag_len"] = lens.pop() if lens else 0
    unknown = {s for g in cfg["groups"].values() for s in g} - set(cfg["samples"])
    if unknown:
        sys.exit(f"groups reference samples not listed in 'samples': {sorted(unknown)}")
    return cfg


def run_group(name, samples, cfg):
    out = os.path.join(cfg["out_dir"], name)
    os.makedirs(out, exist_ok=True)
    tags, excl, T = list(cfg["tag_seqs"].values()), list(cfg["exclude_tag_seqs"].values()), cfg["tag_len"]
    stats = defaultdict(Counter)
    rows = []
    print(f"[{name}] extracting from {len(samples)} samples", flush=True)
    for smp in samples:
        with open(cfg["samples"][smp]) as f:
            for line in f:
                cb, umi, seq = line.rstrip("\n").split("\t")
                bc, anchor, status = extract(seq, cfg, tags, excl, T)
                stats[smp][f"{anchor}:{status}"] += 1
                stats[smp]["kept" if bc else "dropped"] += 1
                if bc:
                    rows.append((smp, cb, umi, bc, anchor, status))
    with open(f"{out}/step1_output.txt", "w") as f:
        f.writelines("\t".join(r) + "\n" for r in rows)
    with open(f"{out}/sc_input.txt", "w") as f:
        f.writelines(r[3] + "\n" for r in rows)
    keys = sorted({k for c in stats.values() for k in c})
    with open(f"{out}/extraction_stats.tsv", "w") as f:
        f.write("sample\t" + "\t".join(keys) + "\n")
        for smp in samples:
            f.write(smp + "\t" + "\t".join(str(stats[smp][k]) for k in keys) + "\n")
    for smp in samples:
        n = stats[smp]["kept"] + stats[smp]["dropped"]
        print(f"  {smp}: kept {stats[smp]['kept']} / {n} anchor reads", flush=True)

    if not rows:
        print(f"[{name}] no barcodes kept - skipping starcode", flush=True)
        return
    print(f"[{name}] starcode -d {cfg['sc_mm']} on {len(rows)} barcodes", flush=True)
    subprocess.run([cfg["starcode"], "-d", str(cfg["sc_mm"]), "-t", str(cfg["threads"]),
                    "-i", f"{out}/sc_input.txt", "-o", f"{out}/sc_output.txt", "--seq-id"], check=True)

    centroid = {}
    with open(f"{out}/sc_output.txt") as f:
        for line in f:
            bc, _, ids = line.rstrip("\n").split("\t")
            for i in map(int, ids.split(",")):
                centroid[i] = bc  # starcode seq-ids are 1-based line numbers of sc_input
    umis, reads = defaultdict(set), Counter()
    with open(f"{out}/EB_output_modified.txt", "w") as f:
        for i, (smp, cb, umi, bc, anchor, status) in enumerate(rows, 1):
            bc = centroid[i]
            f.write("\t".join((smp, cb, umi, bc, anchor, status)) + "\n")
            umis[(smp, cb, bc)].add(umi)
            reads[(smp, cb, bc)] += 1
    with open(f"{out}/combo_counts.txt", "w") as f:
        f.write("sample\tCB\tbarcode\tn_UMI\tn_reads\n")
        for (smp, cb, bc), u in sorted(umis.items()):
            f.write(f"{smp}\t{cb}\t{bc}\t{len(u)}\t{reads[(smp, cb, bc)]}\n")
    print(f"[{name}] done: {len(set(centroid.values()))} clusters, {len(umis)} (sample,CB,barcode) combos -> {out}",
          flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config")
    ap.add_argument("--group", help="run only this group")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    groups = cfg["groups"] if not a.group else {a.group: cfg["groups"][a.group]}
    missing = [f"{s} ({cfg['samples'][s]})" for g in groups.values() for s in g
               if not os.path.exists(cfg["samples"][s])]
    if missing:
        sys.exit("missing step-1 output (run scripts/step1_extract_reads.sh first):\n  " + "\n  ".join(missing))
    if shutil.which(cfg["starcode"]) is None and not os.access(cfg["starcode"], os.X_OK):
        sys.exit(f"starcode not found: {cfg['starcode']} (install via bioconda, see README)")
    for name, samples in groups.items():
        run_group(name, samples, cfg)


if __name__ == "__main__":
    main()
