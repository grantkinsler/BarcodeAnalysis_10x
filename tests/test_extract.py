"""Tests for scripts/extract_barcodes_10x.py. Run from the repo root:  pytest -q"""
import json
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import extract_barcodes_10x as eb  # noqa: E402

U, D = "AAAGGACGAAACACCG", "GTTTTAGAGCTAGAAA"
BC = "ACGTACGTTGCATGCAACGT"          # 20 nt
TAG1, TAG2 = "GCCGAATC", "GATAAGCC"
CFG = dict(eb.DEFAULTS, tag_seqs={"Tag1": TAG1}, exclude_tag_seqs={"Tag2": TAG2})
TAGS, EXCL = [TAG1], [TAG2]


def ex(seq, cfg=CFG, tags=TAGS, excl=EXCL, tag_len=8):
    return eb.extract(seq, cfg, tags, excl, tag_len)


def test_both_anchors():
    assert ex("TTT" + U + BC + TAG1 + D + "AAA") == (BC, "both", "match")


def test_upstream_only():
    assert ex("TT" + U + BC + TAG1 + "GTTT") == (BC, "up", "match")


def test_downstream_only_takes_barcode_before_tag():
    # the upstream script took the 20 nt right before the downstream anchor (= 12 nt of barcode + tag)
    assert ex("CCC" + BC + TAG1 + D + "AA") == (BC, "down", "match")


def test_reverse_complement():
    assert ex(eb.rc("TTT" + U + BC + TAG1 + D + "AAA")) == (BC, "both", "match")


def test_one_mismatch_in_tag_kept():
    assert ex(U + BC + "GCCGAATG" + D)[0] == BC


def test_excluded_tag_dropped():
    assert ex(U + BC + TAG2 + D) == (None, "both", "excluded")


def test_unknown_tag_dropped():
    assert ex(U + BC + "TTTTTTTT" + D) == (None, "both", "mismatch")


def test_wrong_insert_length_dropped():
    assert ex(U + BC + "A" + TAG1 + D) == (None, "both", "insert_len_29")


def test_partial_tag_at_read_end():
    assert ex(U + BC + TAG1[:5]) == (BC, "up", "partial")
    assert ex(U + BC + TAG1[:5], cfg=dict(CFG, keep_partial_tag=False))[0] is None


def test_truncated_barcode_dropped():
    assert ex(U + BC[:12]) == (None, "up", "bc_truncated")


def test_N_in_barcode_dropped():
    assert ex(U + "N" + BC[1:] + TAG1 + D) == (None, "both", "N_in_bc")


def test_untagged_construct():
    cfg = dict(eb.DEFAULTS)
    assert ex(U + BC + D, cfg=cfg, tags=[], excl=[], tag_len=0) == (BC, "both", "untagged")
    assert ex("GG" + BC + D, cfg=cfg, tags=[], excl=[], tag_len=0) == (BC, "down", "untagged")


def test_no_anchor():
    assert ex("ACGT" * 20) == (None, "none", "no_anchor")


@pytest.mark.skipif(shutil.which("starcode") is None, reason="starcode not on PATH")
def test_end_to_end(tmp_path):
    bc2 = "TTGCAACGTACGTACGGATC"
    err = BC[:-1] + "A"  # 1 mismatch -> should collapse onto BC
    reads = {
        "s1": [("CELL1", "UMI1", U + BC + TAG1 + D), ("CELL1", "UMI2", U + err + TAG1 + D),
               ("CELL1", "UMI2", "CC" + BC + TAG1 + D), ("CELL2", "UMI3", U + bc2 + TAG1 + D),
               ("CELL2", "UMI4", U + bc2 + TAG2 + D)],                     # excluded tag
        "s2": [("CELL1", "UMI9", U + BC + TAG1 + D)]                        # same CB, different sample
              + [("CELL3", f"UMI{i}", U + BC + TAG1 + D) for i in range(10, 14)],
    }
    # starcode (message passing) only merges `err` into BC if BC has >= 5x its reads: 7 vs 1 here
    samples = {}
    for s, rs in reads.items():
        p = tmp_path / f"{s}.txt"
        p.write_text("".join("\t".join(r) + "\n" for r in rs))
        samples[s] = p.name
    cfg = dict(out_dir="out", tag_seqs={"Tag1": TAG1}, exclude_tag_seqs={"Tag2": TAG2}, threads=1,
               samples=samples, groups={"g": ["s1", "s2"]})
    (tmp_path / "config.json").write_text(json.dumps(cfg))
    eb.main([str(tmp_path / "config.json")])

    lines = (tmp_path / "out" / "g" / "combo_counts.txt").read_text().splitlines()[1:]
    got = {tuple(l.split("\t")[:3]): tuple(map(int, l.split("\t")[3:])) for l in lines}
    assert got == {("s1", "CELL1", BC): (2, 3), ("s1", "CELL2", bc2): (1, 1), ("s2", "CELL1", BC): (1, 1),
                   ("s2", "CELL3", BC): (4, 4)}
