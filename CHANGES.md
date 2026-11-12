# Changes vs. SydShafferLab/BarcodeAnalysis_v2_10X

Upstream: `step1.sh` + `scripts/ExtractBarcodes_10X.py` at commit `88bed94` (2024-11-27).

## Bugs fixed

1. **Downstream-anchor-only reads gave the wrong barcode.** Upstream takes the `bclen` nt immediately
   before `GTTTTAGAG`. For a tagged construct, that is the last 12 nt of the barcode plus the 8-nt tag.
   The barcode is now the `bclen` nt before `[tag][downstream anchor]`.
2. **Reads with both anchors produced two records**, one from each anchor (and with bug 1, one of them was
   wrong). Now each read produces exactly one record. When both anchors are present, the insert between
   them must be exactly `bclen + tag length`, otherwise the read is dropped (`insert_len_<n>`).
3. **Cells from different samples were merged.** Upstream keys counts by (CB, barcode), but 10x cell
   barcodes are reused across samples, so pooling samples merged unrelated cells. Counts are now keyed by
   (sample, CB, barcode).
4. **Reverse-complement reads were ignored.** Both orientations are now checked. This rarely matters for
   unmapped reads, but it does with `--all-reads`.

## Features added

- **Tag check.** The tag next to the barcode must match one of `tag_seqs` (Hamming ≤ `tag_max_mm`).
  Reads matching `exclude_tag_seqs` (other libraries' tags, e.g. contamination) or no known tag are dropped.
  Reads whose tag runs off the end of the read are kept if the visible part matches (`keep_partial_tag`).
  An empty `tag_seqs` gives the original untagged behaviour.
- **Per-group starcode.** Samples that should share a barcode namespace (e.g. all timepoints of one
  population) are pooled and clustered together, as in the gDNA pipeline. Different groups are clustered
  independently.
- **Step 1 pre-filters on the anchors**, so the step-1 output contains only candidate barcode reads instead
  of every unmapped read, which makes it orders of magnitude smaller. It takes the BAM path directly and runs
  locally (upstream was written for PMACS). `--all-reads` is available as a QC check.
- **`extraction_stats.tsv`** records why each read was kept or dropped, per sample.
- One JSON config with defaults, paths relative to the config, and validation of samples, groups and
  tag lengths.
- Unit tests for every extraction branch, plus an end-to-end test through starcode.
