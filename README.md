# 10x lineage-barcode extraction (tag-aware)

Assigns lineage barcodes to single cells from a 10x Genomics / Cell Ranger run, for lineage-barcode
constructs where the barcode is followed by a short **tag** that identifies the barcode library:

```
AAAGGACGAAACACCG  [20 nt barcode]  GCCGAATC  GTTTTAGAGCTAGAAA
 upstream anchor                     tag     downstream anchor
```

Adapted from [`SydShafferLab/BarcodeAnalysis_v2_10X`](https://github.com/SydShafferLab/BarcodeAnalysis_v2_10X),
with bugs fixed that matter for tagged constructs and for multi-sample runs. See [`CHANGES.md`](CHANGES.md).
Companion to the gDNA pipeline [`grantkinsler/BarcodeAnalysis_gDNA`](https://github.com/grantkinsler/BarcodeAnalysis_gDNA)
and uses the same tag conventions, so 10x and gDNA lineages can be matched by sequence.
Untagged constructs work too: leave `tag_seqs` empty.

## Quick start

```bash
conda env create -f environment.yml        # python, samtools, starcode, pytest
conda activate barcode_analysis_10x

# Step 1: one call per 10x sample. Pull anchor-bearing reads (with cell barcode + UMI) out of the BAM.
scripts/step1_extract_reads.sh \
    /path/to/cellranger/<sample>/outs/per_sample_outs/<sample>/count/sample_alignments.bam \
    runs/step1/<sample>_unmapped.txt

# Step 2: extract + tag-check barcodes, cluster with starcode per group, count UMIs per cell.
cp configs/example_config.json configs/my_run.json      # edit samples / groups / tags
python scripts/extract_barcodes_10x.py configs/my_run.json
```

## Step 1: `scripts/step1_extract_reads.sh`

`step1_extract_reads.sh [--all-reads] [-u UPSTREAM] [-d DOWNSTREAM] [-t THREADS] <bam> <out.txt>`

Streams the BAM through `samtools view` and keeps reads that have both a `CB` and a `UB` tag and contain
either anchor (either orientation). Output is `CB<TAB>UMI<TAB>SEQ`, one line per read.

- By default only **unmapped** reads (`-f 4`) are scanned, because barcode reads don't align to the genome.
  `--all-reads` scans everything. Run it once on a sample as a QC check that `-f 4` isn't missing
  barcode reads that got partially aligned (for example to a transgene in a custom reference).
- Works on `sample_alignments.bam` (cellranger multi) or `possorted_genome_bam.bam` (cellranger count).

## Step 2: `scripts/extract_barcodes_10x.py`

`python scripts/extract_barcodes_10x.py <config.json> [--group NAME]`

For each read: find the anchors, cut out the barcode and the tag next to it, and keep the read only if the
tag matches one of `tag_seqs`. Kept barcodes are then clustered with starcode **per group**, and UMIs are
counted per (sample, cell barcode, consensus barcode).

### Config fields

Relative paths are resolved against the config file's directory. Omitted fields use these defaults.

| field | default | meaning |
|---|---|---|
| `samples` | (required) | `{sample_name: step-1 output file}` |
| `groups` | (required) | `{group_name: [sample_name, ...]}`. Each group is clustered with starcode together, so it shares one barcode namespace. Pool samples whose lineages should be matched (e.g. all timepoints of one population). |
| `out_dir` | (required) | output root; one subfolder per group |
| `starcode` | `starcode` | starcode executable (name on PATH or path) |
| `threads` | 4 | starcode threads |
| `upstream_anchor` / `downstream_anchor` | `AAAGGACGAAACACCG` / `GTTTTAGAGCTAGAAA` | sequences flanking barcode+tag |
| `bclen` | 20 | barcode length |
| `sc_mm` | 4 | starcode `-d` (max Levenshtein distance for clustering) |
| `tag_seqs` | `{}` | tags to **keep** (`{name: seq}`). Empty means an untagged construct. |
| `exclude_tag_seqs` | `{}` | tags of other libraries; reads matching these are dropped and counted as `excluded` |
| `tag_max_mm` | 1 | max Hamming mismatches when matching a tag |
| `keep_partial_tag` | true | keep reads whose tag runs off the end of the read if the visible part matches |

All tags must be the same length.

### Outputs (`<out_dir>/<group>/`)

| file | contents |
|---|---|
| `combo_counts.txt` | **main output**: `sample, CB, barcode, n_UMI, n_reads` |
| `EB_output_modified.txt` | per read: `sample, CB, UMI, consensus barcode, anchor, tag status` |
| `step1_output.txt` | same as above, with the raw (pre-clustering) barcode |
| `sc_input.txt`, `sc_output.txt` | starcode input / output |
| `extraction_stats.tsv` | per sample, the number of reads for every `anchor:status` outcome (`both:match`, `down:excluded`, `up:bc_truncated`, ...) plus `kept` / `dropped` |

Note that cell barcodes (`CB`) are only unique **within** a 10x sample, so always key cells by
(sample, CB).

## Tests

```bash
pytest -q          # the end-to-end test is skipped if starcode isn't on PATH
```

## License

MIT, see [`LICENSE`](LICENSE).
