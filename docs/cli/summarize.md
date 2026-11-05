# poremod summarize

Detection QC without ground truth: for every modification code seen in
the BAM, emit four histograms (the fourth only with a reference FASTA)
of call quality and methylation frequency. The main artifacts are TSV
tables; PNG figures are optional and are always drawn from `bins.tsv`.

## Usage

```bash
# Scan a BAM (writes DIR/bins.tsv and DIR/summary.tsv)
poremod summarize --input in.bam --output DIR [--reference REF] [--figure] [--report]
                 [--sample-reads N] [--seed S]
                 [--region] [--min-probability] [--filter-percentile]
                 [--edge-filter S E] [--skip-mode implicit|explicit]

# Replot existing bins.tsv (does not scan a BAM)
poremod summarize --bins DIR/bins.tsv --output DIR [--figure] [--report]
```

`--input` and `--bins` are mutually exclusive; exactly one is required.

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | one of `--input` / `--bins` | Input modbam. Unmapped reads with MM/ML tags enter the three per-read distributions. |
| `--bins PATH` | one of `--input` / `--bins` | Existing `bins.tsv` to replot. Does **not** open or scan a BAM. Requires `--figure` or `--report`. Rejects `--reference`, `--sample-reads`, `--seed`, and every filter flag. Does not rewrite `bins.tsv` / `summary.tsv`. |
| `--output DIR` | yes | Output directory (created if missing). Scan mode writes `bins.tsv` and `summary.tsv` here; `--figure` writes PNGs here. |
| `--report` | no | Write a self-contained HTML report to `DIR/report.html`: run information, one section per modification code holding its distribution figures with bilingual descriptions, and the `summary.tsv` table with an info box explaining its `metric` column. Implies `--figure`, so the figures exist to be inlined (base64) -- the page has no external references. Also legal with `--bins`. |
| `--figure` | no (`--bins`: yes, or `--report`) | Write one PNG per `(mod_code, dist)` with `n > 0` into `--output`. Filenames are `{mod_code}.{dist}.png` (ChEBI codes use the digits, e.g. `21839.mod_qual.png`). Even when combined with `--input`, figures are drawn by reading the TSV just written — never from the in-memory histogram. |
| `--reference PATH` | no (`--input` only) | Reference FASTA. Enables the `site_freq` distribution via a full pileup (`percent_modified`). Without it there are **no** `site_freq` rows in `bins.tsv` / `summary.tsv` and **no** `*.site_freq.png` files (not zero-filled stubs). A failed pileup exits non-zero and writes no TSV. An input that is not aligned (all-unmapped) or not coordinate-sorted is automatically re-aligned against this reference before the pileup (requires `mappy` and `samtools` on PATH); already-aligned input is used as-is and never re-aligned. |
| `--sample-reads N` | no (`--input` only) | Reservoir-sample `N` eligible reads for the three **per-read** distributions (`mod_qual` / `call_qual` / `read_freq`). **Default: `10000`.** Omitting the flag still samples — this is the opposite of `extract`, where omitting `--sample-reads` means all reads. `0` disables sampling (every eligible read enters the histograms). Negative `N` is a usage error (exit 2). `site_freq` is never sampled. |
| `--seed S` | no (`--input` only) | RNG seed for `--sample-reads` (default `0`). Legal on its own because summarize samples by default. Ignored when `N = 0`. The same `N`, seed, and input give the same summarize sample; that set is **not** promised to match `extract --sample-reads`. |
| `--region STR` | no (`--input` only) | Region filter: contig name (`chr1`) or `chr1:1000-2000` (0-based, half-open). Unmapped reads have no coordinates, so they drop out of a region-filtered run. |
| `--min-probability FLOAT` | no (`--input` only) | Global fail threshold on `call_prob` (`[0, 1]`). Calls below it are `fail=True` and do not enter the per-read histograms. Unset: no probability fail marking. |
| `--filter-percentile FLOAT` | no (`--input` only) | Accepted for Filter compatibility with extract (default: off). summarize does not use it to estimate a threshold; fail marking happens only with `--min-probability`. |
| `--edge-filter S E` | no (`--input` only) | Edge filter (bases to ignore at the start / end of each read). |
| `--skip-mode {implicit,explicit}` | no (`--input` only) | Canonical inference mode. Default parsing is `implicit`. This is not an on/off filter switch: default “no filter” only skips thresholds / edge / region. |
| `--keep-realigned` | no (`--input` only) | With `--reference`, when the input is re-aligned (see `--reference`), keep the result as `--output/realigned.bam` (plus `.bai`) instead of a temporary file deleted after the run. Has no effect on already-aligned input. Rejected with `--bins`. |

By default no `Filter` is built (no fail marking, no edge trim, no region drop). Giving `--min-probability`, `--edge-filter`, `--region`, or `--skip-mode explicit` applies the same filtering as `extract`.

## Sampling vs extract

| | `summarize` | `extract` |
|--|-------------|-----------|
| Omit `--sample-reads` | sample **10000** eligible reads | **all** reads (no sampling) |
| `--sample-reads 0` | disable sampling (all eligible) | (extract has no `0` = all; unset means all) |
| `--sample-reads N` (`N ≥ 1`) | reservoir of `N` eligible reads | reservoir of `N` output reads |
| `--seed` without `--sample-reads` | legal (default sampling is on) | usage error (exit 2) |

Eligible for summarize: a parsed read that still has at least one non-fail call after optional filters, **including unmapped**. Secondary alignments (FLAG 0x100) are always dropped before parsing, matching `extract` / `pileup`; supplementary alignments (0x800) are kept. Reads with no MM/ML tags are skipped. Sampling does not change `site_freq` (full mapped pileup when `--reference` is set). `records_seen` still counts every BAM record (including secondary) as scan progress.

## Output

`--output DIR` receives:

```
DIR/bins.tsv
DIR/summary.tsv
DIR/summarize.log           # scan mode (--input) only
DIR/{mod_code}.{dist}.png   # only with --figure, and only when that dist has n > 0
DIR/report.html             # only with --report
```

### `bins.tsv`

```
mod_code	dist	bin_lo	bin_hi	count
```

- `dist` ∈ {`mod_qual`, `call_qual`, `read_freq`, `site_freq`}
- Dense table: every present `(mod_code, dist)` gets its full bin count of rows (`count=0` is written) — 256 rows for each quality dist (`mod_qual`, `call_qual`), 100 for the frequency dists
- Both quality dists have one bin per ML tag level `0..255`: bin `m` is `[m/256, (m+1)/256)`, so the level's probability `(m+0.5)/256` falls in its own bin. A `call_qual` observation is the winning call's probability, drawn from the same ML level mapping, so it lands one per bin too
- The frequency dists keep 100 equal-width bins on `[0, 1]`; bins `0..98` are half-open `[i/100, (i+1)/100)`; bin `99` is `[0.99, 1.0]` (includes `1.0`)
- `bin_lo` / `bin_hi` are exact decimals with trailing zeros stripped (`0`, `0.01`, `0.00390625`, `1`)
- Row order: `mod_code` in rust_modbam `ModCode` order (ascending ChEBI, then single-character ASCII); within a code, dist order as listed above; within a dist, bin index ascending
- No reference → the `site_freq` block is omitted entirely (no empty rows)
- `--figure` / `--report` reject a quality-dist block that is not 256 rows: a `bins.tsv` written by an earlier version (100-row `mod_qual` / `call_qual`) must be regenerated by re-running the scan

### `summary.tsv`

```
mod_code	mod_base	metric	n	mean	min	max	q25	median	q75
```

- `mod_base` is the display name of the modification code (`m` → `5mC`, `21839` → `4mC`, `a` → `6mA`, ...; same mapping as the figure titles). Codes outside the table ("any" codes `C`/`A`/`G`/`T`, unlisted ChEBI ids) keep the raw code
- `metric` is the distribution name (`mod_qual`, `call_qual`, `read_freq`, `site_freq`) — the same values `bins.tsv` stores under its `dist` column
- A metric with `n == 0` still gets one row; `mean` / `min` / `max` / `q*` are empty fields; no PNG is drawn
- `n` / `mean` / `min` / `max` are exact streaming stats; quartiles are linearly interpolated from the same dense histogram the figure uses
- `--bins` replot does **not** rewrite these two TSV files

### stdout (`--input`)

Two sampling-account lines, then the contents of `summary.tsv`:

```
reads_eligible: <parsed reads that still have a call>
reads_sampled: <reads that entered the three per-read distributions>
mod_code	mod_base	metric	n	mean	min	max	q25	median	q75
...
```

When `N = 0` or eligible ≤ `N`, the two counts are equal. With no modification codes the header is still written and stdout adds `mod_codes: none`; `reads_eligible` / `reads_sampled` still print (they may be 0) and the exit code is 0. Those two account lines are **not** stored in the TSV files (`--bins` does not need them).

`--bins` reprints a sibling `summary.tsv` when one is found; it does not reprint the account lines and it does not scan a BAM.

### `summarize.log`

Scan mode (`--input`) writes a run log to `DIR/summarize.log` and echoes every
line to stderr as it is written (each line is flushed, so `tail -f
DIR/summarize.log` and the terminal both show the run live). `--bins` mode
writes no log.

One line per event, `YYYY-MM-DD HH:MM:SS [phase] message` (local time), with
`phase` one of `run` / `align` / `pileup` / `scan` / `tables` / `figures` /
`report` / `done` / `error`. During the scan the log reports progress about
every 10 seconds (`[scan] records=... eligible=... elapsed=...s`); the
`[align]` / `[pileup]` lines appear only with `--reference`. The file is
overwritten by the next run into the same `--output` -- the log describes one
run. A failing run writes an `[error]` line before the error surfaces and the
lines already written stay.

## Distributions

| `dist` | Each observation | Who enters |
|--------|------------------|------------|
| `mod_qual` | Probability of code `C` in that call's `mod_probs` | Every call that declares `C` (regardless of argmax) |
| `call_qual` | That call's `call_prob` | Only `call_code == C` |
| `read_freq` | Per-read methylation frequency of `C` | A read with at least one call on `C`'s canonical base |
| `site_freq` | Pileup `percent_modified` (`n_mod / n_valid_cov`, a `[0, 1]` fraction) | Only with `--reference`; rows whose `mod_code == C` |

The three per-read distributions include **unmapped** reads and do **not** apply extract's reference-base match. `fail == true` calls are excluded. Quality axes are probabilities `[0, 1]`, not ML bytes `[0, 255]` (both quality dists bin one ML level per bin, so the figures show the tag's discrete structure).

Per-read frequency for code `C` (canonical base `B`): denominator = calls on `B` on that read; numerator = those with `call_code == C`; skip the read when the denominator is 0.

## Figures (`--figure`)

- One PNG per `(mod_code, dist)`, not a 2×2 panel
- Filename `{mod_code}.{dist}.png`; figure size 6×4.5 inches, saved at dpi 300 (1800×1350 px) so the figures stay sharp when the HTML report opens one full screen
- matplotlib Agg; x axis `[0, 1]`, except `call_qual` panels, which start at their first occupied bin's lower edge (a winner's call probability exceeds 0.5 for a single declared code, so the left half would stay empty; the axis follows the data rather than clipping it); y axis density (`count_i / (n * bin_width)`, `bin_width = bin_hi - bin_lo`)
- All plot text is English: panel title `{mod_name} {dist_name}`, no suptitle. `mod_name` is the chemical shorthand of the modification code (`m` → `5mC`, `a` → `6mA`, ChEBI `21839` → `4mC`, `h` → `5hmC`, `o` → `8oxoG`, ...), covering the 17 named codes of rust_modbam; "any" codes (`C`/`A`/`G`/`T`) and unlisted ChEBI ids keep the raw code. `dist_name` is the formal name of the distribution (`mod_qual` → `Modification probability`, `call_qual` → `Call probability`, `read_freq` → `Methylation frequency (per read)`, `site_freq` → `Methylation frequency (per site)`); the TSV files keep the short `dist` values. Axis labels follow the same convention: x axis `Probability` for `mod_qual`/`call_qual`, `Methylation frequency (per read)` for `read_freq`, `Methylation frequency (per site)` for `site_freq`; y axis `Density` on every panel.
- `n == 0` or a missing dist (no FASTA → no `site_freq`): no file, no placeholder PNG
- The plotter only reads `bins.tsv` columns. Missing columns, an illegal `dist`, or a group whose row count is not the one its dist expects (256 for `mod_qual` / `call_qual`, 100 for the frequency dists) is an error; no partial figure set is written for a corrupt TSV

## Input validation

Before any analysis work, the input modbam is gated through the same
structural check as `poremod validate`: the first 10000 records are scanned
for MM/ML structure problems. A passing input — including one whose only
findings are warnings — proceeds silently. A failing input prints the
standard validate report to stderr and exits 1, before any output file is
created.

## Exit codes

| Exit code | Scenario |
|-----------|----------|
| 0 | Tables (and optional figures) written; empty `mod_codes` is success |
| 2 | BAM / reference / `bins.tsv` missing or unreadable |
| 2 | Usage: neither or both of `--input` / `--bins`; `--bins` without `--figure`; `--bins` with scan flags; `--sample-reads` negative |
| 2 | MM/ML syntax error hit while scanning (bad records are not skipped; only “tag not found” is skipped; a bad record inside the scan window (the first 10000 records) is already caught by the input-validation gate (exit code 1); only records outside the window can reach this row) |
| 2 | `--reference` pileup failed (no TSV is left behind) |

Errors print `error: ...` on stderr and do not leave a “looks like an empty match” TSV from a failed BAM scan.

## Examples

```bash
# Default: sample 10000 eligible reads (or all of them if fewer), TSV only
poremod summarize --input in.bam --output qc/

# Keep every eligible read (including unmapped)
poremod summarize --input in.bam --output qc/ --sample-reads 0

# Enable site_freq and write PNGs
poremod summarize --input in.bam --output qc/ --reference ref.fa --figure

# Smaller sample, reproducible
poremod summarize --input in.bam --output qc/ --sample-reads 2000 --seed 1

# Same Filter knobs as extract (fail marking only with --min-probability)
poremod summarize --input in.bam --output qc/ --region chr1:0-2000 --min-probability 0.7

# Replot from TSV — the BAM is not opened
poremod summarize --bins qc/bins.tsv --output qc-replot/ --figure
```

## Notes

- Values come from the core accumulators; the CLI does not re-bin.
- `--figure` with `--input` always writes TSV first, then calls the same `plot_summarize(bins_path, output_dir)` as `--bins`.
- `read_freq` (sampled per-read) and `site_freq` (full mapped pileup) can describe different corpora; stdout prints both `reads_eligible` and `reads_sampled` so the two `n` columns are not mistaken for the same read set.
- `poremod extract` continues to skip unmapped reads. summarize's “include unmapped” contract exists only on this command.
- When `--reference` re-aligns an unaligned input, both the per-read distributions and `site_freq` are computed from the re-aligned BAM (reads that still fail to map drop out of both), rather than from the original unmapped input. The re-alignment uses `align`'s defaults (`preset=map-ont`, 4 mapping threads) and carries each read's MM/ML tags over verbatim; `--keep-realigned` preserves the intermediate BAM.
