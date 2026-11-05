# poremod evaluate

Per-read accuracy evaluation: realigns every dataset against its own reference, scores the filtered modification calls against a five-column BED ground truth, and prints discrete metrics (Acc, confusion, P/R/F1, specificity) plus ROC/PR AUC. `--output-dir` additionally writes four per-read PNG figures, a metrics TSV and the run log.

## Usage

```bash
# single dataset
poremod evaluate --input in.bam --truth gt.bed --reference ref.fa --canonical-base BASE [--json] [filter args]

# the input is realigned against the reference before scoring
poremod evaluate --input reads.bam --reference ref.fa --truth sites.bed \
    --canonical-base C --output-dir out/

# multiple datasets: the Nth --input pairs with the Nth --truth and --reference
poremod evaluate --input a.bam --reference ref.fa --truth a.bed \
                 --input b.bam --reference ref.fa --truth b.bed \
                 --canonical-base C --json
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input modbam path, aligned or not (see "Realignment"). No `.bai` is needed: the input is read sequentially, and `--region` restricts the realigned copy, which carries its own index. Repeatable; the Nth occurrence pairs with the Nth `--truth` and `--reference`. |
| `--truth PATH` | yes | BED ground truth: at least five TAB columns (`chrom`, 0-based `start`, `end`, `canonical_base`, `expected_state`); columns 6+ are ignored. Repeatable. |
| `--reference PATH` | yes | Reference FASTA path. Repeatable; the Nth occurrence pairs with the Nth `--input`. Every dataset is realigned against its reference before scoring, so the scored coordinates always come from this FASTA (an already-aligned input is realigned too). |
| `--json` | no | Emit the full report dict as JSON (default: human-readable summary) |
| `--region STR` | no | Region filter: contig name (`chr1`) or `chr1:1000-2000` (0-based, half-open). It restricts the scored sites; the realignment itself still reads the whole input. |
| `--min-probability FLOAT` | no | Global fail threshold on `call_prob` (argmax probability, `[0, 1]`). Calls below this value are marked `fail=True` and are not scored. Unset: no probability fail marking. |
| `--filter-percentile FLOAT` | no | Accepted for `Filter` compatibility with `extract`/`pileup` (default: off). evaluate does not apply the percentile estimate; probability fail marking happens only with `--min-probability`. |
| `--edge-filter S E` | no | Edge filter (bases to ignore at the start / end of each read) |
| `--skip-mode {implicit,explicit}` | no | Canonical inference mode (default implicit) |
| `--canonical-base BASE` | yes | Target canonical base (A/C/G/T, case-insensitive). GT sites whose BED column 4 differs are skipped -- not an error, not counted as uncovered, excluded from every metric denominator. |
| `--output-dir DIR` | no | Write the four per-read PNG figures, `per_read.metrics.tsv` and the run log `evaluate.log` into `DIR` (created when missing). Combines with `--json`. |
| `--threads N` | no | Mapper processes per dataset, handed to the same `align` call that does the realignment (default 4, align's own default). Datasets are realigned concurrently, up to `max(1, cpu_count // N)` of them at a time; the pooled report is identical to a serial run's. `N < 1` is a usage error (exit code 2). |

The filter arguments apply the same filtering as `extract` (no `--min-coverage`, no `--mod-code`). Each GT site's canonical base comes from its own BED column 4; the reference FASTA is used for the realignment, not for canonical-base resolution.

`expected_state` uses the same tokens as `call_code`: `-` (canonical), a single-character modification code (e.g. `m`), or a decimal ChEBI code. Tokens such as `5mC` are errors.

`canonical_base` and `expected_state` must agree: when `expected_state` is a modification code with a definite canonical base in the mapping table (e.g. `m`→C, `a`→A), that base must equal column 4, or the row is an error carrying its line number; `-` and codes outside the mapping table (including decimal ChEBI codes) are not subject to this check.

## Realignment

Every dataset is realigned against its own `--reference` before scoring, so the scored coordinates always come from that reference and never from whatever the input already carries:

- an already-aligned input is realigned too -- there is no "already aligned, skip" branch;
- the realigned BAM (with its `.bai`) is written into a temporary `poremod-evaluate-*` directory that is removed on success and failure alike; the directory sits under `--output-dir` when that option is given and in the system temporary directory otherwise; the input is only read, never modified;
- the alignment is the same call `poremod align` uses (preset `map-ont`, secondary records (`0x100`) suppressed, matching align's default `--secondary=no`; supplementary and unmapped records kept, `--threads` mapper processes), so every run costs one full alignment pass per dataset;
- the reference must be writable: like `align`, evaluate generates/refreshes the reference's `.fai` at startup (rebuilt even when one already exists), so a read-only reference directory fails with exit code 2 (`failed to index reference ...`).

A dataset that fails to realign reports its 1-based index and keeps the failure's exception type, so the exit code matches the failure kind (see "Exit codes"):

```
error: dataset 1: input not found: reads.bam
```

## Output

### Text mode (default)

A summary, not the full JSON:

```
n: 2
Acc: 1
specificity: 1
P/R/F1:
  -: 1 / 1 / 1
  m: 1 / 1 / 1
AUC roc: 1
AUC pr: 1
confusion: [[1, 0], [0, 1]]
uncovered: 0
```

Printed fields: `n`, `Acc`, `specificity`, per-class (and macro, when present) P/R/F1, AUC (plus the one-vs-rest blocks when they exist), the confusion matrix, and the **count** of uncovered ground-truth sites.

### JSON mode (`--json`)

The full report dict. `None` becomes JSON `null`. Binary reports fill `roc_auc` / `pr_auc` and leave `macro` / `roc_auc_ovr` / `pr_auc_ovr` as `null`; multiclass is the reverse. The branch is anchored to the pooled ground-truth state set (`gt_states`): a prediction-only class does not change the 2class / multiclass split. The `site_*` fields (`site_accuracy` / `site_roc_auc_ovr`) are non-null only with three or more classes.

```json
{
  "n": 2,
  "n_correct": 2,
  "n_incorrect": 0,
  "n_positive": 1,
  "n_negative": 1,
  "accuracy": 1.0,
  "insufficient": false,
  "specificity": 1.0,
  "site_accuracy": null,
  "site_roc_auc_ovr": null,
  "labels": ["-", "m"],
  "gt_states": ["-", "m"],
  "confusion": [[1, 0], [0, 1]],
  "per_class": {
    "-": {"n": 1, "precision": 1.0, "recall": 1.0, "f1": 1.0},
    "m": {"n": 1, "precision": 1.0, "recall": 1.0, "f1": 1.0}
  },
  "macro": null,
  "roc_auc": 1.0,
  "pr_auc": 1.0,
  "roc_auc_ovr": null,
  "pr_auc_ovr": null,
  "uncovered": []
}
```

`gt_states` is the pooled ground-truth state set (sorted) -- the class set -- and is validated to hold two or three states (see "Validation"). It is the anchor for every class-sensitive structure: the `confusion` rows and columns, the `per_class` keys, the `macro` gate, the one-vs-rest keys and the binary / multiclass branch. `labels` is the GT ∪ prediction union (prediction-only classes included) and is kept for `channel_dist` and JSON compatibility. The artifact names' `2class` / `3class` comes from the length of `gt_states`. With `--output-dir` the report additionally carries `roc_curve` / `pr_curve` / `channel_dist` -- the payloads the figures are drawn from, including the decision thresholds on the curves.

Behavior notes:

- One scored observation is a filtered per-read modification call that has `chrom` + `ref_position`, `fail` is false, and that coordinate exists in the BED.
- Calls that miss the BED are ignored. Ground-truth sites with no scored hit go to `uncovered` and are not in the accuracy denominator.
- Unmapped reads and reads without MM/ML tags are skipped. The canonical base of a site comes from BED column 4, never from the reference.
- An empty scored set is not an error: metrics may be `null` with `insufficient=true`.
- The case opposite to the previous bullet: a single class in the class set (`k = 1`, with `k` counted over `gt_states`) with `n > 0` -- `insufficient` is `true` and both AUCs are `null`, while `accuracy` stays defined. With an all-canonical result (`labels == ["-"]`) `specificity` is `1.0` rather than `null`: every site has a defined canonical recall, so the site-level macro average is defined -- a different question from whether the AUC can be computed.

## Artifacts (`--output-dir`)

`--output-dir DIR` writes six files into `DIR` (created when missing) and forces the report payloads on (`roc_curve` / `pr_curve` / `channel_dist`), since the figures need them:

| File | Content |
|------|---------|
| `per_read.mod_prob.dist.{2,3}class.png` | one probability-distribution panel per pooled ground-truth state |
| `per_read.ovr_roc.{2,3}class.png` | ROC curve(s): the single binary curve, or one one-vs-rest curve per ground-truth state |
| `per_read.ovr_pr.{2,3}class.png` | PR curve(s) with the same split, plus a prevalence baseline per class |
| `per_read.cofuse_matrix.{2,3}class.png` | confusion-matrix heatmap, rows = truth, columns = predicted |
| `per_read.metrics.tsv` | the metrics table (format unchanged, see "Metrics TSV") |
| `evaluate.log` | this run's log (see "Parallelism and the run log"); overwritten by a rerun into the same directory |

The GT class count in the names comes from the **pooled ground-truth** state count (`gt_states`), never from `labels`: a prediction-only class adds no panel and no curve. The confusion matrix is anchored to `gt_states` as well: its rows and columns are both `gt_states`, and in a corpus where the two sets disagree an out-of-class prediction lands in no cell.

The figures are produced without a display and need no network access or external asset files; all text on them is English, and PNGs are written at 300 dpi.

Distribution semantics: the x axis is that **channel's** predicted probability -- for a modification class, the probability of its code in that observation's `mod_probs`, **treated as 0 when the code is not declared** (modbam only declares modifications the read calls or infers, so a spike at 0 is input semantics, not a defect); for the canonical class it is `1 - sum of modification probabilities`. Every panel stratifies the observations by their ground-truth class, legend title `true_state`.

**The y axis is per-stratum normalized density** (`counts / (N * bin width)`), **not a count**: a stratum with few observations still fills the panel, which is inherent to density normalization.

A confusion cell's two lines are the count plus the **share of the row total** (`count / row_sum`, where `row_sum` is the sum of that row's in-class-set counts; an empty row prints `0.0%`). **BREAKING** (earlier release): this used to be the share of the row maximum.

Insufficient data renders as a placeholder in that panel (`insufficient data`; `no ground-truth observations` for a GT state with no scored observation), never as an error.

Class display names use the chemical shorthand: the canonical class shows the `--canonical-base` letter, `m` shows `5mC`, `21839` shows `4mC`, `a` shows `6mA`, and so on; a code outside the shorthand table falls back to the raw code.

**BREAKING**: `--report` and the single-file interactive HTML channel are gone -- the threshold slider, in-browser zooming and hovering are **no longer available**, and this version ships no replacement. Callers that need to enumerate thresholds must recompute them in a downstream script.

## Metrics TSV

`per_read.metrics.tsv` (written by `--output-dir`) holds this run's metrics as a single-row TSV: one header row plus one data row, using the same serialization convention as the `extract` output.

- Missing values are written as **empty fields** (not `null`, not `NaN`);
- booleans are written as `true` / `false`;
- numbers keep their own scale (`accuracy` is a `[0, 1]` fraction, never a percentage).

The column set is decided by the report alone; the order is fixed in four blocks:

| Block | Columns | Present when |
|-------|---------|--------------|
| scalars | `n`, `n_correct`, `n_incorrect`, `accuracy`, `insufficient`, `specificity`, `n_positive`, `n_negative`, `site_accuracy`, `site_roc_auc_ovr`, `roc_auc`, `pr_auc` (`roc_auc` / `pr_auc` carry values in a binary report, empty in a multiclass one) | always |
| per-class | `n_<label>` / `precision_<label>` / `recall_<label>` / `f1_<label>` (in `per_class` key order) | at least one scored observation |
| macro | `macro_precision` / `macro_recall` / `macro_f1` | three or more classes |
| per-label | two whole groups: every `roc_auc_ovr_<label>` first, then every `pr_auc_ovr_<label>` -- not interleaved per label. `macro` is a label of both groups; the labels are the class set plus `macro` -- a prediction-only class never appears -- and they follow the report's own key order | multiclass only |

A missing block produces no columns at all rather than empty ones. A multi-dataset run produces **one pooled row** (not one row per dataset) -- the TSV and the report printed to stdout come from the same source.

## Validation

Two structural checks run before any scoring; both fail with a `ValueError` (exit code 1):

- **Pooled 2/3-class check.** The union of the expected states over every dataset, taken after canonical-base filtering, must hold exactly two or three distinct states and must include `-`. The messages are `expected_state pool must contain exactly 2 or 3 distinct states, got N: ...` and `expected_state pool must include '-' (canonical), got: ...`; a pool that canonical-base filtering emptied completely reports `no ground-truth sites remain after canonical-base filtering (0 distinct expected_state)`. A single dataset may hold only one state as long as the pool over all datasets has two or three (an all-`-` dataset beside a `{-, m}` one is fine). These pooled errors carry no `dataset N:` prefix.
- **Per-dataset modification-code coverage.** Every non-`-` state of a dataset's truth must be declared on the target canonical base by at least one read of that dataset. The message names the missing and the declared codes and carries the dataset index:

  ```
  error: dataset 2: expected_state(s) not declared on canonical base C: m (declared: h)
  ```

  Declaring more codes than the truth needs is fine. When no read declared **any** code on the target base (`declared` is empty), the message appends the cause observed by the scan: `scanned 0 records -- the BAM is empty`, `scanned N records, none mapped -- the reads did not align to the reference`, `scanned N mapped reads, none carries MM/ML tags`, or `N mapped reads carry MM/ML tags but none declares a code on canonical base C (declared on other canonical bases: A:a)`. The `[align] ... done` log line also reports `reads_total` and `reads_mapped` for the dataset.

Because the pooled check runs first, a produced report's `gt_states` always holds two or three states.

## Canonical base filtering

`--canonical-base` is always in effect: every GT site carries its own canonical base in BED column 4, and the filtering runs right after the BED is parsed, before any scoring. Sites whose column 4 differs from the target are removed: skipped silently -- they do not raise errors, do not appear in `uncovered`, and do not contribute to any metric denominator. Column 4 accepts only `A`/`C`/`G`/`T` (case-insensitive); anything else, including `N`, is a parse error (exit code 1) carrying the 1-based BED line number of the offending row. The reference FASTA is not consulted here; if you have a legacy 4-column truth BED, upgrade it by inserting the canonical base as its 4th column, between `end` and `expected_state`.

A run that scores no observation is not an error -- the report is then the empty corpus report (`n = 0`, metrics `null`, `insufficient = true`). But a truth that loses **all** of its states to the filtering is a validation failure, not an empty report (see "Validation").

## Multiple datasets

When `--input` / `--truth` / `--reference` are given more than once, the Nth `--input`, the Nth `--truth` and the Nth `--reference` all pair up. Each dataset is realigned against its own reference and scored independently against its own truth, then all observations are pooled into one summary report. The report schema and JSON keys are unchanged.

- Counts are summed and metrics (Acc, confusion, P/R/F1, specificity, AUC) are computed over the pooled corpus -- the same numbers as putting every (GT, prediction) observation into one evaluation set.
- `uncovered` entries are concatenated in dataset order; within each dataset the sites stay internally sorted. A site may appear several times when multiple datasets leave it uncovered (one entry per dataset, `expected_state` from that dataset's truth).
- Contradictory truths across datasets are allowed: each observation is scored against its own dataset's truth.
- Passing the same BAM twice counts its observations twice.
- One global `Filter` applies to every dataset.
- A dataset that fails to realign or score reports its 1-based index: `dataset 2: <error>`. The pooled validation errors are the exception -- they carry no index. With concurrency, the lowest-numbered failure is the one reported -- the same dataset a serial run would have hit first.
- Datasets are realigned and scored concurrently (see "Parallelism and the run log"), and the pooled report is identical to a serial run's.
- A different number of `--input` / `--truth` / `--reference` occurrences is a usage error (exit code 2).
- The index prefix is applied uniformly, so a single-dataset run also reports `dataset 1: <error>`.
- `--output-dir` covers the whole pooled run: the four figures are drawn over the pooled corpus and the metrics TSV holds the pooled row (**BREAKING**: the figure report used to be mutually exclusive with multiple datasets).

## Parallelism and the run log

**Datasets run concurrently.** With more than one dataset, the wrapper realigns up to `max(1, cpu_count // --threads)` of them at a time (never more than the dataset count), and then scores all datasets concurrently in the compiled extension -- one dataset per worker, each with its own reader and accumulator. `--threads 1` therefore does not serialize the datasets: it gives each of them a single in-process mapper and lets them share the box.

The pooled report does not depend on any of this: a serial and a concurrent run produce the same numbers. A failing dataset is reported by its 1-based index -- the lowest-numbered failure when several fail at once, i.e. the one a serial run would have hit first.

**The run log** is an append-only stream of stamped lines (`YYYY-MM-DD HH:MM:SS [phase] message`), flushed per line so `tail -f` follows a long run. Every run echoes it to stderr; `--output-dir DIR` additionally writes `DIR/evaluate.log` with the same lines. Phases:

| Phase | Lines |
|-------|-------|
| `run` | run start: dataset count, canonical base, mapper threads, dataset concurrency, the effective filter (`off`, `default`, or the non-default knobs) |
| `align` | per dataset: start (with both paths) and done (with the elapsed seconds, `reads_total`, and `reads_mapped`) |
| `score` | scan start, periodic per-dataset progress (`dataset N/M records=R elapsed=S`), and the pooled `n` with the scan's elapsed time |
| `outputs` | how many files were written into `--output-dir` (only with that flag) |
| `done` | total elapsed seconds |
| `error` | the exception, written before it is re-raised; the lines already written survive |

Progress lines appear about every 10 seconds and each dataset's counter is monotonic; a run shorter than the interval prints none. The log also names datasets, so it is not a substitute for the single `error: ...` line the CLI prints on failure.

## Breaking changes

Changes in this version that need command lines edited:

- `--figure` and `--metrics-tsv` are gone; `--output-dir DIR` writes the four PNGs and `per_read.metrics.tsv` instead;
- `--reference` is required: every dataset is realigned before scoring, and an already-aligned input is realigned too;
- the report gained the `gt_states` key (the pooled ground-truth states);
- the confusion matrix went from `labels × labels` to `|gt_states| × |gt_states|`; the `per_class` keys and the one-vs-rest maps are keyed by the class set; an out-of-class prediction still counts in `n_incorrect` and accuracy;
- a binary corpus with a prediction-only class takes the binary branch again (`roc_auc` / `pr_auc` carry values, `macro` / `site_*` are `null`), and the TSV column set follows;
- `insufficient` is no longer triggered by a prediction-only class.

Changes from earlier releases that still affect comparability:

- `specificity` is a **macro average over `(chrom, pos)` groups** instead of a global micro average: each site is scored first and the sites are then averaged;
- the top-level JSON key `score_dist` became `channel_dist` (it bins the per-channel probabilities, the same scores the ROC/PR curves use);
- the confusion-cell percentage is the share of the **row total** instead of the share of the row maximum;
- the figure report went from an interactive HTML file to PNG, and is now written by `--output-dir` instead of `--figure`.

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
| 0 | Report printed (and artifacts written) |
| 1 | Invalid ground truth (missing columns, invalid canonical_base, illegal token, conflicting site, empty BED), pooled state-set or per-dataset mod-code validation failure, or any other `ValueError` / `AlignRuntimeError` |
| 2 | Input BAM or reference FASTA missing / unreadable, an unindexable (e.g. read-only) reference, unknown region, `--output-dir` pointing at an existing file, or any other `OSError` / `FileNotFoundError` |
| 2 | `--input` / `--truth` / `--reference` occurrence counts differ (usage error) |
| 2 | Usage errors: `--canonical-base` not given; `--threads` below 1; `--figure` / `--metrics-tsv` removed (unrecognized arguments) |

One known divergence: a reference FASTA with duplicate contig names makes `evaluate` exit 1 (the failure travels as a `ValueError`), while `align` exits 2 for the same file -- one `align()` call, two commands, two codes, pinned by a test.

Errors do not print a partial JSON report.

## Examples

```bash
# Human-readable summary
poremod evaluate --input in.bam --truth gt.bed --reference ref.fa --canonical-base C

# Full JSON report
poremod evaluate --input in.bam --truth gt.bed --reference ref.fa --canonical-base C --json

# Figures, the metrics TSV and the run log (one directory for all six files)
poremod evaluate --input in.bam --truth gt.bed --reference ref.fa --canonical-base C --output-dir out/

# One mapper per dataset, datasets sharing the box
poremod evaluate --input a.bam --truth a.bed --reference a.fa --input b.bam --truth b.bed --reference b.fa --canonical-base C --threads 1 --json

# Same Filter knobs as extract (no min-coverage)
poremod evaluate --input in.bam --truth gt.bed --reference ref.fa --canonical-base C --region chr1:0-2000 --min-probability 0.7

# Pool two datasets into one report; the Nth --reference pairs too
poremod evaluate --input a.bam --truth a.bed --reference a.fa --input b.bam --truth b.bed --reference b.fa --canonical-base C --json
```

## Notes

- Values come from the core accumulator; the CLI does not recompute metrics.
- `--filter-percentile` (default: off) is accepted for Filter compatibility with extract/pileup; evaluate does not apply the percentile estimate, so the output changes only when `--min-probability` is set.
- The pooled state check rejects a degenerate corpus before scoring, so `gt_states` always holds two or three states; `insufficient` is still a per-run property -- a class that never appears in the truth has a `null` one-vs-rest AUC even on a complete run.
