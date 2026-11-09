# poremod convert

Modification-type conversion: drop or merge modification codes, or rewrite the per-segment skip-mode layout, in a modbam and write a new modbam.

Modification codes are described in the [Modification codes](../../README.md#modification-codes) section of the README.

## Usage

```bash
poremod convert --input in.bam --output out.bam [--drop CODE... [--drop-mode {dist,norm,canonical}] | --merge FROM:TO] [--explicit | --compress-threshold FLOAT] [--json]
```

At least one operation argument is required; the two bracketed groups above are independent (see [Arguments](#arguments)).

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input modbam path (no `.bai` index needed) |
| `--output PATH` | yes | Output modbam path (records keep input order; a `.bai` is built when the input has `@HD SO:coordinate`) |
| `--drop CODE...` | one of group A | Modification codes to drop (e.g. `h`). Folding happens **per site** — dropped mass is never spread to other sites — and **locally within each MM segment**: mass is redistributed only among surviving codes of the same segment; survivors in other segments are untouched (see Notes). The default mode `canonical` leaves the remaining codes' probabilities untouched and folds the dropped mass into the implicit canonical class; `dist`/`norm` must not cover every mod code a read declares, `canonical` may (see Notes). Sites left without a mod code are removed from MM |
| `--drop-mode {dist,norm,canonical}` | no | Drop probability mode (default `canonical`). Requires `--drop`; cannot be used with `--merge`. `dist`: per-site redistribute the dropped mass across the remaining codes at that site + the implicit canonical class; `norm`: after deleting the dropped codes, renormalize the remaining codes + canonical at that site to 1; `canonical`: keep survivor probabilities unchanged, fold the dropped mass into the implicit canonical class |
| `--merge FROM:TO` | one of group A | Merge modification codes into target `TO` (e.g. `--merge h:m`): the probability of each code in `FROM` is added to `TO`, and when the sum exceeds 1.0 it is truncated and the remaining surviving codes are renormalized. One `FROM:TO` pair per call; separate multiple `FROM` codes with **commas** (e.g. `--merge h,f:m`), the only form that carries a ChEBI number and a letter code at once (`--merge h,76792:m`) — a per-character form such as `hf:m` is not accepted. Every code in `FROM` must sit on the same canonical base as `TO`: `h,f:m` is valid, `h,a:m` is not |
| `--explicit` | one of group B | Skip-mode op: write the `?` (explicit) skip marker on **every** output segment |
| `--compress-threshold FLOAT` | one of group B | Skip-mode op: keep a position when **any** surviving code probability is `>= FLOAT` (`0`–`1`, inclusive), otherwise remove that position from its segment and recompute the MM deltas and `ML` |
| `--json` | no | Emit conversion stats as JSON (default: line-based `key: value`) |

**Group A** (`--drop` / `--merge`) and **group B** (`--explicit` / `--compress-threshold`) are two independent constraints: at most one argument from each group, and **at least one argument overall**. A group-A and a group-B argument may be combined — the fold runs first, then the skip-mode projection decides which positions are emitted (e.g. `--drop h --compress-threshold 0.8`).

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
| 0 | Conversion succeeded |
| 1 | Invalid argument values (`--merge` not in `FROM:TO` form, an empty `FROM` list or an empty code in it, a `FROM`/`TO` value that is not a single letter code or ChEBI number, illegal `--drop-mode` string, `--compress-threshold` outside `[0, 1]`, etc.), merge codes that do not share the target's canonical base, `dist`/`norm` drop covers all the mod codes a read declares, drop/merge references a code a read does not declare, the input has no MM-tagged read while `--drop`/`--merge` is requested, an input segment carries an explicit `?` marker while a skip-mode operation is requested, etc. Once the output BAM is open, every exit-1 failure removes the output path and its `<output>.bai`. |
| 2 | No operation given (`must specify --drop or --merge, or a skip-mode operation (--explicit/--compress-threshold)`), more than one skip-mode operation given (`--explicit and --compress-threshold are mutually exclusive`), `--compress-threshold` outside `[0, 1]`, `--drop` and `--merge` both given, `--merge` given more than once (`--merge given more than once; one FROM:TO pair per call (e.g. --merge h,f:m)`) or holding a second pair, `--drop-mode` without `--drop`, `--drop-mode` with `--merge`, input missing, output path not writable |

## Examples

```bash
# Drop 5hmC (h); canonical is the default, so m keeps its own probability
poremod convert --input in.bam --output out.bam --drop h

# Drop h with dist instead: h's probability is redistributed per site
poremod convert --input in.bam --output out.bam --drop h --drop-mode dist

# Merge h and f into m (multiple FROM codes are comma separated)
poremod convert --input in.bam --output out.bam --merge h,f:m

# JSON stats
poremod convert --input in.bam --output out.bam --drop h --json

# Make implicit skip positions explicit (write '?' markers on every output segment)
poremod convert --input in.bam --output out.bam --explicit

# Keep only positions where some code reaches probability 0.8
poremod convert --input in.bam --output out.bam --compress-threshold 0.8

# Skip-mode operations combine with drop/merge (fold first, then filter)
poremod convert --input in.bam --output out.bam --drop h --compress-threshold 0.8
```

## Notes

- **BREAKING: `--drop` now defaults to `canonical`** (was `dist`). The fold is **per site** either way — dropped mass is never spread across other sites; only the destination changes. `canonical` folds the dropped mass into the implicit canonical class and leaves the surviving codes' probabilities exactly as they were; `dist` shares the dropped mass among the remaining codes **at that site** plus the canonical class, so their probabilities grow. Pass `--drop-mode dist` explicitly to keep the previous numbers. After drop, a site with no remaining mod code is deleted from MM (implicit canonical); the position set may shrink, never grow.
- **Drop folding is local to each MM segment; mass is never redistributed across segments.** A dropped code only yields its mass to surviving codes **within the same MM segment**. When multiple codes on the same base are written as separate segments (e.g. `C+m,0;C+h,0;`), dropping `h` does not raise `m`'s ML in the other segment even under `dist`/`norm` — the dropped mass silently folds into the implicit canonical class. BAMs merged from different models/software often split codes into separate segments; if you rely on `dist`/`norm` raising survivors, make sure the codes share one segment (e.g. `C+mh,0;`).
- Reads without MM/ML are copied verbatim (QNAME/CIGAR/SEQ/QUAL and all other tags preserved).
- The `MN` tag is removed from the output (the folded canonical implicit probability changed; the old value is stale).
- Index building is skipped when the input is `@HD SO:queryname` / `unsorted` / has no `@HD`.
- `--drop` and `--merge` are mutually exclusive; `--drop` is accepted only once — pass all modification codes in one flag (e.g. `--drop h a`), repeating it is an error (exit code 2); `--drop-mode` requires `--drop` and cannot be used with `--merge`.
- **`--merge` runs at most one fold per call.** Exactly one `FROM:TO` pair is accepted: repeating the flag (`--merge h:m --merge f:m`) exits 2 rather than silently keeping the last pair, and so does a second pair inside one flag. `FROM` holds either a single code (`--merge h:m` — a one-code merge is legal, it still moves that code's mass onto the target) or several comma-separated ones. Commas are what let a ChEBI number ride alongside letter codes; splitting per character would read `76792` as the five codes 7, 6, 7, 9, 2. Every code in `FROM` must sit on the same canonical base as the target, so `h,f:m` is accepted and `h,a:m` fails with `merge code 'a' is on A but the target 'm' is on C; every code in a merge must share one canonical base` (exit code 1); ChEBI codes outside the code table are exempt, since their base is unknown.
- **Skip-mode operations** (`--explicit`, `--compress-threshold FLOAT`) rewrite the per-segment skip layout. They are mutually exclusive with one another (at most one per run) and each may be combined with `--drop` or `--merge`: the drop/merge fold runs first, then the skip-mode projection decides which positions are emitted.
  - `--explicit`: the `?` marker is written on every output segment, so implicit skip positions become explicit.
  - `--compress-threshold FLOAT`: a position is kept when **any** surviving code probability is `>= FLOAT` (inclusive, `0`–`1`), otherwise it is deleted from that segment — its MM delta and its `ML` values are removed and the remaining deltas recomputed.
  - **Volume cost:** `ML` holds exactly one value per (kept position, surviving code) pair, so its length always matches what the output `MM` declares — the same length the input had, minus dropped positions and removed codes.
  - Both operations reject input segments that already carry an explicit `?` marker with `input segment has explicit skip '?' marker; not allowed with skip-mode operations` (exit code 1): a `?` segment declares its unlisted positions as no-calls, and rewriting the segment would lose that distinction.
  - The output keeps the input segment shape — one input segment maps to one output segment; `--compress-threshold` shortens a segment's position list but never adds or splits segments.
- **With no skip-mode operation, the output keeps each segment's own skip marker.** A pure `--drop` / `--merge` does not change skip semantics, so an input segment carrying `.` comes back with `.`, one carrying `?` comes back with `?`, and one with no marker keeps none. No position filtering happens either: only `--drop` / `--merge` take effect, and an input `ML` whose length disagrees with the `MM` deltas is rejected rather than padded.
- **`dist` / `norm` must keep at least one declared mod code in every MM-carrying read.** If the drop list covers every mod type a read declares (including "the BAM only has one mod type and it is dropped"), conversion fails on that read with `drop covers all declared mod codes; nothing to keep` (exit code 1) and no empty MM/ML is written.
- **Drop and merge codes are validated per read.** Every code named by `--drop`, and by `--merge`'s `FROM` set together with its target, must be declared by *each* MM-carrying read. A read that does not declare one fails the conversion at that read — `drop code 'h' not declared in read 'r2'` (exit code 1) — so a file whose reads declare different sets, usually an upstream batch boundary, is rejected instead of silently half-converted. Reads without MM/ML are not subject to this check.
- **An input without any MM-tagged read is refused when a fold operation is requested.** With `--drop` or `--merge`, conversion fails with `input contains no MM-tagged reads; nothing to convert` (exit code 1) instead of writing a copy of the input that looks like a successful conversion. A request with only a skip-mode operation (`--explicit` / `--compress-threshold`) still succeeds, since those operations reference no mod code.
- **A failed conversion leaves no partial output.** Once the output BAM has been opened, any failure — per-read validation, a malformed `ML`, an explicit `?` segment under a skip-mode operation, a write error, or a failed index build — deletes the output path and its `<output>.bai` index (both may already be absent), including a file an earlier run had left there. Failures before the output is opened — a malformed request, `--input` and `--output` naming the same file, an unreadable input, an output that cannot be created — exit non-zero without writing anything and without touching whatever already sits at the output path.
- **`canonical` may drop every declared mod code.** Conversion succeeds and strips MM/ML/MN from reads that had those tags. Since `canonical` is now the default, dropping every declared code is no longer an error unless a redistributing mode is named explicitly.
- **`--fill` was removed in v0.6.1.0.** It only padded absent `ML` slots with `0`; an `ML` whose length disagrees with the `MM` deltas is now rejected as malformed input. Drop the flag.
