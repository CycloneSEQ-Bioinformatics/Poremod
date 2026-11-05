# poremod extract

Per-read modification-call extraction: writes a TSV of sites whose canonical base is in the requested set, listing every modification code declared at that position and its ML byte quality.

## Usage

```bash
poremod extract --input in.bam --output calls.tsv --canonical-base C [--reference ref.fa] [filter args]
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input modbam path (a `.bai` index is required with `--region`) |
| `--output PATH` | yes | Output TSV path. A path ending in exactly `.gz` is written as gzip; decompressed it is byte-identical to the same run without the suffix, and the header carries no timestamp or original filename, so the same input writes byte-identical files twice. Any other suffix -- `.GZ`, `.bgz` and `.gz.bak` included -- is still written as plain text. |
| `--reference PATH` | no | Reference FASTA. Calls are kept only when the genomic base matches the modification (forward: canonical, e.g. `a`→A, `m`/`21839`→C; reverse: complement, `a`→T, C mods→G). A `.fai` is generated automatically before use. Omit it when the input spells its matches with `=`/`X` (e.g. from `poremod align --eqx`): each read's reference bases then come from its CIGAR, and `X` positions are dropped. An input using `M` needs a reference: its reference bases are undetermined, so a read carrying calls exits 2 with an error instead of silently dropping them. |
| `--canonical-base BASE[,BASE...]` | yes | Target canonical base(s), `A`/`C`/`G`/`T`, case-insensitive; comma-separate for several (e.g. `C,A`). One row per call whose `canonical_base` is in the set; the row lists every modification code declared at that position. `--mod-code` is no longer accepted. |
| `--region STR` | no | Region filter: contig name (`chr1`) or `chr1:1000-2000` (0-based, half-open) |
| `--min-probability FLOAT` | no | Global fail threshold on the internal argmax probability (`[0, 1]`; the `call_qual` column is that value's ML byte quantization). Calls below this value are marked `fail=True`; they are still written. Unset: no probability fail marking. |
| `--filter-percentile FLOAT` | no | Opt-in fail-threshold estimation: when `--min-probability` is unset, the input is scanned once and one threshold per mod code is estimated from the calls that reach the threshold stage (lowest `p` fraction dropped, `[0, 1]` clamped; calibrated against the `(2q+1)/512` ML lattice). **The estimate is scoped to `--region`**: with a region only the calls inside it are scanned; without one, the whole file. Calls strictly below their code's threshold get `fail=True` and are still written, so the row count does not change; canonical calls never fail. Default: off. `--min-probability` takes precedence when both are given and the output is then byte-identical to a run without this flag. The estimate costs one extra scan (same scope). |
| `--edge-filter S E` | no | Edge filter (bases to ignore at the start / end of each read) |
| `--skip-mode {implicit,explicit}` | no | How to treat positions the MM tag implies but lists without an ML probability. The tag itself decides which positions can be inferred: a segment whose code part ends in the `?` skip marker marks unlisted positions as no-calls, while a `.` there or no marker at all marks them as implicitly canonical. `implicit` (default) keeps the inferred-canonical calls; `explicit` drops them, so the affected rows disappear entirely. It cannot turn a `?` no-call into a canonical call; only `inferred=true` calls are dropped, and the calls it keeps are unchanged. Read-time filter only -- it does not rewrite the BAM; `convert --explicit` rewrites the MM tags instead. |
| `--exclude-supplementary` | no | Drop supplementary (FLAG 0x800) records; secondary (0x100) records are always dropped regardless of this flag. |
| `--columns COL [COL ...]` | no | Extra TSV columns appended after the default 10 |
| `--sample-reads N` | no | Keep a random sample of N reads (all of them if the input has fewer). The sample is drawn over the reads that reach the output -- those carrying at least one call on a `--canonical-base` base -- so N bounds the reads that reach the output after the other flags, not the reads read. The draw is uniform over the whole input rather than the first N reads, so a coordinate-sorted BAM does not turn the sample into a region bias. Unset: no sampling. |
| `--seed S` | no | Seed for `--sample-reads` (default 0). The same seed and input give the same sample; a different seed gives a different one. Requires `--sample-reads`: a seed on its own exits 2. |

The filter arguments share the same `Filter` pipeline as `pileup` (without `--min-coverage` / `--context`).

## Output (TSV, 10 columns by default)

```
read_id	chrom	read_position	ref_position	ref_mod_strand	canonical_base	mod_code	mod_qual	call_code	call_qual
```

`read_position` is 0-based on the forward query sequence. Coordinates are 0-based.

| Column | Meaning |
|--------|---------|
| `read_id` | BAM QNAME (read identifier) |
| `chrom` | Reference sequence name |
| `read_position` | 0-based position on the forward query (read) sequence |
| `ref_position` | 0-based reference position |
| `ref_mod_strand` | Modification strand on the reference: `mod_strand` XOR-corrected by `ref_strand` (serialized as `+` / `-` / `*`) |
| `canonical_base` | Canonical base the modification sits on (e.g. `m` → C) |
| `mod_code` | Every modification code declared at this position under its canonical base, ascending (numeric ChEBI codes first, then single-char codes), comma-separated |
| `mod_qual` | ML quality byte (`[0, 255]`) of each code in `mod_code`, same order and count, comma-separated; `0` for implicit unmodified |
| `call_code` | Argmax class: a modification code (e.g. `m`) or `-` (canonical) |
| `call_qual` | ML quality byte (`[0, 255]`) of the argmax class; for canonical argmax it is derived from the canonical probability by the same quantization |

`--columns` may append, in the given order, any of: `mod_strand`, `ref_strand`, `fw_soft_clipped_start`, `fw_soft_clipped_end`, `alignment_start`, `alignment_end`, `read_length`, `base_qual`, `modified_primary_base`, `fail`, `inferred`, `within_alignment`, `flag`. Names already in the default 10, or unknown names, are errors.

Optional column meanings:

| Column | Meaning |
|--------|---------|
| `mod_strand` | Strand declared in the MM tag (`+`/`-`), copied verbatim; not flipped by alignment |
| `ref_strand` | Alignment strand (serialized as `+` / `-`) |
| `fw_soft_clipped_start` | Soft-clipped bases at the start of the read (forward coordinates) |
| `fw_soft_clipped_end` | Soft-clipped bases at the end of the read (forward coordinates) |
| `alignment_start` | 0-based alignment start |
| `alignment_end` | 0-based alignment end |
| `read_length` | Read (query) length |
| `base_qual` | Base quality at this position (0-255) |
| `modified_primary_base` | Primary base of the modification; complemented when `mod_strand` is Reverse |
| `fail` | Whether the argmax probability is below `--min-probability` (`true`/`false`); `false` when that flag is unset |
| `inferred` | Implicitly canonical (no probability given; `.` skip marker or unmarked segment) |
| `within_alignment` | Position falls within the alignment region |
| `flag` | BAM FLAG |

`mod_strand` is the second character of each MM segment (the `+` in `C+m`), copied as-is. Alignment orientation does not flip it. When every MM segment is `+` (typical for ONT simplex 5mC, `C+m`), this column is always `+`. Reverse-mapped reads show up in `ref_strand` instead; `ref_mod_strand` is then equal to `ref_strand`. `mod_strand` is `-` only when MM itself declares the complementary strand (e.g. `C-m`).

Serialization rules:

- `None` → empty field
- `bool` → `true` / `false`
- everything else → `str(...)`

Behavior notes:

- One row per call whose canonical base is in `--canonical-base`, including implicit unmodified sites (`mod_qual` is an all-zero list of equal length, `call_code` is `-`; rows are not skipped).
- `call_code` / `call_qual` remain the argmax; they need not equal `mod_code` / `mod_qual`.
- `--min-probability` compares the internal `[0, 1]` probability (not the `call_qual` byte). Calls below the threshold are marked `fail=True` but still written; add `--columns fail` to emit the flag. Unset: no probability fail marking.
- `--filter-percentile` (default: off) estimates that cutoff: the input is scanned once (scoped to `--region` when given -- with a region only the calls inside it are scanned), and each mod code gets the threshold at 0-based rank `floor(clamp(p) * n)` of its modified calls' probabilities (e.g. with `p=0.1` and 10 modified calls for a code, rank `floor(0.1 * 10) = 1` selects the second-lowest value, not the minimum). Calls strictly below their code's threshold get `fail=True` and are still written; canonical calls never fail. `--min-probability` takes precedence when both are given.
- `--skip-mode explicit` removes rows rather than flagging them: only calls the MM tag lists with an ML probability survive, so an unlisted canonical-base position covered by a read emits no row. `--skip-mode implicit` (default) follows the tag, so `?` segments still yield no calls and `.`/unmarked segments yield inferred-canonical rows (`inferred=true` in the optional column). The flag only drops `inferred=true` calls; the calls it keeps are unchanged, including their `fail` marking. It is applied before the threshold stage. Because `--sample-reads` draws from the reads that still reach the output, `explicit` also shrinks that sampling population. This is not the same as `convert --explicit`: that command rewrites the MM tag to `?`, after which both modes select the same rows.
- With no matching rows the header is still written, exit code 0.
- Streaming: rows are written while iterating; no full in-memory accumulation of calls.
- One invocation covers every declared modification type on the requested bases; a site with several codes puts them in one row as a comma-separated list.

### Example output (first 8 of 211 rows)

The command below runs on the fixture BAM committed under `tests/fixtures/extract_ci/`, so this block can be reproduced line for line -- it is that command's real output, verbatim:

```bash
poremod extract --input tests/fixtures/extract_ci/sample.bam --output calls.tsv --canonical-base C,A
```

```
read_id	chrom	read_position	ref_position	ref_mod_strand	canonical_base	mod_code	mod_qual	call_code	call_qual
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	1	1937890	-	A	a	35	-	220
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	7	1937884	-	A	a	49	-	206
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	8	1937883	-	C	21839,m	46,9	-	200
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	10	1937881	-	C	21839,m	247,1	21839	247
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	12	1937879	-	C	21839,m	30,50	-	175
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	15	1937876	-	C	21839,m	9,12	-	234
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	16	1937875	-	A	a	10	-	245
260F300783011_6_3310_3439_102700657_4923_3_18.69	ptg000001c	17	1937874	-	A	a	4	-	251
```

Those eight rows are the first eight calls of a single read: four `A` sites carrying `a`, four `C` sites carrying `21839,m`. At each `C` site the row lists every code declared there, with one `mod_qual` byte per code, while `call_code` -- the argmax -- names a modification on only one of the eight. That is why `mod_code` and `call_code` are separate columns.

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
| 0 | Output written (including no matches, header only) |
| 1 | Structural MM/ML failure in the input modbam (see `Input validation` above) |
| 2 | argparse usage errors such as a missing `--canonical-base` or the removed `--mod-code` |
| 2 | Input missing / cannot be opened, invalid region, output path not writable, invalid filter args, invalid `--columns`, or a read whose `M` CIGAR leaves its reference bases undetermined while no `--reference` was given |

When the input cannot be opened, no header-only empty TSV is left behind.

## Examples

```bash
# Every call on a C or an A, with its declared codes and ML bytes
poremod extract --input in.bam --output calls.tsv --canonical-base C,A --reference ref.fa

# C sites only, with fail/flag appended after the default 10 columns
poremod extract --input in.bam --output calls.tsv --canonical-base C --reference ref.fa \
  --columns fail flag

# Whole contig, or a coordinate interval
poremod extract --input in.bam --output calls.tsv --canonical-base C --reference ref.fa --region chr1
poremod extract --input in.bam --output calls.tsv --canonical-base C --reference ref.fa \
  --region chr1:1000-2000 --min-probability 0.7

# No reference: calls.tsv comes from an =/X rebuild, so each read carries its
# own reference bases (`X` positions are dropped)
poremod extract --input eqx.bam --output calls.tsv --canonical-base C

# A big modbam, capped to a random 2000 reads' worth of C rows. The sample is
# uniform over the whole file and reproducible run to run; --seed takes a
# different one
poremod extract --input big.bam --output sample.tsv --canonical-base C --reference ref.fa \
  --sample-reads 2000
```

## Notes

- Unmapped reads are skipped entirely (their calls have no reference coordinates, so they cannot form rows); empty `chrom` / `ref_position` fields never appear in the output.
- Supplementary (FLAG 0x800) records are kept by default and contribute calls like primary ones; `--exclude-supplementary` drops them. Secondary (0x100) records are always dropped regardless of this flag.
- `--canonical-base` only filters rows; it is not pushed down into the shared `Filter`.
- Calls whose genomic base does not match the modification's expected base, and calls without coordinates, are dropped. The genomic base comes from the reference FASTA, or, without one, from the read's `=` CIGAR positions.
- `--sample-reads` samples **reads**, uniformly over the input; it does not sample rows or genomic positions. A read contributes all of its surviving calls, so a sample of N reads is N reads deep, not N reads' worth of coverage, and the per-site depth of a sampled extract is not uniform across the genome. Use `pileup` with `--min-coverage` for per-site work.
- A `--output` path ending in exactly `.gz` is written as gzip (compression level 6, matching `gzip(1)`). The gzip header carries no timestamp and no original filename, so the same input and arguments write byte-identical `.gz` files across runs. Decompressed, the content is byte-identical to the same run without the `.gz` suffix. The suffix is matched exactly: `.GZ`, `.bgz` and any other spelling still write plain text. Failure semantics are unchanged: an unwritable `.gz` path exits 2 and leaves no file, and an empty result is still a valid gzip stream holding only the header line.
