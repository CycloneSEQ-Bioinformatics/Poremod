# poremod pileup

Per-site modification statistics (bedMethyl counting semantics, identical to rust_modbam), written as TSV. Sites are kept only where the genomic base matches the modification's canonical base; that base comes from `--reference` when it is given and from each read's CIGAR when it is not (which limits the no-reference route to input whose matches are spelled `=`/`X`, i.e. `M`-free). `--context` is optional CpG labeling: it never drops sites.

Modification codes are described in the [Modification codes](../../README.md#modification-codes) section of the README.

## Usage

```bash
poremod pileup --input in.bam --output sites.tsv [--reference ref.fa] [filter args] [--context cpg] [--mod-code m]
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input modbam path (a `.bai` index is required with `--region`) |
| `--output PATH` | yes | Output TSV path. A path ending in exactly `.gz` is written as gzip; decompressed it is byte-identical to the same run without the suffix, and the header carries no timestamp or original filename, so the same input writes byte-identical files twice. Any other suffix -- `.GZ`, `.bgz` and `.gz.bak` included -- is still written as plain text. |
| `--reference PATH` | no | Reference FASTA supplying the genomic base of each call. Sites are kept only when that base matches this row's modification (forward: canonical, e.g. `a`→A, `m`/`21839`→C; reverse: complement, `a`→T, C mods→G). A `.fai` is generated automatically before use. May be omitted when the input spells matches with `=`/`X` (e.g. a `poremod align --eqx` product): the base then comes from the CIGAR, where `=` positions resolve to the read's own base and `X` positions are dropped (exactly for `+`-strand calls, approximately for `-`-strand ones -- see Notes). That rule is about concrete canonical bases: an `N` canonical names no base, so a wildcard call has none to resolve and is kept wherever the CIGAR places it on the reference. A mapped read that still carries calls after filtering while using `M` cannot resolve a base, so it needs the FASTA (an unmapped read has no reference coordinates, so it contributes no site and is skipped; a file whose every read is unmapped is refused as not aligned). Required with `--context cpg`. |
| `--context {all,cpg}` | no | Annotation only; never drops sites. Default `all` does not label CpG (`context` column is always `-`). `cpg` labels C sites as `CpG` / `non-CpG`; non-C sites (including 6mA `a`) stay `-`. Legal with any `--mod-code`. |
| `--region STR` | no | Region filter: contig name (`chr1`) or `chr1:1000-2000` (0-based, half-open) |
| `--min-probability FLOAT` | no | Minimum call probability threshold, on the internal `[0, 1]` probability scale |
| `--filter-percentile FLOAT` | no | Opt-in fail-threshold estimation: when `--min-probability` is unset, the input is scanned once and one threshold per mod code is estimated from the calls that reach the threshold stage (lowest `p` fraction dropped, `[0, 1]` clamped; calibrated against the `(2q+1)/512` ML lattice). **The estimate is scoped to `--region`**: with a region only the calls inside it are scanned; without one, the whole file. Calls strictly below their code's threshold are counted in `n_fail` and dropped from `n_valid_cov`; canonical calls never fail. Default: off. `--min-probability` takes precedence when both are given and the output is then byte-identical to a run without this flag. The estimate costs one extra scan (same scope). |
| `--edge-filter S E` | no | Edge filter (bases to ignore at the start / end of each read) |
| `--skip-mode {implicit,explicit}` | no | How to treat positions the MM tag implies but lists without an ML probability. A segment whose code part ends in the `?` skip marker marks unlisted positions as no-calls; a `.` there or no marker at all marks them as implicitly canonical. `implicit` (default) keeps the inferred-canonical calls, which count towards `n_canonical` and `n_valid_cov`; `explicit` drops them before aggregation, so a site covered only by inferred calls emits no row and a mixed site loses canonical coverage. It cannot turn a `?` no-call into a canonical call. Read-time filter only -- not the same as `convert --explicit`, which rewrites the BAM's MM tags. |
| `--exclude-supplementary` | no | Drop supplementary (FLAG 0x800) records; secondary (0x100) records are always dropped regardless of this flag. |
| `--min-coverage INT` | no | Drop sites below this coverage |
| `--mod-code CODE` | no | Target modification code (e.g. `m`). Omit to emit every modification (the default). |

## Output (TSV, 12 columns)

```
chrom	start	end	mod_code	strand	n_valid_cov	n_mod	n_canonical	n_other_mod	n_fail	percent_modified	context
```

Coordinates are 0-based. `end` is exclusive (a single-base site has `end = start + 1`). `percent_modified` is on a 0–1 scale.

| Column | Meaning |
|--------|---------|
| `chrom` | Reference contig |
| `start` | Site start, 0-based |
| `end` | Site end, 0-based exclusive (single-base site = `start + 1`) |
| `mod_code` | Modification code this row counts (e.g. `m`) |
| `strand` | `+` / `-` / `*` |
| `n_valid_cov` | Valid coverage = `n_mod + n_other_mod + n_canonical` |
| `n_mod` | Calls classified as this row's modification |
| `n_canonical` | Calls classified as canonical |
| `n_other_mod` | Calls with another modification on the same canonical base |
| `n_fail` | Calls below the probability threshold |
| `percent_modified` | `n_mod / n_valid_cov`, on a 0–1 scale (`0` when `n_valid_cov == 0`) |
| `context` | Label column only. Always `-` in `all` mode; in `cpg` mode, C sites are `CpG` / `non-CpG` and non-C sites are `-`. A `-` does not mean the row was dropped. |

Invariant: `n_mod + n_other_mod + n_canonical == n_valid_cov`.

CpG determination rule (on the forward C-G dinucleotide): a forward-strand site is CpG when its `C` is followed by `G`; a reverse-strand site is CpG when its `G` is preceded by `C`; a contig end without a neighbor base is `non-CpG`.

### Example output (first 8 of 17 rows)

The command below runs on the corpus fixture committed under `tests/fixtures/pileup/`, so this block can be reproduced line for line -- it is that command's real output, verbatim:

```bash
poremod pileup --input tests/fixtures/pileup/corpus.bam --output sites.tsv --reference tests/fixtures/pileup/corpus.fa --context cpg --region ptg000001c:36-46
```

```
chrom	start	end	mod_code	strand	n_valid_cov	n_mod	n_canonical	n_other_mod	n_fail	percent_modified	context
ptg000001c	36	37	21839	+	1	0	1	0	0	0.0	CpG
ptg000001c	36	37	m	+	1	0	1	0	0	0.0	CpG
ptg000001c	37	38	21839	-	2	1	1	0	0	0.5	CpG
ptg000001c	37	38	m	-	2	0	1	1	0	0.0	CpG
ptg000001c	38	39	a	-	1	0	1	0	0	0.0	-
ptg000001c	39	40	a	-	2	1	1	0	0	0.5	-
ptg000001c	40	41	21839	-	2	0	2	0	0	0.0	non-CpG
ptg000001c	40	41	m	-	2	0	2	0	0	0.0	non-CpG
```

Those eight rows show the shape of the file: one row per (site, modification code), so one site can appear twice; all three `context` labels; and `percent_modified` following `n_mod / n_valid_cov` (`1/2 = 0.5`). `--region` only narrows the window -- drop it to emit the whole contig.

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
| 0 | Output written |
| 1 | Structural MM/ML failure in the input modbam (see `Input validation` above) |
| 2 | `--reference` omitted where it is needed (`M`-spelled input whose mapped read still carries calls after filtering, or any `--context cpg` run), input not aligned/coordinate-sorted (a non-`coordinate` `@HD` `SO`, or a non-empty file whose every read is unmapped), input/reference missing, output path not writable, invalid filter args or `--mod-code` |

## Examples

```bash
# default --context all: no CpG labels (context column is -)
poremod pileup --input in.bam --output sites.tsv --reference ref.fa

# =/X input without a reference FASTA (same sites as with the matching ref.fa)
poremod pileup --input in.bam --output sites.tsv

# CpG labels only; non-CpG and non-C sites are still written
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --context cpg

# With filters
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --min-coverage 10 --region chr1

# Only 5mC sites
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --mod-code m

# 6mA + CpG labels: succeeds; a rows have context -
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --mod-code a --context cpg
```

## Notes

- Breaking output-contract change: the TSV went from 15 columns to 12. The three always-zero columns `n_delete` / `n_diff` / `n_nocall` are no longer written, so the column order is now `chrom start end mod_code strand n_valid_cov n_mod n_canonical n_other_mod n_fail percent_modified context`. A consumer that parses by position rather than by the header line will silently mis-map every column from `n_fail` onward.
- Supplementary (FLAG 0x800) records are kept by default and contribute sites like primary ones; `--exclude-supplementary` drops them. Secondary (0x100) records are always dropped regardless of this flag.
- Counting semantics are identical to rust_modbam's pileup when the same FASTA is passed. Calls whose reference base does not match the modification's expected canonical are not counted; a wildcard `N` canonical names no base and matches every base, so that rule never drops a call on that account.
- `--skip-mode explicit` drops inferred-canonical calls before aggregation. A site covered only by inferred calls emits no row; at a mixed site `n_canonical` and `n_valid_cov` fall while `n_mod` is unchanged, so `percent_modified = n_mod / n_valid_cov` can rise. The stage runs before `--min-probability` / `--filter-percentile` and before `--min-coverage`: skip mode decides which sites and calls exist, the threshold decides which calls fail, and the reduced coverage can then drop further sites. Inferred-canonical calls sit at probability 1.0 on the canonical state, so within the documented `[0, 1]` threshold range they never fail. `implicit` follows the tag: `?` segments still produce no calls and cannot be resurrected. This is not the same as `convert --explicit`: that command rewrites the MM tag to `?`, after which both modes select the same sites.
- Omitting `--reference` on `=`/`X` input matches passing the alignment's own reference FASTA on `=` positions, which is where every `+`-strand call and most `-`-strand calls sit. The one divergence is a `-`-strand call on an `X` position: dropped without the FASTA, kept with it when the base there turns out to be the call's expected (complementary) base. A wildcard (`N`) call is not part of that divergence: it names no base, so no CIGAR operator can contradict it, and it is kept on an `X` position (as anywhere else the CIGAR places it on the reference) both with and without the FASTA -- the FASTA route needs a resolvable base (A/C/G/T/N) there, so an IUPAC ambiguity code such as `R` drops the call where the CIGAR route keeps it. `M`-spelled input whose mapped read still carries calls after filtering cannot be resolved this way at all and exits 2 rather than being silently thinned; an unmapped read is skipped outright (its calls could not have become a site), so its tags are never parsed. That skip is per read: a file whose *every* read is unmapped is refused up front as not aligned and coordinate-sorted, since no site could exist in it.
- `--context` does not filter. The site set and counting columns match the same command with `--context all`; `cpg` only fills the `context` column. To keep CpG rows only, filter that column afterwards.
- `--context cpg --mod-code a` does not error; A sites are written with `context` = `-`.
- CpG annotation needs the reference `.fai` (generated automatically by the command, no manual preparation); the reference FASTA and its `.fai` must stay consistent.
- Annotation loads reference sequence lazily per contig; memory usage in whole-genome runs is roughly the total length of the contigs touched.
- A `--output` path ending in exactly `.gz` is written as gzip (compression level 6, matching `gzip(1)`). The gzip header carries no timestamp and no original filename, so the same input and arguments write byte-identical `.gz` files across runs. Decompressed, the content is byte-identical to the same run without the `.gz` suffix. The suffix is matched exactly: `.GZ`, `.bgz` and any other spelling still write plain text. Failure semantics are unchanged: an unwritable `.gz` path exits 2 and leaves no file, and an empty result is still a valid gzip stream holding only the header line.
