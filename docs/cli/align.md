# poremod align

Re-aligns modification-carrying sequences (FASTA / FASTQ / modbam inputs) to a reference FASTA with mappy, producing a modbam that preserves the modification information (MM/ML tags).

## Usage

```bash
poremod align --input INPUT --reference ref.fa --output out.bam \
  [--preset map-ont] [--threads N[M]] [--eqx] [--secondary {yes,no}] \
  [-k N] [-w N] [-r N[,N]] [-A N] [-B N] [-O N[,N]] [-E N[,N]] [--score-N N]
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input path: FASTA / FASTQ / modbam (.bam/.sam, detected by extension) |
| `--reference PATH` | yes | Reference FASTA (a `.fai` is generated automatically at startup) |
| `--output PATH` | yes | Output modbam path (a `.bai` is built alongside) |
| `--preset STR` | no | mappy alignment preset (default `map-ont`) |
| `--threads N[M]` | no | Mapping worker processes, optionally followed by the `samtools sort` thread count: `--threads N` maps with N workers and sorts with the default 4; `--threads N M` maps with N workers and sorts with M (default `4` mappers and `4` sort threads; `--threads 1` maps serially) |
| `-k N` | no | k-mer size, 1-28 (default: preset's k) |
| `-w N` | no | minimizer window size, >= 1 (default: preset's w) |
| `-r N[,N]` | no | bandwidth `bw[,bw_long]` (default: preset's bandwidth) |
| `-A N` | no | matching score (other scoring slots follow the preset) |
| `-B N` | no | mismatch penalty (other scoring slots follow the preset) |
| `-O N[,N]` | no | gap open penalty `q[,q2]` |
| `-E N[,N]` | no | gap extension penalty `e[,e2]` |
| `--score-N N` | no | bonus score on ambiguous bases |
| `--eqx` | no | Emit `=`/`X` CIGAR operators instead of `M` (default off) |
| `--secondary {yes,no}` | no | Write secondary (`0x100`) records (default `no`); `no` suppresses them, keeping supplementary (`0x800`) and unmapped records, matching minimap2 `--secondary=no` |

## Input requirements

- **FASTA/FASTQ**: the readid description must carry SAMtag-style modification info:

  ```
  >read1 MM:Z:C+m,0,1; ML:B:C,128,200
  ```

  `MM:Z:` is the BAM `MM` tag string verbatim; `ML:B:C` is a comma-separated list of 0–255 probabilities, byte-identical to the output BAM's `ML` tag (zero translation). The ML value count must equal the probability count declared by MM.
- **modbam**: extracted and re-aligned from the primary records (when the input is an already-aligned BAM): only primary records (flag carrying 0x000/0x010/0x004) are taken; secondary (0x100) and supplementary (0x800) records are skipped wholesale; a read present only in those two forms is dropped, so `reads_total` counts reads, not records. From each primary, the full forward sequence is restored (a reverse-strand record's stored SEQ is flipped back to the forward orientation) and the full raw `MM`/`ML` is passed through verbatim (no re-derivation, no trimming to the aligned segment); the original qname is preserved. A read with MM but without ML raises an error naming the read (no silent drop).
- Missing MM/ML, parse failures, or count mismatches: an error naming the read (no silent drops).

## Output

- `@HD SO:coordinate` + reference `@SQ` lines + an `@PG` line (with mappy version, preset, threads; passed tuning flags and `--eqx` are recorded as a `tune:` tail, e.g. `tune: -k17 -r500,10000 -A3 --eqx`; the default `--secondary=no` appends its own `; secondary=no` segment after that tail)
- Records sorted by coordinate; a read with multiple hits produces one record per hit, classified by minimap2 flag semantics: the first non-secondary hit is primary (`0x000`/`0x010`), later non-secondary hits are supplementary (`0x800`), every remaining hit is secondary (`0x100`) and is written only under `--secondary=yes` (the default `no` drops it); unmapped reads are written with `flag=4` (no CIGAR, MM/ML preserved); reverse-mapped reads get `flag=16` (sequence stored reverse-complemented, MM/ML values unchanged)
- Every mapped record — primary or non-primary — clips with `S` and writes the whole read's `SEQ`/`QUAL`, stored in that record's own orientation (forward keeps the read, reverse holds its reverse complement). Flag, coordinates, CIGAR and mapq are all that separate a non-primary record from a primary one; no record clips with `H`. `minimap2 -ax map-ont --eqx -Y` produces the same shape, so its records and these can be compared as-is
- `MM`/`ML` are written on every record that carries them, non-primary included, verbatim — never re-derived or trimmed to the alignment segment. A record's stored `SEQ` is flipped back to the read's original orientation before the tags are decoded, so each record contributes only the calls inside its own alignment segment: the soft-clipped remainder maps to no reference position and is not counted. `-Y` itself writes no modification tags, so the shape matches it while the tags are this entry's addition — ignore `MM`/`ML` when diffing against a `-Y` product
- Downstream counts follow from that: a read whose hits span several sites now contributes at each of them, so a site's `n_valid_cov` can exceed the number of reads covering it
- `QUAL` is preserved: every record carries the input quality values, base for base with `SEQ`. A reverse (`0x10`) record stores them reversed, matching its reverse-complemented `SEQ` — reversed only, never complemented. Input without quality (FASTA, or a modbam record whose `QUAL` is `*`) writes `*`, and a FASTQ quality byte outside the printable Phred+33 range (33..126) aborts the run with exit 2
- `--secondary` defaults to `no`: `0x100` records are suppressed, supplementary and unmapped records stay, and every surviving record is field-identical to the `--secondary=yes` product (the switch never changes the alignment parameters, so primary mapq is untouched). Since non-primary records carry `MM`/`ML`, the switch is not transparent to `pileup`: the calls made by those records are not counted when they are suppressed. Note that this is *not* equivalent to the pre-0.4.0 output, which emitted exactly one record per read — filter by flag (e.g. `samtools view -F 0x900`) if you need that semantics
- A `.bai` index is built; the `MN` tag is not kept (the canonical probability is optional)
- With `--eqx` the CIGAR carries `=` (identical base) and `X` (differing base) instead of `M`; alignment coordinates, mapq, flag and `MM`/`ML` passthrough are unchanged

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
| 0 | Alignment finished |
| 1 | Runtime alignment failures: mappy failed to load the reference, samtools sort failed, etc. |
| 2 | Missing/unparsable header, invalid `--threads`, input or reference missing, reference contains duplicate contig names, output path not writable, invalid tuning values (`-k` out of range, three-part `-r`, non-integers), a FASTQ quality byte outside the printable Phred+33 range |

## Examples

```bash
# FASTQ (readid carries MM/ML); default 4-process parallel mapping
poremod align --input reads.fastq --reference ref.fa --output out.bam

# Serial mapping (lower reference-index memory)
poremod align --input reads.fastq --reference ref.fa --output out.bam --threads 1

# modbam re-alignment (modification calls survive the new alignment)
poremod align --input in.bam --reference ref.fa --output reout.bam

# Custom preset
poremod align --input reads.fastq --reference ref.fa --output out.bam --preset map-hifi

# Tuning knobs share minimap2 CLI names; partial-override semantics
poremod align --input reads.fastq --reference ref.fa --output out.bam -k 17 -r 500,10000 -A 3

# --eqx: =/X CIGAR operators instead of M (minimap2 name, no value)
poremod align --input reads.fastq --reference ref.fa --output out.bam --eqx

# --secondary=no: drop 0x100 (secondary) records; 0x800 (supplementary) stay
poremod align --input reads.fastq --reference ref.fa --output out.bam --secondary=no
```

## Notes

- Contig names in the reference FASTA must be unique (duplicates are rejected).
- `.cram` input is untested (depends on htslib's local reference resolution configuration).
- With the default `--threads 4`, each worker process loads its own mappy reference index; peak memory scales roughly with the process count. Use `--threads 1` for large genomes.
- Alignment results are streamed through an OS pipe to `samtools sort` (`samtools` >= 1.11 must be on PATH), which writes the sorted BAM and its `.bai` index in the same pass (`--write-index`); no full in-memory sorting in the CLI, and no second read of the sorted BAM to build the index.
- `chunk_size` (internal default 1000) is not exposed.
- Tuning flags use minimap2 CLI names with partial-override semantics: only the slots you pass are overridden; everything else follows `--preset`.
- `-z` / `--end-bonus` / `-q` have no mappy constructor entry and are not supported.
- `=` means literal base identity, not "agrees with the reference": an `N` aligned against an `N` is reported as `=` (minimap2 semantics, kept as-is).
- Downstream tooling that concatenates `M` runs must understand the `=`/`X` operators (SAM ops 7/8) before consuming an `--eqx` product.
- `--eqx` is a bare flag with no value; `extra_flags` and the other minimap2 `MM_F_*` bits are not exposed.
- Secondary/supplementary records are classified from mappy's hit order (the first non-secondary hit is treated as `sam_pri`); the read's own query name is passed to `mappy.Aligner.map`, so this reproduces minimap2's `mm_set_sam_pri` including when hits tie on score -- minimap2 folds a hash of the query name into its hit order, and poremod now hands it the same name. Upstream, the secondary flag itself is never suppressed: `MM_F_NO_PRINT_2ND` is a no-op through mappy's map path and is not passed, so `--secondary` is implemented in poremod's own mapping layer, not through mappy flags.
- Non-primary records match `minimap2 -Y`, which soft-clips both kinds and writes the whole read's `SEQ`/`QUAL`, so a `-Y` product and this one can be compared directly. Without `-Y` minimap2 is inconsistent with itself — `S` + `SEQ=*` on secondary, `H` + the aligned segment on supplementary — so a plain `minimap2 -ax map-ont` product necessarily differs on those records. The `MM`/`ML` this tool writes on them are an addition: ignore those two tags when comparing.
- `-x sr` and `-x splice:sr` set `NO_PRINT_2ND` in the stock minimap2 CLI (2.31), so those presets write no secondary records by default; plain `-x splice` does not. mappy does not consume that bit, so this tool's `--secondary` switch decides on its own: the default `no` writes no secondary records under any preset; `--secondary=yes` emits them even under those presets.
- The output side's primary/secondary/supplementary classification is regenerated by mappy from the new alignment and does not inherit the input record's flag. Each output record contributes only the modification calls inside its own aligned segment (the soft-clipped bases map to no reference position and are skipped).
- Re-aligning align's own output is content-idempotent field-by-field: with the same input modbam, the same reference, and the same alignment parameters (preset/tune/secondary/threads), the output records' qname, flag, reference position, CIGAR, SEQ, QUAL, mapq and MM/ML match the first run field-for-field. If any round uses `secondary=True` and the next uses the default, the next round drops the secondary records and the record count changes -- that is not idempotence.
- The output is never byte-equal to the first run: the samtools `@PG` `CL` line carries the output path (different each round), and align rebuilds the header from scratch every round (the `@PG` chain does not accumulate; poremod's own `DS` line is unchanged). To assert idempotence, compare the record multiset (all fields) -- not the file hash or record order.
- The relative order of records at identical coordinates can shift once between round 1 and round 2: round 2 re-extracts reads in the coordinate order of their primary records, so the input order to the sort changes. samtools sort preserves input order within identical coordinates (verified in isolation); round 2 onward converges to a fixed point.
