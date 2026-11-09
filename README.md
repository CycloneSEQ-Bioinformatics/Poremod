# poremod

A high-level base-modification analysis toolkit built on `rust_modbam`:
a unified command-line interface for the complete modbam workflow —
modification-type conversion (`convert`), MM/ML structural validation
(`validate`), mod-tagged re-alignment (`align`), per-site pileup with
optional CpG context annotation (`pileup`), per-read extraction
(`extract`), per-read accuracy evaluation (`evaluate`), and detection QC
histograms (`summarize`).

## Installation

Distribution is via a pre-built wheel. Download the wheel for your platform
from the [releases page](https://github.com/CycloneSEQ-Bioinformatics/Poremod/releases)
and install it:

```bash
pip install poremod-*.whl
```

Requirements:

- Python >= 3.9
- Linux x86_64, Linux aarch64, or macOS on Apple Silicon (macOS 11+)

The runtime dependencies (`mappy`, `pysam`, `matplotlib`, `jinja2`, `pyyaml`,
`markdown`, `colorama`) are installed automatically, and the `poremod`
command becomes available on `PATH`.

`align` additionally needs `samtools` (>= 1.11) on `PATH` -- it streams
the alignment through `samtools sort` to write the sorted BAM and its
`.bai` index in one pass. `samtools` is not a Python package, so install
it separately:

```bash
# Debian / Ubuntu
sudo apt-get install -y samtools

# conda
conda install -c bioconda samtools

# macOS
brew install samtools
```

The other commands (`convert`, `validate`, `pileup`, `extract`,
`evaluate`, `summarize`) need no external tools.

**macOS on Apple Silicon with Python 3.9-3.13:** `mappy` ships no macOS
arm64 wheel for those interpreters, so pip builds it from source and its
universal build no longer compiles with the current clang. Build it
arm64-only:

```bash
ARCHFLAGS="-arch arm64" pip install poremod-*.whl
```

Python 3.14 needs no workaround (`mappy` ships an arm64 wheel there).

On Linux the wheel requires glibc 2.28 or newer.

Do not install the standalone `rust-modbam` package alongside poremod --
the wheel already bundles the `rust_modbam` extension module, and both
packages ship files under the `rust_modbam/` path.

Check the installation with:

```bash
poremod -V
poremod -h
```

## Usage

### Modification codes

A modbam is a SAM/BAM file that follows the official format definition
(<https://samtools.github.io/hts-specs/>) and additionally carries two
standard tags, defined in SAMtags §1.7
(<https://samtools.github.io/hts-specs/SAMtags.pdf>), that store per-read
modification calls such as methylation: `MM` (`Z`, the call positions) and
`ML` (`B:C`, one probability byte per call). A tag pair looks like this:

```text
MM:Z:C+m?,8,8,4,...   ML:B:C,183,171,173,...
```

Each `MM` segment has the layout `<base><strand><code>[,<delta>...]`;
segments are separated by `;` (the tag normally ends with one; `poremod` also
accepts it without the final `;`):

- `base` is the canonical base the modification can sit on: `C`, `A`, `G` or
  `T` (`U` is accepted as `T`)
- `strand` is `+` or `-`
- `code` is one or more single-character modification codes, or one ChEBI
  numeric identifier
- each `delta` is the number of bases of that canonical type to skip before
  the next modified position: the first delta counts from the start of the
  read (delta `0` is the first such base), and each following delta counts
  from the position after the previous one
- a trailing `.` or `?` is the skip marker for the segment's unlisted
  positions of that canonical base: `?` declares them unknown (no call),
  while `.` or no marker declares them implicitly unmodified (probability 0)

When one base carries more than one modification type at the same positions,
both encodings below mean the same thing:

```text
MM:Z:C+mh.,0;         ML:B:C,200,180    two codes in one segment
MM:Z:C+m.,0;C+h.,0;   ML:B:C,200,180    one segment per code
```

`C+mh` lists two codes sharing one delta list. `ML` values are consumed in
`MM` order -- segment by segment, within a segment position by position, and
within a position code by code -- so in both encodings `200` is `m`'s
probability and `180` is `h`'s at the same base.

A modification code names a chemical modification of a base. `poremod` reads
codes from each read's `MM` tag and reports them in the `mod_code` columns of
its TSVs (and in `validate`'s `mod_code_counts`). The code part is either one
or more single-character codes -- several may be combined in one segment
(`C+mh`) -- or a ChEBI numeric identifier, which must stand alone in its
segment. ChEBI codes use their decimal identifier (`21839`, not
`CHEBI:21839`). The segment's canonical base fixes the positions it covers.
`validate` reports a known code used on a different base as
`mm_code_base_mismatch`, and a code outside the table below as
`mm_code_unknown`; both are warnings, not errors.

The 21 supported codes, grouped by canonical base:

| Code | Base | Modification |
|------|------|--------------|
| `m` | C | 5-methylcytosine (5mC) |
| `h` | C | 5-hydroxymethylcytosine (5hmC) |
| `f` | C | 5-formylcytosine (5fC) |
| `c` | C | 5-carboxylcytosine (5caC) |
| `C` | C | any modification of C |
| `21839` | C | 4-methylcytosine (4mC) |
| `19228` | C | 2'-O-methylcytosine |
| `a` | A | 6-methyladenine (6mA) |
| `A` | A | any modification of A |
| `17596` | A | inosine |
| `69426` | A | 2'-O-methyladenosine |
| `g` | T | 5-hydroxymethyluracil (5hmU) |
| `e` | T | 5-formyluracil (5fU) |
| `b` | T | 5-carboxyuracil (5caU) |
| `T` | T | any modification of T/U |
| `17802` | T | pseudouridine |
| `16450` | T | deoxyuridine |
| `19227` | T | 2'-O-methyluracil |
| `o` | G | 8-oxoguanine (8oxoG) |
| `G` | G | any modification of G |
| `19229` | G | 2'-O-methylguanine |

Codes are accepted by `pileup --mod-code` (restrict output to one code),
`convert --drop` / `--merge` (fold one code into another), and the
`expected_state` column of `evaluate`'s BED truth. `summarize` and `extract`
report every code they see; `summarize`'s `mod_base` column and figures add
the chemical display name for codes in the table (`m` -> `5mC`, `21839` ->
`4mC`), while "any" and unlisted codes stay raw.

### convert

Convert modification codes (drop / merge; `--drop-mode dist|norm|canonical`,
default `canonical`):

```bash
poremod convert --input in.bam --output out.bam --drop h
poremod convert --input in.bam --output out.bam --drop h --drop-mode dist
poremod convert --input in.bam --output out.bam --merge mh:m
poremod convert --input in.bam --output out.bam --drop h --json
```

Skip-mode rewrite: at most one of the two flags, each combinable with
`--drop` / `--merge`. `--explicit` writes `?` markers on every output
segment; `--compress-threshold` drops positions whose every code
probability is below the threshold. Without any of them the output is
unchanged apart from the folds.

```bash
poremod convert --input in.bam --output out.bam --explicit
poremod convert --input in.bam --output out.bam --compress-threshold 0.8
```

### validate

Validate MM/ML structure. Exit 0 for clean input or warnings only, 1 when
errors are found. By default the first 10000 records are checked; `--all`
scans the entire file.

```bash
poremod validate --input in.bam
poremod validate --input in.bam --json
poremod validate --input in.bam --all
```

### align

Build a modbam by re-aligning mod-tagged sequences against a reference.
Input may be FASTA, FASTQ, or modbam. Exit 2 on input errors, 1 on mappy
failures.

```bash
poremod align --input reads.fastq --reference ref.fa --output out.bam

# Tuning knobs mirror minimap2 CLI names
poremod align --input reads.fastq --reference ref.fa --output out.bam -k 17 -r 500,10000

# --eqx: =/X CIGAR operators instead of M
poremod align --input reads.fastq --reference ref.fa --output out.bam --eqx

# --secondary=no: suppress 0x100 (secondary) records; 0x800 stay
poremod align --input reads.fastq --reference ref.fa --output out.bam --secondary=no
```

### pileup

Per-site pileup with bedMethyl semantics. `--reference` is required for
`M`-input on a mapped read and for `--context cpg`; it is optional for
`=`/`X` input. An unmapped read is skipped because its calls cannot
become a site.

```bash
poremod pileup --input in.bam --output sites.tsv --reference ref.fa
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --context cpg
poremod pileup --input in.bam --output sites.tsv --reference ref.fa --mod-code m
poremod pileup --input eqx.bam --output sites.tsv
```

### extract

Per-read modification calls as TSV. `--canonical-base` is required;
`--reference` is needed unless the CIGAR spells matches with `=`/`X`
(e.g. `align --eqx` output).

```bash
poremod extract --input in.bam --output calls.tsv --canonical-base C,A --reference ref.fa
poremod extract --input eqx.bam --output calls.tsv --canonical-base C

# Cap a large modbam to a random 2000 reads' worth of rows
poremod extract --input in.bam --output sample.tsv --canonical-base C --reference ref.fa \
  --sample-reads 2000
```

### evaluate

Score per-read modification calls against a five-column BED ground truth
(`chrom`, 0-based `start`, `end`, `canonical_base`, `expected_state`).
`--canonical-base` is required, and GT sites whose column 4 differs are
skipped.

```bash
poremod evaluate --input in.bam --truth gt.bed --canonical-base C
poremod evaluate --input in.bam --truth gt.bed --canonical-base C --json
```

### summarize

Per-mod-code QC histograms. By default 10000 eligible reads are sampled;
`--sample-reads 0` keeps all eligible reads, including unmapped ones.
`--reference` adds site-frequency output.

```bash
poremod summarize --input in.bam --output qc/
poremod summarize --input in.bam --output qc/ --reference ref.fa --figure
poremod summarize --bins qc/bins.tsv --output qc-replot/ --figure
```

### Read-side `--skip-mode`

`--skip-mode` decides how positions the `MM` tag implies but lists
without an `ML` probability are treated. It is a read-side setting and
moves the numbers `pileup` reports; it is unrelated to
`convert --explicit` and `--compress-threshold`, which rewrite `MM`/`ML`
in the output BAM.

- `implicit` (default) keeps them as inferred canonical calls with
  probability 0. They enter `n_canonical` and therefore
  `n_valid_cov = n_mod + n_canonical + n_other_mod` and the
  `percent_modified` denominator; a site whose only calls are inferred
  still produces a row, with `n_mod = 0`.
- `explicit` drops them. They no longer count towards `n_canonical` or
  `n_valid_cov`, and a site left without any call disappears from the
  output entirely.

`pileup`, `extract`, and `evaluate` share the same filter, so use one
value consistently when comparing runs. `extract` marks the affected rows
with its optional `inferred` column and omits them altogether under
`explicit`.

FASTA/FASTQ read IDs must carry modification information in SAMtag style:
`MM:Z:<mm string> ML:B:C,<comma-separated probabilities>`.

## Command docs

Detailed per-command documentation (arguments, output columns, exit codes,
edge cases):

- [convert](docs/cli/convert.md) -- drop / merge modification codes
- [validate](docs/cli/validate.md) -- MM/ML structural validation report
- [align](docs/cli/align.md) -- mappy alignment to build a modbam
- [pileup](docs/cli/pileup.md) -- per-site pileup with CpG context annotation
- [extract](docs/cli/extract.md) -- per-read modification calls as TSV
- [evaluate](docs/cli/evaluate.md) -- per-read accuracy vs BED ground truth
- [summarize](docs/cli/summarize.md) -- per-mod-code QC histograms (no ground truth)

## Authors

Ming Lian 连明 (lianming@genomics.cn)

Jiayuan Zhang 张嘉远 (zhangjiayuan@genomics.cn)

## License and Usage Restrictions
**Research Use Only**  

This software is provided strictly for individual research purposes. Commercial use is strictly prohibited. This means:  
**Allowed:** Personal academic research, personal learning, and non-commercial experimentation.  
**Not Allowed:** Any form of commercial application, distribution, or use that generates revenue directly or indirectly. This includes, but is not limited to, integration into commercial products, offering this software as a service, or using it for commercial gain.  

For commercial licensing or permissions, please contact us.
