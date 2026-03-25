# poremod

A Python package to analyze CycloneSEQ methylation data and produce HTML report

## Installation

Install with a clean environment:

```shell
$ python -m pip install --upgrade pip && pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
$ conda env create -n poremod_env -f environment.yaml
$ conda activate poremod_env

# Install from source code
$ pip install .
# Or install from compiled wheel file
$ pip install poremod-x.x.x-py3-none-any.whl
```

## Usage

To list all available subcommands, run:
 
```
$ poremod -h
usage: poremod [-h] {align,make_annoDb,bedMethyl,stats,report,run} ...

A Python package to analyze CycloneSEQ methylation data and produce HTML report
Version: x.x.x

optional arguments:
  -h, --help            show this help message and exit

Commands:
  {align,make_annoDb,bedMethyl,stats,report,run}
    align               Align modbam to the reference genome with minimap2
    make_annoDb         Construct annotation files for key genomic functional regions
    bedMethyl           Summarise one BAM with modified base tags to bedMethyl
    stats               Based on the bedMethyl generated in the previous step, the statistics of genome-wide methylation detection results are obtained
    report              Based on the statistical results obtained from the previous step, generate an HTML report
    run                 One command to run poremod
```

### align

Align modbam to the reference genome with minimap2

```
$ poremod align -h
usage: poremod align [-h] --bam BAM --ref REF --threads THREADS --prefix PREFIX

optional arguments:
  -h, --help         show this help message and exit
  --bam BAM          the modbam file
  --ref REF          the reference genome file
  --threads THREADS  set the threads for minimap2
  --prefix PREFIX    set the output prefix, the output filename would be <prefix>.align.sort.bam
```

for example:

```shell
poremod align \
 --bam <bam file> \
 --ref <ref fasta file> \
 --threads <thread number> \
 --prefix <output prefix>
```

### make_annoDb

Construct annotation files for key genomic functional regions

```
$ poremod make_annoDb -h
usage: poremod make_annoDb [-h] --ref REF --gtf GTF --cgi CGI [--upstream UPSTREAM] [--downstream DOWNSTREAM] --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --ref REF             the reference genome file
  --gtf GTF             the gtf file
  --cgi CGI             the annotation file for CpG islands (BED file), which needs to be downloaded from the UCSC Table Browser.
  --upstream UPSTREAM   set the upstream from TSS to define promoter (default: 20000)
  --downstream DOWNSTREAM
                        set the downstream from TSS to define promoter (default: 20000)
  --prefix PREFIX       set the prefix of output files, the output filenames would be <prefix>.promoters.bed.gz, <prefix>.exons.bed.gz ...
```

for example:

```shell
$ poremod make_annoDb \
 --ref <ref fasta file> \
 --gtf <gtf file> \
 --cgi <CpG Island bed file> \
 --prefix <output prefix>
```

Those interval annotation files can be downloaded from the corresponding databases:
- `gtf file`: encode or ensembel database
- `CpG Island bed file`: [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables), select the genome of interest, then select "Group" as "Expression and Regulation" or "Regulation" and "Track" as "CpG Islands", finally get outfile file called "\<genome\>.cpgIsland.bed" by setting "Output format" as "BED - browser extensible data" and "Output filename" as "\<genome\>.cpgIsland.bed"

### bedMethyl

Summarise one BAM with modified base tags to bedMethyl

```
$ poremod bedMethyl -h
usage: poremod bedMethyl [-h] --bam BAM --ref REF [--chrom CHROM] [--start START] [--end END] [--threads THREADS] [--filt] [--combine_strands]
                         [--force] [--temp] --prefix PREFIX

optional arguments:
  -h, --help         show this help message and exit
  --bam BAM          the aligned modbam file
  --ref REF          the reference genome file
  --chrom CHROM      chromosome for which to fetch read. If not specified, it includes all primary chromosomes
  --start START      reference start coordinate. If not specified, the entire chromosome is considered
  --end END          reference end coordinate. If not specified, the entire chromosome is considered
  --threads THREADS  set the threads to process multiple chromosomes
  --filt             filter bedMethyl records without coverage
  --combine_strands  output additonal bedMethy file, which sum the counts from the positive and negative strands into the counts for the positive
                     strand position
  --force            force overwrite if the output file exists
  --temp             keep temporary files
  --prefix PREFIX    set the prefix of output filename(s). The output filename(s) would be <prefix>.bedMethyl.gz and
                     <preffix>.combine_strands.bedMethyl.gz, if add the '--combine_strands' option
```

for example:

```shell
poremod bedMethyl \
 --bam <bam file> \
 --ref <ref fasta file> \
 --threads <thread number> \
 --combine_strands \
 --force \
 --prefix <output prefix>
```

### stats

Based on the bedMethyl generated in the previous step, the statistics of genome-wide methylation detection results are

```shell
$ poremod stats -h
usage: poremod stats [-h] --ref REF --infile INFILE --indir INDIR --outdir OUTDIR [--annoDb {hg38,mm39}] [--annoPrefix ANNOPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --ref REF             the reference genome file
  --infile INFILE       the input bedMethyl file
  --indir INDIR         the input directory to store bedMethyl files for each chromosomes, which should be the <bedMethy_temp> subfolder under the output directory from the previous step (the bedMethyl
                        subcommand)
  --outdir OUTDIR       the output directory
  --annoDb {hg38,mm39}
                        using the built-in annotation database. Please choose to set either the --annoDb or --annoPrefix parameter
  --annoPrefix ANNOPREFIX
                        set the prefix of annotation database. Please choose to set either the --annoDb or --annoPrefix parameter
```

for example:

```shell
$ poremod stats \
 --ref <ref fasta file> \
 --infile <bedMethyl file> \
 --indir <dir> \
 --outdir <outdir> \
 --annoDb <hg38 or mm39>

# or
$ poremod stats \
 --ref <ref fasta file> \
 --infile <bedMethyl file> \
 --indir <dir> \
 --outdir <outdir> \
 --annoPrefix <annotation prefix>
```

### run

One command to run poremod

```
$ poremod run -h
usage: poremod run [-h] --name NAME [--bam BAM] [--dir DIR] --ref REF [--align] [--threads THREADS] --prefix PREFIX [--annoDb {hg38,mm39}] [--annoPrefix ANNOPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --name NAME           the sample name
  --bam BAM             the unaligned or aligned modbam file
  --dir DIR             the directory to store unaligned bam files
  --ref REF             the reference genome file
  --align               set the alignment status of the provided BAM file(s), defaulting to unaligned (if this parameter is not set), otherwise please set this parameter
  --threads THREADS     set the thread number for minimap2 and subcommand bedMethyl (default: 16)
  --prefix PREFIX       set the prefix of output filename(s)
  --annoDb {hg38,mm39}
                        using the built-in annotation database. Please choose to set either the --annoDb or --annoPrefix parameter
  --annoPrefix ANNOPREFIX
                        set the prefix of annotation database. Please choose to set either the --annoDb or --annoPrefix parameter
```

for example:

```shell
# input as unaligned bam files
poremod run \
 --name <sample name> \
 --dir <bam dir> \
 --ref <ref fasta file> \
 --threads <thread number> \
 --annoDb <hg38 or mm39> \
 --prefix <output prefix>

# input as an aligned bam file
poremod run \
 --name <sample name> \
 --bam <bam file> \
 --ref <ref fasta file> \
 --align \
 --threads <thread number> \
 --annoPrefix <annotation prefix> \
 --prefix <output prefix>
```

use by docker

```shell
# input as unaligned bam files
sample_name=<sample_name>
bam_dir=<bam_dir>
ref_file=<ref_file>
threads=<thread number>
out_dir=<out_dir>
docker run -d --rm \
-v $bam_dir:/data/modbams \
-v $(dirname $ref_file):/data/ref \
-v $out_dir:/data/analysis_out \
poremod:latest \
poremod run \
 --name $sample_name \
 --dir /data/modbams \
 --ref /data/ref/$(basename $ref_file) \
 --threads $threads \
 --annoDb <hg38 or mm39> \
 --prefix /data/analysis_out/$sample_name

# input as an aligned bam file
sample_name=<sample_name>
bam=<bam file>
ref_file=<ref_file>
out_dir=<out_dir>
docker run -d --rm \
-v $(dirname $bam):/data/modbams \
-v $(dirname $ref_file):/data/ref \
-v $out_dir:/data/analysis_out \
poremod:latest \
poremod run \
 --name $sample_name \
 --bam $bam \
 --ref /data/ref/$(basename $ref_file) \
 --align \
 --threads $threads \
 --annoDb <hg38 or mm39> \
 --prefix /data/analysis_out/$sample_name
```
