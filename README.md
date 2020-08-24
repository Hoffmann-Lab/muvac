# MUVAC
---

MUVAC implements multiple variant calling options from Exome-Seq/WG-Seq or RNA-Seq data. It offers free GATK best-prectices in an optimized, parallelized fashion.

MUVAC leverages on bashbone, which is a bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses. MUVAC makes use of bashbones best-practice parameterized and run-time tweaked software wrappers and compiles them into a multi-threaded pipeline for analyses of model AND non-model organisms. 

## Features

- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks

- For paired-end and single-end derived raw sequencing or prior mapped read data
- Data preprocessing (quality check, adapter clipping, quality trimming, error correction, artificial rRNA depletion)
- Read alignment and post-processing
  - knapsack problem based slicing of alignment files for parallel task execution
  - sorting, filtering, unique alignment extraction, removal of optical duplicates
  - read group modification, split N-cigar reads, left-alignment and base quality score recalibration
- Germline and somatic variant detection from DNA or RNA sequencing experiments plus VCF normalization

# License
---

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download
---

```bash
git clone --recursive https://github.com/koriege/muvac.git
git checkout $(git describe --tags)
```

# Quick start (without installation)
Please see installation section to get all third-party tools set-upped and subsequently all functions to work properly.

Load the library and list available functions. Each function comes with a usage. Or check out the MUVAC usage.

```bash
source activate.sh
bashbone -h
muvac.sh -h
```

# Installation
---

```bash
setup -i all -d <path/to/installation>
source <path/of/installation/activate.sh>
muvac.sh -h
```

## Update to a newer release

```bash
setup -i upgrade -d <path/of/installation>
source <path/of/installation/activate.sh>
bashbone -h
```

# Usage
---

## Retrieve SRA datasets

Use the enclosed script to fetch sequencing data from SRA

```bash
source <path/of/installation/activate.sh> -c true
sra-dump.sh -h
```

## Genome setup

Human genome chromosomes must follow GATK order and naming schema: chrM,chr1..chr22,chrX,chrY. This requirement needs to be fulfilled in all associated VCF files, too. To obtain properly configured genomes and dbSNP files, see blow.

### Retrieve genomes

Use the enclosed script to fetch human hg19/hg38 or mouse mm9/mm10 genomes, annotations and dbSNP.

```bash
source <path/of/installation/activate.sh> -c true
dlgenome.sh -h
```

### Retrieve GATK resources

To obtain panel of normals, common somatic variants and population variants with allele frequencies visit

- HG38: <https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/>
- HG19: <https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37/>

Afterwards, place the files next to your genome fasta file with equal name plus extension suffix as shown below.

- genome.fa.somatic_common.vcf.gz.tbi
- genome.fa.pon.vcf.tbi
- genome.fa.af_only_gnomad.vcf.gz.tbi

Analogously, place a custom dbSNP file next to the genome as genome.fa.vcf. Possible resources are

- HG38: ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
- HG19: ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/	

## Examples

This section showcases some usages without explaining each parameter in a broader detail. Check out the MUVAC help page for more configuration options. Most of them will be opt-out settings.

### Data pre-processing and mapping

Data pre-processing without mapping by segemehl or STAR and without SortMeRNA for artificial rRNA depletion.

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-rrm -no-sege -no-star
```

Data pre-processing, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing, removal of duplicons, clipping of overlapping mate pairs).

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>]
```

Data pre-processing with Illumina universal adapter removal, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing, removal of duplicons, clipping of overlapping mate pairs). More sequences can be found at the resource of Trimmomatic (<https://github.com/timflutre/trimmomatic/tree/master/adapters>).

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -a1 AGATCGGAAGAG [-a2 AGATCGGAAGAG]
```

Data pre-processing, mapping by segemehl and STAR and disabled post-processing (i.e. unique read extraction, sorting, indexing, removal of duplicons, clipping of overlapping mate pairs).

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-uniq -no-sort -no-idx -no-rmd -no-cmo
```

Multiple inputs can be submitted as comma separated list.

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>]
```

Tweak the amount of allowed mismatches in % during mapping

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -d 10
```

### Variant calling

Perform pre-processing, mapping and post-processing with subsequent germline variant calling.

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>]
```

Perform pre-processing, mapping and post-processing with subsequent somatic variant calling with known sites and a panel of normals provided by GATK resources (stored as genome.fa.pon.tar.gz, see additional information section of MUVAC usage).

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>]
```

Perform generation of a custom panel of normals for subsequent somatic variant calling.

```bash
source <path/of/installation/activate.sh>
muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -pon

muvac.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>] -mypon
```

### Start, redo or resume

List all possible break points and keywords to control Rippchen.

```bash
source <path/of/installation/activate.sh>
muvac.sh -dev
```

Use comma separated lists to e.g. skip md5 check and quality analysis.

```bash
source <path/of/installation/activate.sh>
muvac.sh [...] -skip md5,qual 
```

Example how to resume from the segemehl mapping break point after previous data pre-processing.

```bash
source <path/of/installation/activate.sh>
muvac.sh [...] -resume sege
```

Single tasks can be re-computed with the `redo` parameter and a comma separated list of arguments.

```bash
source <path/of/installation/activate.sh>
muvac.sh [...] -redo bqsr,idx,hc
```

# Third-party software
---

## In production

| Tool | Source | DOI |
| ---  | ---    | --- |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BCFtools      | <http://www.htslib.org/doc/bcftools.html>                           | 10.1093/bioinformatics/btr509 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| GATK          | <https://github.com/broadinstitute/gatk>                            | 10.1101/gr.107524.110 <br> 10.1038/ng.806 |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| ReSeqC        | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| vcflib        | <https://github.com/vcflib/vcflib>                                  | NA |
| Vt            | <https://genome.sph.umich.edu/wiki/Vt>                              | 10.1093/bioinformatics/btv112 |

## In preparation

| Tool | Source | DOI |
| ---  | ---    | --- |
| BWA             | <https://github.com/lh3/bwa>                               | 10.1093/bioinformatics/btp324 |
| freebayes       | <https://github.com/ekg/freebayes>                         | arXiv:1207.3907 |
| Platypus        | <https://rahmanteamdevelopment.github.io/Platypus>         | 10.1038/ng.3036 |
| SnpEff          | <https://pcingola.github.io/SnpEff>                        | 10.4161/fly.19695 |
| STAR-Fusion     | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>          | 10.1101/120295 |
| VarDict         | <https://github.com/AstraZeneca-NGS/VarDict>               | 10.1093/nar/gkw227 |
| VarScan         | <http://dkoboldt.github.io/varscan>                        | 10.1101/gr.129684.111 |

# Supplementary information

MUVAC can be executed in parallel instances and thus is able to be submitted as a job into a queuing system like a Sun Grid Engine (SGE). This could be easily done by amending the following code snipped.

```bash
for i in *R1.fastq.gz; do
	j=${i/R1/R2}
	sh=job_$(basename $i .R1.fastq.gz)
	cat <<- EOF > $sh.sh
		#!/usr/bin/env bash
		source <path/of/installation/activate.sh>
		muvac.sh -v 2 -t <threads> -g <fasta> -o <outdir> -l <logfile> -tmp <tmpdir> -1 $i -2 $j
	EOF
	echo "rm -f $sh.+(log|err) && qsub -pe <env> <threads> -l 'h=<hostname>|<hostname>' -S /bin/bash -e $sh.err -o $sh.log -V -cwd $sh.sh"
done
```

In some cases a glibc pthreads bug (<https://sourceware.org/bugzilla/show_bug.cgi?id=23275>) may cause pigz failures (`internal threads error`) and premature termination of toola leveraging on it e.g. Cutadapt. One can circumvent this by upgrading the operating system or making use of an alternative pthreads library and `LD_PRELOAD`

```bash
source <path/of/installation/activate.sh>
LD_PRELOAD=/lib64/noelision/libpthread.so.0 muvac.sh [...]
```