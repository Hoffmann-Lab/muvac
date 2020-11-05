#! /usr/bin/env bash
# (c) Konstantin Riege

options::usage() {
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		MUVAC is an ultra fast germline and somatic variant calling pipeline for model and non-model organisms.
		It implements GATK best practices in an optimized, parallelized fashion.

		VERSION
		$VERSION
		utilizing bashbone $BASHBONE_VERSION

		SYNOPSIS INDEXING
		muvac.sh -x -g genome.fa -gtf genome.gtf

		SYNOPSIS GERMLINE
		muvac.sh -1 R1.fq -2 R2.fq -g genome.fa

		SYNOPSIS SOMATIC
		muvac.sh -n1 normal.fq,normal.fq -t1 tumor1.fq,tumor2.fq -g genome.fa

		BASIC OPTIONS
		-h       | --help                     : prints this message
		-v       | --verbosity [value]        : set level of verbosity. default: 0
		                                        0 - get simple status updates
		                                        1 - get status updates and commands
		                                        2 - get full output
		-o       | --out [path]               : output directory. default: $OUTDIR
		-l       | --log [path]               : output directory. default: $OUTDIR/run.log
		-tmp     | --tmp                      : temporary directory. default: $TMPDIR/rippchen.XXXXXXXXXX
		-r       | --remove                   : remove temporary and unnecessary files upon succesful termination
		-t       | --threads [value]          : number of threads. default: $THREADS
		-mem     | --memory [value]           : amout of memory for creating bam slices and processing them in parallel instances
		                                        available: $MAXMEMORY
		                                        default: 30000 (allows for $MTHREADS instances)
		                                        NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors

		ADVANCED OPTIONS
		-dev     | --devel                    : prints list of keywords in processing order for advanced pipeline control
		-resume  | --resume-from [string]     : resume from a specific pipeline step (see -dev)
		-skip    | --skip [string,..]         : skip specific pipeline step(s). comma seperated (see -dev)
		-redo    | --redo [string,..]         : just rerun specific pipeline step(s). comma seperated (see -dev)

		GENOME OPTIONS
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)
		-gtf     | --gtf [path]               : annotation gtf input. default: [-g].gtf (see dlgenome.sh)
		-s       | --snp [path]               : genome dbSNP input - optional, default: [-g].vcf (see dlgenome.sh)
		-no-dbsnp| --no-dbsnp                 : disbale dbSNP usage for BQSRecalibration and variant calling
		-x       | --index                    : create all requiered genome indices and md5 sums and exit. otherwise create necessary indices on the fly
		-no-sege | --no-segemehl              : disables indexing for segemehl when used with -x
		-no-star | --no-star                  : disables indexing for STAR when used with -x. use when indexing is applied on plug-n-play CTAT resource
		                                        NOTE: md5sum of [-g].star.idx/SA file needs to be manually added to [-g].md5.sh file
		-no-bwa  | --no-bwa                   : disables indexing for BWA when used with -x

		PREPROCESSING OPTIONS
		-no-qual | --no-qualityanalysis       : disables read quality analysis
		-no-trim | --no-trimming              : disables quality trimming
		-no-clip | --no-clipping              : disables removal of poly N, mono- and di-nucleotide ends as well as adapter sequences when used with -a
		-a1      | --adapter1 [string,..]     : adapter sequence(s). single or first pair. comma seperated (e.g. Illumina universal adapter AGATCGGAAGAGC)
		-a2      | --adapter2 [string,..]     : adapter sequence(s). second pair. comma seperated (can be the same as -a1. no revere complement required)
		-cor     | --correction               : enable majority based raw read error correction
		-rrm     | --rrnafilter               : enable rRNA filter

		ALIGNMENT OPTIONS
		-d       | --distance                 : maximum read alignment edit distance in % - default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs - default: 200000
		-split   | --split                    : enable split read mapping e.g. to call variants from RNA-Seq data
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-no-rmd  | --no-removeduplicates      : disables removing duplicates
		-rx      | --regex                    : regex of read name identifier with grouped tile information - default: ^\S+:(\d+):(\d+):(\d+)\s*.*
		                                        NOTE: necessary for sucessful deduplication, if unavailable set to 'null'
		-no-cmo  | --no-clipmateoverlaps      : disables clipping of read mate overlaps
		-no-adgrp| --no-addreadgroup          : disables proper read group modification by Picard
		-no-reo  | --no-reordering            : disables reordering according to genome file by Picard
		-no-laln | --no-leftalign             : disables left alignment by GATK
		-no-realn| --no-realign               : disables indel realignment by GATK
		-no-bqsr | --no-qualrecalibration     : disables any base quality score recalibration (BQSR)
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA, when -no-split is used. default: no BWA mapping
		-no-stats| --no-statistics            : disables preprocessing statistics

		GERMLINE/PON OPTIONS
		-1       | --fq1 [path,..]            : fastq input - single or first pair, comma seperated
		-2       | --fq2 [path,..]            : fastq input - optional. second pair, comma seperated
		-m       | --mapped [path,..]         : SAM/BAM input - comma seperated (replaces -1 and -2)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input - optional. default: custom
		-rgn     | --readgroup-name [string]  : sets custom read group name - use TUMOR or NORMAL for subsequent somatic calls - default: 'SAMPLE'
		-pon     | --panelofnormals           : disables variant calling and instead prepares a panel of normals for subsequent somatic calls
		-no-pondb| --no-pondatabase           : disables creation of panel of normals database

		SOMATIC OPTIONS
		-n1      | --normalfq1 [path,..]      : normal fastq input - single or first pair, comma seperated
		-n2      | --normalfq2 [path,..]      : normal fastq input - optional. second pair, comma seperated
		-t1      | --tumorfq1 [path,..]       : tumor fastq input - single or first pair, comma seperated
		-t2      | --tumorfq2 [path,..]       : tumor fastq input - optional. second pair, comma seperated
		-nm      | --normalmapped [path,..]   : normal SAM/BAM input - comma seperated (replaces -n1 -n2 -t1 -t2)
		-tm      | --tumormapped [path,..]    : tumor SAM/BAM input - comma seperated (replaces -n1 -n2 -t1 -t2)
		-mypon   | --my-panelofnormals        : priorize own panel of normals database over [-g].pon.vcf.gz

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de

		ADDITIONAL INFORMATION
		Human genome chromosomes must follow GATK order and naming schema: chrM,chr1..chr22,chrX,chrY
		This requierement needs to be fulfilled in all additional VCF files, too - see below.

		To obtain panel of normals, common somatic variants and population variants with allele frequencies visit
		HG38: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/
		HG19: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37/
		After download, place files next to your genome fasta file with equal name plus extension suffix as shown
		genome.fa.somatic_common.vcf.gz, genome.fa.somatic_common.vcf.gz.tbi
		genome.fa.pon.vcf.gz, genome.fa.pon.vcf.tbi
		genome.fa.af_only_gnomad.vcf.gz, genome.fa.af_only_gnomad.vcf.gz.tbi

		Analogously, obtain a dbSNP file, extract and re-name it: genome.fa.vcf
		HG38: ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
		HG19: ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/
	EOF
	exit 1
}

options::developer() {
	cat <<- EOF
		DESCRIPTION
		In case of restarting or to resume an analysis use the identifiers below, listed in processing order

		DEVELOPER OPTIONS
		md5   : check for md5sums and if necessary trigger genome indexing
		qual  : quality analysis for input and trim, clip, cor, rrm
		trim  : trimming
		clip  : adapter clipping
		cor   : raw read correction
		rrm   : rRNA filtering
		sege  : Segemehl mapping
		star  : STAR mapping
		uniq  : extraction of properly paired and uniquely mapped reads
		sort  : sorting and indexing of sam/bam files
		slice : better dont touch! slicing of bams for parallelization, needs -prevtmp | --previoustmp [path]
		rg    : read group modification
		rmd   : removing duplicates
		cmo   : clipping mate overlaps
		stats : proprocessing and mapping statistics
		nsplit: splitting split-read alignments
		reo   : bam reordering according to genome
		laln  : left alignment
		bqsr  : BQSRecalibration
		idx   : intermediate and final bam indexing
		pon   : panel of normals
		pondb : panel of normals database
		hc    : haplotypecaller
		mu    : mutect
	EOF
	exit 0
}
		# bt    : bcftools
		# fb    : freebayes
		# pp    : platypus
		# vs    : varscan
		# vd    : vardict

options::checkopt (){
	local arg=false
	case $1 in
		-h   | --help) (options::usage); exit 0;;
		-dev | --devel) options::developer;;

		-r   | --remove) CLEANUP=true;;
		-v   | --verbosity) arg=true; VERBOSITY=$2;;
		-t   | --threads) arg=true; THREADS=$2;;
		-mem | --memory) arg=true; MEMORY=$2;;
		-g   | --genome) arg=true; GENOME=$2;;
		-s   | --snp) arg=true; DBSNP=$2;;
		-gtf | --gtf) arg=true; GTF=$2;;
		-o   | --out) arg=true; OUTDIR=$2;;
		-l   | --log) arg=true; LOG=$2;;
		-tmp | --tmp) arg=true; TMPDIR=$2;;
		-prevtmp | --previoustmp) arg=true; PREVIOUSTMPDIR=$2;;

		-1   | --fq1 | -n1 | --normalfq1) arg=true; mapfile -t -d ',' NFASTQ1 < <(printf '%s' "$2");;
		-2   | --fq2 | -n2 | --normalfq2) arg=true; mapfile -t -d ',' NFASTQ2 < <(printf '%s' "$2");;
		-t1  | --tumorfq1) arg=true; mapfile -t -d ',' TFASTQ1 < <(printf '%s' "$2");;
		-t2  | --tumorfq2) arg=true; mapfile -t -d ',' TFASTQ2 < <(printf '%s' "$2");;
		-m   | --mapped | -nm | --normalmapped) arg=true; mapfile -t -d ',' NMAPPED < <(printf '%s' "$2");;
		-mn  | --mapper-name) arg=true; MAPNAME=$2;;
		-tm  | --tumormapped) arg=true; mapfile -t -d ',' TMAPPED < <(printf '%s' "$2");;
		-rx  | --regex) arg=true; REGEX=$2;;
		-a1  | --adapter1) arg=true; mapfile -t -d ',' ADAPTER1 < <(printf '%s' "$2");;
		-a2  | --adapter2) arg=true; mapfile -t -d ',' ADAPTER2 < <(printf '%s' "$2");;
		-d   | --distance) arg=true; DISTANCE=$2;;
		-i   | --insertsize) arg=true; INSERTSIZE=$2;;
		-rgn | --readgroup-name) arg=true; RGPREFIX=$2;;
		-pon | --panelofnormals) PON=true;;
		-mypon | --my-panelofnormals) MYPON=true;;

		-resume | --resume-from) arg=true; options::resume "$2";;
		-skip | --skip) arg=true; options::skip "$2";;
		-redo | --redo) arg=true; options::redo "$2";;

		-x  | --index) INDEX=true;;

		-cor      | --correction) NOcor=false;;
		-rrm      | --rrnafilter) NOrrm=false;;
		-split    | --split) NOsplitreads=false; NOnsplit=false;;

		-no-qual  | --no-qualityanalysis) NOqual=true;;
		-no-clip  | --no-clipping) NOclip=true;;
		-no-trim  | --no-trimming) NOtrim=true;;
		-no-stats | --no-statistics) NOstats=true;;

		-no-sege  | --no-segemehl) NOsege=true;;
		-no-star  | --no-star) NOstar=true;;
		-no-bwa   | --no-bwa) NObwa=true;;
		-no-uniq  | --no-uniqify) NOuniq=true;;
		-no-sort  | --no-sort) NOsort=true;;
		-no-rmd   | --no-removeduplicates) NOrmd=true;;
		-no-cmo   | --no-clipmateoverlaps) NOcmo=true;;
		-no-addrg | --no-addreadgroup) NOaddrg=true;;
		-no-reo   | --no-reordering) NOreo=true;;
		-no-laln  | --no-leftalign) NOlaln=true;;
		-no-realn | --no-realign) NOrealn=true;;
		-no-bqsr  | --no-qualrecalibration) NObqsr=true;;
		-no-dbsnp | --no-dbsnp) NOdbsnp=true;;
		-no-pondb | --no-pondatabase) NOpondb=true;;

		-no-hc    | --no-haplotypecaller) NOhc=true;;
		-no-mu    | --no-mutect) NOmu=true;;
		-no-fb    | --no-freebayes) NOfb=true;;
		-no-bt    | --no-bcftools) NObt=true;;
		-no-pp    | --no-platypus) NOpp=true;;
		-no-vs    | --no-varscan) NOvs=true;;
		-no-vd    | --no-vardict) NOvd=true;;

		-*) commander::printerr "illegal option $1"; return 1;;
		*) commander::printerr "illegal option $2"; return 1;;
	esac

	$arg && {
		[[ ! $2 ]] && commander::printerr "argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	}
}

options::resume(){
	local s enable=false
	# don't Smd5, Sslice !
	for s in qual trim clip cor rrm sege star bwa uniq sort rg rmd cmo stats nsplit reo laln bqsr idx pon pondb hc mu bt fb pp vs vd; do
		eval "\${SKIP$s:=true}" # unless SKIP$s already set to false by -redo, do skip
		$enable || [[ "$1" == "$s" ]] && {
			enable=true
			eval "SKIP$s=false"
		}
	done
}

options::skip(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for x in "${mapdata[@]}"; do
		for s in md5 qual trim clip cor rrm sege star bwa uniq sort slice rg rmd cmo stats nsplit reo laln bqsr idx pon pondb hc mu bt fb pp vs vd; do
			[[ "$x" == "$s" ]] && eval "SKIP$s=true"
		done
	done
}

options::redo(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for s in qual trim clip cor rrm sege star bwa uniq sort rg rmd cmo stats nsplit reo laln bqsr idx pon pondb hc mu bt fb pp vs vd; do
		eval "\${SKIP$s:=true}" # unless SKIP$s alredy set to false by -resume, do skip
	done
	for x in "${mapdata[@]}"; do
		for s in qual trim clip cor rrm sege star bwa uniq sort rg rmd cmo stats nsplit reo laln bqsr idx pon pondb hc mu bt fb pp vs vd; do
			[[ "$x" == "$s" ]] && eval "SKIP$s=false"
		done
	done
}
