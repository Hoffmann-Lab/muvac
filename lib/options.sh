#! /usr/bin/env bash
# (c) Konstantin Riege

function options::usage(){
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		MUVAC is a multiple, universal germline and somatic variant caller pipeline.
		It implements GATK best practices and other variant callers in an optimized, parallelized fashion.


		VERSION
		$VERSION
		utilizing bashbone $BASHBONE_VERSION


		BASIC OPTIONS
		-h       | --help                     : prints this message
		-dev     | --devel                    : prints list of keywords in processing order for advanced pipeline control
		-e       | --env                      : list tools and versions in setupped environment
		-v       | --verbosity [value]        : set level of verbosity. default: 0
		                                        0 - get simple status updates
		                                        1 - get status updates and commands
		                                        2 - get full output
		-o       | --out [path]               : output directory. default: $OUTDIR
		-l       | --log [path]               : output directory. default: $OUTDIR/run.log
		-tmp     | --tmp                      : temporary directory. default: ${TMPDIR:-/tmp}/muvac.XXXXXXXXXX
		                                        NOTE: respects TMPDIR environment variable
		-r       | --remove                   : remove temporary and unnecessary files upon successful termination
		-rr      | --remove-remove            : remove temporary and unnecessary files upon termination
		-t       | --threads [value]          : number of threads. default: $THREADS
		-xmem    | --max-memory [value]       : fraction or total amount of allocatable memory in MB. default: $MAXMEMORY MB i.e. currently available memory
		-mem     | --memory [value]           : allocatable memory per instance of memory greedy tools in MB. defines internal number of parallel instances
		                                        default: $MEMORY which allows for $MTHREADS instances and $MTHREADS SAM/BAM slices according to -xmem
		                                        NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors
		-resume  | --resume-from [string]     : resume from a specific pipeline step (see -dev)
		                                        NOTE: define entry point before skip and redo
		-skip    | --skip [string,..]         : skip specific pipeline step(s). comma separated (see -dev)
		-redo    | --redo [string,..]         : redo specific pipeline step(s). comma separated (see -dev)


		INDEXING OPTIONS
		-x       | --index                    : triggers creation of all requiered genome and annotation indices plus md5 sums
		-g       | --genome [path]            : genome fasta input
		-no-sege | --no-segemehl              : disables indexing for segemehl
		-no-star | --no-star                  : disables indexing for STAR
		-no-bwa  | --no-bwa                   : disables indexing for BWA


		GERMLINE OPTIONS
		-do-pon  | --do-panelofnormals        : triggers panel of normals variant calling
		                                        NOTE: implies -no-call except for Mutect2
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
		                                        NOTE: no fasta file implies -no-map
		-s       | --snp [path]               : genome matching dbSNP, compressed and tabix indexed
		-1       | --fq1 [path,..]            : fastq input. single or first mate. comma separated
		-2       | --fq2 [path,..]            : fastq input. mate pair. comma separated
		-3       | --fq3 [path,..]            : fastq input. UMI sequences. comma separated
		-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
		                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' quality trimming
		-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
		                                      : NOTE: clipping also includes simple 3' quality trimming
		-pclip   | --polyntclipping           : enables removal of trailing mono- and di-nucleotide sequences i.e. poly-(A|C|G|T) and poly-(AC|AG..|GT)
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
		-cor     | --correction               : enables majority based raw read error correction
		-rrm     | --rrnafilter               : enables rRNA filter
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-d       | --distance                 : maximum read alignment edit distance in %. default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 200000
		-split   | --split                    : enables split read mapping and afterwards split alignments by N-cigar strings to call variants from RNA-Seq data
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA
		-m       | --mapped [path,..]         : SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-no-addrg| --no-addreadgroup          : disables proper read group modification by Picard
		-no-rmd  | --no-removeduplicates      : disables removing duplicates - not recommended unless reads were mapped on a transcriptome
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+).*
		                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
		-no-cmo  | --no-clipmateoverlaps      : disables clipping of read mate overlaps
		-reo     | --reorder                  : enables reordering according to genome file by Picard
		-rgn     | --readgroup-name [string]  : sets custom read group name - use TUMOR or NORMAL for subsequent somatic calls. default: SAMPLE
		-no-laln | --no-leftalign             : disables left alignment by GATK
		-no-bqsr | --no-qualrecalibration     : disables any base quality score recalibration (BQSR)
		-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
		-no-pondb| --no-pondatabase           : disables creation of a panel of normals database from panel of normals variant calling
		-no-call | --no-call                  : disables variant calling and downstream analyses
		-no-gatk | --no-gatk                  : disables variant calling by HaplotypeCaller (Mutect2 in case of -do-pon)
		-no-bt   | --no-bcftools              : disables variant calling by BCFtools
		-no-fb   | --no-freebayes             : disables variant calling by freebayes
		-no-vs   | --no-varscan               : disables variant calling by VarScan
		-no-vd   | --no-vardict               : disables variant calling by VarDict
		-no-pp   | --no-platypus              : disables variant calling by Platypus


		SOMATIC OPTIONS
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
		                                        NOTE: no fasta file implies -no-map
		-s       | --snp [path]               : genome matching dbSNP, compressed and tabix indexed
		-p       | --pon [path]               : genome panel of normals input
		-n1      | --normalfq1 [path,..]      : normal fastq input. single or first mate. comma separated
		-n2      | --normalfq2 [path,..]      : normal fastq input. mate pair. comma separated
		-n3      | --normalfq3 [path,..]      : normal fastq input. UMI sequences. comma separated
		-t1      | --tumorfq1 [path,..]       : tumor fastq input. single or first mate. comma separated
		-t2      | --tumorfq2 [path,..]       : tumor fastq input. mate pair. comma separated
		-t3      | --tumorfq3 [path,..]       : tumor fastq input. UMI sequences. comma separated
		-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
		                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' quality trimming
		-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
		                                      : NOTE: clipping also includes simple 3' quality trimming
		-pclip   | --polyntclipping           : enables removal of trailing mono- and di-nucleotide sequences i.e. poly-(A|C|G|T) and poly-(AC|AG..|GT)
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
		-cor     | --correction               : enables majority based raw read error correction
		-rrm     | --rrnafilter               : enables rRNA filter
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-d       | --distance                 : maximum read alignment edit distance in %. default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 200000
		-split   | --split                    : enables split read mapping and afterwards split alignments by N-cigar strings to call variants from RNA-Seq data
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA
		-nm      | --normalmapped [path,..]   : normal SAM/BAM input. comma separated (replaces fastq input)
		-tm      | --tumormapped [path,..]    : tumor SAM/BAM input. comma separated (replaces fastq input)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-no-addrg| --no-addreadgroup          : disables proper read group modification by Picard
		-no-rmd  | --no-removeduplicates      : disables removing duplicates - not recommended unless reads were mapped on a transcriptome
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+).*
		                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
		-no-cmo  | --no-clipmateoverlaps      : disables clipping of read mate overlaps
		-reo     | --reorder                  : enables reordering according to genome file by Picard
		-rgn     | --readgroup-name [string]  : sets custom read group name - use TUMOR or NORMAL for subsequent somatic calls. default: SAMPLE
		-no-laln | --no-leftalign             : disables left alignment by GATK
		-no-bqsr | --no-qualrecalibration     : disables any base quality score recalibration (BQSR)
		-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
		-no-call | --no-call                  : disables variant calling and downstream analyses
		-no-gatk | --no-gatk                  : disables variant calling by HaplotypeCaller/Mutect2
		-no-bt   | --no-bcftools              : disables variant calling by BCFtools
		-no-fb   | --no-freebayes             : disables variant calling by freebayes
		-no-vs   | --no-varscan               : disables variant calling by VarScan
		-no-vd   | --no-vardict               : disables variant calling by VarDict
		-no-pp   | --no-platypus              : disables variant calling by Platypus


		ADDITIONAL INFORMATION
		Chromosome order for all input files (genome, annotation, dbsnp, pon, gnomad, ..) must be identical.
	    If possible, pleas provide them in karyotypic order and following naming schema: chrM,chr1,chr2,..,chrX,chrY

		To obtain comprehensive panel of normals, small common and full gnomAD/Exac population variants with allele frequencies visit
		https://gatk.broadinstitute.org/hc/en-us/articles/360035890811
		-> for hg38: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/
		-> for hg19: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37/
		After download, place the gnomAD/Exac vcf and tbi index files next to your genome fasta file.
	    Rename gnomad/exac according to your genome fasta file and add the following suffixes
		<genome.fa>.small_common.vcf.gz
	    <genome.fa>.small_common.vcf.gz.tbi
		<genome.fa>.af_only_gnomad.vcf.gz
	    <genome.fa>.af_only_gnomad.vcf.gz.tbi

		To obtain dbSNP common variants vcf and tbi index along with a matching genome and annotation use
	    either the supplied dlgenome.sh script
	    or visit
		https://ftp.ncbi.nlm.nih.gov/snp/organisms
		-> for hg38: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF
		or download 1000genomes common variants, dbSNP by extraction from chromosomal vcf files respectively via
	    http://ftp.ensembl.org/pub/current_variation/vcf
		-> for hg38: ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens
		-> for hg19: ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens


		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

function options::developer(){
	cat <<- EOF
		DESCRIPTION
		In case of restarting or to resume an analysis use the identifiers below, listed in processing order

		DEVELOPER OPTIONS
		md5   : check for md5sums and if necessary trigger genome indexing
		fqual : input quality metrics
		trim  : trimming
		clip  : adapter clipping (& simple trimming)
		pclip : poly- mono-and di-nucleotide clipping
		cor   : raw read correction
		rrm   : rRNA filtering
		sege  : Segemehl mapping
		star  : STAR mapping
		mqual : mapping/input quality metrics
		uniq  : extraction of properly paired and uniquely mapped reads
		sort  : sorting and indexing of sam/bam files
		slice : better dont touch! slicing of bams for parallelization, needs -prevtmp | --previoustmp [path to rippchen.XXXXXXXXXX]
		rg    : read group modification
		rmd   : removing duplicates
		cmo   : clipping mate overlaps
		stats : proprocessing and mapping statistics
		reo   : bam reordering according to genome
		nsplit: splitting split-read alignments
		laln  : left alignment
		bqsr  : BQSRecalibration
		idx   : intermediate and final bam indexing
		pon   : panel of normals
		pondb : panel of normals database
		gatk  : haplotypecaller/mutect
		bt    : bcftools
		fb    : freebayes
		vs    : varscan
		vd    : vardict
		pp    : platypus

	EOF
	exit 0
}

function options::checkopt(){
	local arg=false skipredo=false
	declare -a mapdata

	case $1 in
		-h        | --help) options::usage || exit 0;;
		-e        | --env) bashbone -e; exit 0;;
		-dev      | --devel) options::developer;;
		-prevtmp  | --previoustmp) arg=true; PREVIOUSTMPDIR="$2";;
		-resume   | --resume-from) $skipredo && commander::printerr "define entry point via $1 before skip and redo" && return 1; arg=true; options::resume "$2";;
		-skip     | --skip) skipredo=true; arg=true; options::skip "$2";;
		-redo     | --redo) skipredo=true; arg=true; options::redo "$2";;

		-tmp      | --tmp) arg=true; TMPDIR="$2";;
		-r        | --remove) CLEANUP=true;;
		-rr       | --remove-remove) FORCECLEANUP=true;;
		-v        | --verbosity) arg=true; VERBOSITY=$2;;
		-t        | --threads) arg=true; THREADS=$2;;
		-mem      | --memory) arg=true; MEMORY=$2;;
		-xmem     | --max-memory) arg=true; [[ ${2%.*} -ge 1 ]] && MAXMEMORY=${2%.*} || MAXMEMORY=$(grep -F MemTotal /proc/meminfo | awk -v i=$2 '{printf("%d",$2/1024*0.95*i)}');;

		-x        | --index) INDEX=true;;
		-g        | --genome) arg=true; GENOME="$2";;
		-o        | --out) arg=true; OUTDIR="$2";;
		-l        | --log) arg=true; LOG="$2";;

		-s        | --snp) arg=true; DBSNP="$2";;
		-p        | --pon) arg=true; PONDB="$2";;
		-do-pon   | --do-panelofnormals) NOpon=false; NOfb=true; NObt=true; NOpp=true; NOvs=true; NOvd=true;;
		-no-pondb | --no-pondatabase) NOpondb=true;;

		-1        | --fq1 | -n1 | --normalfq1) arg=true; mapfile -t -d ',' NFASTQ1 < <(printf '%s' "$2");;
		-2        | --fq2 | -n2 | --normalfq2) arg=true; mapfile -t -d ',' NFASTQ2 < <(printf '%s' "$2");;
		-3        | --fq3 | -n3 | --normalfq3) arg=true; mapfile -t -d ',' NFASTQ3 < <(printf '%s' "$2");;
		-t1       | --tumorfq1) arg=true; mapfile -t -d ',' TFASTQ1 < <(printf '%s' "$2");;
		-t2       | --tumorfq2) arg=true; mapfile -t -d ',' TFASTQ2 < <(printf '%s' "$2");;
		-t3       | --tumorfq3) arg=true; mapfile -t -d ',' TFASTQ3 < <(printf '%s' "$2");;

		-a1       | --adapter1) arg=true; mapfile -t -d ',' ADAPTER1 < <(printf '%s' "$2");;
		-a2       | --adapter2) arg=true; mapfile -t -d ',' ADAPTER2 < <(printf '%s' "$2");;
		-d        | --distance) arg=true; DISTANCE=$2;;
		-i        | --insertsize) arg=true; INSERTSIZE=$2;;
		-split    | --split) NOsplitreads=false; NOnsplit=false;;
		-no-qual  | --no-qualityanalysis) NOqual=true;;
		-no-trim  | --no-trimming) NOtrim=true;;
		-no-clip  | --no-clipping) NOclip=true;;
		-pclip    | --polyntclipping) nopclip=false;;
		-cor      | --correction) NOcor=false;;
		-rrm      | --rrnafilter) NOrrm=false;;
		-no-sege  | --no-segemehl) NOsege=true;;
		-no-star  | --no-star) NOstar=true;;
		-no-bwa   | --no-bwa) NObwa=true;;
		-no-map   | --no-mapping) NOsege=true; NOstar=true; NObwa=true;;

		-m        | --mapped | -nm | --normalmapped) arg=true; mapfile -t -d ',' NMAPPED < <(printf '%s' "$2");;
		-mn       | --mapper-name) arg=true; MAPNAME=$2;;
		-tm       | --tumormapped) arg=true; mapfile -t -d ',' TMAPPED < <(printf '%s' "$2");;
		-rx       | --regex) arg=true; REGEX="$2";;

		-no-uniq  | --no-uniqify) NOuniq=true;;
		-no-sort  | --no-sort) NOsort=true;;
		-no-rmd   | --no-removeduplicates) NOrmd=true;;
		-no-cmo   | --no-clipmateoverlaps) NOcmo=true;;
		-rgn      | --readgroup-name) arg=true; RGPREFIX="$2";;
		-no-addrg | --no-addreadgroup) NOaddrg=true;;
		-reo      | --reorder) NOreo=false;;
		-no-laln  | --no-leftalign) NOlaln=true;;
		-no-bqsr  | --no-qualrecalibration) NObqsr=true;;
		-no-idx   | --no-index) NOidx=true;;
		-no-stats | --no-statistics) NOstats=true;;

		-no-call  | --no-call) NOgatk=true; NOfb=true; NObt=true; NOpp=true; NOvs=true; NOvd=true;;
		-no-gatk  | --no-gatk) NOgatk=true;;
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

function options::resume(){
	local s enable=false
	# don't Smd5, Sslice !
	for s in fqual trim clip pclip cor rrm sege star bwa mqual uniq sort addrg rmd cmo stats reo nsplit laln bqsr idx pon pondb gatk bt fb vs vd pp; do
		eval "\${SKIP$s:=true}" # unless SKIP$s already set to false by -redo, do skip
		$enable || [[ "$1" == "$s" ]] && {
			enable=true
			eval "SKIP$s=false"
		}
	done
}

function options::skip(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for x in "${mapdata[@]}"; do
		for s in md5 fqual trim clip pclip cor rrm sege star bwa mqual uniq sort slice addrg rmd cmo stats reo nsplit laln bqsr idx pon pondb gatk bt fb vs vd pp; do
			[[ "$x" == "$s" ]] && eval "SKIP$s=true"
		done
	done
}

function options::redo(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for s in fqual trim clip cor pclip rrm sege star bwa mqual uniq sort addrg rmd cmo stats reo nsplit laln bqsr idx pon pondb gatk bt fb vs vd pp; do
		eval "\${SKIP$s:=true}" # unless SKIP$s already set to false by -resume, do skip
	done
	for x in "${mapdata[@]}"; do
		for s in fqual trim clip cor pclip rrm sege star bwa mqual uniq sort addrg rmd cmo stats reo nsplit laln bqsr idx pon pondb gatk bt fb vs vd pp; do
			[[ "$x" == "$s" ]] && eval "SKIP$s=false"
		done
	done
}
