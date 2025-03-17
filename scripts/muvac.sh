#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(dirname "$(readlink -e "$0")")")/activate.sh" -l true -c false -r true -x cleanup -a "$@" || exit 1

cleanup() {
	[[ -e "$LOG" ]] && {
		echo "date: $(date)" | tee -ia "$LOG"
		[[ $1 -eq 0 ]] && echo "success" | tee -ia "$LOG" || echo "failed" | tee -ia "$LOG"
	}
	[[ -e "$CLEANUP_TMPDIR" ]] && {
		find -L "$CLEANUP_TMPDIR" -type f -name "cleanup.*" -exec rm -f "{}" \; &> /dev/null || true
		find -L "$CLEANUP_TMPDIR" -depth -type d -name "cleanup.*" -exec rm -rf "{}" \; &> /dev/null || true
	}
	${CLEANUP:=true} && {
		[[ -e "$CLEANUP_TMPDIR" ]] && {
			find -L "$CLEANUP_TMPDIR" -type f -exec rm -f "{}" \; &> /dev/null || true
			find -L "$CLEANUP_TMPDIR" -depth -type d -exec rm -rf "{}" \; &> /dev/null || true
			rm -rf "$CLEANUP_TMPDIR"
		}
		[[ -e "$OUTDIR" ]] && {
			local b f
			for f in "${FASTQ1[@]}"; do
				readlink -e "$f" | file -b --mime-type -f - | grep -qF -e 'gzip' -e 'bzip2' && b=$(basename "$f" | rev | cut -d '.' -f 3- | rev) || b=$(basename "$f" | rev | cut -d '.' -f 2- | rev)
				find -L "$OUTDIR" -depth -type d -name "$b*._STAR*" -exec rm -rf {} \; &> /dev/null || true
				# find -L "$OUTDIR" -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s "$1" ]] && rm -f "$(dirname "$1")/$(basename "$1" .sorted.bam).bam"' bash {} \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.*.gz" -exec bash -c '[[ -s "$1" ]] && rm -f "$(dirname "$1")/$(basename "$1" .gz)"' bash {} \; &> /dev/null || true
			done
			for f in "${MAPPED[@]}"; do
				b=$(basename "$f" | rev | cut -d '.' -f 2- | rev)
				find -L "$OUTDIR" -depth -type d -name "$b*._STAR*" -exec rm -rf "{}" \; &> /dev/null || true
				# find -L "$OUTDIR" -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s "$1" ]] && rm -f "$(dirname "$1")/$(basename "$1" .sorted.bam).bam"' bash {} \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.*.gz" -exec bash -c '[[ -s "$1" ]] && rm -f "$(dirname "$1")/$(basename "$1" .gz)"' bash {} \; &> /dev/null || true
			done
		}
	}
	return 0
}

VERSION=$version
CMD="$(basename "$0") $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
MAXMEMORY=$(grep -F MemTotal /proc/meminfo | awk '{printf("%d",$2/1024*0.95)}')
MEMORY=20000
MTHREADS=$((MAXMEMORY/MEMORY))
VERBOSITY=0
OUTDIR="$PWD/results"
TMPDIR="${TMPDIR:-$OUTDIR}"
DISTANCE=5

BASHBONE_ERROR="parameterization issue"
options::parse "$@"
bashbone -c

MTHREADS=$((MAXMEMORY/MEMORY))
[[ $MTHREADS -gt $THREADS ]] && MTHREADS=$THREADS
BASHBONE_ERROR="too less memory available ($MAXMEMORY)"
[[ $MTHREADS -ne 0 ]]
# check for varscan and vardict
BASHBONE_ERROR="too less memory per thread available ($((MAXMEMORY/THREADS))) < 3072Mb"
[[ $((THREADS*3072)) -lt $MAXMEMORY ]]

[[ ! $LOG ]] && LOG="$OUTDIR/run.log"
BASHBONE_ERROR="cannot access $LOG"
mkdir -p "$(dirname "$LOG")"
commander::printinfo "muvac $VERSION utilizing bashbone $BASHBONE_VERSION started with command: $CMD" | tee -i "$LOG"
LOG="$(realpath -se "$LOG")"

BASHBONE_ERROR="cannot access $OUTDIR"
mkdir -p "$OUTDIR"
OUTDIR="$(realpath -se "$OUTDIR")"

BASHBONE_ERROR="cannot access $TMPDIR"
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR="$(realpath -se "$PREVIOUSTMPDIR")"
else
	mkdir -p "$TMPDIR"
	TMPDIR="$(realpath -se "$TMPDIR")"
	TMPDIR="$(mktemp -d -p "$TMPDIR" muvac.XXXXXXXXXX)"
fi
CLEANUP_TMPDIR="$TMPDIR"

${INDEX:=false} || {
	BASHBONE_ERROR="fastq or sam/bam file input missing"
	[[ ! $NFASTQ1 ]] && [[ ! $TFASTQ1 ]] && [[ ! $NMAPPED ]] && [[ ! $TMAPPED ]] && false
}

[[ ! $NFASTQ2 && "$NOcmo" == "false" ]] && {
	commander::warn "second mate fastq file missing. proceeding without mate overlap clipping"
	NOcmo=true
}

[[ $GENOME ]] && {
	BASHBONE_ERROR="genome file does not exists or is compressed $GENOME"
	readlink -e "$GENOME" | file -b --mime-type -f - | grep -qF 'text'
	[[ ! -s "$GENOME.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$GENOME.md5.sh"
	source "$GENOME.md5.sh"
} || {
	BASHBONE_ERROR="genome file missing"
	! ${INDEX:=false}
	commander::warn "genome file missing. proceeding without mapping"
	NOsege=true
	NOstar=true
	NObwa=true
	Smd5=true
}

if [[ $DBSNP ]]; then
	BASHBONE_ERROR="dbSNP file does not exists $DBSNP"
	[[ -s "$DBSNP" ]]
else
	if [[ -e "$GENOME.dbSNP.vcf.gz" ]]; then
		DBSNP="$GENOME.dbSNP.vcf.gz"
		commander::warn "using dbSNP file $DBSNP"
	else
		commander::warn "dbSNP file missing. proceeding without dbSNP file"
		NOdbsnp=true
	fi
fi

if [[ $PONDB ]]; then
	BASHBONE_ERROR="pon file does not exists $PONDB"
	[[ -s "$PONDB" ]]
else
	if [[ -e "$GENOME.PON.vcf.gz" ]]; then
		PONDB="$GENOME.PON.vcf.gz"
		commander::warn "using pon file $DBSNP"
	else
		commander::warn "pon file missing. proceeding without pon file"
		NOpon=true
	fi
fi

declare -a FASTQ1 FASTQ2 FASTQ3 MAPPED NIDX TIDX
helper::addmemberfunctions -v FASTQ1 -v FASTQ2 -v FASTQ3 -v MAPPED -v NIDX -v TIDX
helper::addmemberfunctions -v NFASTQ1 -v NFASTQ2 -v NFASTQ3 -v NMAPPED
helper::addmemberfunctions -v TFASTQ1 -v TFASTQ2 -v TFASTQ3 -v TMAPPED

if [[ $NFASTQ1 ]]; then
	FASTQ1.join NFASTQ1 TFASTQ1
	FASTQ2.join NFASTQ2 TFASTQ2
	FASTQ3.join NFASTQ3 TFASTQ3
	NIDX.push $(NFASTQ1.idxs)
	TIDX.push $(seq $(NFASTQ1.length) $(($(NFASTQ1.length)+$(TFASTQ1.length)-1)))
else
	MAPPED.join NMAPPED TMAPPED
	NIDX.push $(NMAPPED.idxs)
	TIDX.push $(seq $(NMAPPED.length) $(($(NMAPPED.length)+$(TMAPPED.length)-1)))
fi

if [[ $FASTQ3 ]]; then
	NOrmd=false
fi

commander::printinfo "temporary files go to: $HOSTNAME:$TMPDIR" | tee -ia "$LOG"
commander::printinfo "date: $(date)" | tee -ia "$LOG"
x=$(ulimit -Hn)
[[ $((x/100)) -lt $THREADS ]] && {
	commander::warn "detected a low user limit of open file descriptors (ulimit -Hn : $x) for too many threads ($THREADS)"
	commander::warn "in case of memory allocation errors, you may decrease the number of threads to $((x/100))." | tee -ia "$LOG"
	commander::warn "possible memory allocation errors are 'bash: fork: Cannot allocate memory', 'Failed to read from standard input', 'Failed to open -', 'Too many open files'" | tee -ia "$LOG"
}

if ${INDEX:=false}; then
	BASHBONE_ERROR="indexing failed"
	progress::log -v $VERBOSITY -o "$LOG" -f pipeline::index
else
	if [[ $TFASTQ1 || $TMAPPED ]]; then
		BASHBONE_ERROR="somatic variant calling pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::somatic
	else
		BASHBONE_ERROR="germline variant calling pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::germline
	fi
fi
unset BASHBONE_ERROR

exit 0
