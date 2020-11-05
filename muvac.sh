#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(readlink -e "$0")")/activate.sh" -c true -x cleanup || exit 1

cleanup() {
	[[ -e $TMPDIR ]] && {
		find $TMPDIR -type f -name "cleanup.*" -exec rm -f {} \;
		find $TMPDIR -depth -type d -name "cleanup.*" -exec rm -rf {} \;
	}
	[[ $1 -eq 0 ]] && ${CLEANUP:=false} && {
		[[ -e $TMPDIR ]] && {
			find $TMPDIR -type f -exec rm -f {} \;
			find $TMPDIR -depth -type d -exec rm -rf {} \;
			rm -rf $TMPDIR
		}
		[[ -e $OUTDIR ]] && {
			for f in "${FASTQ1[@]}"; do
				readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)' && b=$(basename $f | rev | cut -d '.' -f 3- | rev) || b=$(basename $f | rev | cut -d '.' -f 2- | rev)
				find $OUTDIR -depth -type d -name "$b*._STAR*" -exec rm -rf {} \;
				find $OUTDIR -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
				find $OUTDIR -type f -name "$b*.*.gz" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .gz)' \;
			done
			local b
			for f in "${MAPPED[@]}"; do
				b=$(basename $f | rev | cut -d '.' -f 2- | rev)
				find $OUTDIR -depth -type d -name "$b*._STAR*" -exec rm -rf {} \;
				find $OUTDIR -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
			done
		}
	}
}

VERSION=$version
CMD="$(basename $0) $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
MAXMEMORY=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
MEMORY=30000
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
[[ $MTHREADS -eq 0 ]] && die "too less memory available ($MAXMEMORY)"
VERBOSITY=0
OUTDIR=$PWD/results
TMPDIR=$OUTDIR
REGEX='\S+:(\d+):(\d+):(\d+)\s*.*'
DISTANCE=5

BASHBONE_ERROR="parameterization issue"
options::parse "$@"

BASHBONE_ERROR="cannot access $OUTDIR"
mkdir -p $OUTDIR
OUTDIR=$(readlink -e $OUTDIR)
[[ ! $LOG ]] && LOG=$OUTDIR/run.log
BASHBONE_ERROR="cannot access $LOG"
mkdir -p "$(dirname "$LOG")"

BASHBONE_ERROR="cannot access $TMPDIR"
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR=$PREVIOUSTMPDIR
	mkdir -p $TMPDIR
	TMPDIR=$(readlink -e $TMPDIR)
else
	mkdir -p $TMPDIR
	TMPDIR=$(readlink -e $TMPDIR)
	TMPDIR=$(mktemp -d -p $TMPDIR muvac.XXXXXXXXXX)
fi

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
	readlink -e $GENOME | file -f - | grep -qF ASCII
	[[ ! -s $GENOME.md5.sh ]] && cp $(dirname $(readlink -e $0))/bashbone/lib/md5.sh $GENOME.md5.sh
	source $GENOME.md5.sh
} || {
	BASHBONE_ERROR="genome file missing"
	! ${INDEX:=false}
	commander::warn "genome file missing. proceeding without mapping"
	SKIPmd5=true
	NOsege=true
	NOstar=true
}

if [[ $GTF ]]; then
	BASHBONE_ERROR="annotation file does not exists or is compressed $GTF"
	readlink -e $GTF | file -f - | grep -qF ASCII
else
	readlink -e $GENOME.gtf | file -f - | grep -qF ASCII && {
		GTF=$GENOME.gtf
	} || {
		BASHBONE_ERROR="annotation file missing"
		! ${INDEX:=false}
		commander::warn "proceeding without gtf file"
	}
fi

[[ $DBSNP ]] && {
	BASHBONE_ERROR="dbSNP file does not exists $DBSNP"
	readlink -e $DBSNP &> /dev/null
} || {
	[[ $TFASTQ1 ]] && commander::warn "proceeding without dbSNP file"
}

declare -a FASTQ1 FASTQ2 MAPPED NIDX TIDX
helper::addmemberfunctions -v FASTQ1 -v FASTQ2 -v MAPPED -v NIDX -v TIDX
helper::addmemberfunctions -v NFASTQ1 -v NFASTQ2 -v NMAPPED
helper::addmemberfunctions -v TFASTQ1 -v TFASTQ2 -v TMAPPED

if [[ $NFASTQ1 ]]; then
	FASTQ1.join NFASTQ1 TFASTQ1
	FASTQ2.join NFASTQ2 TFASTQ2
	NIDX.push $(NFASTQ1.idxs)
	TIDX.push $(seq $(NFASTQ1.length) $(($(NFASTQ1.length)+$(TFASTQ1.length)-1)))
else
	MAPPED.join NMAPPED TMAPPED
	NIDX.push $(NMAPPED.idxs)
	TIDX.push $(seq $(NMAPPED.length) $(($(NMAPPED.length)+$(TMAPPED.length)-1)))
fi
[[ ! $FASTQ2 ]] && NOcmo=true

progress::log -v $VERBOSITY -o $LOG
commander::printinfo "muvac $VERSION utilizing bashbone $BASHBONE_VERSION started with command: $CMD" | tee -ai "$LOG"
commander::printinfo "temporary files go to: $HOSTNAME:$TMPDIR" | tee -ia "$LOG"

if ${INDEX:=false}; then
	BASHBONE_ERROR="indexing failed"
	progress::observe -v $VERBOSITY -o "$LOG" -f pipeline::index
else
	if [[ $TFASTQ1 ]]; then
		BASHBONE_ERROR="somatic variant calling pipeline failed"
		progress::observe -v $VERBOSITY -o "$LOG" -f pipeline::somatic
	else
		BASHBONE_ERROR="germline variant calling pipeline failed"
		progress::observe -v $VERBOSITY -o "$LOG" -f pipeline::germline
	fi
fi
unset BASHBONE_ERROR

${Smd5:=false} || {
	commander::printinfo "finally updating genome and annotation md5 sums" >> $LOG
	thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
	[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	thismd5gtf=$(md5sum $GTF | cut -d ' ' -f 1)
	[[ "$md5gtf" != "$thismd5gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" $GENOME.md5.sh
}

commander::printinfo "success" >> $LOG
exit 0
