#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
shopt -s extglob
shopt -s expand_aliases

die() {
	unset CLEANUP
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

cleanup() {
	if [[ $CLEANUP ]]; then
		local b e
		for f in "${FASTQ1[@]}"; do
			helper::basename -f "$f" -o b -e e
			f=$b
			[[ -e $TMPDIR ]] && find $TMPDIR -type f -name "$f*" -exec rm -f {} \;
			if [[ -e $OUTDIR ]]; then
				find $OUTDIR -type f -name "$f*.all" -exec rm -f {} \;
				find $OUTDIR -type f -name "$f*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
				find $OUTDIR -type f -name "$f*.*.gz" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .gz)' \;
			fi
		done
	fi
}

[[ ! $MUVAC ]] && die "cannot find installation. please run setup and/or do: export MUVAC=/path/to/install/dir"
INSDIR=$MUVAC
source $INSDIR/latest/bashbone/activate.sh || die "install directory cannot be found"
unset JAVA_HOME #activate.sh loads java 12 into path which lets gatk4 fail with IncompatibleClassChangeError
for f in $INSDIR/latest/muvac/lib/*.sh; do
	source $f
done

CMD="$(basename $0) $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
MAXMEMORY=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
MEMORY=30000
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
VERBOSITY=0
OUTDIR=$PWD/results
TMPDIR=$OUTDIR
REGEX='\S+:(\d+):(\d+):(\d+)\s*.*'
DISTANCE=5

options::parse "$@" || die "parameterization issue"

mkdir -p $OUTDIR || die "cannot access $OUTDIR"
OUTDIR=$(readlink -e $OUTDIR)
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR=$PREVIOUSTMPDIR
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
else
	SKIPslice=false
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
	TMPDIR=$(mktemp -p $TMPDIR -d --suffix=.muvac) || die "cannot access $TMPDIR"
fi

[[ ! $LOG ]] && LOG=$OUTDIR/run.log
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
[[ $MTHREADS -eq 0 ]] && die "too less memory available ($MAXMEMORY)"
[[ ! $NFASTQ1 ]] && [[ ! $TFASTQ1 ]] && [[ ! $NMAPPED ]] && [[ ! $TMAPPED ]] && die "fastq or sam/bam file input missing - call "$(basename $0)" -h for help"
if [[ $GENOME ]]; then
	readlink -e $GENOME | file -f - | grep -qF ASCII || die "genome file does not exists or is compressed $GENOME"
else
	commander::warn "proceeding without genome file"
	SKIPmd5=true
	NOsege=true
	NOstar=true
	NObwa=true
fi
if [[ $DBSNP ]]; then
	readlink -e $DBSNP | file -f - || die "dbSNP file does not exists $DBSNP"
else
	[[ $TFASTQ1 ]] && commander::warn "proceeding without dbSNP file"
fi

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

commander::print "muvac v$version started with command: $CMD" > $LOG || die "cannot access $LOG"
commander::print "temporary files go to $HOST:$TMPDIR" >> $LOG
progress::log -v $VERBOSITY -o $LOG

if [[ $TFASTQ1 ]]; then
	pipeline::somatic >> $LOG 2> >(tee -a $LOG >&2) || die
else
	pipeline::germline >> $LOG 2> >(tee -a $LOG >&2) || die
fi

commander::print "success" >> $LOG
exit 0
