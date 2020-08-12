#! /usr/bin/env bash
# (c) Konstantin Riege

die() {
	echo ":ERROR: $*" >&2
	exit 1
}

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

# defines INSDIR and by sourcing bashbone it defines BASHBONEVERSION variable as well
source $(dirname $(readlink -e $0))/activate.sh -c true || die

trap 'configure::exit -p $$ -f cleanup $?' EXIT
trap 'die "killed"' INT TERM

VERSION=$version
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
	TMPDIR=$(mktemp -d -p $TMPDIR muvac.XXXXXXXXXX) || die "cannot access $TMPDIR"
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

commander::printinfo "muvac $VERSION utilizing bashbone $BASHBONEVERSION started with command: $CMD" > $LOG || die "cannot access $LOG"
commander::printinfo "temporary files go to $HOSTNAME:$TMPDIR" >> $LOG
progress::log -v $VERBOSITY -o $LOG

${Smd5:=false} || {
	[[ ! -s $GENOME.md5.sh ]] && cp $(dirname $(readlink -e $0))/bashbone/lib/md5.sh $GENOME.md5.sh
	source $GENOME.md5.sh
}
[[ ! $FASTQ2 ]] && NOcmo=true
if [[ $TFASTQ1 ]]; then
	pipeline::somatic 2> >(tee -ai $LOG >&2) >> $LOG || die
else
	pipeline::germline 2> >(tee -ai $LOG >&2) >> $LOG || die
fi
${Smd5:=false} || {
	commander::printinfo "finally updating genome and annotation md5 sums" >> $LOG
	thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
	[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	thismd5gtf=$(md5sum $GTF | cut -d ' ' -f 1)
	[[ "$md5gtf" != "$thismd5gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" $GENOME.md5.sh
}

commander::printinfo "success" >> $LOG
exit 0
