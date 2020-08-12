#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::muvac -i "$insdir" -t $threads && \
		compile::conda -i "$insdir" -t $threads && \
		compile::trimmomatic -i "$insdir" -t $threads && \
		compile::sortmerna -i "$insdir" -t $threads && \
		compile::segemehl -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::muvac() {
	local insdir threads src=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
	compile::_parse -r insdir -s threads "$@" || return 1

	local version bashboneversion
	source $src/bashbone/lib/version.sh
	bashboneversion=$version
	source $src/lib/version.sh
	shopt -s extglob

	commander::printinfo "installing muvac"
	{	rm -rf $insdir/muvac-$version && \
		mkdir -p $insdir/muvac-$version && \
		cp -r $src/!(bashbone|setup*) $insdir/muvac-$version && \
		mkdir -p $insdir/latest && \
		ln -sfn $insdir/muvac-$version $insdir/latest/muvac && \
		ln -sfn $insdir/bashbone-$bashboneversion $insdir/muvac-$version/bashbone
	} || return 1

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::muvac -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::conda() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing conda and tools"
	{	url='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && \
		wget -q -O $insdir/miniconda.sh $url && \
		version=$(bash $insdir/miniconda.sh -h | grep -F Installs | cut -d ' ' -f 3) && \
		rm -rf $insdir/conda && \
		mkdir -p $insdir/conda && \
		bash $insdir/miniconda.sh -b -f -p $insdir/conda && \
		rm $insdir/miniconda.sh && \
		source $insdir/conda/bin/activate && \

		conda env remove -y -n py2 && \
		conda env remove -y -n py3 && \
		conda create -y -n py2 python=2 && \
		conda create -y -n py2r python=2 && \
		conda create -y -n py3 python=3 && \
		
		# macs2, tophat2/hisat2 and R stuff needs python2 whereas cutadapt,idr,rseqc need python3 env
		
		conda install -n py2 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			perl perl-threaded perl-db-file perl-dbi perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-set-intervaltree perl-uri \
			numpy scipy pysam cython matplotlib \
			datamash \
			fastqc rcorrector \
			star bwa hisat2 \
			samtools picard bamutil bedtools vcflib vt \
			bcftools gatk4 freebayes varscan platypus-variant vardict vardict-java \
			snpeff snpsift && \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda list -n py2 -f "fastqc|rcorrector|star|star-fusion|bwa|hisat2|macs2|samtols|picard" | grep -v '^#' > $insdir/condatools.txt && \

		conda install -n py3 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			cutadapt rseqc && \
		conda list -n py3 -f "cutadapt|rseqc" | grep -v '^#' >> $insdir/condatools.txt && \

		conda install -n py2r -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			r-devtools bioconductor-biocinstaller bioconductor-biocparallel \
			bioconductor-genomicfeatures bioconductor-genefilter \
			r-dplyr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-ggpubr r-treemap r-rngtools && \
		conda list -n py2r -f "r-treemap" | grep -v '^#' >> $insdir/condatools.txt && \

		conda clean -y -a
	} || return 1

	return 0
}

