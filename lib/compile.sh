#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::muvac -i "$insdir" -t $threads && \
		compile::bashbone -i "$insdir" -t $threads && \
		compile::conda -i "$insdir" -t $threads && \
		compile::java -i "$insdir" -t $threads && \
		compile::perlmodules -i "$insdir" -t $threads && \
		compile::sortmerna -i "$insdir" -t $threads && \
		compile::segemehl -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads
		# compile::annovar -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::muvac() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	local version
	source $SRC/lib/version.sh
	shopt -s extglob

	commander::print "installing muvac"
	{	mkdir -p $insdir/muvac-$version && \
		cp -r $SRC/!(bashbone|setup*) $insdir/muvac-$version && \
		mkdir -p $insdir/latest && \
		ln -sfn $insdir/muvac-$version $insdir/latest/muvac
	} || return 1

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	compile::muvac -i "$insdir" -t $threads || return 1
	return 0
}

compile::conda() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	shift $# # necessary for conda activate

	commander::print "installing conda and tools"
	{	mkdir -p $insdir/conda && \
		url='https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh' && \
		wget -q -O $insdir/miniconda.sh $url && \
		bash $insdir/miniconda.sh -b -f -p $insdir/conda && \
		rm $insdir/miniconda.sh && \
		source $insdir/conda/bin/activate && \
		conda create -y -n py2 python=2 && \
		conda create -y -n py2r python=2 && \
		conda create -y -n py3 python=3 && \
		
		# readline 7 causes library version number to be lower on the shared object warnings
		# use quantstack gcc for r and perl module installation - defaults gcc has a weird usage
		# under perl > 5.22 List::MoreUtils installation fails
		# macs2, tophat2/hisat2 and R stuff needs python2 whereas cutadapt,idr,rseqc need python3 env
		conda activate py2 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript htslib readline=6 perl=5.22 perl-threaded=5.22 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			fastqc trimmomatic rcorrector \
			star bwa hisat2 \
			samtools picard bedtools vcflib \
			bcftools gatk4 freebayes varscan platypus-variant vardict vardict-java \
			snpeff snpsift && \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		pip install numpy scipy pysam cython && \
		conda clean -y -a && \

		conda activate py3 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript readline=6 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			cutadapt rseqc && \
		conda clean -y -a && \

		conda activate py2r && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript readline=6 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			r-devtools bioconductor-biocinstaller bioconductor-biocparallel \
			bioconductor-genomicfeatures bioconductor-genefilter && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			r-dplyr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-ggpubr r-treemap r-rngtools && \
		conda clean -y -a && \

		conda deactivate
	} || return 1

	return 0
}

compile::annovar() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::print "installing annovar"
	{	url="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz" && \
		wget -q $url -O $insdir/annovar.tar.gz && \
		tar -xzf $insdir/annovar.tar.gz -C $insdir && \
		rm $insdir/annovar.tar.gz && \
		cd $insdir/annovar && \
		url='http://www.openbioinformatics.org/annovar/download/table_annovar.pl' && \
		wget -q $url -O table_annovar.pl && \
		chmod 755 table_annovar.pl && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/annovar
	} || return 1

	return 0
}

compile::_setup_annovar() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::print "configuring annovar databases"
	{	source $insdir/conda/bin/activate py2 && \
		cd -P $insdir/latest/annovar && \
		#refSeq
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/ && \
		# ./annotate_variation.pl -buildver hg19 -downdb ensGene humandb/ #UCSC - needs seq to get mRNA
		# ./annotate_variation.pl -buildver hg19 -downdb knownGene humandb/ #UCSC - needs seq to get mRNA
		# ./annotate_variation.pl -buildver hg19 -downdb seq humandb/hg19_seq
		# ./retrieve_seq_from_fasta.pl humandb/hg19_ensGene.txt -seqdir humandb/hg19_seq -format ensGene -outfile humandb/hg19_ensGeneMrna.fa
		# ./retrieve_seq_from_fasta.pl humandb/hg19_knownGene.txt -seqdir humandb/hg19_seq -format knownGene -outfile humandb/hg19_knownGeneMrna.fa
		url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensemblToGeneName.txt.gz' && \
		wget -q $url -O humandb/ensemblToGeneName.txt.gz && \
		gzip -d humandb/ensemblToGeneName.txt.gz  && \
		# ./annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
		# ./annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/ 
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp33a humandb/ && \
		# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp142 humandb
		./annotate_variation.pl -buildver hg19 -downdb tfbsConsSites humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb targetScanS humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb wgRna humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb gwasCatalog humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20170130 humandb/
	} || return 1

	return 0
}

compile::_setup_snpeff() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::print "configuring snpeff databases"
	{	source $insdir/conda/bin/activate py2 && \
		java -jar snpEff.jar download -v GRCh37.75 && \
		#java -jar snpEff.jar download -v hg19 #hg19: UCSC, hg19kg: UCSC knownGenes, GRCh37.75: Ensembl 
		url='http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/promoter_predictions/master_known.bed' && \
		wget -q $url -O data/promoter.bed && \
		url='http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/mirna_tss/miRNA_promoters_hg19_edited_data.bed' && \
		wget -q $url -O data/miRNApromoter.bed && \
		url='ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz' && \
		wget -q $url -O data/clinvar.vcf.gz && \
		url='ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi' && \
		wget -q $url -O data/clinvar.vcf.gz.tbi && \
		url='ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz' && \
		wget -q $url -O data/exac.vcf.gz && \
		url='ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi' && \
		wget -q $url -O data/exac.vcf.gz.tbi && \
		# url='https://drive.google.com/open?id=0B60wROKy6OqceTNZRkZnaERWREk'
		#(see https://sites.google.com/site/jpopgen/dbNSFP)
		url='ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.9.3.zip' && \
		wget -q $url -O data/dbnsfp.zip && \
		unzip -q -o -d data dbnsfp.zip && \
		rm data/dbnsfp.zip && \
		head -n 1 data/dbNSFP*_variant.chr1 > data/dbNSFP.txt && \
		cat dbNSFP*_variant.chr* | grep -v "^#" >> data/dbNSFP.txt && \
		rm dbNSFP*_variant.chr* && \
		$MUVAC/bin/samtools/bgzip -f -@ $THREADS < data/dbNSFP.txt > data/dbNSFP.txt.gz && \
		$MUVAC/bin/samtools/tabix -f -s 1 -b 2 -e 2 data/dbNSFP.txt.gz && \
		# url='http://www.genome.gov/admin/gwascatalog.txt'
		url='ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv' && \
		wget -q $url -O data/gwas.txt
	} || return 1

	return 0
}