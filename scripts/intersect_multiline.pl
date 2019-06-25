#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;
use Data::Dumper qw(Dumper);
use List::Util qw(max);
use List::MoreUtils qw(firstidx any);

use Bio::DB::SeqFeature::Store;

if ($#ARGV < 2 || any {/^\s*-+(h|help|man)\s*$/} @ARGV){
    say "usage intersect.pl [regex-for-samplename] [annotation.gtf|.gff] [dbSNP.vcf] [vcf] [vcf] ...";
    say "example: intersect.pl 'sample\\d+' genome.(gff|gtf) dbsnp.vcf sample1.caller1.vcf sample1.caller2.vcf ... sample2.caller1.vcf sample2.caller2.vcf ...";
    exit;
}

my $r = shift @ARGV;
my $regex = qr($r);

my $exons = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
$exons->init_database([1]);
my $genes = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
$genes->init_database([1]);

my $gtfgff = shift @ARGV;
open F,"<".$gtfgff or die $!;
if ($gtfgff=~/\.gtf$/){
    while(<F>) {
        next if $_=~/^\s*(#+|$)/;
        chomp $_;
        my @l = split /\t/,$_;
        if ($l[2] eq "exon") {
            $l[-1]=~/gene_id\s+"(.+?)"/;
            my $geneid = $1;
            $l[-1]=~/gene_name\s+"(.+?)"/;
            my $genename = $1 ? $1 : '.';
            $exons->new_feature(
                -seq_id => $l[0],
                -start => $l[3],
                -stop => $l[4],
                -type => $geneid,
                -name => $genename,
                -index => 1
                );
        } elsif ($l[2] eq "gene"){
            $l[-1]=~/gene_id\s+"(.+?)"/;
            my $geneid = $1;
            $l[-1]=~/gene_name\s+"(.+?)"/;
            my $genename = $1 ? $1 : '.';
            $genes->new_feature(
                -seq_id => $l[0],
                -start => $l[3],
                -stop => $l[4],
                -type => $geneid,
                -name => $genename,
                -index => 1
                );
        }
        # last unless $l[0] eq 'chrM';
    }
} else {
    my %parentid;
    my %parentname;
    while(<F>) {
        next if $_=~/^\s*(#+|$)/;
        chomp $_;
        my @l = split /\t/,$_;
        if ($l[2] eq "mRNA"){
            $l[-1]=~/parent=([^;]+)/i;
            my $geneid = $1;
            $l[-1]=~/id=([^;]+)/i;
            $parentid{$1} = $geneid;
            $l[-1]=~/name=([^;]+)/i;
            my $genename = $1 ? $1 : '.';
            $parentname{$1} = $genename;
            $genes->new_feature(
                -seq_id => $l[0],
                -start => $l[3],
                -stop => $l[4],
                -type => $geneid,
                -name => $genename,
                -index => 1
            );
        }
	# last unless $l[0] eq 'chrM';
    }
    seek F, 0, 0;
    while(<F>) {
        next if $_=~/^\s*(#+|$)/;
        chomp $_;
        my @l = split /\t/,$_;
        if ($l[2] eq "exon") {
            $l[-1]=~/parent=([^;]+)/i;
            my $geneid = $parentid{$1};
            my $genename = $parentname{$1};
            $exons->new_feature(
                -seq_id => $l[0],
                -start => $l[3],
                -stop => $l[4],
                -type => $geneid,
                -name => $genename,
                -index => 1
            );
        }
    }
}
close F;   

my %dbsnp;
open F,"<".(shift @ARGV) or die $!;
while(<F>) {
	next if $_=~/^\s*(#+|$)/;
	chomp $_;
	my @l = split /\t/,$_;
	$dbsnp{$l[0]}{$l[1]} = 1;
	# last;
}
close F;

my %maprefalt;
my %mapgene;
my %mapgeneid;
my %mapexonintron;

my %mapcaller;
my %mapmaxmaf;
my %mapmaxcov;

my %mapsamplemaf;
my %mapsamplecov;
my %mapsamplecaller;

for my $file (@ARGV){
	$file=~/recalibrated\.(.+)\.splitpoly/;
	my $caller = $1;
	$file=~/($regex)/;
	my $sample = $1;
	my $mapsamplemaf = \%mapsamplemaf;
	my $mapsamplecov = \%mapsamplecov;
	my $mapsamplecaller = \%mapsamplecaller;
	my @t;
	my $mafi;
	my $covi;

	open F,"<$file" or die $!;
	while(<F>){
		next if $_=~/^\s*(#+|$)/;
		chomp $_;
		my @l = split /\s+/,$_;
		if ($#t == -1) {
			@t = split /:/,$l[-3];
			$mafi = firstidx { $_ eq "MAF" } @t;
			$covi = firstidx { $_ eq "COV" } @t;
		}
		# say $_;
		my @v = split /:/,$l[-1]; #tumor
		unless (exists $mapcaller{$l[0]}{$l[1]}) {
			my @f = $genes->features(-seq_id => $l[0], -start => $l[1], -stop => $l[1], -range_type => 'overlaps');
			if ($#f > -1) {
				for (@f) {
					push @{$mapgene{$l[0]}{$l[1]}} , $_->name;
					push @{$mapgeneid{$l[0]}{$l[1]}} , $_->type;
					my @e = $exons->features(-type => $_->type, -seq_id => $l[0], -start => $l[1], -stop => $l[1], -range_type => 'overlaps');
					push @{$mapexonintron{$l[0]}{$l[1]}} , $#e > -1 ? 'e' : 'i';
				}
			} else {
				push @{$mapexonintron{$l[0]}{$l[1]}} , '.';
				push @{$mapgene{$l[0]}{$l[1]}} , '.';
				push @{$mapgeneid{$l[0]}{$l[1]}} , '.';
			}

			$maprefalt{$l[0]}{$l[1]} = $l[3]."\t".$l[4];
			$mapmaxmaf{$l[0]}{$l[1]} = $v[$mafi];
			$mapmaxcov{$l[0]}{$l[1]} = $v[$covi];
		} else {
			$mapmaxmaf{$l[0]}{$l[1]} = max($mapmaxmaf{$l[0]}{$l[1]},$v[$mafi]);
			$mapmaxcov{$l[0]}{$l[1]} = max($mapmaxcov{$l[0]}{$l[1]},$v[$covi]);
		}
		$mapcaller{$l[0]}{$l[1]}{$caller} = 1;
		$mapsamplecaller->{$l[0]}{$l[1]}{$sample}{$caller} = 1;
		if (exists $mapsamplemaf->{$l[0]}{$l[1]}{$sample}){
			$mapsamplemaf->{$l[0]}{$l[1]}{$sample} = max($mapsamplemaf->{$l[0]}{$l[1]}{$sample},$v[$mafi]);
			$mapsamplecov->{$l[0]}{$l[1]}{$sample} = max($mapsamplecov->{$l[0]}{$l[1]}{$sample},$v[$covi]);
		} else {
			$mapsamplemaf->{$l[0]}{$l[1]}{$sample} = $v[$mafi];
			$mapsamplecov->{$l[0]}{$l[1]}{$sample} = $v[$covi];
		}
		# last;
	}
	close F;
}



say "#chr\tpos\tref\talt\tdbsnp\tgene\tgeneid\tintron/exon\tsample\tCOV\tMAF\tcaller\t#caller";
for my $chr (keys %mapcaller){
	for my $pos (keys %{$mapcaller{$chr}}){
		my @allcaller = sort {$a cmp $b} keys %{$mapcaller{$chr}{$pos}};
		my @samplecaller;
		my $samplecallermax = 0;
		my @samplemaf;
		my @samplecov;
		my @samplecovmafcaller;

		my @samples = sort {$a cmp $b} keys %{$mapsamplemaf{$chr}{$pos}};
		my $sc = 0;
		if ($#samples == -1){
			push @samples , '.';
		} else {
			$sc = $#samples + 1;
			for my $s (@samples){
				push @samplemaf, $mapsamplemaf{$chr}{$pos}{$s};
				push @samplecov, $mapsamplecov{$chr}{$pos}{$s};
				my @caller = sort {$a cmp $b} keys %{$mapsamplecaller{$chr}{$pos}{$s}};
				push @samplecaller, join(",",@caller);
				$samplecallermax = max($samplecallermax,$#caller+1);
				push @samplecovmafcaller, join(",",($samplecov[-1],$samplemaf[-1],$s,$samplecaller[-1]));
			}
		}
        
        for (0..$#samples) {
            say "$chr\t$pos\t".$maprefalt{$chr}{$pos}."\t".
                (exists $dbsnp{$chr}{$pos} ? "known" : "novel")."\t".
                join(',',@{$mapgene{$chr}{$pos}})."\t".
                join(',',@{$mapgeneid{$chr}{$pos}})."\t".
                join(',',@{$mapexonintron{$chr}{$pos}})."\t".
                $samples[$_]."\t".
                $samplecov[$_]."\t".
                $samplemaf[$_]."\t".
                $samplecaller[$_]."\t".
                scalar(split/,/,$samplecaller[$_])
        }
	}
}
