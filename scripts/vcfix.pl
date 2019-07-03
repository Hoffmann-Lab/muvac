#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use List::MoreUtils qw(indexes);
use List::Util qw(min max sum);
use Getopt::Long;
die "option -i|--in is missing" if $#ARGV == -1;

my %format = (
    GQ => '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Phred-scaled genotype quality">',
    PL => '##FORMAT=<ID=PL,Number=.,Type=Float,Description="List of Phred-scaled genotype likelihoods">', #Number=G dont work for bcftools
    GL => '##FORMAT=<ID=GL,Number=.,Type=Float,Description="List of genotype likelihoods">', #Number=G dont work for bcftools
    MAF => '##FORMAT=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequence">',
    COV => '##FORMAT=<ID=COV,Number=1,Type=Integer,Description="Position read coverage">',
    AD => '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
    DP => '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth at the locus">',
    DP4 => '##FORMAT=<ID=DP4,Number=1,Type=Integer,Description="Strand specific ref and alt read counts: ref-fwd, ref-rev, alt-fwd, alt-rev">',
    ASF => '##FORMAT=<ID=ASF,Number=1,Type=Integer,Description="Alt strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
    RSF => '##FORMAT=<ID=ASF,Number=1,Type=Integer,Description="Ref strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
);
my $caller='';

(Getopt::Long::Parser->new)->getoptions(
    'c|caller=s' => \$caller,
    'i|in' => sub {
        while(<>){
            my $l = $_;
            chomp $l;
            my @l = split /\s+/,$l;
            if ($l[0]=~/^#/){
                if ($l[0]=~/^##FORMAT=<ID=([^,]+)/) {
                    say exists $format{$1} ? $format{$1} : $l;
                    delete $format{$1};
                } elsif ($l[0]=~/^#CHROM/) {
                    delete $format{PL};
                    delete $format{GL};
                    say $_ for values %format;
                    say "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
                } else {
                    say $l;
                }
                next;
            }

            # skip undefined GT
            next if $l[-1] =~ /^\./;

            my @t = split /:/,$l[-2];
            my @v = split /:/,$l[-1];
            my $i;

            #bcftools norm adapts DP, AD, DP4, PL, GL during slipt of multiallelic sites

            my ($dp4) = indexes { $_=~/^DP4$/ } @t;
            if(! defined $dp4){
                my @dp4;
                if(($i) = indexes { $_=~/^SB$/ } @t){ #gatk4 haplotypecaller fix (requieres -A StrandBiasBySample)
                    @dp4 = split/,/,$v[$i];
                } elsif($l[-3]=~/[\s;]SRF=([^\s;]+)/){ #freebayes fix
                    @dp4 = ($1);
                    $l[-3]=~/[\s;]SRR=([^\s;]+)/;
                    push @dp4 , $1;
                    $l[-3]=~/[\s;]SAF=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);
                    $l[-3]=~/[\s;]SAR=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);
                } elsif($l[-4]=~/[\s;]TCF=([^\s;]+)/){ #platypus fix
                    my $allf = $1;
                    $l[-4]=~/[\s;]TCR=([^\s;]+)/;
                    my $allr=$1;
                    $l[-4]=~/[\s;]NF=([^\s;]+)/;
                    my $altf=sum(split/,/,$1);
                    $l[-4]=~/[\s;]NR=([^\s;]+)/;
                    my $altr=sum(split/,/,$1);
                    @dp4 = ($allf-$altf,$allr-$altr,$altf,$altr);
                } elsif(($i) = indexes { $_=~/^ALD$/ } @t){ #vardict fix
                    @dp4 = ($v[$i]);
                    ($i) = indexes { $_=~/^RD$/ } @t;
                    push @dp4 , $v[$i];
                } else {
                    next;
                }

                if(defined $dp4){
                    $v[$dp4] = join(",",@dp4);
                } else {
                    splice @t, 1, 0, 'DP4';
                    splice @v, 1, 0, join(",",@dp4);
                }
            }

            ($i) = indexes { $_=~/^AD$/ } @t;
            if (defined $i && $v[$i]!~/,/){ #varscanfix (only correct for 1 or 2 alleles)
                if (my ($j) = indexes { $_=~/^RD$/ } @t){
                    my @alleles = split /,/,$l[4];
                    pop @alleles;
                    my @ad;
                    push @ad, $v[$j];
                    push @ad, $v[$j]-$v[$i] for @alleles;
                    push @ad, $v[$i];
                    $v[$i] = join ",",@ad;
                }
            }

            unless (grep { $_ =~ /^DP$/ } @t){
                if (($i) = indexes { $_=~/^NR$/ } @t) { #platypus fix
                    my ($j) = indexes { $_=~/^NV$/ } @t;
                    my $dp = max(split /,/,$v[$i]);
                    my @ad = split /,/,$v[$j];
                    unshift @ad, $dp-sum(@ad);
                    splice @t, 1, 0, 'DP';
                    splice @t, 1, 0, 'AD';
                    splice @v, 1, 0, $dp;
                    splice @v, 1, 0, join(",",@ad);
                } elsif (($i) = indexes { $_=~/^DP4$/ } @t) {
                    my $dp = sum split /,/,$v[$i];
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'DP';
                } elsif (($i) = indexes { $_=~/^AD$/ } @t) {
                    my $dp = sum split /,/,$v[$i];
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'DP';
                } else {
                    next;
                }
            }

            # add or recalculate GQ as difference from best and second best observed (pred scaled)-genotype likelihood            
            my $gq=0;
            if(($i) = indexes { $_ =~ /^PL$/ } @t){
                my @q = sort {$a <=> $b} split /,/,$v[$i];
                $gq = sprintf("%.4f",$q[1] - $q[0]); 
            } elsif(($i) = indexes { $_ =~ /^GL$/ } @t){
                my @pl = split /,/,$v[$i];
                $_ = ($_*-1)/10 for @pl;
                my @q = sort {$a <=> $b} @pl;
            } elsif($l[-3]=~/[\s;]TLOD=([^\s;]+)/){ #gatk3 mutect2 fix
                $gq = $1;
            } elsif (($i) = indexes { $_ =~ /^GQ$/ } @t){
                $gq = $v[$i];
            }
            $gq = sprintf("%.4f",$gq); #fix for vt normalize which cannot handle E-10 like representation
            $gq=~s/0+$//g;
            $gq=~s/\.$//g;
            ($i) = indexes { $_ =~ /^GQ$/ } @t;
            if (defined $i){
                $v[$i]=$gq;
            } else {
                splice @t, 1, 0, 'GQ';
                splice @v, 1, 0, $gq;
            }

            # finally add COV, MAV and strand usage fraction info tags
            ($i) = indexes { $_=~/^DP$/ } @t;
            my ($j) = indexes { $_=~/^DP4$/ } @t;
            my @dp4 = split/,/,$v[$j];
            ($j) = indexes { $_=~/^AD$/ } @t;
            my @ad = split /,/,$v[$j]; 

            $v[$i] = min(sum(@ad), $v[$i]); #varscan, vardict and freebayes fix DP >= sum(AD|DP4) - (== total COV, not filtered depth), whereas AD values are based on filtered reads
            my $cov = max(sum(@ad), sum(@dp4), $v[$i]); #try to find real COV
            shift @ad;
            my @maf;
            push @maf, $cov == 0 ? 0 : sprintf("%.4f",$_/$cov) for @ad;
            $_=~s/(\.[^0]*)0+$/\1/ for @maf;
            $_=~s/\.$// for @maf;
            my $rsf = max($dp4[0],$dp4[1]) == 0 ? 0 : sprintf("%.4f",min($dp4[0],$dp4[1])/max($dp4[0],$dp4[1]));
            $rsf=~s/(\.[^0]*)0+$/\1/;
            $rsf=~s/\.$//g;
            my $asf = max($dp4[2],$dp4[3]) == 0 ? 0 : sprintf("%.4f",min($dp4[2],$dp4[3])/max($dp4[2],$dp4[3]));
            $asf=~s/(\.[^0]*)0+$/\1/;
            $asf=~s/\.$//;

            ($i) = indexes { $_=~/^COV$/ } @t;
            if(defined $i){
                $v[$i]=$cov
            } else {
                splice @t, 1, 0, 'COV';
                splice @v, 1, 0, $cov;
            }
            ($i) = indexes { $_=~/^MAF$/ } @t;
            if(defined $i){
                $v[$i]=join(",",@maf);
            } else {
                splice @t, 1, 0, 'MAF';
                splice @v, 1, 0, join(",",@maf);
            }
            ($i) = indexes { $_=~/^ASF$/ } @t;
            if(defined $i){
                $v[$i]=join(",",@maf);
            } else {
                splice @t, 1, 0, 'ASF';
                splice @v, 1, 0, $asf;
            }
            ($i) = indexes { $_=~/^RSF$/ } @t;
            if(defined $i){
                $v[$i]=join(",",@maf);
            } else {
                splice @t, 1, 0, 'RSF';
                splice @v, 1, 0, $rsf;
            }

            # print
            $l[-1] = join ':' , @v;
            $l[-2] = join ':' , @t;
            $l = join "\t" , @l;
            say $l;
        }
    }
) or die "Unkown option";
