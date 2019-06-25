#! /usr/bin/env perl
use strict;
use warnings;
use v5.10;
use List::MoreUtils qw(indexes);
use List::Util qw(min max sum);
use Getopt::Long;
die "option -i|--in is missing" if $#ARGV == -1;

(Getopt::Long::Parser->new)->getoptions(
    'i|in' => sub {
        while(<>){
            my $l = $_;
            chomp $l;
            my @l = split /\s+/,$l;
            if ($l[0]=~/^#/){
                if ($l[0]=~/^#CHROM/) {
                    say '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Phred-scaled genotype quality">';
                    say '##FORMAT=<ID=PL,Number=.,Type=Float,Description="List of Phred-scaled genotype likelihoods">'; #Number=G dont work for bcftools
                    say '##FORMAT=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequence">';
                    say '##FORMAT=<ID=COV,Number=1,Type=Integer,Description="Position read coverage">';
                    say '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">';
                    say '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth at the locus">';
                    say '##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Strand specific alt read counts: alt-fwd, alt-rev">';
                    say "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
                } elsif ($l[0]=~/^##FORMAT=<ID=(PL|GQ),/) {
                    next;
                } else {
                    say $l;
                }
                next;
            }

            next if $l[-1] =~ /^\./;

            my @t = split /:/,$l[-2];
            my @v = split /:/,$l[-1];
            my $i;

            my ($dp) = grep { $_ =~ /^DP$/ } @t;
            unless ($dp){ 
                if (($i) = indexes { $_=~/^NR$/ } @t) { #platypus fix
                    my ($j) = indexes { $_=~/^NV$/ } @t;

                    $dp = max(split /,/,$v[$i]);
                    my @ad = split /,/,$v[$j];
                    unshift @ad, $dp-sum(@ad);

                    splice @v, 1, 0, join(",",@ad);
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'AD';
                    splice @t, 1, 0, 'DP';
                } elsif($l[-3]=~/[\s;]DP=([^\s;]+)/){ #haplotypecaller fix
                    $dp = $1;
                    splice @t, 1, 0, 'DP';
                    splice @v, 1, 0, $dp;
                } elsif($l[-3]=~/[\s;]TLOD=([^\s;]+)/){ #mutect2 fix
                    my $gq = $1;
                    ($i) = indexes { $_ =~ /^AD$/ } @t;
                    $dp = sum split /,/,$v[$i];
                    splice @t, 1, 0, 'DP';
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'GQ';
                    splice @v, 1, 0, $gq;
                } elsif (($i) = indexes { $_=~/^DP4$/ } @t) {
                    $dp = sum split /,/,$v[$i];
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'DP';
                } elsif (($i) = indexes { $_=~/^AD$/ } @t) {
                    $dp = sum split /,/,$v[$i];
                    splice @v, 1, 0, $dp;
                    splice @t, 1, 0, 'DP';
                } else {
                    next;
                }
            }

            my ($dp4) = grep { $_ =~ /^DP4$/ } @t;
            unless($dp4){
                if($l[-3]=~/[\s;]SRF=([^\s;]+)/){ #freebayes fix
                    my @dp4 = ($1);
                    $l[-3]=~/[\s;]SRR=([^\s;]+)/;
                    push @dp4 , $1;
                    $l[-3]=~/[\s;]SAF=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);
                    $l[-3]=~/[\s;]SAR=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);

                    splice @t, 1, 0, 'DP2';
                    splice @v, 1, 0, "$dp4[2],$dp4[3]";
                }elsif($l[-4]=~/[\s;]TCF=([^\s;]+)/){ #platypus fix
                    my $allf = $1;
                    $l[-4]=~/[\s;]TCR=([^\s;]+)/;
                    my $allr=$1;
                    $l[-4]=~/[\s;]NF=([^\s;]+)/;
                    my $altf=sum(split/,/,$1);
                    $l[-4]=~/[\s;]NR=([^\s;]+)/;
                    my $altr=sum(split/,/,$1);
                    my @dp4 = ($allf-$altf,$allr-$altr,$altf,$altr);

                    splice @t, 1, 0, 'DP2';
                    splice @v, 1, 0, "$dp4[2],$dp4[3]";
                }elsif(($i) = indexes { $_=~/^REF_F1R2$/ } @t){ #mutect2 fix
                    my @dp4 = ($v[$i]);
                    ($i) = indexes { $_=~/^REF_F2R1$/ } @t;
                    push @dp4 , $v[$i];
                    ($i) = indexes { $_=~/^ALT_F1R2$/ } @t;
                    push @dp4 , $v[$i];
                    ($i) = indexes { $_=~/^ALT_F2R1$/ } @t;
                    push @dp4 , $v[$i];

                    splice @t, 1, 0, 'DP2';
                    splice @v, 1, 0, "$dp4[2],$dp4[3]";
                } elsif(($i) = indexes { $_=~/^ALD$/ } @t){ #vardict fix
                    my @dp4 = ($v[$i]);
                    ($i) = indexes { $_=~/^RD$/ } @t;
                    push @dp4 , $v[$i];

                    splice @t, 1, 0, 'DP2';
                    splice @v, 1, 0, "$dp4[2],$dp4[3]";
                } elsif(($i) = indexes { $_=~/^DP4$/ } @t){ #varscan fix
                    my @dp4 = split/,/,$v[$i];
                    splice @t, 1, 0, 'DP2';
                    splice @v, 1, 0, "$dp4[2],$dp4[3]";
                } else {
                    next;
                }
            }

            # add MAF, COV and varscan fix 
            if (($i) = indexes { $_=~/^AD$/ } @t) {
                my @ad;
                if ($v[$i]!~/,/){ #varscanfix (only correct for 1 or 2 alleles)
                    if (my ($j) = indexes { $_=~/^RD$/ } @t){
                        my @alleles = split /,/,$l[4];
                        pop @alleles;
                        push @ad, $v[$j];
                        push @ad, $v[$j]-$v[$i] for @alleles;
                        push @ad, $v[$i];
                        $v[$i] = join ",",@ad;
                    }
                }

                my ($j) = indexes { $_=~/^DP4$/ } @t;
                my @dp4 = split/,/,$v[$j];
                ($j) = indexes { $_=~/^DP$/ } @t;
                @ad = split /,/,$v[$i]; 

                $v[$j] = min(sum(@ad), $v[$j]); #varscan, vardict and freebayes fix DP >= sum(AD|DP4) - (== total COV, not filtered depth), whereas AD values are based on filtered reads
                my $cov = max(sum(@ad), sum(@dp4), $v[$j]); #try to find real COV

                shift @ad;
                my @maf;
                push @maf, $cov == 0 ? 0 : sprintf("%.4f",$_/$cov) for @ad;

                if (($i) = indexes { $_ =~ /^MAF$/ } @t){
                    $v[$i] = join(",",@maf);
                } else {
                    splice @t, 1, 0, 'MAF';
                    splice @v, 1, 0, join(",",@maf);
                }
                if (($i) = indexes { $_ =~ /^COV$/ } @t){
                    $v[$i] = $cov;
                } else {
                    splice @t, 1, 0, 'COV';
                    splice @v, 1, 0, $cov;
                }
            }

            # add or recalculate (bctools norm fix) GQ as difference from best and second best observed (pred scaled)-genotype likelihood
            ($i) = indexes { $_ =~ /^PL$/ } @t;
            if ($i) {
                my @q = sort {$a <=> $b} split /,/,$v[$i];
                my $gq = sprintf("%.4f",$q[1] - $q[0]); #vt normalize fix to avoid e+X/e-X representation
                if (($i) = indexes { $_ =~ /^GQ$/ } @t){
                    $v[$i] = $gq;
                } else {
                    splice @t, 1, 0, 'GQ';
                    splice @v, 1, 0, $gq;
                }
            } else {
                ($i) = indexes { $_ =~ /^GL$/ } @t;
                if ($i) {
                    my @pl = split /,/,$v[$i];
                    $_ = ($_*-1)/10 for @pl;
                    my @q = sort {$a <=> $b} @pl;
                    my $gq = sprintf("%.4f",$q[1] - $q[0]);
                    splice @t, 1, 0, 'PL';
                    splice @v, 1, 0, join(",",@pl);
                    if (($i) = indexes { $_ =~ /^GQ$/ } @t){
                        $v[$i] = $gq;
                    } else {
                        splice @t, 1, 0, 'GQ';
                        splice @v, 1, 0, $gq;
                    }
                }
            }

            ($i) = indexes { $_=~/^GQ$/ } @t;
            if ($i){
                $v[$i] = 0 if $v[$i] eq ".";
            } else {
                splice @t, 1, 0, 'GQ';
                splice @v, 1, 0, 0;
            }  

            $l[-1] = join ':' , @v;
            $l[-2] = join ':' , @t;

            # $l[-1]=~/^(\d)\/(\d)/;
            # $l[-1]=~s/^(\d)\/(\d)/$2\/$1/ if $2 < $1;

            $l = join "\t" , @l;
            say $l;
        }
    }
) or die "Unkown option";
