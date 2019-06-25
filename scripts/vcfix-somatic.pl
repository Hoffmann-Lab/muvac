#! /usr/bin/env perl
use strict;
use warnings;
use v5.10;
use List::MoreUtils qw(indexes uniq);
use List::Util qw(min max sum);
use Getopt::Long;
die "option -i|--in and/or -s|--switch is missing" if $#ARGV == -1;
my $switch='';

(Getopt::Long::Parser->new)->getoptions(
    's|switch' => \$switch,
    'i|in' => sub {
        
        while(<>){
            my $l = $_;
            chomp $l;
            my @l = split /\s+/,$l;
            if ($l[0]=~/^#/){
                if ($l[0]=~/^#CHROM/) {
                    say '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Phred-scaled genotype quality">';
                    say '##FORMAT=<ID=PL,Number=.,Type=Float,Description="List of Phred-scaled genotype likelihoods">'; #Number=G dont work for bcftools
                    say '##FORMAT=<ID=MAF,Number=.,Type=Float,Description="Minor allele frequency">';
                    say '##FORMAT=<ID=COV,Number=1,Type=Integer,Description="Position read coverage">';
                    say '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Filtered allelic depths for the ref and alt alleles">';
                    say '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth at the locus">';
                    say '##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Strand specific alt read counts: alt-fwd, alt-rev">';
                    say "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";
                } elsif ($l[0]=~/^##FORMAT=<ID=(PL|GQ|AD|DP),/) {
                    next;
                } else {
                    say $l;
                }
                next;
            }

            next if $l[-1]=~/^\./ || $l[-2]=~/^\./;

            my @t = split /:/,$l[-3];
            my @vt;
            my @vn;
            my $i;

            if($switch){ #needed for mutect, vardict, deepsnv, freebayes/vcfsamplediff
                @vt = split /:/,$l[-2];
                @vn = split /:/,$l[-1];
            } else {
                @vt = split /:/,$l[-1];
                @vn = split /:/,$l[-2];
            }

            my ($dp) = grep { $_ =~ /^DP$/ } @t;
            unless($dp){
                if($l[-4]=~/[\s;]TLOD=([^\s;]+)/){ #mutect2 fix
                    my $gq = $1;
                    ($i) = indexes { $_ =~ /^AD$/ } @t;
                    $dp = sum split /,/,$vt[$i];
                    splice @vt, 1, 0, $dp;
                    splice @vt, 1, 0, $gq;

                    $l[-4]=~/[\s;]NLOD=([^\s;]+)/;
                    $gq = $1;
                    $dp = sum split /,/,$vn[$i];
                    splice @vn, 1, 0, $dp;
                    splice @vn, 1, 0, $gq;

                    splice @t, 1, 0, 'DP';
                    splice @t, 1, 0, 'GQ';
                } elsif (($i) = indexes { $_=~/^NR$/ } @t) { #platypus fix
                    my ($j) = indexes { $_=~/^NV$/ } @t;

                    $dp = max(split /,/,$vt[$i]);
                    my @ad = split /,/,$vt[$j];
                    unshift @ad, $dp-sum(@ad);

                    splice @vt, 1, 0, join(",",@ad);
                    splice @vt, 1, 0, $dp;

                    $dp = max(split /,/,$vn[$i]);
                    @ad = split /,/,$vn[$j];
                    unshift @ad, $dp-sum(@ad);
                    
                    splice @vn, 1, 0, join(",",@ad);
                    splice @vn, 1, 0, $dp;

                    splice @t, 1, 0, 'AD';
                    splice @t, 1, 0, 'DP';
                } elsif (($i) = indexes { $_=~/^DP4$/ } @t) {
                    $dp = sum split /,/,$vt[$i];
                    splice @vt, 1, 0, $dp;

                    $dp = sum split /,/,$vn[$i];
                    splice @vn, 1, 0, $dp;

                    splice @t, 1, 0, 'DP';
                } elsif (($i) = indexes { $_=~/^AD$/ } @t) {
                    $dp = sum split /,/,$vt[$i];
                    splice @vt, 1, 0, $dp;

                    $dp = sum split /,/,$vn[$i];
                    splice @vn, 1, 0, $dp;

                    splice @t, 1, 0, 'DP';
                } else {
                    next;
                }
            }

            #add dp2 (dp4) - better move to vcfixuniq to handle mutliallelic sites better than just using sum like varscan (ok if order in file corresponds to comma sep. list i.e. e.g. SAF=10,20)
            my ($dp4) = grep { $_ =~ /^DP2$/ } @t;
            unless($dp4){ #dp4 hard to predict from caller using paired data and vcfsamplediff afterwards due to mixing up normal and tumor ref + alt counts
                if($l[-4]=~/[\s;]SAF=([^\s;]+)/){ #freebayes fix
                    my @dp4 = ($1);
                    $l[-4]=~/[\s;]SRR=([^\s;]+)/;
                    push @dp4 , $1;
                    $l[-4]=~/[\s;]SAF=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);
                    $l[-4]=~/[\s;]SAR=([^\s;]+)/;
                    push @dp4 , sum(split/,/,$1);

                    splice @t, 1, 0, 'DP2';
                    splice @vt, 1, 0, "$dp4[2],$dp4[3]";
                    splice @vn, 1, 0, '.,.';
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
                    splice @vt, 1, 0, "$dp4[2],$dp4[3]";
                    splice @vn, 1, 0, '.,.';
                }elsif(($i) = indexes { $_=~/^REF_F1R2$/ } @t){ #mutect2 fix
                    my @dp4 = ($vt[$i]);
                    ($i) = indexes { $_=~/^REF_F2R1$/ } @t;
                    push @dp4 , $vt[$i];
                    ($i) = indexes { $_=~/^ALT_F1R2$/ } @t; #only caller which reports first pair strandness aka fragment strandness to be used in DP4 instead of read strandness + in some cases count error (+1)
                    push @dp4 , $vt[$i]; #does not report multiple alt
                    ($i) = indexes { $_=~/^ALT_F2R1$/ } @t;
                    push @dp4 , $vt[$i];

                    splice @t, 1, 0, 'DP2';
                    splice @vt, 1, 0, "$dp4[2],$dp4[3]";
                    splice @vn, 1, 0, '.,.';
                } elsif(($i) = indexes { $_=~/^ALD$/ } @t){ #vardict fix
                    my @dp4 = ($vt[$i]); #des not report mult alt
                    ($i) = indexes { $_=~/^RD$/ } @t;
                    unshift @dp4 , $vt[$i];

                    splice @t, 1, 0, 'DP2';
                    splice @vt, 1, 0, $dp4[1];
                    splice @vn, 1, 0, '.,.';
                } elsif(($i) = indexes { $_=~/^DP4$/ } @t){ #varscan fix
                    my @dp4 = split/,/,$vt[$i];
                    splice @t, 1, 0, 'DP2';
                    splice @vt, 1, 0, "$dp4[2],$dp4[3]";
                    splice @vn, 1, 0, '.,.';
                } else {
                    next;
                }
            }

            # add MAF, COV and AD varscan fix
            if (($i) = indexes { $_=~/^AD$/ } @t) {
                my @ad;
                if ($vt[$i]!~/,/){ #varscanfix (only correct for 1 or 2 alleles)
                    if (my ($j) = indexes { $_=~/^RD$/ } @t){
                        my @alleles = split /,/,$l[4];
                        pop @alleles;
                        push @ad, $vt[$j];
                        push @ad, $vt[$j]-$vt[$i] for @alleles;
                        push @ad, $vt[$i];
                        $vt[$i] = join ",",@ad;

                        @ad = ();
                        push @ad, $vn[$j];
                        push @ad, $vn[$j]-$vn[$i] for @alleles;
                        push @ad, $vn[$i];
                        $vn[$i] = join ",",@ad;
                    }
                }

                my ($j) = indexes { $_=~/^DP$/ } @t;
                my @dp4 = split/,/,$vt[$j];
                ($j) = indexes { $_=~/^DP$/ } @t;
                @ad = split /,/,$vt[$i];
                $vt[$j] = min(sum(@ad), $vt[$j]); #varscan, vardict and freebayes fix DP >= sum(AD|DP4) - (== total COV, not filtered depth), whereas AD values are based on filtered reads
                my $covt = max(sum(@ad), sum(@dp4), $vt[$j]); #try to find real COV
                
                shift @ad;
                my @maft;
                push @maft, $covt == 0 ? 0 : sprintf("%.4f",$_/$covt) for @ad;

                @ad = split /,/,$vn[$i];
                $vn[$j] = min(sum(@ad), $vn[$j]);
                my $covn = max(sum(@ad), $vn[$j]);

                shift @ad;
                my @mafn;
                push @mafn, $covn == 0 ? 0 : sprintf("%.4f",$_/$covn) for @ad;

                if (($i) = indexes { $_ =~ /^MAF$/ } @t){
                    $vt[$i] = join(",",@maft);
                    $vn[$i] = join(",",@mafn);
                } else {
                    splice @t, 1, 0, 'MAF';
                    splice @t, 1, 0, 'COV';
                    splice @vt, 1, 0, join(",",@maft);
                    splice @vn, 1, 0, join(",",@mafn);
                    splice @vt, 1, 0, $covt;
                    splice @vn, 1, 0, $covn;
                }
                if (($i) = indexes { $_ =~ /^COV$/ } @t){
                    $vt[$i] = $covt;
                    $vn[$i] = $covn;
                } else {
                    splice @t, 1, 0, 'COV';
                    splice @vt, 1, 0, $covt;
                    splice @vn, 1, 0, $covn;
                }
            } else {
                next;
            }

            # add or recalculate (bctools norm fix) GQ as difference from best and second best observed (pred scaled)-genotype likelihood
            if (($i) = indexes { $_ =~ /^PL$/ } @t) {
                my @qt = sort {$a <=> $b} split /,/,$vt[$i];
                my @qn = sort {$a <=> $b} split /,/,$vn[$i];
                my $gqt = sprintf("%.4f",$qt[1] - $qt[0]); #vt normalize fix to avoid e+X/e-X representation
                my $gqn = sprintf("%.4f",$qn[1] - $qn[0]);
                if (($i) = indexes { $_ =~ /^GQ$/ } @t){
                    $vt[$i] = $gqt;
                    $vn[$i] = $gqn;
                } else {
                    splice @t, 1, 0, 'GQ';
                    splice @vt, 1, 0, $gqt;
                    splice @vn, 1, 0, $gqn;
                }
            } elsif(($i) = indexes { $_ =~ /^GL$/ } @t) {
                my @plt = split /,/,$vt[$i];
                my @pln = split /,/,$vn[$i];
                $_ = ($_*-1)/10 for @plt;
                $_ = ($_*-1)/10 for @pln;
                my @qt = sort {$a <=> $b} @plt;
                my @qn = sort {$a <=> $b} @pln;
                my $gqt = sprintf("%.4f",$qt[1] - $qt[0]);
                my $gqn = sprintf("%.4f",$qn[1] - $qn[0]);
                splice @t, 1, 0, 'PL';
                splice @vt, 1, 0, join(",",@plt);
                splice @vn, 1, 0, join(",",@pln);
                if (($i) = indexes { $_ =~ /^GQ$/ } @t){
                    $vt[$i] = $gqt;
                    $vn[$i] = $gqn;
                } else {
                    splice @t, 1, 0, 'GQ';
                    splice @vt, 1, 0, $gqt;
                    splice @vn, 1, 0, $gqn;
                }
            }

            ($i) = indexes { $_=~/^GQ$/ } @t; #vardict fix
            if ($i){
                $vt[$i] = 0 if $vt[$i] eq ".";
                $vn[$i] = 0 if $vn[$i] eq ".";
            } else {
                splice @t, 1, 0, 'GQ';
                splice @vt, 1, 0, 0;
                splice @vn, 1, 0, 0;
            }

            $l[-1] = join ':' , @vt;
            $l[-2] = join ':' , @vn;
            $l[-3] = join ':' , @t;

            # $l[-1]=~/^(\d)\/(\d)/;
            # $l[-1]=~s/^(\d)\/(\d)/$2\/$1/ if $2 < $1;
            # $l[-2]=~/^(\d)\/(\d)/;
            # $l[-2]=~s/^(\d)\/(\d)/$2\/$1/ if $2 < $1;

            $l = join "\t" , @l;
            say $l;
        }
    }
) or die "Unkown option";
