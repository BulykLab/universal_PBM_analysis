#!/usr/bin/perl
#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#
use strict;
use warnings;

#use lib './';
use lib '/data/bulyk/pipelines/universal_PBM_analysis/PBM_analysis_suite';
use seed_and_wobble_modules;

#########################################################
### Takes as input a two-column list of all "combinatorial"
###   probe intensities and sequences as well as a
###   probability-based PWM.
### The PWM should be formatted exactly as direct output from
###   seed_and_wobble (labeled with 'Probability Matrix').
### Outputs a list of sequences re-ranked by the obs/exp
###   intensity, based on the sequence score of the PWM.
### ***To score sequences, uses the probability matrix
###   and the GOMER framework adapted by Chen and
###   Morris in Bioinformatics 2007 (RankMotif).
###
### M. Berger
#########################################################

if($#ARGV !=2){die "Usage is:\n\tPBM data (two column list of probe intensity and sequence, ranked from highest intensity)\n\tposition weight matrix (output from seed_and_wobble.pl)\n\toutput file name\nFor example, \"perl rerank.pl proteinA_combinatorial.txt proteinA_pwm.txt proteinA_reranked.txt\".\n\n";}

my $intensity_file=shift; #sequences and intensities
my $pwm_file=shift; #file containing PWM
my $output_file=shift; # beginning of output file

my $spotlength=36; # total number of positions to consider

##########################################################
### Read in intensities and create array of sequences
##########################################################

my @sequences;

open (FH, "<$intensity_file") or die "Cannot open intensity file.\n";

my $spotnumber = 0;

while (my $text = <FH>) {
    chomp $text;
    my @line = split ("\t", $text);
    $sequences[$spotnumber][0] = $line[0]; ### spot intensity
    $sequences[$spotnumber][1] = $line[1]; ### spot sequence
    $spotnumber++;
}

close FH;

###########################################################
### Read in PWM and create hash
###########################################################

my %areapwm;
my $k;
my @check;
for (my $x=0; $x<=3; $x++) {$check[$x]=0;}

open (FH, "<$pwm_file") or die "Cannot open PWM file.\n";

my $marker=0;
my $linenumber=0;

while (my $text = <FH>) {
    chomp $text;
    if ($text =~ "Probability") {$marker=1;}
    if ($marker==1 && $text){
	my @line = split ("\t", $text);
	if ($line[0] eq "A\:") {
	    for ($k=0; $k<$#line; $k++) {$areapwm{A}[$k] = $line[$k+1];}
	    $check[0]=1;
	}
	elsif ($line[0] eq "C\:") {
	    for ($k=0; $k<$#line; $k++) {$areapwm{C}[$k] = $line[$k+1];}
	    $check[1]=1;
	}
	elsif ($line[0] eq "G\:") {
	    for ($k=0; $k<$#line; $k++) {$areapwm{G}[$k] = $line[$k+1];}
	    $check[2]=1;
	}
	elsif ($line[0] eq "T\:") {
	    for ($k=0; $k<$#line; $k++) {$areapwm{T}[$k] = $line[$k+1];}
	    $check[3]=1;
	}
    }
    $linenumber++;
    if ($check[0]+$check[1]+$check[2]+$check[3]>=4) {last;}
}

close (FH);

if ($marker==0) {
    die "Matrix must contain the heading 'Probability Matrix', and columns must add to 1.\n";
}

my $pwm_length = $#{$areapwm{A}};
if ($#{$areapwm{C}} ne $pwm_length || $#{$areapwm{G}} ne $pwm_length || $#{$areapwm{T}} ne $pwm_length) {
    die "PWM file is in wrong format. Each row must have same number of columns.\n";
}

#######################################################
### Re-rank spots according to sequence score for PWM
#######################################################

my @reranked_sequences; # array to store re-ranked output

rerank_spots_unexpected_signals_gomer(\@sequences, \@reranked_sequences, \%areapwm, $spotlength);

open (OUTPUT1, ">$output_file") or die "Can't open output file\n";

for ($k=0; $k<=$#reranked_sequences; $k++) {
    my $value = 1/$reranked_sequences[$k][0];
    print OUTPUT1 "$value\t$reranked_sequences[$k][1]\n";
}

close (OUTPUT1);
