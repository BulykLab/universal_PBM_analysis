#!/usr/bin/perl
#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#
use strict;
use warnings;

#use lib './';
use lib '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite';
use seed_and_wobble_modules;

#########################################################
### Takes as input a single seed and TWO lists of all
###   "deBruijn" intensities from an Agilent 4x44K array.
###   MUST BE PRE-ORDERED FROM BRIGHTEST TO DIMMEST.
### Also takes as input a file with list of gapped patterns
###   covered for that "k".
### Outputs the PWM **combined by averaging the PWM from the
###   two arrays** generated from "wobbling" each position
###   within and outside the seed.
###
### M. Berger 11/18/07 (modified from 04/27/07)
#########################################################

if($#ARGV !=3){die "Usage is:\n\tPBM data file #1 (sorted by intensities)\n\tPBM data file #2 (sorted by intensities)\n\tseed\n\tlist of all covered gapped patterns (e.g., 111..1.1.11)\nFor example: \"perl wobble_single_seed_twoarray.pl proteinA_v1_combinatorial.txt proteinA_v2_combinatorial.txt TGA.GTCAT all_patterns.txt > proteinA_pwm.txt\"\n\n";}

my $intensity_file_v1 = shift; #spot intensities, v1
my $intensity_file_v2 = shift; #spot intensities, v2
my $seed = shift; #k-mer seed, potentially gapped
my $pattern_file = shift; #file with patterns covered on array

### PARAMETERS for the width of the motif
my $spotlength=36; # total number of positions after common primer
my $startposition=2; # position from end of strand to consider, starting at 1

my @seed_characters = split//,$seed;
my $seed_numinfopos = 0;
for (my $pos=0; $pos<=$#seed_characters; $pos++) {
	if ($seed_characters[$pos] ne ".") {$seed_numinfopos++;}
}

###########################################################################
### Read in list of patterns covered in universal PBM design
###########################################################################

my @patterns;

open (FHH, "<$pattern_file") or die "Cannot open patterns file.\n";

while (my $text = <FHH>) {
	my $numinfopos=0;
	chomp $text;
	push @patterns, $text;	
	my @characters = split//,$text;
	for (my $p=0; $p<=$#characters; $p++) {
		if ($characters[$p] ne ".") {$numinfopos++;}
	}
	if ($numinfopos != $seed_numinfopos) {
		die "Number of positions in seed $seed does not agree with $text in $pattern_file.\n";
	}
}	

close FHH;

###########################################################################
### Seed-n-Wobble PWM construction
###########################################################################

my $key;

###

my %kmerranks_v1;
my %areapwm_v1;
my $array_v1 = $intensity_file_v1;
print "$seed\n\n";
$seed = ".......".$seed.".......";

wobble_seed_rerank(\%kmerranks_v1,$seed,\%{$areapwm_v1{$seed}},$array_v1,$spotlength,$startposition,1);

my $minimuminfopos_v1=find_minimum_info_pos(\%{$areapwm_v1{$seed}},$seed,log(10000));

extend_seed_allpatterns_rerank(\%kmerranks_v1,$seed,$minimuminfopos_v1,\%{$areapwm_v1{$seed}},$array_v1,$spotlength,$startposition,1,\@patterns);

###

my %kmerranks_v2;
my %areapwm_v2;
my $array_v2 = $intensity_file_v2;

wobble_seed_rerank(\%kmerranks_v2,$seed,\%{$areapwm_v2{$seed}},$array_v2,$spotlength,$startposition,1);

my $minimuminfopos_v2=find_minimum_info_pos(\%{$areapwm_v2{$seed}},$seed,log(10000));

extend_seed_allpatterns_rerank(\%kmerranks_v2,$seed,$minimuminfopos_v2,\%{$areapwm_v2{$seed}},$array_v2,$spotlength,$startposition,1,\@patterns);

###

#Must have data for both PWM-v1 and PWM-v2 before averaging

while (($areapwm_v1{$seed}{A}[0]==0) || ($areapwm_v2{$seed}{A}[0]==0)) {
    shift @{$areapwm_v1{$seed}{A}}; shift @{$areapwm_v1{$seed}{C}}; shift @{$areapwm_v1{$seed}{G}}; shift @{$areapwm_v1{$seed}{T}};
    shift @{$areapwm_v2{$seed}{A}}; shift @{$areapwm_v2{$seed}{C}}; shift @{$areapwm_v2{$seed}{G}}; shift @{$areapwm_v2{$seed}{T}};
}

while (($areapwm_v1{$seed}{A}[-1]==0) || ($areapwm_v2{$seed}{A}[-1]==0)) {
    pop @{$areapwm_v1{$seed}{A}}; pop @{$areapwm_v1{$seed}{C}}; pop @{$areapwm_v1{$seed}{G}}; pop @{$areapwm_v1{$seed}{T}};
    pop @{$areapwm_v2{$seed}{A}}; pop @{$areapwm_v2{$seed}{C}}; pop @{$areapwm_v2{$seed}{G}}; pop @{$areapwm_v2{$seed}{T}};
}

#Average matrix elements

my %areapwm_combined;

print "Enrichment Score Matrix\n\n";
foreach $key (sort keys %{$areapwm_v1{$seed}}) {
    print "$key:";
    for (my $y=0; $y<=$#{$areapwm_v1{$seed}{$key}}; $y++) {
	$areapwm_combined{$seed}{$key}[$y] = ($areapwm_v1{$seed}{$key}[$y] + $areapwm_v2{$seed}{$key}[$y]) / 2;
	printf "\t$areapwm_combined{$seed}{$key}[$y]";
    }
    print "\n";
}

print "\nEnergy matrix for enoLOGOS\n\n";
print "PO";
for (my $counter=1; $counter<=($#{$areapwm_combined{$seed}{A}}+1); $counter++) {
	print "\t$counter";
}
print "\n";
foreach $key (sort keys %{$areapwm_combined{$seed}}) {
    print "$key:";
    for (my $y=0; $y<=$#{$areapwm_combined{$seed}{$key}}; $y++) {
	my $logscaled = $areapwm_combined{$seed}{$key}[$y]*(-log(10000));
	print "\t$logscaled";
    }
    print "\n";
}
#print "\nReverse complement for enoLOGOS\n\n";
#print "PO";
#for (my $counter=1; $counter<=($#{$areapwm_combined{$seed}{A}}+1); $counter++) {
#	print "\t$counter";
#}
#print "\n";
#foreach $key (sort keys %{$areapwm_combined{$seed}}) {
#    my $compkey;
#    if ($key eq "A") {$compkey="T";}
#    if ($key eq "C") {$compkey="G";}
#    if ($key eq "G") {$compkey="C";}
#    if ($key eq "T") {$compkey="A";}
#    print "$compkey:";
#    for (my $y=$#{$areapwm_combined{$seed}{$key}}; $y>=0; $y--) {
#	my $logscaled = $areapwm_combined{$seed}{$key}[$y]*(-log(10000));
#	print "\t$logscaled";
#    }
#    print "\n";
#}
print "\nProbability matrix\n\n";
foreach $key (sort keys %{$areapwm_combined{$seed}}) {
	print "$key:";
	for (my $y=0; $y<=$#{$areapwm_combined{$seed}{$key}}; $y++) {
      	my $numerator = exp(log(10000)*$areapwm_combined{$seed}{$key}[$y]);
		my $denominator = exp(log(10000)*$areapwm_combined{$seed}{A}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{C}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{G}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{T}[$y]);
		my $probability = $numerator/$denominator;
		print "\t$probability";
	}
	print "\n";
}
print "\n";
