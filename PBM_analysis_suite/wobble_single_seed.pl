#!/usr/bin/perl
#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#
use strict;
use warnings;

#use lib './';
use lib '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite';
use seed_and_wobble_modules;

##################################################################
### Takes as input a two-column list of all "combinatorial" probe
###   intensities and sequences from an Agilent 'all 10-mer' 4x44K
###   array.  MUST BE PRE-ORDERED FROM BRIGHTEST TO DIMMEST.
### Also takes as input a single k-mer (with "." for gaps).
### Also takes as input a file with list of gapped patterns
###   covered for that "k" (e.g., 111..11.1.11).
### Outputs the PWM generated from "wobbling" each position
###   within and outside the seed.
###
### M. Berger 04/27/07
##################################################################

if($#ARGV !=2){die "Usage is:\n\tPBM data (sorted by intensities)\n\tseed\n\tlist of all covered gapped patterns\nFor example: \"perl wobble_single_seed.pl proteinA_combinatorial.txt TGA.GTCAT all_patterns.txt > proteinA_pwm.txt\"\n\n";}

my $intensity_file = shift; #spot intensities
my $seed = shift; #k-mer seed, potentially gapped
my $pattern_file = shift; #file with patterns covered on array

my @data_matrix;

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

##########################################################
### Read in intensities and create array of sequences
##########################################################

open (FH, "<$intensity_file") or die "Cannot open intensity file.\n";

my $spot = 1;

while (my $text = <FH>) {
    chomp $text;
    my @line = split ("\t", $text);
    $data_matrix[$spot][0] = $line[0]; ### spot intensity
    $data_matrix[$spot][1] = $line[1]; ### spot sequence
    if ($spot>1) {
	if ($data_matrix[$spot][0] > $data_matrix[$spot-1][0]) {
	    die "Probes in input file are not sorted from brightest to dimmest.\n";
	}
    }
    $spot++;
}

close FH;

my %kmerranks; #data structure that stores ranks of every k-mer
my %kmerareas; #stores area for each gapped k-mer
my $numberspotsarray=$#data_matrix; #number of spots on array

my %observedpatterns;
my $key;


###########################################################################
### Seed-n-Wobble PWM construction
###########################################################################

my %areapwm;
my $array=$intensity_file;
print "$seed\n\n";
$seed = ".......".$seed.".......";

wobble_seed_rerank(\%kmerranks,$seed,\%{$areapwm{$seed}},$array,$spotlength,$startposition,1);

my $minimuminfopos=find_minimum_info_pos(\%{$areapwm{$seed}},$seed,log(10000));

extend_seed_allpatterns_rerank(\%kmerranks,$seed,$minimuminfopos,\%{$areapwm{$seed}},$array,$spotlength,$startposition,1,\@patterns);

while (($areapwm{$seed}{A}[0]==0) && ($areapwm{$seed}{C}[0]==0) && ($areapwm{$seed}{G}[0]==0) && ($areapwm{$seed}{T}[0]==0)) {
	shift @{$areapwm{$seed}{A}}; shift @{$areapwm{$seed}{C}}; shift @{$areapwm{$seed}{G}}; shift @{$areapwm{$seed}{T}};
}

while (($areapwm{$seed}{A}[-1]==0) && ($areapwm{$seed}{C}[-1]==0) && ($areapwm{$seed}{G}[-1]==0) && ($areapwm{$seed}{T}[-1]==0)) {
	pop @{$areapwm{$seed}{A}}; pop @{$areapwm{$seed}{C}}; pop @{$areapwm{$seed}{G}}; pop @{$areapwm{$seed}{T}};
}

print "Enrichment Score Matrix\n\n";
foreach $key (sort keys %{$areapwm{$seed}}) {
	print "$key:";
	for (my $y=0; $y<=$#{$areapwm{$seed}{$key}}; $y++) {
		print "\t$areapwm{$seed}{$key}[$y]";
	}
	print "\n";
}
print "\nEnergy matrix for enoLOGOS\n\n";
print "PO";
for (my $counter=1; $counter<=($#{$areapwm{$seed}{A}}+1); $counter++) {
	print "\t$counter";
}
print "\n";
foreach $key (sort keys %{$areapwm{$seed}}) {
	print "$key:";
        for (my $y=0; $y<=$#{$areapwm{$seed}{$key}}; $y++) {
	    my $logscaled = $areapwm{$seed}{$key}[$y]*(-log(10000));
	    print "\t$logscaled";
	}
	print "\n";
}
#print "\nReverse complement for enoLOGOS\n\n";
#print "PO";
#for (my $counter=1; $counter<=($#{$areapwm{$seed}{A}}+1); $counter++) {
#	print "\t$counter";
#}
#print "\n";
#foreach $key (sort keys %{$areapwm{$seed}}) {
#	my $compkey;
#	if ($key eq "A") {$compkey="T";}
#	if ($key eq "C") {$compkey="G";}
#	if ($key eq "G") {$compkey="C";}
#	if ($key eq "T") {$compkey="A";}
#	print "$compkey:";
#	for (my $y=$#{$areapwm{$seed}{$key}}; $y>=0; $y--) {
#	    my $logscaled = $areapwm{$seed}{$key}[$y]*(-log(10000));
#	    print "\t$logscaled";
#	}
#	print "\n";
#}
print "\nProbability matrix\n\n";
foreach $key (sort keys %{$areapwm{$seed}}) {
	print "$key:";
	for (my $y=0; $y<=$#{$areapwm{$seed}{$key}}; $y++) {
		my $numerator = exp(log(10000)*$areapwm{$seed}{$key}[$y]);
		my $denominator = exp(log(10000)*$areapwm{$seed}{A}[$y]) + exp(log(10000)*$areapwm{$seed}{C}[$y]) + exp(log(10000)*$areapwm{$seed}{G}[$y]) + exp(log(10000)*$areapwm{$seed}{T}[$y]);
		my $probability = $numerator/$denominator;
		print "\t$probability";
	}
	print "\n";
}
print "\n";
