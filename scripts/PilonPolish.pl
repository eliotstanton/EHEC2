#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;

# ---------------------------------------------------------------------------- #

# File Name:		PilonPolish
# Date Created:		13 February. 2020
# Last Modified:	13 February. 2020
# Created By:		Eliot Stanton (estanton@Wisc.Edu)

# Description:		This script iteratively polishes a FASTA file using Bowtie
#					and Pilon.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
my $file_FASTA	= $ARGV[0];
my $file_reads1	= $ARGV[1];
my $file_reads2	= $ARGV[2];
my $var_prefix	= "pilon";
my $file_SAM	= "pilon/$var_prefix\.SAM";
my $file_BAM	= "pilon/$var_prefix\.BAM";
my $file_sorted	= "pilon/$var_prefix\.sorted\.BAM";

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nPilonPolish.pl [FASTA] [FASTQ1] [FASTQ2]\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_FASTA, $var_reads1, or $var_reads2 are missing, print $var_help and 
# stop script:
unless ( scalar @ARGV == 3 ) {

	print "$var_help\n";
	
	print "\tFASTA and FASTQ files are required!\n" and exit;

}

# Check if directory called $var_prefix exists and create it if needed:
Check::Directory ( "pilon", \%hash_variables );

# ---------------------------------------------------------------------------- #

my $var_counter	= 0;

Align ( $file_FASTA, $file_reads1, $file_reads2, $var_prefix, $file_SAM, $file_BAM, $file_sorted, $var_counter );

Pilon ( $file_FASTA, $file_sorted, $var_prefix );

my $var_changes	= "$var_prefix/pilon.changes";

$file_FASTA		= "$var_prefix/pilon.fasta";

while ( -s $var_changes ) {

	$var_counter++;

	Align ( $file_FASTA, $file_reads1, $file_reads2, $var_prefix, $file_SAM, $file_BAM, $file_sorted, $var_counter );

	Pilon ( $file_FASTA, $file_sorted, $var_prefix );

	last if $var_counter > 3;

}

# Move and rename fasta file:

# ---------------------------------------------------------------------------- #

sub Pilon {

	# Arguments:
	my ( $file_FASTA, $file_sorted, $var_prefix )	= @_;

	# ----------------------------------

	print "\tpilon --genome $file_FASTA --frags $file_sorted --outdir $var_prefix --changes\n";
	system ( "pilon --genome $file_FASTA --frags $file_sorted --outdir $var_prefix --changes" );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Align
# Description:		Align Illumina data to FASTA file using Bowtie.

# ----------------------------------------------------------

sub Align {

	# Arguments:
	my ( $file_FASTA, $file_reads1, $file_reads2, $var_prefix, $var_SAM, $var_BAM, $var_sorted, $var_counter )	= @_;

	# ----------------------------------

	# Build Bowtie index:
	print "\t$var_counter: bowtie-build --quiet --threads 8 $file_FASTA $var_prefix/$var_prefix\n";
	system ( "bowtie-build --quiet --threads 8 $file_FASTA $var_prefix/$var_prefix" );

	# Aligned paired-end reads to contigs using Bowtie:
	print "\t$var_counter: bowtie $var_prefix/$var_prefix --threads 8 -1 $file_reads1 -2 $file_reads2 -X 1000 -S $var_SAM\n";
	system ( "bowtie $var_prefix/$var_prefix --threads 8 -1 $file_reads1 -2 $file_reads2 -X 1000 -S $var_SAM" );

	# ----------------------------------

	# Convert SAM files to sorted BAM files
	print "\t$var_counter: samtools view -@ 8 -u $var_SAM > $var_BAM\n";
	system ( "samtools view -@ 8 -u $var_SAM > $var_BAM" );

	# Sort the BAM file:
 	print "\t$var_counter: samtools sort -@ 6 $var_BAM -o $var_sorted\n";
	system ( "samtools sort -@ 4 $var_BAM -o $var_sorted" );

	# Index the BAM file:
	print "\t$var_counter: samtools index -@ 8 $var_sorted\n";
	system ( "samtools index -@ 8 $var_sorted" );

	# ----------------------------------

	# End subroutine:
	return;


}

# ---------------------------------------------------------------------------- #
