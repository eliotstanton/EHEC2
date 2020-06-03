#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use Check;
use General;

# ---------------------------------------------------------------------------- #

# File name:		AlignReads.pl
# Date created:		21 November, 2018
# Last modified:	20 January, 2020
# Created by:		Eliot Stanton

# Description:		This is a script for aligning reads to a FASTA sequence or
#					sequences using Bowtie.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('c:o:s:t:u:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_FASTQ		= $ARGV[1];
my $var_scale		= $hash_variables{c} || 10;
my $var_output		= $hash_variables{o};
my $var_threads		= $hash_variables{p} || 4;
my $var_prefix		= $hash_variables{s} || "default";
my $var_start		= $hash_variables{t} || 1;
my $var_end			= $hash_variables{u};
my $file_SAM		= "$var_output/$var_prefix.sam";
my $file_aligned	= "$var_output/$var_prefix.aligned.fastq";
my $var_length_all	= 0;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nAlignReads.pl [FASTA] [FASTQ1],[FASTQ2]

	-c Scaling factor (default: 10)
	-o Output directory (required)
	-p Number of threads (default: 4)
	-s Output files prefix (default: default)
	-t Start coordinate (default: 1)
	-u End coordinate (default: end of sequence)\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output or $var_FASTA is missing, print $var_help and stop script:
unless ( $var_FASTA && $hash_variables{o} ) {

	print "$var_help\n";

	print "\tOutput directory required!\n" unless $hash_variables{o};

	print "\tFASTA file required!\n" unless $var_FASTA;

	print "\tFASTQ file(s) required!\n" unless $var_FASTQ;

	exit;

}

# ---------------------------------------------------------------------------- #

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Store $var_prefix in %hash_variables:
$hash_variables { var_prefix }	= $var_prefix;

# Store $var_start in %hash_variables:
$hash_variables { var_start }	= $var_start;

# Store $var_scale in %hash_variables:
$hash_variables { var_scale }	= $var_scale;

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;
my @array_FASTQ		= split /\,/, $var_FASTQ;

# Exit if coordinates are specified and there is more than one FASTA file:
if ( scalar @array_FASTA > 1 && $var_start && $var_end ) {

	print "$var_help\n";

	print "Start and end coordinates can only be used with one FASTA file\n";

	exit;

}

# Check if files in @array_FASTA and @array_FASTQ are present:
@array_FASTA		= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };
@array_FASTQ		= @{ Check::Files ( \@array_FASTQ, \%hash_variables ) };

# Merge FASTA files back into a comma-separated list:
$var_FASTA	= join ",", @array_FASTA;

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# TODO: Check if bowtie is present:
Check::Bowtie ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import each FASTA sequence to @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# Define overall length of sequence:
$var_length_all				= length $array_FASTA[0][2];

# Store default length for $var_end if required:
$var_end					= $var_length_all unless $var_end;

# Store $var_end in %hash_variables:
$hash_variables { var_end }	= $var_end;

# ---------------------------------------------------------------------------- #

# Create Bowtie index unless it already exists:
unless ( -e "$var_output/$var_prefix.1.ebwt" ) {

	print "bowtie-build \\\n\t$var_FASTA \\\n\t$var_output/$var_prefix \\\n";
	print "\t--threads $var_threads \\\n\t--quiet\n";
	system ( "bowtie-build $var_FASTA $var_output/$var_prefix --threads $var_threads --quiet" );

}

# ---------------------------------------------------------------------------- #

# Run bowtie:
print "bowtie \\\n\t$var_output/$var_prefix \\\n";
print "\t$array_FASTQ[0] \\\n" if scalar @array_FASTQ == 1;
print "\t-1 $array_FASTQ[0] \\\n" if scalar @array_FASTQ == 2;
print "\t-2 $array_FASTQ[1] \\\n" if scalar @array_FASTQ == 2;
print "\t-n 0 \\\n\t-I 500 \\\n\t-X 1000 \\\n\t" if scalar @array_FASTQ == 2;
print "-k 1 \\\n\t-S $file_SAM \\\n\t--al $file_aligned \\\n\t-p 4\n";

if ( scalar @array_FASTQ == 1 ) {

	system ( "bowtie $var_output/$var_prefix $array_FASTQ[0] -n 0 -k 1 -S $file_SAM --al $file_aligned -p $var_threads" );

}

else {

	system ( "bowtie $var_output/$var_prefix -1 $array_FASTQ[0] -2 $array_FASTQ[1] -n 0 -I 500 -X 1000 -k 1 -S $file_SAM --al $file_aligned -p $var_threads" );

}

# ---------------------------------------------------------------------------- #

# Process SAM file:
my ( $array_coverage, $var_reads, $var_read_length, $var_total_reads )	= ProcessSAM ( $file_SAM );

# Generate coverage information:
$array_coverage		= CovInfo ( $array_coverage );

# Trim and scale coverage data
Scale ( $array_coverage, \%hash_variables );

# ---------------------------------------------------------------------------- #

# Evaluate read coverage and report:
my $var_length		= $var_end - $var_start + 1;
my $var_coverage	= $var_reads/$var_length*$var_read_length;
$var_coverage		= sprintf ( "%.2f", $var_coverage );

# Print number of reads, total length of sequence, and total read coverage:
print "\n* Total aligned reads:    $var_total_reads reads\n";
print "* Coverage aligned reads: $var_reads\n";
print "* Total coverage length:  $var_length bp\n";
print "* Average coverage:       $var_coverage reads/bp\n\n";

# ---------------------------------------------------------------------------- #

# Subroutine name:	ProcessSAM
# Description:		This subroutine stores extracts read data and stores
#					coverage information in a hash.

# ----------------------------------------------------------

sub ProcessSAM {

	# Arguments:
	my ( $file_SAM )	= @_;

	# Data structures:
#	my %hash_coverage;
	my @array_coverage;

	# Variables:
	my $var_reads		= 0;
	my $var_counter		= 3;
	my $var_read_length	= 0;
	my $var_total_reads	= 0;

	# ----------------------------------

	# Open the $file_SAM:
	open ( my $file_read, "<", "$file_SAM" ) or die "Unable to open $file_SAM!\n";

	# Iterate through each line:
	while (<$file_read>) {

		# Skip header:
		if ( $var_counter ) {

			$var_counter--;

			next;

		}

		my $var_line	= $_;

		chomp $var_line;

		# Split $var_line into an array:
		my @array_SAM	= split " ", $var_line;

		# Define variables holding information needed:
		my $var_ref		= $array_SAM[2];
		my $var_loc0	= $array_SAM[3];
		my $var_loc1	= $array_SAM[7];
		my $var_length	= $array_SAM[8];

		# Define $var_read_length:
		$var_read_length	= length $array_SAM[9] unless $var_read_length;

		# Move on if there is no alignment or it's the second read:
		next unless $var_loc0;
		next unless $var_length > 0;

		# Iterate $var_total_reads by one:
		$var_total_reads++;

		# Move on if region is outside of $var_start and $var_end:
		next unless $var_loc0 >= $var_start - $var_read_length;
		next unless $var_loc1 <= $var_end + $var_read_length;

		# Iterate $var_reads by one:
		$var_reads++;

		my @array_temp	= ( $var_loc0, $var_loc1, $var_length );

		push @array_coverage, \@array_temp;

	}

	# Close $file_SAM:
	close $file_read or die "Unable to close $file_SAM!\n";

	# ----------------------------------

	# Return reference to %hash_coverage and the variables $var_reads,
	# $var_read_length, and $var_total_reads:
	return \@array_coverage, $var_reads, $var_read_length, $var_total_reads;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	CovInfo
# Description:		This subroutine create coverage information using 

# ----------------------------------------------------------

sub CovInfo {

	# Arguments:
	my ( $array_coverage )	= @_;

	# Data structures:
	my @array_coverage	= @$array_coverage;
	my @array_out;

	# ----------------------------------

	# Iterate through alignments and create coverage map:
	for ( my $i = 0; $i < scalar @array_coverage; $i++ ) {

		# Define start and end locations:
		my $var_loc0	= $array_coverage[$i][0];
		my $var_loc1	= $array_coverage[$i][1];

		# Add coverage for each location:
		for ( my $j = $var_loc0; $j <= $var_loc1; $j++ ) { $array_out[$j]++; }

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Scale
# Description:		Scales down coverage data by a given factor.

# ----------------------------------------------------------

sub Scale {

	# Arguments:
	my ( $array_coverage, $hash_variables )	= @_;

	# Data structures:
	my @array_coverage	= @$array_coverage;
	my %hash_variables	= %$hash_variables;

	# Variables:
	my $var_output	= $hash_variables { var_output };
	my $var_prefix	= $hash_variables { var_prefix };
	my $var_start	= $hash_variables { var_start };
	my $var_end		= $hash_variables { var_end };
	my $var_scale	= $hash_variables { var_scale };
	my $file_out	= "$var_output/$var_prefix.out";

	# ----------------------------------

	# Trim coverage data into a temporary array:
	my @array_temp	= @array_coverage[$var_start..$var_end];

	# Define a second temporary array to hold scaled down sequence data:
	my @array_temp2;

	# Add file name to @array_temp2:
	push @array_temp2, "# $file_out";

	# Iterate through first temporary array by the number of elements defined by
	# $var_scale:
	for ( my $i = 0; $i < scalar @array_temp; $i += $var_scale ) {

		# Define the value of the element at hand:
		my $var_value	= $array_temp[$i];

		# Set $var_value equal to zero if empty:
		$var_value		= 0 unless $var_value;

		# Push $var_value to second temporary array:
		push @array_temp2, $var_value;

	}

	# Print output file name to user:
	print "$file_out\n";

	# Print data to $file_out:
	General::ArrayToFile ( \@array_temp2, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #
