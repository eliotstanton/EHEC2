#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
use Homology;
use Quantify;
use Circos;

# ---------------------------------------------------------------------------- #

# File name:		HomologyAnalyzer.pl
# Date created:		27 October, 2018
# Last modified:	12 January, 2020
# Created by:		Eliot Stanton

# Description:		This is a master script for calculating and visualising
#					homology within a circular bacterial genome.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('cdim:n:o:qs:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_circos		= $hash_variables{c};
my $var_min			= $hash_variables{m} || 1000;
my $var_nmer		= $hash_variables{n} || 100;
my $var_output		= $hash_variables{o} || "homology";
my $var_quantify	= $hash_variables{q};
my $var_prefix		= $hash_variables{s} || "default";

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nHomologyAnalyzer.pl [FASTA] [Features]
	-d Prohibit direct links from being drawn
	-c Visualise data using Circos
	-i Prohibit inverted links from being drawn
	-m Minimum merged repeat length (default: 1000)
	-n nmer length for establishing homology (default: 100 bp)
	-o OUTPUT DIRECTORY (default: homology)
	-q Quantify homology data using genmic features
	-s Output files prefix (default: default)

	Feature file format:
	seqID	feature_type	feature_name	start	end

	example:
	0	prophage	prophage0	103894	163432\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Store $var_nmer in %hash_variables:
$hash_variables { var_nmer }	= $var_nmer;

# Store $var_output in %hash_variables:
$hash_variables { var_minimum }	= $var_min;

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# --------------------------------------

# Print out variables in %hash_variables:
foreach my $var_variable ( keys %hash_variables ) {

#	print "$var_variable: $hash_variables{$var_variable}\n";

}

# ---------------------------------------------------------------------------- #

# If $var_FASTA or $var_features are missing, print $var_help and stop script:
unless ( scalar @ARGV == 2 ) {

	print "$var_help\n";
	
	print "\tFASTA and feature files are required!\n" and exit;

}

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

if ( scalar @array_FASTA > 1 || scalar @array_features > 1 ) {

	print "$var_help\n";
	
	print "\tOnly one FASTA and Feature file allowed!\n" and exit;

}

# ---------------------------------------------------------------------------- #

# Check if Circos will run:
Check::Circos ( \%hash_variables ) if $var_circos;

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };

# Check if files in @array_features are present:
@array_features		= @{ Check::Files ( \@array_features, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import data from @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# Import data from @array_features:
@array_features		= @{ General::FileToArray ( $array_features[0] ) };

# Add size information to @array_features:
@array_features		= @{ Homology::ProcessFeatures ( \@array_features ) };

# --------------------------------------

# Store length of the FASTA sequence in %hash_variables:
my $var_chr_length					= length $array_FASTA[0][2];
$hash_variables { var_chr_length }	= $var_chr_length;

# ---------------------------------------------------------------------------- #

# Generate a hash of each nmer present in all the FASTA sequences:
my $hash_nmer		= Homology::Nmer ( \%hash_variables, \@array_FASTA, $var_nmer );

if ( $var_quantify ) {

	# Quantify distribution of repetitiveness:
	Quantify::CopyNumber ( \%hash_variables, $hash_nmer );

	# Quantify overlap between raw repeats and genomic features:
	Quantify::OverlapRaw ( \%hash_variables, $hash_nmer, \@array_features );

}

# Create and save formatted link data:
Homology::Links ( \%hash_variables, $hash_nmer, $var_nmer );

# Add genomic feature data to each link:
Homology::FeatureData ( \%hash_variables, "direct", \@array_features );
Homology::FeatureData ( \%hash_variables, "inverted", \@array_features );

if ( $var_quantify ) {

	# Quantify distance between links:
#	RawDistance ( \%hash_variables, \@array_FASTA, "direct" );

	# Quantify overlap between raw repeats and genomic features:
	Quantify::RawLinks ( \%hash_variables, "direct" );
	Quantify::RawLinks ( \%hash_variables, "inverted" );

}

# Merge direct and inverted links whose locations are off by one:
Homology::MergeLinks (  \%hash_variables, "direct", $var_min );
Homology::MergeLinks ( \%hash_variables, "inverted", $var_min );

# Flatten direct and inverted links into repeat regions:
Homology::FlattenLinks ( \%hash_variables, "direct" );
Homology::FlattenLinks ( \%hash_variables, "inverted" );
Homology::FlattenLinks ( \%hash_variables, "combined" );

# ---------------------------------------------------------------------------- #

# Quantify data if $var_quantify is defined from command-line:
if ( $var_quantify ) {

#	Distance ( \%hash_variables, \@array_FASTA, "direct" );

	# Quantify distribution of merged direct and inverted repeats:
	Quantify::Distribution ( \%hash_variables, "direct" );
	Quantify::Distribution ( \%hash_variables, "inverted" );

	print "\n    + Quantification of repetitive regions:\n";

	# Quantify overlap of genomic features with flattened repeat regions:
	Quantify::Regions ( \%hash_variables, "direct", \@array_features );
	Quantify::Regions ( \%hash_variables, "inverted", \@array_features );
	Quantify::Regions ( \%hash_variables, "combined", \@array_features );

}

# ---------------------------------------------------------------------------- #

# Visualise data using Circos if flagged:
if ( $var_circos ) {

	# Create Circos config file:
	Homology::CircosConfig ( \%hash_variables );

	# Create Circos links files:
	Homology::CircosLinks ( "direct", \@array_features, \%hash_variables );
	Homology::CircosLinks ( "inverted", \@array_features, \%hash_variables );

	# Create Circos Karyotype file:
	Homology::CircosKaryotype ( \@array_FASTA, \%hash_variables );

	# Create Circos Ideogram file:
	Circos::Ideogram ( \%hash_variables );

	# Create Circos highlights features:
	Homology::CircosFeatures ( \@array_features, \@array_FASTA, \%hash_variables );

	# Create ticks file:
	Homology::CircosTicks ( \@array_FASTA, \%hash_variables );

	# Create Circos Image file:
	Circos::Image ( \%hash_variables );

	# Run Circos:
	Circos::Run ( \%hash_variables );

}

# ---------------------------------------------------------------------------- #

sub Distance {

	# Arguments:
	my ( $hash_in, $array_FASTA, $var_variable )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_FASTA	= @$array_FASTA;
	my @array_in;
	my @array_out;

	# Variables:
	my $file_in;
	my $file_out;

	# ----------------------------------

	# Define $file_in based upon $var_variable:
	if ( $var_variable eq "direct" ) {

		$file_in	= $hash_out { file_direct_merged };

	}

	if ( $var_variable eq "inverted" ) {

		$file_in	= $hash_out { file_inverted_merged };

	}

	# Import $file_in to @array_in:
	@array_in	= @{ General::FileToArray( $file_in ) };

	# ----------------------------------

	# Store length of sequences:
	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

		# Define variable to hold sequence:
		my $var_sequence	= $array_FASTA[$i][2];

		# Define variable to hold length of sequence:
		my $var_length		= length $var_sequence;

		$array_FASTA[$i][3]	= $var_length;

	}

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define locations of each link:
		my $var_seq0	= $array_in[$i][0];
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_seq1	= $array_in[$i][3];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];	
		my $var_length	= $array_in[$i][6];		

		my $var_distance	= $var_loc2 - $var_loc1;

		# Define variable to hold length of sequence:
		my $var_chr_length		= $array_FASTA[$var_seq0][3];

		# Calculate the shorter distance if needed:
		if ( $var_distance > $var_chr_length/2 ) {

			my $var_distance	= $var_chr_length - $var_loc3 + $var_loc0;

			print "$i: @{$array_in[$i]} - $var_length $var_distance\n";

		}

		$array_in[$i][7]	= $var_distance;

	}

	# ----------------------------------

	# Sort @array_in by distance:
#	@array_in	= sort { $a -> [6] <=> $b -> [6] } @array_in;

	# Save @array_in to $file_in:
	General::ArrayToFile ( \@array_in, $file_in );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

sub RawDistance {

	# Arguments:
	my ( $hash_in, $array_FASTA, $var_variable )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_FASTA	= @$array_FASTA;
	my @array_in;
	my @array_out;

	# Variables:
	my $file_in;
	my $file_out;

	# ----------------------------------

	# Define $file_in based upon $var_variable:
	if ( $var_variable eq "direct" ) {

		$file_in	= $hash_out { file_direct_features };

	}

	if ( $var_variable eq "inverted" ) {

		$file_in	= $hash_out { file_inverted_features };

	}

	# Import $file_in to @array_in:
	@array_in	= @{ General::FileToArray( $file_in ) };

	# ----------------------------------

	# Store length of sequences:
	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

		# Define variable to hold sequence:
		my $var_sequence	= $array_FASTA[$i][2];

		# Define variable to hold length of sequence:
		my $var_length		= length $var_sequence;

		$array_FASTA[$i][3]	= $var_length;

	}

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define locations of each link:
		my $var_seq0	= $array_in[$i][0];
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_seq1	= $array_in[$i][3];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];		

		my $var_distance	= $var_loc2 - $var_loc1;

		# Define variable to hold length of sequence:
		my $var_length		= $array_FASTA[$var_seq0][3];

		if ( $var_distance > $var_length/2 ) {

			my $var_distance	= $var_length - $var_loc3 + $var_loc0;

#			print "$i: @{$array_in[$i]} - $var_distance - $var_length\n";

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}
