#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
use Homology;
use SVG;

# ---------------------------------------------------------------------------- #

# File name:		HomologyLite.pl
# Date created:		21 November, 2018
# Last modified:	20 January, 2020
# Created by:		Eliot Stanton

# Description:		This is a script for determining homology shared between
#					one or more short FASTA sequences.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('m:n:o:s:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_min			= $hash_variables{m} || 100;
my $var_nmer		= $hash_variables{n} || 20;
my $var_output		= $hash_variables{o} || "lite";
my $var_prefix		= $hash_variables{s} || "default";
my $var_scale		= 10;
my $var_height		= 100;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nHomologyLite.pl [OPTIONS] [FASTA]
	-m Minimum repeat length (default: 100)
	-n nmer length for homology (default: 20 bp)
	-o OUTPUT DIRECTORY (default: lite)
	-s Output files prefix (default: default)\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output or $var_FASTA is missing, print $var_help and stop script:
unless ( $var_FASTA ) {

	print "$var_help\n";

	print "\tFASTA file(s) required!\n";

	exit;

}

# ---------------------------------------------------------------------------- #

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Store $var_scale in %hash_variables:
$hash_variables { var_scale }	= $var_scale;

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# Import data from @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# ---------------------------------------------------------------------------- #

# Generate a hash of each nmer present in all the FASTA sequences:
my $hash_nmer		= Homology::Nmer ( \%hash_variables, \@array_FASTA, $var_nmer );

# Create and save formatted link data:
Homology::Links ( \%hash_variables, $hash_nmer, $var_nmer );

# Merge direct and inverted links whose locations are off by one:
my @array_direct	= @{Homology::MergeLinks ( \%hash_variables, "direct", $var_min )};
my @array_inverted	= @{Homology::MergeLinks ( \%hash_variables, "inverted", $var_min )};

# ---------------------------------------------------------------------------- #

if ( @array_direct ) {

	print "Direct Links:\n";

	# Print data to command line:
	PrintData ( \@array_direct );

	# Format data for SVG:
	@array_FASTA	= @{FormatData ( \@array_direct, \@array_FASTA, \%hash_variables, "Direct" )};

}

if ( @array_inverted ) {

	print "Inverted Links:\n";

	# Print data to command line:
	PrintData ( \@array_inverted );

	# Format data for SVG:
#	@array_FASTA	= @{FormatData ( \@array_inverted, \@array_FASTA, \%hash_variables, "Inverted" )};

}

# Complete formatting for SVG file:
FormatSVG ( \@array_FASTA );

# ---------------------------------------------------------------------------- #

# Subroutine name:	PrintData
# Description:		This subroutine prints data to the command line:

# ----------------------------------------------------------

sub PrintData {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;

	# ----------------------------------

	# Iterate through regions in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		if ( ref $array_in[$i] ) {

			my $var_loc0	= $array_in[$i][1];
			my $var_loc1	= $array_in[$i][2];

			my $var_length	= $var_loc1 - $var_loc0 + 1;

			print "\t$i: @{$array_in[$i]} - $var_length bp\n" 

		}

		else {

			my @array_temp	= split " ", $array_in[$i];

			my $var_loc0	= $array_temp[1];
			my $var_loc1	= $array_temp[2];

			my $var_length	= $var_loc1 - $var_loc0 + 1;

			print "\t$i: @array_temp - $var_length bp\n";

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatData
# Description:		This subroutine formats data for export as SVG object.

# ----------------------------------------------------------

sub FormatData {

	# Arguments:
	my ( $array_in, $array_out, $hash_in, $var_variable )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out	= @$array_out;
	my %hash_in		= %$hash_in;

	# Variables:
	my $var_scale	= $hash_in { var_scale };

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my ( $var_seq0, $var_loc0, $var_loc1 );
		my ( $var_seq1, $var_loc2, $var_loc3 );

		if ( ref $array_in[$i] ) {

			$var_seq0	= $array_in[$i][0];
			$var_loc0	= $array_in[$i][1];
			$var_loc1	= $array_in[$i][2];
			$var_seq1	= $array_in[$i][3];
			$var_loc2	= $array_in[$i][4];
			$var_loc3	= $array_in[$i][5];

		}

		else {

			my @array_temp	= split " ", $array_in[$i];

			$var_seq0	= $array_temp[0];
			$var_loc0	= $array_temp[1];
			$var_loc1	= $array_temp[2];
			$var_seq1	= $array_temp[3];
			$var_loc2	= $array_temp[4];
			$var_loc3	= $array_temp[5];

		}

		# Define length of region:
		my $var_length	= $var_loc1 - $var_loc0;

		# Scale variables:
		$var_loc0		/= $var_scale;
		$var_loc2		/= $var_scale;
		$var_length		/= $var_scale;

		# Create SVG rectangles:
		my $var_string0;
		my $var_string1;

		if ( $var_variable eq "Direct" ) {

			$var_string0	= SVG::Rectangle ( $var_loc0, 0, $var_length, 50, "(100,100,100)", 0 );
			$var_string1	= SVG::Rectangle ( $var_loc2, 0, $var_length, 50, "(100,100,100)", 0 );

			push @{$array_out[$var_seq0][3]}, $var_string0;
			push @{$array_out[$var_seq1][3]}, $var_string1;

		}

		if ( $var_variable eq "Inverted" ) {

			$var_string0	= SVG::Rectangle ( $var_loc0, 50, $var_length, 50, "(200,200,200)", 0 );
			$var_string1	= SVG::Rectangle ( $var_loc2, 50, $var_length, 50, "(200,200,200)", 0 );

			push @{$array_out[$var_seq0][4]}, $var_string0;
			push @{$array_out[$var_seq1][4]}, $var_string1;

		}


	}
	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatSVG
# Description:		This subroutine complete formatting for SVG file.

# ----------------------------------------------------------

sub FormatSVG {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Iterate through @array_FASTA and print SVG data to file:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define file path of FASTA sequence:
		my $var_file	= $array_in[$i][0];

		# Define length of FASTA sequence:
		my $var_length	= length $array_in[$i][2];

		# Define path for file to be printed:
		my $file_out	= "$var_output/$var_prefix.$i.svg";

		print "$i: $var_file -> $file_out\n";

		# Define an array to hold information to be printed:
		my @array_out;

		# Define temporary arrays holding direct and inverted repeats:
		my @array_direct 	= @{$array_in[$i][3]} if $array_in[$i][3];
		my @array_inverted	= @{$array_in[$i][4]} if $array_in[$i][4];

		# Define the width of image using length of FASTA sequence:
		my $var_width		= $var_length / $var_scale;

		# Add header information:
		push @array_out, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

		# Add background:
		push @array_out, "<g id=\"bg\">";
		my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
		$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";
		push @array_out, $var_string;
		push @array_out, "</g>";

		# ------------------------------

		# Print direct repeats:
		for my $var_string ( @array_direct ) {

			push @array_out, $var_string;

		}

		# print inverted repeats
		for my $var_string ( @array_inverted ) {

			push @array_out, $var_string;

		}

		# ----------------------------------

		push @array_out, "</svg>";

		# Print data to file:
		General::ArrayToFile ( \@array_out, $file_out );

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #
