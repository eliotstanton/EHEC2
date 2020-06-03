#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
use Synteny;
use Synteny2;
use Quantify;
use SVG;
use Circos;

# ---------------------------------------------------------------------------- #

# File name:		Synteny.pl
# Date created:		31 October, 2018
# Last modified:	14 January, 2020
# Created by:		Eliot Stanton

# Description:		This is a master script for calculating and visualising
#					synteny between closely related strains using Mauve and
#					Circos.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('cfm:o:pqs:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_circos		= $hash_variables{c};
my $var_min			= $hash_variables{m} || 100;
my $var_output		= $hash_variables{o} || "synteny";
my $var_pM			= $hash_variables{p};
my $var_quantify	= $hash_variables{q};

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nSynteny.pl [OPTIONS] [FASTA1,FASTA2,etc] [FEATURES1,FEATURES2,etc]
	-c Visualise data using Circos
	-m Minimum length for region alignment (default: 100)
	-o Output directory (default: synteny)
	-p Force progressiveMauve to run
	-q Quantify genomic feature data and synteny
	-s Output files prefix (default: default)

	Feature file format:
	seqID	feature_type	feature_name	start	end

	example:
	0	prophage	prophage0	103894	163432

\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Store $var_output in %hash_variables:
$hash_variables { var_minimum }	= $var_min;

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# If $var_FASTA or $var_features are missing, print $var_help and stop script:
unless ( scalar @ARGV == 2 ) {

	print "$var_help\n";

	print "\tFASTA and Feature files are required!\n" and exit;

}

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

# ---------------------------------------------------------------------------- #

# Check if Circos will run:
Check::Circos ( \%hash_variables ) if $var_circos;

# Check if ProgressiveMauve will run:
Check::ProgressiveMauve ( \%hash_variables);

# Check if files in @array_FASTA are present:
@array_FASTA	= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };

# Check if files in @array_features are present:
@array_features	= @{ Check::Files ( \@array_features, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import data from @array_features:
@array_features	= @{ General::ImportFiles ( \@array_features ) };

# Run progressiveMauve if $var_pM is flagged:
Synteny::progressiveMauve ( \%hash_variables, \@array_FASTA ) if $var_pM;

# Process $file_backbone created by progressiveMauve and split backbone and
# insert alignments into separate arrays:
%hash_variables	= %{Synteny::FormatAlignment ( \%hash_variables, \@array_FASTA )};

# Organise backbone regions and split transposed backbone regions into a
# separate array:
%hash_variables	= %{Synteny::OrganiseBackbone ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Add bracketing reference numbers to inserts and transposed regions, this
# subroutine all returns @array_backbone rearranged by sequence:
%hash_variables	= %{Synteny::BracketRefs ( \%hash_variables )};

# Split insert regions sitting within junctions of backbone regions into a
# separate array:
%hash_variables	= %{Synteny::JunctionInserts ( \%hash_variables)};

# Check for transposed regions in @array_inserts:
%hash_variables	= %{Synteny::TransposedInserts ( \%hash_variables )};

# Organise insert regions:
%hash_variables	= %{Synteny::OrganiseInserts ( \%hash_variables )};

%hash_variables = %{ FudgeInserts ( \%hash_variables )} if $hash_variables{f};

# Merge regions all together:
%hash_variables	= %{Synteny::MergeRegions ( \%hash_variables )};

# Check that all regions are present in array_merged:
Synteny::CheckRegions ( \%hash_variables );

# Assign adjusted global locations to all regions:
%hash_variables	= %{ Synteny::GlobalLocations ( \%hash_variables ) };

# ---------------------------------------------------------------------------- #

# If $var_quantify is flagged determine overlap between synteny and various
# genomic features:
if ( $var_quantify ) {

	# Calculate overlap between different genomic features and 
	# conserved/non-conserved regions:
	Quantify::Overlap ( \@array_features, \@array_FASTA, \%hash_variables );

	# Import FASTA files to @array_FASTA:
	@array_FASTA	= @{General::ImportFASTA ( \@array_FASTA )};

	# Write conserved and non-conserved sequences to file in FASTA format:
	Synteny::WriteSeq ( \%hash_variables, \@array_FASTA );

}

# Format genomic features:
%hash_variables	= %{ Synteny::FormatFeatures ( \@array_features, \%hash_variables ) };

# Store final alignment data as files:
Synteny::StoreFinal ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Create SVG data for syntenic regions:
%hash_variables	= %{SVG::RegionsSVG ( \%hash_variables )};

# Create SVG data for features:
%hash_variables	= %{SVG::FeaturesSVG ( \%hash_variables )};

# Create SVG data for ticks:
%hash_variables	= %{SVG::TicksSVG ( \%hash_variables )};

# Merge SVG data and print to file:
SVG::MergeSVG ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Run Circos visualise data if flagged:
if ( $var_circos ) {

	# Create Circos config file:
	Synteny2::Config ( \%hash_variables );

	# Create Circos Karyotype file:
	Synteny2::Karyotype ( \%hash_variables );

	# Create highlights files:
	Synteny2::Highlights ( \%hash_variables );

	# Create Circos Ideogram file:
	Circos::Ideogram ( \%hash_variables );

	# Create Circos Image file:
	Circos::Image ( \%hash_variables );

	# Create ticks files:
	Synteny2::Ticks ( \%hash_variables );

	# Run Circos:
	Circos::Run ( \%hash_variables );

}

# ---------------------------------------------------------------------------- #

sub FudgeInserts {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my %hash_inserts		= %{$hash_out{hash_inserts}};
	my %hash_junctions		= %{ $hash_out { hash_junctions }};

	my %hash_regions;

	# Variables:
	my $var_scalar			= scalar @array_backbone;

	# ----------------------------------

	# Iterate through %hash_inserts looking for regions in junctions:
	foreach my $var_bracketIDs	( sort { $a cmp $b } keys %hash_inserts ) {

#		print "$var_bracketIDs\n";

		if ( $var_bracketIDs eq "124 126" || $var_bracketIDs eq "114 126" ||
			$var_bracketIDs eq "113 114" || $var_bracketIDs eq "113 124" ||
			$var_bracketIDs eq "150 156" || $var_bracketIDs eq "96 111" ||
			$var_bracketIDs eq "119 124" || $var_bracketIDs eq "141 143" ||
			$var_bracketIDs eq "143 144" ) {

#			print "$var_bracketIDs\n";

			# Define @array_regions from %hash_inserts:
			my @array_regions	= @{ $hash_inserts { $var_bracketIDs } };

			$hash_junctions	{ "$var_bracketIDs" } = "1";

			# Transfer regions from %hash_inserts to %hash_regions:		
			for ( my $i = 0; $i < scalar @array_regions; $i++ ) {

#			print "\t$i\n";

				for ( my $j = 0; $j < scalar @{$array_regions[$i]}; $j++ ) {

					# Define initial location of each region:
					my $var_loc0		= $array_regions[$i][$j][1];

					# Move on if handling a blank region:
					next unless $var_loc0;

					# Define bracketIDs:
					my $var_bracket0	= $array_regions[$i][$j][9];
					my $var_bracket1	= $array_regions[$i][$j][10];

					next unless "$var_bracket0 $var_bracket1" eq "$var_bracketIDs";

					$array_regions[$i][$j][5]	= $array_regions[$i][$j][6];

#					print "\t\t@{$array_regions[$i][$j]}\n";

					# Add regions to %hash_regions:
					push( @{ $hash_regions { "$var_bracket0 $var_bracket1" } }, $array_regions[$i][$j] ); 

				}

			}

			undef $hash_inserts { $var_bracketIDs };

		}

#		unless ( $var_bracketIDs eq "141 143" ) {

#			print "$var_bracketIDs\n";

#			undef $hash_inserts { $var_bracketIDs };

#		}

	}

	# ----------------------------------

	# Iterate through @array_backbone and identify junctions for reintegrating
	# regions:
	for ( my $i = 0; $i < scalar @{$array_backbone[0]}; $i++ ) {

#		print "\t$i: \n";

		# Define temporary hash and array:
		my %hash_temp;
		my @array_temp;

		# Define array to hold organised regions:
		my @array_regions;

		# Iterate across regions looking at bracket IDs:
		for ( my $j = 0; $j < scalar @array_backbone; $j++ ) {

			# Define bracket IDs:
			my $var_bracket0	= $array_backbone[$j][$i][3];
			my $var_bracket1	= $array_backbone[$j][$i][10];

#			print "\t\t$j: @{$array_backbone[$j][$i]} - $var_bracket0 $var_bracket1\n";

			if ( $hash_junctions { "$var_bracket0 $var_bracket1" } ) {

#				print "\t\t$j: @{$array_backbone[$j][$i]} - $var_bracket0 $var_bracket1\n";

				$hash_temp { "$var_bracket0 $var_bracket1" }	= "";

			}

		}

		next unless %hash_temp;

		# Iterate through bracketIDs stored in %hash_temp and grab associated
		# regions:
		foreach my $var_bracketIDs ( keys %hash_temp ) {

#			print "\t* $var_bracketIDs\n";

			if ( $hash_regions { "$var_bracketIDs" } ) {

			# Define temporary array holding regions associated with bracketIDs:
			my @array_temp2	= @{ $hash_regions { "$var_bracketIDs" }};

				for ( my $j = 0; $j < scalar @array_temp2; $j++ ) {

#					print "\t\t$j: @{$array_temp2[$j]}\n";

					# Merge regions into @array_temp:
					push @array_temp, $array_temp2[$j];

				}

			}

		}

		# Sort regions in @array_temp by location:
		@array_temp	= sort { $a -> [1] <=> $b -> [1] } @array_temp;

		# Iterate through @array_temp and sort regions by sequence into
		# @array_regions:
		for ( my $j = 0; $j < scalar @array_temp; $j++ ) {

			# Define sequence ID for each region:
			my $var_seqID	= $array_temp[$j][0];

			push @{$array_regions[$var_seqID]}, $array_temp[$j];

#			print "\t\t$j: @{$array_temp[$j]}\n";

		}

		# Define variable for maximum length:
		my $var_length0		= 0;

		# Define variable for maximum number of regions:
		my $var_regions0	= 0;

		# Identify maximum length of regions in each sequence:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_length1	= 0;

			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

				my $var_length_temp	= $array_regions[$j][$k][5];

#				print "\t\t\t$k: @{$array_regions[$j][$k]}\n";

				$var_length1	+= $var_length_temp;

			}

			$var_length0 = $var_length1 if $var_length1 > $var_length0;

		}

#		print "$var_length0 - $var_scalar\n";

		# Add blank regions with extra length as required:
		for ( my $j = 0; $j < $var_scalar; $j++ ) {

#			print "\t\t$j:\n";

			unless ( $array_regions[$j][0] ) {

				my @array_temp	= ( $j, 0, 0, 0, 0, $var_length0, 0, 0, 0, 0, 0, 0, 0 );

				push @{$array_regions[$j]}, \@array_temp;

				next;

			}

			my $var_length1	= 0;

			# Determine length of regions in sequence:
			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

				my $var_length_temp	= $array_regions[$j][$k][5];

#				print "\t\t\t$k: @{$array_regions[$j][$k]}\n";

				$var_length1	+= $var_length_temp;

			}

#			print "\t- $var_length1\n";

			# Add region with extra length if required:
			if ( $var_length0 > $var_length1 ) {

				my $var_length_temp	= $var_length0 - $var_length1;

				my @array_temp	= ( $j, 0, 0, 0, 0, $var_length_temp, 0, 0, 0, 0, 0, 0, 0 );

#				print "\t@array_temp\n";

				push @{$array_regions[$j]}, \@array_temp;

			}

		}

		# Identify maximum number of regions in each sequence:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_regions1	= scalar @{$array_regions[$j]};

			$var_regions0		= $var_regions1 if $var_regions1 > $var_regions0;

		}

		# Add regions to each sequence as needed:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

#			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_regions1	= scalar @{$array_regions[$j]};

			if ( $var_regions1 < $var_regions0 ) {

				for ( my $k = scalar @{$array_regions[$j]}; $k <= $var_regions0; $k++ ) {

					my @array_temp	= ( $j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

					push @{$array_regions[$j]}, \@array_temp;

				}

			}

		}

		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

#			print "\t$j\n";

			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

#				print "\t\t$k: @{$array_regions[$j][$k]}\n";

			}

		}

		# Add data back to %hash_inserts:
		foreach my $var_bracketIDs ( keys %hash_temp ) {

#			print "$var_bracketIDs\n";

			$hash_inserts { "$var_bracketIDs" } = \@array_regions;

		}

	}

	# ----------------------------------

	$hash_out { hash_inserts }	= \%hash_inserts;

	# ----------------------------------

	# Retrun reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #


# ----------------------------------------------------------

sub FudgeInsertsOld {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my %hash_inserts		= %{$hash_out{hash_inserts}};
	my %hash_junctions		= %{ $hash_out { hash_junctions }};

	my %hash_regions;

	# Variables:
	my $var_scalar			= scalar @array_backbone;

	# ----------------------------------

	# Iterate through %hash_inserts looking for regions in junctions:
	foreach my $var_bracketIDs	( sort { $a cmp $b } keys %hash_inserts ) {

# 			$var_bracketIDs eq "119 124" ||
 #			$var_bracketIDs eq "150 156" ||
 #			$var_bracketIDs eq "170 173"
#						$var_bracketIDs eq "113 114" || $var_bracketIDs eq "113 124"

		if (

			$var_bracketIDs eq "124 126" || $var_bracketIDs eq "114 126"

		) {

			print "$var_bracketIDs\n";

			# Define @array_regions from %hash_inserts:
			my @array_regions	= @{ $hash_inserts { $var_bracketIDs } };

			# Transfer regions from %hash_inserts to %hash_regions:		
			for ( my $i = 0; $i < scalar @array_regions; $i++ ) {

#			print "\t$i\n";

				for ( my $j = 0; $j < scalar @{$array_regions[$i]}; $j++ ) {

					# Define initial location of each region:
					my $var_loc0		= $array_regions[$i][$j][1];

					# Move on if handling a blank region:
					next unless $var_loc0;

					# Define bracketIDs:
					my $var_bracket0	= $array_regions[$i][$j][9];
					my $var_bracket1	= $array_regions[$i][$j][10];

					next unless "$var_bracket0 $var_bracket1" eq "$var_bracketIDs";

					$array_regions[$i][$j][5]	= $array_regions[$i][$j][6];

#					print "\t\t@{$array_regions[$i][$j]}\n";

					# Add regions to %hash_regions:
					push( @{ $hash_regions { "$var_bracket0 $var_bracket1" } }, $array_regions[$i][$j] ); 

				}

			}

			undef $hash_inserts { $var_bracketIDs };

		}

		else {

#			unless ( 
#					$var_bracketIDs eq "209 212" ||
#					$var_bracketIDs eq "207 208" ||
#					$var_bracketIDs eq "206 207" ||
#					$var_bracketIDs eq "200 202" ||
#					$var_bracketIDs eq "195 197" ||
#					$var_bracketIDs eq "193 195" ||
#					$var_bracketIDs eq "189 191" ||
#					$var_bracketIDs eq "185 187" ||
#					$var_bracketIDs eq "181 183" ||
#					$var_bracketIDs eq "179 181" ||
#					$var_bracketIDs eq "168 170" ||
#					$var_bracketIDs eq "166 168" ||
#					$var_bracketIDs eq "158 160" ||
#					$var_bracketIDs eq "147 150" ||
#					$var_bracketIDs eq "145 147" ||
#					$var_bracketIDs eq "143 144" ||
#					$var_bracketIDs eq "141 143" ||
#					$var_bracketIDs eq "140 141" ||
#					$var_bracketIDs eq "138 140" ||
#					$var_bracketIDs eq "134 136" ||
#					$var_bracketIDs eq "132 134" ||
#					$var_bracketIDs eq "130 132" ||
#					$var_bracketIDs eq "128 130" ||
#					$var_bracketIDs eq "126 128" ||

#					$var_bracketIDs eq "116 118" ||
#					$var_bracketIDs eq "115 116" ||
#					$var_bracketIDs eq "114 115" ||

#					$var_bracketIDs eq "111 113" ||


#					$var_bracketIDs eq "96 111" ||
#					$var_bracketIDs eq "91 95" ||
#					$var_bracketIDs eq "89 91" ||
#					$var_bracketIDs eq "87 89" ||
#					$var_bracketIDs eq "82 87" ||
#					$var_bracketIDs eq "81 82" ||
#					$var_bracketIDs eq "79 81" ||
#					$var_bracketIDs eq "77 79" ||
#					$var_bracketIDs eq "75 77" ||
#					$var_bracketIDs eq "73 75" ||
#					$var_bracketIDs eq "71 73" ||
#					$var_bracketIDs eq "69 71" ||
#					$var_bracketIDs eq "54 67" ||
#					$var_bracketIDs eq "50 52" ||
#					$var_bracketIDs eq "54 67" ||
#					$var_bracketIDs eq "47 49" ||
#					$var_bracketIDs eq "45 47" ||
#					$var_bracketIDs eq "38 41" ||
#					$var_bracketIDs eq "35 38" ||
#					$var_bracketIDs eq "54 67"
#					) {

#				print "$var_bracketIDs\n";

				undef $hash_inserts { $var_bracketIDs };

#			}

		}

#		undef $hash_inserts { $var_bracketIDs } if $var_bracketIDs eq "114 126";
#		undef $hash_inserts { $var_bracketIDs } if $var_bracketIDs eq "119 126";

	}

	# ----------------------------------

	# Iterate through @array_backbone and identify junctions for reintegrating
	# regions:
	for ( my $i = 0; $i < scalar @{$array_backbone[0]}; $i++ ) {

		print "\t$i: \n";

		# Define temporary hash and array:
		my %hash_temp;
		my @array_temp;

		# Define array to hold organised regions:
		my @array_regions;

		# Iterate across regions looking at bracket IDs:
		for ( my $j = 0; $j < scalar @array_backbone; $j++ ) {

			# Define bracket IDs:
			my $var_bracket0	= $array_backbone[$j][$i][3];
			my $var_bracket1	= $array_backbone[$j][$i][10];

			if ( $hash_junctions { "$var_bracket0 $var_bracket1" } ) {

				print "\t\t$j: @{$array_backbone[$j][$i]} - $var_bracket0 $var_bracket1\n";

				$hash_temp { "$var_bracket0 $var_bracket1" }	= "";

			}

		}

		next unless %hash_temp;

		# Iterate through bracketIDs stored in %hash_temp and grab associated
		# regions:
		foreach my $var_bracketIDs ( keys %hash_temp ) {

			print "\t* $var_bracketIDs\n";

			# Define temporary array holding regions associated with bracketIDs:
			my @array_temp2	= @{ $hash_regions { "$var_bracketIDs" }};

			for ( my $j = 0; $j < scalar @array_temp2; $j++ ) {

#				print "\t\t$j: @{$array_temp2[$j]}\n";

				# Merge regions into @array_temp:
				push @array_temp, $array_temp2[$j];

			}

		}

		# Sort regions in @array_temp by location:
		@array_temp	= sort { $a -> [1] <=> $b -> [1] } @array_temp;

		# Iterate through @array_temp and sort regions by sequence into
		# @array_regions:
		for ( my $j = 0; $j < scalar @array_temp; $j++ ) {

			# Define sequence ID for each region:
			my $var_seqID	= $array_temp[$j][0];

			push @{$array_regions[$var_seqID]}, $array_temp[$j];

#			print "\t\t$j: @{$array_temp[$j]}\n";

		}

		# Define variable for maximum length:
		my $var_length0		= 0;

		# Define variable for maximum number of regions:
		my $var_regions0	= 0;

		# Identify maximum length of regions in each sequence:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_length1	= 0;

			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

				my $var_length_temp	= $array_regions[$j][$k][5];

#				print "\t\t\t$k: @{$array_regions[$j][$k]}\n";

				$var_length1	+= $var_length_temp;

			}

			$var_length0 = $var_length1 if $var_length1 > $var_length0;

		}

#		print "$var_length0 - $var_scalar\n";

		# Add blank regions with extra length as required:
		for ( my $j = 0; $j < $var_scalar; $j++ ) {

#			print "\t\t$j:\n";

			unless ( $array_regions[$j][0] ) {

				my @array_temp	= ( $j, 0, 0, 0, 0, $var_length0, 0, 0, 0, 0, 0, 0, 0 );

				push @{$array_regions[$j]}, \@array_temp;

				next;

			}

			my $var_length1	= 0;

			# Determine length of regions in sequence:
			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

				my $var_length_temp	= $array_regions[$j][$k][5];

#				print "\t\t\t$k: @{$array_regions[$j][$k]}\n";

				$var_length1	+= $var_length_temp;

			}

#			print "\t- $var_length1\n";

			# Add region with extra length if required:
			if ( $var_length0 > $var_length1 ) {

				my $var_length_temp	= $var_length0 - $var_length1;

				my @array_temp	= ( $j, 0, 0, 0, 0, $var_length_temp, 0, 0, 0, 0, 0, 0, 0 );

#				print "\t@array_temp\n";

				push @{$array_regions[$j]}, \@array_temp;

			}

		}

		# Identify maximum number of regions in each sequence:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_regions1	= scalar @{$array_regions[$j]};

			$var_regions0		= $var_regions1 if $var_regions1 > $var_regions0;

		}

		# Add regions to each sequence as needed:
		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

#			next unless $array_regions[$j];

#			print "\t\t$j:\n";

			my $var_regions1	= scalar @{$array_regions[$j]};

			if ( $var_regions1 < $var_regions0 ) {

				for ( my $k = scalar @{$array_regions[$j]}; $k <= $var_regions0; $k++ ) {

					my @array_temp	= ( $j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

					push @{$array_regions[$j]}, \@array_temp;

				}

			}

		}

		for ( my $j = 0; $j < scalar @array_regions; $j++ ) {

			print "\t$j\n";

			for ( my $k = 0; $k < scalar @{$array_regions[$j]}; $k++ ) {

				print "\t\t$k: @{$array_regions[$j][$k]}\n";

			}

		}


		# Add data back to %hash_inserts:
		foreach my $var_bracketIDs ( keys %hash_temp ) {

#			print "$var_bracketIDs\n";

			my ( $var_bracket0, $var_bracket1 )	= split " ", $var_bracketIDs;

			print "$var_bracket0 $var_bracket1\n";

			$hash_inserts { "$var_bracket0 $var_bracket1" } = \@array_regions;
#			$hash_inserts { "$var_bracket1 $var_bracket0" } = \@array_regions;

		}

	}

exit;

	# ----------------------------------

	$hash_out { hash_inserts }	= \%hash_inserts;

	# ----------------------------------

	# Retrun reference to %hash_out and end subroutine:
	return \%hash_out;

}

