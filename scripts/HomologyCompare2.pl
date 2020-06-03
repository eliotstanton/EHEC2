#!/usr/bin/perl

use strict;
use warnings;

use lib '/etc/perl/modules';
use JSON::XS;

#use Getopt::Std;
#use lib '/etc/perl/modules';

use General;
use Synteny;
use Check;
use Homology;
#use Quantify;
#use Circos;

# ---------------------------------------------------------------------------- #

# File name:		HomologyCompare2.pl
# Date created:		02 June, 2020
# Last modified:	03 June, 2020
# Created by:		Eliot Stanton

# Description:		This is a script for comparing repeat abundance stored in 
#					two hashes.

# ---------------------------------------------------------------------------- #

# Variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_nmer		= "20";
my $var_output		= "output";
my $var_min			= "20";
my $var_help		= "";

# Data structures used by this script:
my %hash_variables;
my @array_backbone;

# ---------------------------------------------------------------------------- #

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Store $var_nmer in %hash_variables:
$hash_variables { n }	= $var_nmer;

# Split $var_FASTA into @array_features:
my @array_FASTA		= split /\,/, $var_FASTA;
my @array_FASTA2	= split /\,/, $var_FASTA;

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

# ---------------------------------------------------------------------------- #

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# Add variables from %hash_variables:
%hash_variables		= %{ General::Variables ( \%hash_variables )};

# Perform alignment using ProgressiveMauve:
#Synteny::progressiveMauve ( \%hash_variables, \@array_FASTA );

my $file_backbone	= $hash_variables { file_backbone };

# Import data from $file_backbone to @array_backbone
@array_backbone		= @{ General::FileToArray ( $file_backbone ) };

splice @array_backbone, 0, 1;

# Remove contiguous or non-contiguous blocks from the array:
for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

	my $var_loc0	= $array_backbone[$i][0];
	my $var_loc1	= $array_backbone[$i][1];
	my $var_loc2	= $array_backbone[$i][2];
	my $var_loc3	= $array_backbone[$i][3];

	# Remove regions that are present in both genomes:
	if ( $var_loc0 && $var_loc1 && $var_loc2 && $var_loc3 ) {

		splice @array_backbone, $i, 1;

		$i--;

		next;

	}

	$var_loc0 *= -1 if $var_loc0 < 0;
	$var_loc1 *= -1 if $var_loc1 < 0;
	$var_loc2 *= -1 if $var_loc2 < 0;		
	$var_loc3 *= -1 if $var_loc3 < 0;

	$array_backbone[$i][0]	= $var_loc0;
	$array_backbone[$i][1]	= $var_loc1;
	$array_backbone[$i][2]	= $var_loc2;
	$array_backbone[$i][3]	= $var_loc3;

#	print "$i: @{$array_backbone[$i]}\n";

}

# ---------------------------------------------------------------------------- #

# Import data from @array_FASTA:
@array_FASTA					= @{ General::ImportFASTA ( \@array_FASTA, "" ) };

my @array_temp					= ( $array_FASTA[0] );

my $file_hash					= $hash_variables { file_hash };

# Generate a hash of each nmer present in all the FASTA sequences:
my %hash_nmer0		= %{ Homology::Nmer ( \%hash_variables, \@array_temp, $var_nmer )};

my $var_scalar = scalar keys %hash_nmer0;
print "$var_scalar\n";

%hash_nmer0			= %{ Remove ( \%hash_nmer0, \@array_backbone, "0" ) };

$var_scalar = scalar keys %hash_nmer0;
print "$var_scalar\n";

SaveHash (\%hash_nmer0, "$file_hash" );

system ( "mkdir $var_output/0" ) unless -e "$var_output/0";

system ( "mv $file_hash $var_output/0" );

system ( "~/Desktop/Tools/scripts/HomologyAnalyzer.pl -q -c -o $var_output/0 -n $var_nmer -m $var_min $array_FASTA2[0] $array_features[0]" );

# ---------------------------------------------------------------------------- #

@array_temp						= ( $array_FASTA[1] );

$hash_variables { file_hash }	= "$file_hash";

my %hash_nmer1					= %{Homology::Nmer ( \%hash_variables, \@array_temp, $var_nmer )};

$var_scalar 					= scalar keys %hash_nmer1;
print "$var_scalar\n";

%hash_nmer1						= %{ Remove ( \%hash_nmer1, \@array_backbone, "2" ) };

$var_scalar 					= scalar keys %hash_nmer1;
print "$var_scalar\n";

SaveHash (\%hash_nmer1, "$file_hash" );

system ( "mkdir $var_output/1" ) unless -e "$var_output/1";

system ( "mv $file_hash $var_output/1" );

system ( "~/Desktop/Tools/scripts/HomologyAnalyzer.pl -q -c -o $var_output/1 -n $var_nmer -m $var_min $array_FASTA2[1] $array_features[1]" );

# ---------------------------------------------------------------------------- #

sub Remove {

	# Arguments:
	my ( $hash_in, $array_backbone, $var_increment )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_backbone	= @$array_backbone;

	# ----------------------------------

	@array_backbone		= sort { $a -> [0+$var_increment] <=> $b -> [0+$var_increment] } @array_backbone;

	# Iterate through hash and remove locations that are non-conserved:
	foreach my $var_key ( keys %hash_in ) {

#		print "$var_key\n";

		my @array_temp	= @{$hash_in { $var_key }};

		# Iterate through locations in @array_temp:
		for ( my $i = 0; $i < scalar @array_temp; $i++ ) {

			my $var_loc0 	= $array_temp[$i][1];
			my $var_loc1	= $var_loc0 + $var_nmer;

#			print "\t$i: $var_loc0 $var_loc1\n";

			# Identify and remove locations in non-conserved areas:
			for ( my $j = 1; $j < scalar @array_backbone; $j++ ) {

				my $var_loc2	= $array_backbone[$j][0+$var_increment];
				my $var_loc3	= $array_backbone[$j][1+$var_increment];

				if ( $var_loc0 >= $var_loc2 && $var_loc1 <= $var_loc3 ) {

					splice @array_temp, $i, 1;

					$i--;

#					print "\t\t$j: $var_loc2 $var_loc3\n";

					last;

				}

				last if $var_loc2 > $var_loc1;

			}

		}

		# Store @array_temp back in %hash_in:
		$hash_in { $var_key } = \@array_temp if scalar @array_temp > 0;

		delete $hash_in { $var_key } if scalar @array_temp == 0;

	}

	# ----------------------------------

	# End subroutine:
	return \%hash_in;

}

# ---------------------------------------------------------------------------- #

sub SaveHash {

	# Arguments:
	my ( $hash_out, $file_out )	= @_;

	# VData structures:
	my %hash_out	= %$hash_out;

	# ----------------------------------

	# Encode %hash_out in JSON format:
	my $hash_json = encode_json( \%hash_out );

	# Open $file_hash to write data:
	open ( my $file_write, '>', $file_out ) or die " ERROR: UNABLE TO OPEN $file_out!\n";

	# Print the hash to file.
	print $file_write $hash_json;

	# Close $file_out:
	close $file_write or die " ERROR: UNABLE TO CLOSE $file_out!\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

#sub LoadHash {

	# Argument:
#	my ( $file_in )	= @_;

	# Data-structures:
#	my %hash_out;

	# ----------------------------------

#	print "$file_in\n";

	# Open file containing stored hash.
#	open ( my $file_read, '<', $file_in ) or die " ERROR: UNABLE TO OPEN $file_in!\n";

#	my $hash_json = <$file_read>;

#	%hash_out = %{decode_json( $hash_json )};

#	close $file_read or die " ERROR: UNABLE TO CLOSE $file_in!\n";

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
#	return \%hash_out;

#}

# ---------------------------------------------------------------------------- #
