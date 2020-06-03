#!/usr/bin/perl

use strict;
use warnings;
use JSON::XS;

#use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
use Homology;
use Quantify;
#use Circos;

# ---------------------------------------------------------------------------- #

# File name:		HomologyCompare.pl
# Date created:		02 June, 2020
# Last modified:	02 June, 2020
# Created by:		Eliot Stanton

# Description:		This is a script for comparing repeat abundance stored in 
#					two hashes.

# ---------------------------------------------------------------------------- #

# Define variables used by this script:
my $var_hash		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_nmer		= "20";
my $var_output		= "output";
my $var_stem		= "default";
my $var_min			= "20";

# Define data-structures used by this script:
my %hash_0;
my %hash_1;
my %hash_variables;
my %hash_new_0;
my %hash_new_1;

# Store variables in %hash_variables:
$hash_variables { var_nmer }				= $var_nmer;
$hash_variables { var_min }					= $var_min;
$hash_variables { file_direct }				= "$var_output/$var_stem.0.$var_nmer.direct";
$hash_variables { file_inverted }			= "$var_output/$var_stem.0.$var_nmer.inverted";
$hash_variables { file_direct_links }		= "$var_output/$var_stem.0.direct.$var_nmer.links";
$hash_variables { file_inverted_links }		= "$var_output/$var_stem.0.inverted.$var_nmer.links";
$hash_variables { file_direct_features }	= "$var_output/$var_stem.direct.$var_nmer.0.features";
$hash_variables { file_inverted_features }	= "$var_output/$var_stem.inverted.$var_nmer.0.features";
$hash_variables { file_direct_merged }		= "$var_output/$var_stem.direct.$var_nmer.0.merged";
$hash_variables { file_inverted_merged }	= "$var_output/$var_stem.inverted.$var_nmer.0.merged";
$hash_variables { file_direct_flat }		= "$var_output/$var_stem.direct.$var_nmer.0.flat";
$hash_variables { file_inverted_flat }		= "$var_output/$var_stem.inverted.$var_nmer.0.flat";
$hash_variables { file_flat }				= "$var_output/$var_stem.$var_nmer.0.flat";

# ---------------------------------------------------------------------------- #

# Split $var_hash into @array_hash:
my @array_hash		= split /\,/, $var_hash;

# Define filepaths for hashes:
my $file_hash0		= $array_hash[0];
my $file_hash1		= $array_hash[1];

# Load hashes into memory:
%hash_0	= %{ LoadHash ( $file_hash0 )};
%hash_1	= %{ LoadHash ( $file_hash1 )};

# ---------------------------------------------------------------------------- #

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

# Import data from @array_features:
@array_features		= @{ General::FileToArray ( $array_features[0] ) };

# Add size information to @array_features:
@array_features		= @{ Homology::ProcessFeatures ( \@array_features ) };

# --------------------------------------

print "\n- $file_hash0\n";

# Iterate through %hash_0 and store repeats also present in %hash_1 in 
# %hash_new_0:
foreach my $var_key ( keys %hash_0 ) {

	if ( $hash_1 { $var_key } ) {

		my $var_scalar	= scalar @{ $hash_1 { $var_key } };

		$hash_new_0 { $var_key }	= $hash_0 { $var_key };

	}

}

my $var_scalar	= scalar keys %hash_0;
print "\t- Pre-processed unique repeats: $var_scalar\n";

$var_scalar	= scalar keys %hash_new_0;
print "\t- Post-processed unique repeats: $var_scalar\n";

# --------------------------------------

# Quantify overlap between raw repeats and genomic features:
Quantify::OverlapRaw ( \%hash_variables, \%hash_new_0, \@array_features );

# ---------------------------------------------------------------------------- #

# Split $var_features into @array_features:
@array_features	= split /\,/, $var_features;

# Import data from @array_features:
@array_features		= @{ General::FileToArray ( $array_features[1] ) };

# Add size information to @array_features:
@array_features		= @{ Homology::ProcessFeatures ( \@array_features ) };

# --------------------------------------

# Create and save formatted link data:
Homology::Links ( \%hash_variables, \%hash_new_0, $var_nmer );

# Add genomic feature data to each link:
Homology::FeatureData ( \%hash_variables, "direct", \@array_features );
Homology::FeatureData ( \%hash_variables, "inverted", \@array_features );

Quantify::RawLinks ( \%hash_variables, "direct" );
Quantify::RawLinks ( \%hash_variables, "inverted" );

# Merge direct and inverted links whose locations are off by one:
Homology::MergeLinks (  \%hash_variables, "direct", $var_min );
Homology::MergeLinks ( \%hash_variables, "inverted", $var_min );

# Flatten direct and inverted links into repeat regions:
Homology::FlattenLinks ( \%hash_variables, "direct" );
Homology::FlattenLinks ( \%hash_variables, "inverted" );
Homology::FlattenLinks ( \%hash_variables, "combined" );

# Quantify distribution of merged direct and inverted repeats:
Quantify::Distribution ( \%hash_variables, "direct" );
Quantify::Distribution ( \%hash_variables, "inverted" );

print "\n    + Quantification of repetitive regions:\n";

# Quantify overlap of genomic features with flattened repeat regions:
Quantify::Regions ( \%hash_variables, "direct", \@array_features );
Quantify::Regions ( \%hash_variables, "inverted", \@array_features );
Quantify::Regions ( \%hash_variables, "combined", \@array_features );

# ---------------------------------------------------------------------------- #

print "\n- $file_hash1";

foreach my $var_key ( keys %hash_1 ) {

	if ( $hash_0 { $var_key } ) {

		my $var_scalar	= scalar @{ $hash_0 { $var_key } };

		$hash_new_1 { $var_key }	= $hash_1 { $var_key };

	}

}

# Quantify overlap between raw repeats and genomic features:
Quantify::OverlapRaw ( \%hash_variables, \%hash_new_1, \@array_features );

# --------------------------------------

$hash_variables { file_direct }				= "$var_output/$var_stem.1.$var_nmer.direct";
$hash_variables { file_inverted }			= "$var_output/$var_stem.1.$var_nmer.inverted";
$hash_variables { file_direct_links }		= "$var_output/$var_stem.1.direct.$var_nmer.links";
$hash_variables { file_inverted_links }		= "$var_output/$var_stem.1.inverted.$var_nmer.links";
$hash_variables { file_direct_features }	= "$var_output/$var_stem.direct.$var_nmer.1.features";
$hash_variables { file_inverted_features }	= "$var_output/$var_stem.inverted.$var_nmer.1.features";
$hash_variables { file_direct_merged }		= "$var_output/$var_stem.direct.$var_nmer.1.merged";
$hash_variables { file_inverted_merged }	= "$var_output/$var_stem.inverted.$var_nmer.1.0merged";
$hash_variables { file_direct_flat }		= "$var_output/$var_stem.direct.$var_nmer.1.flat";
$hash_variables { file_inverted_flat }		= "$var_output/$var_stem.inverted.$var_nmer.1.flat";
$hash_variables { file_flat }				= "$var_output/$var_stem.$var_nmer.1.flat";

# Create and save formatted link data:
Homology::Links ( \%hash_variables, \%hash_new_1, $var_nmer );

# Add genomic feature data to each link:
Homology::FeatureData ( \%hash_variables, "direct", \@array_features );
Homology::FeatureData ( \%hash_variables, "inverted", \@array_features );

Quantify::RawLinks ( \%hash_variables, "direct" );
Quantify::RawLinks ( \%hash_variables, "inverted" );

# Merge direct and inverted links whose locations are off by one:
Homology::MergeLinks (  \%hash_variables, "direct", $var_min );
Homology::MergeLinks ( \%hash_variables, "inverted", $var_min );

# Flatten direct and inverted links into repeat regions:
Homology::FlattenLinks ( \%hash_variables, "direct" );
Homology::FlattenLinks ( \%hash_variables, "inverted" );
Homology::FlattenLinks ( \%hash_variables, "combined" );

# Quantify distribution of merged direct and inverted repeats:
Quantify::Distribution ( \%hash_variables, "direct" );
Quantify::Distribution ( \%hash_variables, "inverted" );

print "\n    + Quantification of repetitive regions:\n";

# Quantify overlap of genomic features with flattened repeat regions:
Quantify::Regions ( \%hash_variables, "direct", \@array_features );
Quantify::Regions ( \%hash_variables, "inverted", \@array_features );
Quantify::Regions ( \%hash_variables, "combined", \@array_features );

# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #

sub LoadHash {

	# Argument:
	my ( $file_in )	= @_;

	# Data-structures:
	my %hash_out;

	# ----------------------------------

	print "$file_in\n";

	# Open file containing stored hash.
	open ( my $file_read, '<', $file_in ) or die " ERROR: UNABLE TO OPEN $file_in!\n";

	my $hash_json = <$file_read>;

	%hash_out = %{decode_json( $hash_json )};

	close $file_read or die " ERROR: UNABLE TO CLOSE $file_in!\n";

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}
