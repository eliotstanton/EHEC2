#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
#use Homology;
#use SVG;

# ---------------------------------------------------------------------------- #

# File name:		HomologyLite.pl
# Date created:		09 February, 2020
# Last modified:	09 February, 2020
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This is a script for generating an in silico restriction
#					site map from a FASTA file.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_out;

# Import variables from the command-line provided by the user:
my $var_FASTA		= $ARGV[0];
my $var_sequence	= $ARGV[1];
my $file_out		= $ARGV[2];
my $var_minimum		= $ARGV[3] || 1000;

# Define variables used by this script:
my $var_help		= "";

# ---------------------------------------------------------------------------- #

# Define help file:

# If $var_output or $var_FASTA is missing, print $var_help and stop script:
unless ( $var_FASTA && $var_sequence && $file_out ) {

	print "$var_help\n";

	exit;

}

# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };

# Import data from @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# ---------------------------------------------------------------------------- #

# Redefine $var_FASTA for nucleotide sequence:
$var_FASTA			= $array_FASTA[0][2];

# Define reverse complement of $var_sequence:
#my $var_reverse		= reverse($var_sequence);

# Search for $var_sequence and $var_reverse in $var_FASTA and return their
# locations:
while ( $var_FASTA =~ m/$var_sequence/g )  { 
 
    my $var_location = pos ( $var_FASTA ); 

#    print "$var_location\n";

	push @array_out, $var_location;

} 

#while ( $var_FASTA =~ m/$var_reverse/g )  { 
 
#    my $var_location = pos ( $var_FASTA ); 

#    print "$var_location\n"; 

#	push @array_out, $var_location;

#} 

# ---------------------------------------------------------------------------- #

# Sort locations in @array_out:
@array_out	= sort { $a <=> $b } @array_out;

# Define initial location:
my $var_loc0	= 1;

# Define the size of each restriction fragment and store:
for ( my $i =0; $i < scalar @array_out; $i++ ) {

	my $var_loc1	= $array_out[$i];

	my $var_size	= $var_loc1 - $var_loc0;

	$var_loc0		= $var_loc1 + 1;

#	print "$i: $var_loc0 $var_loc1 - $var_size\n";

	if ( $var_size >= $var_minimum ) {

		$array_out[$i]	= $var_size;

#		print "$i: $var_loc0 $var_loc1 - $var_size\n";

	}

	else {

		splice @array_out, $i, 1;

		$i--;

	}



}


# Print @array_out to file:
General::ArrayToFile ( \@array_out, $file_out );
