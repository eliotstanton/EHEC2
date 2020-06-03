#!/usr/bin/perl

use strict;
use warnings;
#use Getopt::Std;
use lib '/etc/perl/modules';
use General;
use Check;
#use Homology;
#use SVG;

# ---------------------------------------------------------------------------- #

# File name:		ReverseContig.pl
# Date created:		06 February, 2020
# Last modified:	06 February, 2020
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This is script reverse complements contigs in a FASTA file

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_out;

# Import variables from the command-line provided by the user: 
my $var_FASTA	= $ARGV[0];
my $file_out	= $ARGV[1];
my $var_length	= 80;

# Define variables used by this script:

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "";


# If $var_FASTA or $var_output is missing, print $var_help and stop script:
unless ( $var_FASTA && $file_out ) {

	print "$var_help\n";

	exit;

}


# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ Check::Files ( \@array_FASTA, \%hash_variables ) };

# ---------------------------------------------------------------------------- #

# Import FASTA file to @array_FASTA:
@array_FASTA	= @{General::ImportFASTA ( \@array_FASTA, $var_help )};

# Iterate through @array_FASTA:
for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

	# Define FASTA header:
	my $var_header	= $array_FASTA[$i][1];

	# Print header to command-line:
	print "$var_header\n";

	# Define nucleotide sequence as a string:
	my $var_string	= $array_FASTA[$i][2];

	# Create reverse complement:
	my $var_reverse	= General::ReverseComplement ( $var_string );

	# Add line breaks:
	$var_reverse =~ s/(.{1,$var_length})/$1\n/gs;

	# Store $var_header in @array_out:
	push @array_out, $var_header;

	# Store $var_reverse in @array_out:
	push @array_out, $var_reverse;

#	print "$var_reverse\n";

}

# Print @array_out to $file_out:
General::ArrayToFile ( \@array_out, $file_out );

# ---------------------------------------------------------------------------- #
