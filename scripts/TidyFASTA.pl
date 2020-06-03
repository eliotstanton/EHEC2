#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/Perl/Modules';
use General;

# ---------------------------------------------------------------------------- #

# File Name:		TidayFasta
# Date Created:		25 November, 2019
# Last Modified:	25 November, 2019
# Created By:		Eliot Stanton (Estanton@Wisc.Edu)

# Description:		This is a script for re-Editing FASTA files so they are eeat
#					and tidy.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_in;
#my @array_out;

# Import variables from the command-line provided by the user:
getopts('l:', \%hash_variables);
my $var_length	= $hash_variables{l} || 80;
my $file_in		= $ARGV[0];
my $file_out	= $ARGV[1];

# ---------------------------------------------------------------------------- #

# Define string containing helpful information:
my $var_help	= "TidyFASTA.pl [OPTIONS] [FASTA IN]\n";
$var_help		.= "    -l    LINE LENGTH (default: 80)\n";

unless ( $file_in ) {

	print "$var_help\n";

}

# ---------------------------------------------------------------------------- #

# Temporarily store file name in @array_in for subroutine:
$array_in[0]	= $file_in;

# Import FASTA file to @array_in:
@array_in	= @{General::ImportFASTA ( \@array_in, $var_help )};

# Reformat length of sequence lines:
for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#	print "$i: $array_in[$i][0]\n";

	my $var_header		= "$array_in[$i][1]";

	my $var_sequence	= $array_in[$i][2];

	$var_sequence		=~ s/(.{1,$var_length})/$1\n/gs;

	print "$var_header\n";

	print "$var_sequence\n";

}

# ---------------------------------------------------------------------------- #
