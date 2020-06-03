package Check;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# File name:		Check.pm
# Date created:		08 Januar, 2020
# Last modified:	08 Januar, 2020
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This package contains subroutines used checking if other
#					required elements are present.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Bowtie
# Description:		This subroutine checks if bowtie is present, if it isn't it
#					prints an error message and kills the script.

# ----------------------------------------------------------

sub Bowtie {

	# Arguments:
	( my $hash_in )	= @_;

	# Data structures:
	my %hash_in		= %$hash_in;

	# Variables:
	my $var_help	= $hash_in { var_help };
	my $var_return	= 0;

	# ----------------------------------

	$var_return			= system ( "which bowtie > /dev/null" );

	if ( $var_return != "0" ) {

		print "$var_help\n";

		print "\tERROR: Bowtie not found!\n";

		exit;

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Circos
# Description:		This subroutine checks if Circos is present, if it isn't
#					it prints an error message and kills the script.

# ----------------------------------------------------------

sub Circos {

	# Arguments:
	( my $hash_in )	= @_;

	# Data structures:
	my %hash_in		= %$hash_in;

	# Variables:
	my $var_help	= $hash_in { var_help };
	my $var_return	= 0;

	# ----------------------------------

	$var_return			= system ( "which circos > /dev/null" );

	if ( $var_return != "0" ) {

		print "$var_help\n";

		print "\tERROR: Circos not found!\n";

		exit;

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Directory
# Description:		This subroutine checks if a directory exists and if it
#					doesn't, it creates it.

# ----------------------------------------------------------

sub Directory {

	# Arguments:
	my ( $var_directory, $hash_in )	= @_;

	# Data-structures:
	my %hash_in			= %$hash_in;

	# Variables:
	my $var_help		= $hash_in { var_help };
	my $var_path;

	# ----------------------------------

	unless ( $var_directory ) {

		print "$var_help\n";

		print "\tERROR: Output directory required!\n" and exit;

	}

	# ----------------------------------

	# If $var_directory doesn't exist split the path into an array:
	unless ( -e $var_directory ) {

		my @array_temp	= split "\/", $var_directory;

		# Iterate through @array_temp:
		foreach my $var_temp ( @array_temp ) {

			$var_path .= "/" if $var_path;

			$var_path	.= "$var_temp";

			mkdir ( $var_path ) unless -e $var_path;

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Files
# Description:		This subroutine checks if a file exists and contains data
#					if it doesn't exist, or doesn't contain data, it returns
#					null.

# ----------------------------------------------------------

sub Files {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data-structures:
	my @array_in	= @$array_in;
	my %hash_in		= %$hash_in;

	# Variables:
	my $var_help	= $hash_in { var_help };

	# ----------------------------------

	# Iterate through each file and check if it is present and contains data:
	foreach my $var_file ( @array_in ) {

		my $var_return	= $var_file if -s "$var_file";

		unless ( $var_return ) {

			print "$var_help\n";

			print "\tERROR: $var_file not found or contains no data!\n";

		}

	}

	# ----------------------------------

	return \@array_in;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	ProgressiveMauve:
# Description:		This subroutine checks if progressiveMauve is present, if it
#					isn't prints an error message and kills the script.

# ----------------------------------------------------------

sub ProgressiveMauve {

	# Arguments:
	( my $hash_in )	= @_;

	# Data structures:
	my %hash_in		= %$hash_in;

	# Variables:
	my $var_help	= $hash_in { var_help };
	my $var_return	= 0;

	# ----------------------------------

	$var_return			= system ( "which progressiveMauve > /dev/null" );

	if ( $var_return != "0" ) {

		print "$var_help\n";

		print "\tERROR: progressiveMauve not found!\n" and exit;

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
