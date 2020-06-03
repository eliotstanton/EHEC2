package Quantify2;

use strict;
use warnings;
use JSON::XS;

# ---------------------------------------------------------------------------- #

# File name:		Quantify.pm
# Date created:		08 January, 2020
# Last modified:	02 June, 2020
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This package contains subroutines used for quantifing the
#					abundance and location of repeats.

# ---------------------------------------------------------------------------- #

# Subroutine name:	CopyNumber
# Description:		This subroutine quantifies the distribution of repeats by
#					their copy number and link number:

# ----------------------------------------------------------

sub CopyNumber {

	# Arguments:
	my ( $hash_in, $hash_nmer )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_nmer		= %$hash_nmer if $hash_nmer;
	my %hash_copies;
	my %hash_direct;
	my %hash_inverted;
	my %hash_both;

	# Variables:
	my $var_nmer	= $hash_out { var_nmer };

	# ----------------------------------

	# If %hash_nmer is undefined, load it:
	unless ( %hash_nmer ) {

		my $file_hash	= $hash_out { file_hash };

		print "LOADING HASH from $file_hash\n";

		# Open file containing stored hash.
		open ( my $file_read, '<', $file_hash ) or die " ERROR: UNABLE TO OPEN $file_hash!\n";

		my $hash_json = <$file_read>;

		%hash_nmer = %{decode_json( $hash_json )};

		close $file_read or die " ERROR: UNABLE TO CLOSE $file_hash!\n";

	}

	# ----------------------------------

	# Iterate through %hash_nmer:
	foreach my $var_key ( keys %hash_nmer ) {

		# Define array reference mapped by key:
		my $var_value	= $hash_nmer { $var_key };

		next unless $var_value;

		# Define number of locations in array:
		my $var_scalar	= scalar @{$var_value};

		# Define variables for holding copy number and link number:
		my $var_copies	= $var_scalar;
		my $var_both	= 0;

		# If direct repeats are present:
		if ( $var_scalar > 1 ) {

			# Define number of links:
			my $var_scalar2	= ($var_scalar*($var_scalar-1))/2;

			# Iterate up by one for number of linkes:
			$hash_direct { $var_scalar2 }	+= 1;

			# Add number of links to $var_both:
			$var_both	+= $var_scalar2;

		}

		# Define reverse complement of nucleotide sequence:
		my $var_complement	= General::ReverseComplement ( $var_key );

		# If $var_complement is present in %hash_nmer:
		if ( $hash_nmer { $var_complement } ) {

			# Define number of locations in the array:
			my $var_scalar2	= scalar @{$hash_nmer {$var_complement} };

			# Add numbe rof locations to copy number:
			$var_copies	+= $var_scalar2;

			# Iterate up by one for number of inverted links present:
			$hash_inverted	{ $var_scalar2 }	+= 1;

			# Add number of links to $var_both:
			$var_both		+= $var_scalar2;

			# Remove $var_complement from %hash_nmer:
			delete $hash_nmer { $var_complement };

		}

		$hash_copies { $var_copies }	+= 1;

		$hash_both { $var_both }	+= 1;

	}

	# ----------------------------------

	print "\n    + Quantification of repeat sequence abundance:\n";

	print "\t- $var_nmer.mer repeat copy numbers (direct and inverted):\n";
	foreach my $var_key ( sort { $a <=> $b } keys %hash_copies ) {

		my $var_value	= $hash_copies { $var_key };

		print "\t    $var_key - $var_value\n";

	}

	# ----------------------------------

#	print "\t- Direct raw repeat links:\n";

	# Iterate through %hash_direct:
	foreach my $var_key ( sort { $a <=> $b } keys %hash_direct ) {

		my $var_value	= $hash_direct { $var_key };

#		print "\t    $var_key - $var_value\n";

	}

#	print "\t- Inverted raw repeat links:\n";

	# Iterate through %hash_inverted:
	foreach my $var_key ( sort { $a <=> $b } keys %hash_inverted ) {

		my $var_value	= $hash_inverted { $var_key };

#		print "\t    $var_key - $var_value\n";

	}

#	print "\t- Combined raw repeat links:\n";

	foreach my $var_key ( sort { $a <=> $b } keys %hash_both ) {

		my $var_value	= $hash_both { $var_key };

#		print "\t    $var_key - $var_value\n";

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	OverlapRaw
# Description:		Quantifies overlap between raw repeats and genomic features.

# ----------------------------------------------------------

sub OverlapRaw {

	# Arguments:
	my ( $hash_in, $hash_nmer, $array_features )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_nmer		= %$hash_nmer;
	my @array_features	= @$array_features;
	my @array_nmer;
	my %hash_class;

	# Variables:
	my $var_nmer		= $hash_out { var_nmer };
	my $var_total		= 0;
	
	# ----------------------------------

	# Export locations in %hash_nmer to @array_nmer:
	foreach my $var_key ( keys %hash_nmer ) {

		my $var_value	= $hash_nmer { $var_key };

		push( @array_nmer, @{$var_value} );

	}

	$var_total	= scalar @array_nmer;

	# Sort @array_features by size:
	@array_features	= sort { $a -> [5] <=> $b -> [5] } @array_features;

	# Sort @array_nmer by location:
	@array_nmer		= sort { $a -> [1] <=> $b -> [1] } @array_nmer;

	# ----------------------------------

	# Reduce the amount of feature information:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		my $var_type	= $array_features[$i][1];
		my $var_loc0	= $array_features[$i][3];
		my $var_loc1	= $array_features[$i][4];

#		print "$i: @{$array_features[$i]} - $var_type $var_loc0 $var_loc1\n";

		my @array_temp	= ( $var_type, $var_loc0, $var_loc1 );

		$array_features[$i]	= \@array_temp;

	}

	# ----------------------------------

	# Iterate through each repeat and define feature:
	for ( my $i = 0; $i < scalar @array_nmer; $i++ ) {

		my $var_loc0	= $array_nmer[$i][1];
		my $var_loc1	= $var_loc0 + $var_nmer;

#		print "$i: @{$array_nmer[$i]} - $var_loc0 $var_loc1\n";

		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			my $var_type	= $array_features[$j][0];
			my $var_loc2	= $array_features[$j][1];
			my $var_loc3	= $array_features[$j][2];

#			last if $var_loc2 > $var_loc1;

			if ( $var_loc0 >= $var_loc2 && $var_loc1 <= $var_loc3 ) {

#				print "\t$j: @{$array_features[$j]}\n";

				$hash_class { $var_type } += 1;

				last;

			}

			elsif ( $j == scalar @array_features - 1 ) {

				$hash_class { "chromosome" } += 1;

			}

		}

	}

	# ----------------------------------

	# Print data to command-line:
	print "\n    + Breakdown of repeat sequences by genomic feature:\n";

	foreach my $var_class ( sort keys %hash_class ) {

		my $var_count	= $hash_class{$var_class};

		print "\t- $var_class - $var_count\n";

		$var_total -= $var_count;

	}

	print "\t- Total: $var_total\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #

1;
