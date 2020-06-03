package Quantify;

use strict;
use warnings;
use JSON::XS;

# ---------------------------------------------------------------------------- #

# File name:		Quantify.pm
# Date created:		08 January, 2020
# Last modified:	23 January, 2020
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This package contains subroutines used for quantifing the
#					abundance and location of repeats.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Distribution
# Description:		Prints distribution of length of raw links:

# ----------------------------------------------------------

sub Distribution {

	# Arguments:
	my ( $hash_in, $var_variable )	= @_;

	# Data-structures:
	my %hash_out	= %$hash_in;
	my @array_in;

	# Variables:
	my $file_in;

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

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define locations of inital link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		# Calulate length of link:
		my $var_length	= $var_loc1 - $var_loc0 + 1;

		# Store information in array:
		$array_in[$i][6]	= $var_length;

	}

	# ----------------------------------

	# Sort @array_in by length:
	@array_in	= sort { $a -> [6] <=> $b -> [6] } @array_in;

	# Save @array_in to $file_in:
	General::ArrayToFile ( \@array_in, $file_in );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	RawLinks
# Description:		This subroutine quantifies the abundance of repeats in
#					different genomic features.

# ----------------------------------------------------------

sub RawLinks {

	# Arguments:
	my ( $hash_in, $var_variable )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_class;
	my %hash_split;
	my @array_in;

	# Variables:
	my $file_in;
	my $var_total	= 0;
	my $var_nmer	= $hash_out { var_nmer };

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

	# Print message to command-line:
	print "\n    + Quantification of $var_variable pairs of $var_nmer.mer $var_variable links:\n";

	# ----------------------------------

	# Iterate through @array_in grabbing each genomic classification:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define locations of each link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Define $var_class0 and $var_class1:
		my $var_class0	= $array_in[$i][6];
		my $var_class1	= $array_in[$i][7];

		# If both classes are the same add information to %hash_class:
		if ( $var_class0 eq $var_class1 ) {

			$hash_class { "$var_class0" } += 1;

		}

		# If the two classes aren't the same and one of the classes is
		# chromosome, store in %hash_class under chromosome key:
		elsif ( $var_class0 eq "chromosome" || $var_class1 eq "chromosome" ) {

			$hash_class { "$var_class0" } += 1 if $var_class1 eq "chromosome";
			$hash_class { "$var_class1" } += 1 if $var_class0 eq "chromosome";

		}

		# If the two classes aren't the same and one of the classes is IS
		# store in %hash_class under IS key:
		elsif ( $var_class0 eq "IS" || $var_class1 eq "IS" ) {

			$hash_class { "IS" } += 1;

		}

		# If the two classes aren't the same and one of the classes is rhs,
		# store in %hash_class under rhs key:
		elsif ( $var_class0 eq "rhs" || $var_class1 eq "rhs" ) {

			$hash_class { "rhs" } += 1;

		}

		# If the two classes aren't the same and one of the classes is
		# prophage, store in %hash_class under prophage/PLE key:
		elsif ( $var_class0 eq "prophage" || $var_class1 eq "prophage" ) {

			$hash_class { "prophage/PLE" } += 1;

		}

		# If the two classes aren't the same and one of the classes is PLE
		# store in %hash_class under prophage/PLE key:
		elsif ( $var_class0 eq "PLE" || $var_class1 eq "PLE" ) {

			$hash_class { "prophage/PLE" } += 1;

		}

		else {

			print "$i: @{$array_in[$i]}\n";

		}

		# Store locations and class in %hash_split:
		$hash_split { "$var_loc0 $var_loc1" }	= $var_class0;
		$hash_split { "$var_loc2 $var_loc3" }	= $var_class1;

	}

	# ----------------------------------

	# Print data to command-line:
	print "\t- Breakdown of $var_variable links by shared genomic feature:\n";

	foreach my $var_class ( sort keys %hash_class ) {

		my $var_count	= $hash_class{$var_class};

		print "\t    * $var_class - $var_count\n";

		$var_total += $var_count;

	}

	print "\t    * Total: $var_total\n";

	# ----------------------------------

	# Reset %hash_class and $var_total for reuse:
	undef %hash_class;
	$var_total	= 0;

	# Iterate each individual pair of locations in %hash_split and enumerate:
	while ( my ( $var_locations, $var_class ) = each %hash_split ) {

		$hash_class { $var_class }	+= 1;

		$var_total	+= 1;

	}

	# ----------------------------------

	# Print breakdown:
	print "\t- Breakdown of $var_variable links by individual sequence location:\n";

	foreach my $var_class ( sort keys %hash_class ) {

		my $var_count	= $hash_class { $var_class };

		print "\t    * $var_class - $var_count\n";

	}

	print "\t    * Total: $var_total\n";

	# ----------------------------------

	# End subroutine:
	return;

}

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

# Subroutine name:	Overlap
# Description:		This subroutine briefly calculates overlap between genomic
# 					features in each sequence and conserved/non-conserved
#					regions.

# ----------------------------------------------------------

sub Overlap {

	# Arguments:
	my ( $array_features, $array_FASTA, $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_features	= @$array_features;
	my @array_FASTA		= @$array_FASTA;
	my @array_out;

	# Variables:
	my $var_scalar	= scalar @array_merged;

	# ----------------------------------

	# Iterate through each collection of genomic features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# Define temporary hashes to store data:
		my %hash_temp0;
		my %hash_temp1;

		# Iterate down through features:
		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			my $var_type	= $array_features[$i][$j][1];
			my $var_loc0	= $array_features[$i][$j][3];
			my $var_loc1	= $array_features[$i][$j][4];

			# Iterate down through @array_merged identifying and tallying
			# overlapping regions:
			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				my $var_loc2	= $array_merged[$i][$k][1];
				my $var_loc3	= $array_merged[$i][$k][2];
				my $var_copy	= $array_merged[$i][$k][4];

				next if $var_loc2 == 0;

				my $var_length	= 0;

				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc0 && $var_loc3 < $var_loc1 ) {

					$var_length	= $var_loc3 - $var_loc0 + 1;

				}

				if ( $var_loc2 > $var_loc0 && $var_loc2 < $var_loc1 && $var_loc3 > $var_loc1 ) {

					$var_length	= $var_loc1 - $var_loc2 + 1;

				}

				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc1 ) {

					$var_length	= $var_loc1 - $var_loc0 + 1;

				}

				if ( $var_loc2 > $var_loc0 && $var_loc3 < $var_loc1 ) {

					$var_length	= $var_loc3 - $var_loc2 + 1;

				}

				# Add conserved regions to %hash_temp0:
				if ( $var_copy == $var_scalar ) {

					$hash_temp0 { $var_type }	+= $var_length;

				}

				# Add non-conserved regions to %hash_temp1:
				else {

					$hash_temp1 { $var_type }	+= $var_length;

				}

				last if $var_loc2 > $var_loc0;

			}

		}

		# Store temporary hashes in @array_out:
		$array_out[$i][0]	= \%hash_temp0;
		$array_out[$i][1]	= \%hash_temp1;

	}

	# ----------------------------------

	# Iterate through @array_out reporting data to user:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		# Define temporary array holding regions:
		my @array_regions	= @{$array_merged[$i]};

		my $var_seq_length	= $array_regions[$#array_regions][2];

		print "$i: $array_FASTA[$i] - $var_seq_length bp\n";

		# Define temporary arrays holding length data:
		my %hash_temp0	= %{$array_out[$i][0]};
		my %hash_temp1	= %{$array_out[$i][1]};

		print "    Conserved:\n";

		# Print conserved data:
		for my $var_temp ( keys %hash_temp0 ) {

			my $var_length	= $hash_temp0{$var_temp};

			print "\t* $var_temp -\t$var_length bp\n";

		}

		print "    Non-conserved:\n";

		# Print non-conserved data:
		for my $var_temp ( keys %hash_temp1 ) {

			my $var_length	= $hash_temp0{$var_temp};

			print "\t* $var_temp -\t$var_length bp\n";

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:		Regions
# Description:			This subroutine quantifies overlap between flattened
#						repeat regions and genomic features.

# --------------------------------------

sub Regions {

	# Arguments:
	my ( $hash_in, $var_variable, $array_features )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @$array_features;
	my @array_in;
	my @array_out;
	my %hash_features;
	my %hash_class;

	# Variables:
	my $file_in;
	my $var_chromosome	= $hash_out { var_chr_length };

	# ----------------------------------

	# Set $file_in:
	if ( $var_variable eq "direct" ) {

		$file_in		= $hash_out { file_direct_flat };

	}

	if ( $var_variable eq "inverted" ) {

		$file_in		= $hash_out { file_inverted_flat };

	}

	if ( $var_variable eq "combined" ) {

		$file_in		= $hash_out { file_flat };

	}

	# Import data in $file_in to @array_in:
	@array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	# Sort regions in @array_feature by size:
	@array_features	= sort { $a -> [5] <=> $b -> [5] } @array_features;

	# Iterate through @array_features and store total length of each feature
	# in %hash_features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# Define variables holding feature classification and length:
		my $var_class	= $array_features[$i][1];
		my $var_length	= $array_features[$i][5];

		# Add length information to %hash_features using classification as key:
		$hash_features { $var_class } += $var_length;

		$var_chromosome	-= $var_length;

	}

	# Add chromosome length to %hash_features:
	$hash_features { "chromosome" }	= $var_chromosome;

	# ----------------------------------

	# Iterate through each link and determine which genomic features overlap:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations:
		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		# Iterate down through @array_features looking for a match:
		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			# Define element class and start/end locations:
			my $var_class	= $array_features[$j][1];
			my $var_loc2	= $array_features[$j][3];
			my $var_loc3	= $array_features[$j][4];

			# If repeat is entirely within the genomic region:
			if ( $var_loc0 >= $var_loc2 && $var_loc1 <= $var_loc3 ) {

				$array_in[$i][2]	= $var_class;

				last;

			}

			# If genomic feature is entirely within the repeat:
			elsif ( $var_loc0 < $var_loc2 && $var_loc1 > $var_loc3 ) {

				# Remove the element from @array_in;
				splice @array_in, $i, 1;

				# Create temporary arrays with updated locations for the repeat:
				my @array_temp0	=	( $var_loc0, $var_loc2-1 );
				my @array_temp1	=	( $var_loc2, $var_loc3 );
				my @array_temp2	=	( $var_loc3+1, $var_loc1 );

				# Add references to temporary arrays to @array_in:
				splice @array_in, $i, 0, \@array_temp2;
				splice @array_in, $i, 0, \@array_temp1;
				splice @array_in, $i, 0, \@array_temp0;

				$i = -1;

				last;

			}

			# If repeat overlaps the left side of the feature:
			elsif ( $var_loc0 < $var_loc2 && $var_loc1 >= $var_loc2 ) {

				# Remove the element from @array_in;
				splice @array_in, $i, 1;

				# Create temporary arrays with updated locations for the repeat:
				my @array_temp0	=	( $var_loc0, $var_loc2-1 );
				my @array_temp1	=	( $var_loc2, $var_loc1 );

				# Add references to temporary arrays to @array_in:
				splice @array_in, $i, 0, \@array_temp1;
				splice @array_in, $i, 0, \@array_temp0;

				$i = -1;

				last;

			}

			# If repeat overlaps the right side of the feature:
			elsif ( $var_loc0 <= $var_loc3 && $var_loc1 > $var_loc3 ) {

				# Remove the element from @array_in;
				splice @array_in, $i, 1;

				# Create temporary arrays with updated locations for the repeat:
				my @array_temp0	=	( $var_loc0, $var_loc3 );
				my @array_temp1	=	( $var_loc3+1, $var_loc1 );

				# Add references to temporary arrays to @array_in:
				splice @array_in, $i, 0, \@array_temp1;
				splice @array_in, $i, 0, \@array_temp0;

				$i = -1;

				last;

			}

		}

	}

	# ----------------------------------

	# Iterate through @array_in and add length of each region to appropiate
	# key in %hash_class:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations of link:
		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		# Define the links genomic feature classification:
		my $var_class	= $array_in[$i][2] || "chromosome";

		# Define the length of the genomic region:
		my $var_length	= $var_loc1 - $var_loc0 + 1;

		# Store length in %hash_class using genomic feature as the key:
		$hash_class { $var_class } 	+= $var_length;

	}

	# ----------------------------------

	# Print data to command-line:
	print "\t- Breakdown of $var_variable flattened repetitive region by genomic feature:\n";

	# Iterate through %hash_class and print data to command-line:
	foreach my $var_class ( sort keys %hash_class ) {

		my $var_length			= $hash_class { $var_class };

		my $var_class_length	= $hash_features { $var_class };

		print "\t    * $var_class - $var_length ($var_class_length total)\n";

	}

	# ----------------------------------

	# Return @array_out and end subroutine:
	return \@array_out;

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

	# Quantify number of elements in @array_nmer:
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

	for ( my $i = 0; $i < scalar @array_nmer; $i++ ) {

		my $var_loc0	= $array_nmer[$i][1];
		my $var_loc1	= $var_loc0 + $var_nmer;

#		print "$i: @{$array_nmer[$i]} - $var_loc0 $var_loc1\n";

		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			my $var_type	= $array_features[$j][0];
			my $var_loc2	= $array_features[$j][1];
			my $var_loc3	= $array_features[$j][2];

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

#		$var_total += $var_count;

	}

	print "\t- Total: $var_total\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
