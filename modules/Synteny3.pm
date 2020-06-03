package Synteny3;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatAlignment
# Description:		This subroutine takes output from Circos and formats data
#					for use by other subroutines.
#
# Format for each region is as follows:
# [0: SeqID][1: Loc0][2: Loc1][3: RefID][4: Copy][5: Length0][6: Length1]
# [7: Inverted0][8: Inverted1][9: Bracket0][10: Bracket1]
# [11: Global0][12: Global1]

# ----------------------------------------------------------

sub FormatAlignment {

	# Arguments:
	my ( $hash_in, $array_FASTA )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_FASTA		= @$array_FASTA;
	my @array_blocks;
	my @array_in;
	my @array_inserts;
	my @array_backbone;
	my %hash_all;

	# Variables:
	my $file_backbone	= $hash_out { file_backbone };
	my $file_align		= $hash_out { file_align };
	my $file_bbcols		= $hash_out { file_bbcols };
	my $var_minimum		= $hash_out { var_minimum };

	# ----------------------------------

	# If $file_backbone does not exist, run progressiveMauve:
	Synteny::progressiveMauve ( $hash_in, $array_FASTA ) unless -e $file_backbone;

	# ----------------------------------

	# Load $file_bbcols into @array_in:
#	@array_in		= @{ General::FileToArray ( "$file_bbcols" ) };

	# Define the number of blocks:
#	my $var_blocks	= $array_in[$#array_in][0];

	# Populate @array_blocks by the number of sequences and number of blocks:
#	for ( my $i = 0; $i <= $var_blocks; $i++ ) {

#		my @array_temp	= ( ("0") x scalar @array_FASTA );

#		push @array_blocks, \@array_temp;

#	}

	# Identify which blocks are occupied for each sequence:
#	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define block ID number:
#		my $var_blockID	= $array_in[$i][0];

		# Add conserved regions for the block:
#		for ( my $j = 0; $j < scalar @array_FASTA; $j++ ) {

#			my $var_genomeID	= $array_in[$i][$j+3];

#			next unless defined $var_genomeID;

#			$array_blocks[$var_blockID][$var_genomeID]	= 1;

#		}

#	}

	# ----------------------------------

	# Load $file_align into @array_in:
#	@array_in	= @{ General::FileToArray ( "$file_align" ) };

	# Remove extraneous data from @array_in:
#	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#		my $var_line	= join " ", @{$array_in[$i]};

#		unless ( $var_line =~ /\>/  ) {

#			splice @array_in, $i, 1 and $i--;

#		}

#	}

	# Add sequence ID, block ID, and start/end locations for each block:
#	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

		# Iterate down through each sequence:
#		for ( my $j = 0; $j <= $var_blocks; $j++ ) {

			# If present search and add data:
#			if ( $array_blocks[$j][$i] ) {

				# Iterate down through @array_in to find the first hit:
#				for ( my $k = 0; $k < scalar @array_in; $k++ ) {

					# Define sequence and locations:
#					my ( $var_seq, $var_loc0, $var_loc1 )	= split ( /[\:,\-]/, $array_in[$k][1]);

					# If there is a match add it to @array_blocks and remove it
					# from @array_in:
#					if ( $i == $var_seq - 1 ) {

#						my @array_temp	= ( $i, $j, $var_loc0, $var_loc1 );

#						$array_blocks[$j][$i]	= \@array_temp;

#						splice @array_in, $k, 1;

#						last;

#					}

#				}

#			}

#		}

#	}

	# ----------------------------------

	# TODO: Identify block IDs that are transposed:
#	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

#		print "$i:\n";

#		my $var_loc0			= $array_blocks[0][$i][2];
#		$array_blocks[0][$i][4]	= 0;

#		for ( my $j = 1; $j < scalar @array_blocks; $j++ ) {

#			next unless $array_blocks[$j][$i];

#			my $var_loc1	= $array_blocks[$j][$i][2];

#			if ( $var_loc1 < $var_loc0 ) {

#				$array_blocks[$j][$i][4]	= 1;

#				print "\t$j: @{$array_blocks[$j][$i]} - $var_loc0 $var_loc1\n";

#			}

#			else {

#				$array_blocks[$j][$i][4]	= 0;

#				print "\t$j: @{$array_blocks[$j][$i]} - $var_loc0 $var_loc1\n";

#			}

#			$var_loc0 = $var_loc1;

#		}

#	}

#	exit;

	# ----------------------------------

	# Load $file_backbone into @array_in:
	@array_in	= @{ General::FileToArray ( $file_backbone ) };

	# Iterate through each alignment in @array_in:
	for ( my $i = 1; $i < scalar @array_in; $i++ ) {

#		print "$i: @{$array_in[$i]}\n";

		# Define variable to hold copy number:
		my $var_copy		= 0;

		# Define variable to hold global length:
		my $var_length0		= 0;

		# Define array to hold alignment data:
		my @array_alignment;

		# Define variable to hold inverted status:
		my $var_inverted0	= 0;

		# Define variable to hold sequence ID:
		my $var_seq			= 0;

		# Define variable to hold transposed status:
		my $var_transposed	= 0;

		# Iterate through each group of locations:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j+=2 ) {

			# Define locations:
			my $var_loc0		= $array_in[$i][$j];
			my $var_loc1		= $array_in[$i][$j+1];

			# Define variable to hold inverted status of region:
			my $var_inverted1	= 0;

			# If region is inverted straighten it up and flag both inversion
			# variables:
			if ( $var_loc0 < 0 ) {

				$var_loc0 *= -1;
				$var_loc1 *= -1;

				$var_inverted1	= 1;
				$var_inverted0	= 1;

			}

			# Increase copy number by one if region isn't empty:
			$var_copy++ if $var_loc0 > 0;

			# Define length of region:
			my $var_length1	= $var_loc1 - $var_loc0 + 1;

			# If $var_length1 is greater than $var_length0, reassign 
			# $var_length0 to equal $var_length0:
			$var_length0 = $var_length1 if $var_length1 > $var_length0;

			# Create a temporary array to hold region data:
			my @array_temp	= (0)x13;

			# Add informaiton to temporary array:
			$array_temp[0]	= $var_seq;
			$array_temp[1]	= $var_loc0;
			$array_temp[2]	= $var_loc1;
			$array_temp[3]	= $i;
			$array_temp[6]	= $var_length1;
			$array_temp[8]	= $var_inverted1;

			# Determine if region is located in transposed block:
#			for ( my $k = 0; $k < scalar @array_blocks; $k++ ) {

#				next unless $array_blocks[$k][$var_seq];

#				my $var_loc2		= $array_blocks[$k][$var_seq][2];
#				my $var_loc3		= $array_blocks[$k][$var_seq][3];
#				my $var_transposed0	= $array_blocks[$k][$var_seq][4];

#				if ( $var_loc0 >= $var_loc2 && $var_loc1 <= $var_loc3 ) {

#					print "\t$var_seq $k: @{$array_blocks[$k][$var_seq]} - $var_loc2 $var_loc3 - $var_transposed0\n";

#					print "\t$var_seq $k: @{$array_blocks[$k][$var_seq]} - $var_loc2 $var_loc3 - $var_transposed\n";

#					last;

#				}

#				print "\t$k:\n";

#			}

			# Add reference to temporary array to @array_alignment:
			push @array_alignment, \@array_temp;

			# Increment sequence ID for next alignment:
			$var_seq++;

		}
	 
		# Move on if $var_length0 is less than $var_minimum:
		next if $var_length0 < $var_minimum;

		# Store $var_length0, $var_copy and $var_inverted0 in each alignment in
		# %hash_all:
		for ( my $j = 0; $j < scalar @array_alignment; $j++ ) {

			my $var_seqID			= $array_alignment[$j][0];
			my $var_loc0			= $array_alignment[$j][1];
			my $var_refID			= $array_alignment[$j][3];

			$array_alignment[$j][4] = $var_copy;
			$array_alignment[$j][5]	= $var_length0;
			$array_alignment[$j][7]	= $var_inverted0;

			next if $var_loc0 == 0;

			$hash_all{"$var_seqID $var_refID" }	= $array_alignment[$j];

		}

		# Add regions to @array_backbone if all the regions are present:
		if ( $var_copy == $var_seq ) {

			push @array_backbone, \@array_alignment;

		}

		# Otherwise add in to @array_inserts if there are regions absent:
		else {

			push @array_inserts, \@array_alignment;

		}

	}

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store %hash_all in %hash_out:
	$hash_out { hash_all }			= \%hash_all;

	# ----------------------------------

	$hash_out { array_merged }			= \@array_backbone;

	# ----------------------------------

	# Return referene to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	SortBySequence
# Description:		This subroutine sorts regions by their sequence.

# ----------------------------------------------------------

sub SortBySequence {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Sort regions by sequence:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			push @{$array_out[$j]}, $array_in[$i][$j];

		}

	}

	# Sort regions in each sequence:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		@{$array_out[$i]}	= sort { $a -> [1] <=> $b -> [1] } @{$array_out[$i]};

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	GlobalLocations
# Description		This subroutine assigns final locations for Circos.

# ----------------------------------------------------------

sub GlobalLocations {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_in	= @{ $hash_out { array_merged }};

	# Variables:
	my $var_start	= 1;
	my $var_end		= 0;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variables to hold start and end locations:
		my $var_start	= 1;
		my $var_end		= 0;

		# Iterate down through regions assigning locations:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define variables holding global and local location region lengths:
			my $var_length0	= $array_in[$i][$j][5];
			my $var_length1	= $array_in[$i][$j][6];

			$var_end	+= $var_length0;

			unless ( $var_length1 == 0 ) {

				my $var_diff	= int( ( $var_length0 - $var_length1 )/2 );

				my $var_loc0	= $var_start + $var_diff;
				my $var_loc1	= $var_end - $var_diff;

				$array_in[$i][$j][11]	= $var_loc0;
				$array_in[$i][$j][12]	= $var_loc1;

			}

			$var_start	= $var_end + 1;

		}

	}

	# ----------------------------------

	# Store @array_in back in %hash_out:
	$hash_out { array_merged }	= \@array_in;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	JunctionInserts
# Description:		This subroutine sequesters insert regions that sit within
#					a junction of backbone regions into a separate array.

# ----------------------------------------------------------

sub JunctionInserts {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my @array_inserts		= @{$hash_out{array_inserts}};
	my %hash_junctions;
	my @array_junctions;

	# ----------------------------------

	# Identify $var_refIDs associated with junctions:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

#		print "$i:\n";

		# Define variable to hold previous reference ID:
		my $var_refID0		= $array_backbone[$i][-1][3];

		# Define variable to hold previous alignment inverted status:
		my $var_inverted0	= $array_backbone[$i][-1][7];

		# Iterate down through regions in @array_backbone:
		for ( my $j = 0; $j < scalar @{$array_backbone[$i]}; $j++ ) {

			# Define current alignment inverted status:
			my $var_inverted1	= $array_backbone[$i][$j][7];

			# Define current reference ID:
			my $var_refID1		= $array_backbone[$i][$j][3];

			# If current alignment and previous alignment don't agree on
			# inverted status store reference ID in %hash_junctionIDs:
			if ( $var_inverted0 != $var_inverted1 ) {

#				print "\t$j: @{$array_backbone[$i][$j-1]}\n";
#				print "\t$j: @{$array_backbone[$i][$j]}\n";

				# Store current and previous reference IDs in %hash_junctions:
				$hash_junctions { "$var_refID0 $var_refID1" } = 1;

#				print "\t\t$var_refID0 $var_refID1\n";

			}

			# Store current alignments inverted status for next iteration:
			$var_inverted0		= $var_inverted1;

			# Store current reference ID for next iteration:
			$var_refID0			= $var_refID1;

		}

	}

	# ----------------------------------

	# Sequester inserts in @array_inserts associated with junctions and store
	# them in @array_junctions:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

#		print "$i\n";

		# Iterate through regions in alignment:
		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			# Define variables to hold bracketing reference IDs:
			my $var_bracket0	= $array_inserts[$i][$j][9];
			my $var_bracket1	= $array_inserts[$i][$j][10];

			next if $var_bracket0	== 0;

			my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

			if ( $hash_junctions {"$var_bracket0 $var_bracket1"} ) {

#				print "\t$j: @{$array_inserts[$i][$j]}\n";

				push @array_junctions, $array_inserts[$i];

				splice @array_inserts, $i, 1;

				$i--;

#				last;

			}

		}

	}

	# ----------------------------------

	# Store new version of @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store @array_junctions in %hash_out:
	$hash_out { array_junctions }	= \@array_junctions;

	# Store %hash_junctionIDs in %hash_out:
	$hash_out { hash_junctions }	= \%hash_junctions;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	IdentifyTransposed
# Description:		This subroutine identifies and stores transposed regions
#					present in all sequences in a new array.

# ----------------------------------------------------------

sub IdentifyTransposed {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_IDs;
	my @array_backbone	= @{$hash_out{array_backbone}};
	my @array_out		= @{$hash_out{array_backbone}};
	my @array_IDs;
	my @array_transposed;

	# ----------------------------------

	# Sort regions in @array_backbone by sequence and by location:
	@array_backbone		= @{ SortBySequence ( \@array_backbone )};

	# Store reference IDs from each region as elements in @array_IDs:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

#		print "$i\n";

		for ( my $j = 0; $j < scalar @{$array_backbone[$i]}; $j++ ) {

			my $var_refID	= $array_backbone[$i][$j][3];

#			print "\t$j: $var_refID\n";

			push @{$array_IDs[$i]}, $var_refID;

		}

	}

	# ----------------------------------

	# Perform gapped alignment:
	for ( my $i = 0; $i < scalar @array_IDs - 1; $i++ ) {

		my %hash_ID	= %{ GappedAlignment ( \%hash_IDs, $array_IDs[$i], $array_IDs[$i+1] )};

		%hash_IDs	= ( %hash_IDs, %hash_ID );

#		foreach my $var_keys ( keys %hash_IDs ) {

#			print "\t$var_keys\n";

#		}

#		print "\n";

	}

	# ----------------------------------

	# Remove transposed regions from @array_out:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my $var_refID	= $array_out[$i][0][3];

#		print "$i: $var_refID\n";

		if ( $hash_IDs { $var_refID } ) {

#			print "***\n";

			push @array_transposed, $array_out[$i];

			splice @array_out, $i, 1;

			$i--;

		}

	}

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_out;

	# Store @array_transposed in %hash_out:
	$hash_out { array_transposed }	= \@array_transposed;

	# ----------------------------------

	$hash_out { array_merged }		= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	GappedAlignment
# Description:		This subroutine uses Needleman-Wunsch algorithm to identify
#					transposed regions.

# ----------------------------------------------------------

sub GappedAlignment {

	# Arguments:
	my ( $hash_in, $array_in0, $array_in1 )	= @_;

	# Data-structures:
	my %hash_out	= %$hash_in;
	my @array_in0	= @{$array_in0};
	my @array_in1	= @{$array_in1};
	my @array_out0;
	my @array_out1;

	# ----------------------------------

	# Scoring function:
	sub score_sub {

		# Gap penalty:
	  	if (!@_) { return -2; }

		# Mismatch scores -1, match +1
		return ($_[0] eq $_[1]) ? 1 : -1;

	} 

	# ----------------------------------

	my $var_matcher	= Algorithm::NeedlemanWunsch->new(\&score_sub);
	my $var_score	= $var_matcher->align( \@array_in0, \@array_in1, {
	   
		align => sub { 

			unshift @array_out0, $array_in0[$_[0]];
			unshift @array_out1, $array_in1[$_[1]];

		},

		shift_a => sub { 

			unshift @array_out0, $array_in0[$_[0]];
			unshift @array_out1, " ";
			$hash_out{$array_in0[$_[0]]}	= "1";

		},

		shift_b => sub { 

			unshift @array_out0, " ";
			unshift @array_out1, $array_in1[$_[0]];
			$hash_out{$array_in1[$_[0]]}	= "1";

		}
	} );

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketRefs
# Description:		This subroutine is a wrapper for adding bracketing reference
#					IDs to insert and transposed regions.

# ----------------------------------------------------------

sub BracketRefs {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my @array_inserts		= @{$hash_out{array_inserts}};
	my @array_transposed	= @{$hash_out{array_transposed}};

	# ----------------------------------

	# Sort regions in @array_backbone by sequence:
	@array_backbone		= @{ Synteny3::SortBySequence ( \@array_backbone) };

	# Add bracketing reference IDs to @array_inserts:
	@array_inserts		= @{ BracketInserts(\@array_inserts,\@array_backbone) };

	# Add bracketing reference IDs to @array_transposed:
	@array_transposed	= @{ BracketInserts(\@array_transposed,\@array_backbone) };

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store @array_inserts in %hash_out:
	$hash_out { array_transposed }	= \@array_transposed;

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_merged }	= \@array_backbone;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketInserts
# Description:		This subroutine assigns bracketing reference IDs to each
#					insert region.

# ----------------------------------------------------------

sub BracketInserts {

	# Arguments:
	my ( $array_inserts, $array_backbone )	= @_;

	# Data structures:
	my @array_backbone		= @$array_backbone;
	my @array_inserts		= @$array_inserts;

	# ----------------------------------

	# Iterate through @array_inserts:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

#		print "$i: @{$array_inserts[$i]}\n";

		# Iterate through regions:
		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			# Define sequence ID and initial location:
			my $var_seq		= $array_inserts[$i][$j][0];
			my $var_loc0	= $array_inserts[$i][$j][1];

			# Move on if region is blank:
			next if $var_loc0 == 0;

			# Iterate down through backbone regions:
			for ( my $k = 0; $k < scalar @{$array_backbone[$var_seq]}; $k++ ) {

				# Define initial backbone region location:
				my $var_loc1	= $array_backbone[$var_seq][$k][1];

				# If backbone region is larger define bracketing reference IDs
				# and end loop:
				if ( $var_loc1 > $var_loc0 ) {

					# Define reference IDs of the current and former backbone
					# regions:
					my $var_refID0		= $array_backbone[$var_seq][$k-1][3];
					my $var_refID1		= $array_backbone[$var_seq][$k][3];

					# Define inverted status of both regions:
					my $var_inverted0	= $array_backbone[$var_seq][$k-1][8];
					my $var_inverted1	= $array_backbone[$var_seq][$k][8];

					if ( $var_refID0 < $var_refID1 ) {

						# Store bracketing data in region:
						$array_inserts[$i][$j][9]	= $var_refID0;
						$array_inserts[$i][$j][10]	= $var_refID1;

					}

					else {

						# Store bracketing data in region:
						$array_inserts[$i][$j][9]	= $var_refID1;
						$array_inserts[$i][$j][10]	= $var_refID0;

					}

					# If $var_inverted0 and $var_inverted1 both are positive,
					# redefine region as inverted as well:
					if ( $var_inverted0 == 1 && $var_inverted1 == 1 ) {

						$array_inserts[$i][$j][8]	= 1;

					}

					last;

				}

			}

		}

	}

	# ----------------------------------

	# Return reference to \@array_inserts and end subroutine:
	return \@array_inserts;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	OrganiseTransposed
# Description:		This subroutine organises convserved transposed regions and
#					integrates them back in to the backbone:

# ----------------------------------------------------------

sub OrganiseTransposed {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_transposed	= @{$hash_out{array_transposed}};
	my @array_backbone		= @{$hash_out{array_backbone}};

	# Variables:
	my $var_scalar			= scalar @array_backbone;

	# ----------------------------------

	# Organise transposed regions into a hash based upon bracketing reference IDs:
	my %hash_transposed		= %{ BracketHash($hash_out{array_transposed}) };

	# ----------------------------------

	# Iterate through each pair of bracket IDs in %hash_transposed:
	foreach my $var_key ( keys %hash_transposed ) {

		# Define bracket IDs:
		my ( $var_bracket0, $var_bracket1 )	= split " ", $var_key;

		# Define temporary arrays and hash:
		my @array_temp;
		my @array_temp2;
		my %hash_temp;

		# Find regions that have the same bracketing IDs:
		for ( my $i = 0; $i < scalar @array_transposed; $i++ ) {

			for ( my $j = 0; $j < scalar @{$array_transposed[$i]}; $j++ ) {

				my $var_bracket2	= $array_transposed[$i][$j][9];
				my $var_bracket3	= $array_transposed[$i][$j][10];

				# Add regions to @array_temp:
				if ( $var_bracket0 == $var_bracket2 && $var_bracket1 == $var_bracket3 ) {

					push @array_temp, $array_transposed[$i][$j];

				}

			}

		}

		# Rearrange regions by sequence and store reference IDs in %hash_temp:
		for ( my $i = 0; $i < scalar @array_temp; $i++ ) {

			my $var_seqID	= $array_temp[$i][0];
			my $var_refID	= $array_temp[$i][3];

			push @{$array_temp2[$var_seqID]}, $array_temp[$i];

			$hash_temp { $var_refID }	= "";

		}

		# Define the number of reference IDs in %hash_temp:
		my $var_scalar2	= scalar keys %hash_temp;

		# Define variable to store maximum lenght of regions:
		my $var_length0	= 0;

		# Determine the maximum length of regions::
		for ( my $i = 0; $i < $var_scalar; $i++ ) {

			my $var_length1	= 0;

			next unless $array_temp2[$i][0];

			for ( my $j = 0; $j < scalar @{$array_temp2[$i]}; $j++ ) {

				$var_length1	+= $array_temp2[$i][$j][5];

			}

			$var_length0	= $var_length1 if $var_length1 > $var_length0;

		}

		# Add empty regions to unoccupied sequences:
		for ( my $i = 0; $i < $var_scalar; $i++ ) {

			unless ( $array_temp2[$i][0] ) {

				for ( my $j = $var_scalar2; $j > 0; $j-- ) {

					my @array_temp3	= ( $i, 0, 0, 0, 0, 0, 0, 0, 0, $var_bracket0, $var_bracket1, 0, 0 );

					push @{$array_temp2[$i]}, \@array_temp3;

				}

				$array_temp2[$i][0][5]	= $var_length0;

			}
	
		}

		# Store @array_temp back in %hash_transposed:
		$hash_transposed { $var_key } = \@array_temp2;

	}

	# ----------------------------------

	# Itegrate regions into @array_backbone:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

		# Iterate down through regions:
		for ( my $j = 0; $j < scalar @{$array_backbone[$i]}; $j++ ) {

			# Define bracketing reference IDs:
			my $var_bracket0	= $array_backbone[$i][$j-1][3];
			my $var_bracket1	= $array_backbone[$i][$j][3];

			if ( $hash_transposed { "$var_bracket0 $var_bracket1" } ) {

				my @array_temp	= @{$hash_transposed { "$var_bracket0 $var_bracket1" }};

				for ( my $k = scalar @{$array_temp[$i]} - 1; $k >= 0; $k-- ) {

					splice @{$array_backbone[$i]}, $j, 0, $array_temp[$i][$k]; 

				}

			}

		}

	}

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store %hash_transposed in %hash_out:
	$hash_out { array_transposed }	= \@array_backbone;

	$hash_out { array_merged }		= \@array_backbone;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketHash
# Description:		This subroutine stores regions in a hash of arrays using
#					pairs of bracketing reference IDs as keys.

# ----------------------------------------------------------

sub BracketHash {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my %hash_out;

	# ----------------------------------

	# Iterate through @array_in and store regions in %hash_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate down through each region:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define bracketing reference IDs:
			my $var_bracket0	= $array_in[$i][$j][9];
			my $var_bracket1	= $array_in[$i][$j][10];

			# Move on if region is a blank:
			next if $var_bracket0 == 0 && $var_bracket1 == 0;

			# Add region from @array_in to %hash_merged:
			push ( @{$hash_out{"$var_bracket0 $var_bracket1"}}, $array_in[$i][$j]);

		}

	}

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

1;
