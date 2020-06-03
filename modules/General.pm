package General;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# File Name:		General.Pm
# Date Created:		25 November, 2019
# Last Modified:	12 January, 2020
# Created By:		Eliot Stanton (estanton@Wisc.Edu)

# Description:		This Package Contains Routine And Commonly Used Subroutines.

# ---------------------------------------------------------------------------- #

# Subroutine name:	ArrayToFile
# Description:		This subroutine takes an array and writes it to file.

# ----------------------------------------------------------

sub ArrayToFile {

	# Arguments:
	my ( $array_in, $file_out )	= @_;

	# Data structures:
	my @array_in		= @$array_in;

	# Variables:
	my $var_string;

	# ----------------------------------

	# Iterate through each element in @array_in adding each to $var_string:
	foreach my $var_element ( @array_in ) {

		# If the element contains an element unpack it:
		if ( ref $var_element ) {

			$var_string	.= "@{$var_element}\n";

		}

		# Otherwise concatenate directly:
		else {

			$var_string	.= "$var_element\n";

		}

	}

	# Open the file to be read:
	open ( my $file_write, ">", "$file_out" ) or die "Unable to open $file_out!\n";

	# Print $var_string to file:
	print $file_write "$var_string";

	# Close file being read:
	close $file_write or die "Unable to close $file_out!\n";

	# ----------------------------------

	# End subroutine;
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name: 	FileToArray
# Description:		This subroutine imports data from a file and converts each
#					line to an element in an array.

# ----------------------------------------------------------

sub FileToArray {

	# Arguments:
	my ( $file_in )	= @_;

	# Data structures:
	my @array_out;

	# ----------------------------------

	# Open the file to be read:
	open ( my $file_read, "<", "$file_in" ) or die "Unable to open $file_in!\n";

	# Iterate through each line and store as appropiate:
	while (<$file_read>) {

		my $var_line	= $_;

		chomp $var_line;

		# Skip the line if it is commented out:
		next if $var_line =~ /^\#/;

		my @array_temp	= split " ", $var_line;

		push @array_out, \@array_temp;

	}

	# Close file being read:
	close $file_read or die "Unable to close $file_in!\n";

	# ----------------------------------

	# Return @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	ImportFASTA
# Description:		This subroutine opens a FASTA file and stores each entry in
#					the file as a string in an array.

# ----------------------------------------------------------

sub ImportFASTA {

	# Arguments:
	my ( $array_in, $var_help )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# Variables:
	my $var_counter	= -1;

	# ----------------------------------

	# Iterate through each FASTA file in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variable holding file path:
		my $file_in	= $array_in[$i];

		# Split file path into a temporary array to check file extension:
		my @array_temp	= split /\./, $file_in;

		# Print help file and exit if file doesn't have the correct extension:
		unless ( ( $array_temp[$#array_temp] eq "fa" ) || 
			( $array_temp[$#array_temp] eq "fasta" ) ) {

			print "$var_help\n";

			print "$file_in not recognised as FASTA file!\n";

			exit;

		}

		# Redefine temporary array initially holding just the file path:
		@array_temp	= ( $array_in[$i] );

		# Store temporary array as a reference in @array_out:
		push @array_out, \@array_temp;

		# Open the file to be read:
		open ( my $file_read, "<", "$file_in" ) or die "Unable to open $file_in!\n";

		# Iterate through each line and store as appropiate:
		while (<$file_read>) {

			# Define each line of file:
			my $var_line	= $_;

			# Remove line-break:
			chomp $var_line;

			# If entry header, add to second element:
			if ( $var_line =~ />/ ) {

				$var_counter++;

				$array_out[$var_counter][0]	= $file_in;

				$array_out[$var_counter][1]	= $var_line;

			}

			# Otherwise, add to third element:
			else {

				$array_out[$var_counter][2]	.= $var_line;

			}

		}

		# Close the file after reading:
		close $file_read or die "Unable to close $file_in!\n";

	}

	# ----------------------------------

	return \@array_out;

}


# ---------------------------------------------------------------------------- #

# Subroutine name:	ImportFiles
# Description:		This subroutine takes an array of file names and imports
#					the data into an array:

# ----------------------------------------------------------

sub ImportFiles {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Iterate through each file in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $file_in	= $array_in[$i];

		# Import data into a temporary array:
		my @array_temp	= @{General::FileToArray ( $file_in )};

		# Store reference to the temporary array in @array_out:
		push @array_out, \@array_temp;

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	ReverseComplement
# Description:		This subroutine takes a nucleotide sequence and returns the
#					reverse complement sequence.

# ----------------------------------------------------------

sub ReverseComplement {

	# Arguments:
	my ( $var_string )	= @_;

	# Variables:
	my $var_reverse;
	my $var_temp;

	# ----------------------------------

	$var_reverse = reverse($var_string);

	$var_reverse =~ tr/ACGTacgt/TGCAtgca/;

	# ----------------------------------

	# End subroutine and return $var_reverse:
	return $var_reverse;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Variables
# Description:		This subroutine stores standardised variables in a hash.

# ----------------------------------------------------------

sub Variables {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;

	# ----------------------------------

	# Variables:
	my $var_minimum		= $hash_out{m} || "1000";
	my $var_nmer		= $hash_out{n} || "100";
#	my $var_output		= $hash_out{ va} || "homology";
	my $var_stem		= $hash_out{s} || "default";
	my $var_output		= $hash_out { var_output };
	my $var_radius		= 1500;
#	my $var_link_width	= 100;
	my $var_radius0		= 0.99;
	my $var_radius1		= 0.92;
	my $var_span		= ($var_radius0 - $var_radius1)/2;
	my $var_factor		= 10;

	# ----------------------------------

	# Store defaults:
#	$hash_out { var_stem }				= $var_stem;
	$hash_out { var_minimum }			= $var_minimum;

	# ----------------------------------

	# Store Mauve filepaths in %hash_out:
	$hash_out { file_align }			= "$var_output/$var_stem.align";
	$hash_out { file_backbone }			= "$var_output/$var_stem.backbone";
	$hash_out { file_guide }			= "$var_output/$var_stem.guide";
	$hash_out { file_bbcols }			= "$var_output/$var_stem.align.bbcols";

	# ----------------------------------

	# Store output filepaths in %hash_out:
	$hash_out { file_config }				= "$var_output/$var_stem.config";
	$hash_out { file_conserved }			= "$var_output/$var_stem.conserved";
	$hash_out { file_direct }				= "$var_output/$var_stem.$var_nmer.direct";
	$hash_out { file_direct_links }			= "$var_output/$var_stem.direct.$var_nmer.links";
	$hash_out { file_direct_features }		= "$var_output/$var_stem.direct.$var_nmer.features";
	$hash_out { file_direct_merged }		= "$var_output/$var_stem.direct.$var_nmer.merged";
	$hash_out { file_direct_flat }			= "$var_output/$var_stem.direct.$var_nmer.flat";
#	$hash_out { file_direct_fasta }			= "$var_output/$var_stem.direct.fasta";
	$hash_out { file_features }				= "$var_output/$var_stem.features";
	$hash_out { file_flat }					= "$var_output/$var_stem.$var_nmer.flat";
	$hash_out { file_final }				= "$var_output/$var_stem.final";
	$hash_out { file_hash }					= "$var_output/$var_stem.$var_nmer.hash";
	$hash_out { file_highlights }			= "$var_output/$var_stem.highlights";
	$hash_out { file_ideogram }				= "$var_output/$var_stem.ideogram";
	$hash_out { file_image }				= "$var_output/$var_stem.image";
	$hash_out { file_inverted }				= "$var_output/$var_stem.$var_nmer.inverted";
	$hash_out { file_inverted_links }		= "$var_output/$var_stem.inverted.$var_nmer.links";
	$hash_out { file_inverted_features }	= "$var_output/$var_stem.inverted.$var_nmer.features";
	$hash_out { file_inverted_merged }		= "$var_output/$var_stem.inverted.$var_nmer.merged";
	$hash_out { file_inverted_flat }		= "$var_output/$var_stem.inverted.$var_nmer.flat";
#	$hash_out { file_inverted_fasta }		= "$var_output/$var_stem.inverted.fasta";
	$hash_out { file_karyotype }			= "$var_output/$var_stem.karyotype";
#	$hash_out { file_matrix }				= "$var_output/$var_stem.matrix";
	$hash_out { file_nonconserved }			= "$var_output/$var_stem.nonconserved";
	$hash_out { file_png }					= "$var_output/$var_stem.png";
	$hash_out { file_svg }					= "$var_output/$var_stem.linear.svg";
	$hash_out { file_ticks }				= "$var_output/$var_stem.ticks";

	# ----------------------------------

	$hash_out { var_image_radius }			= "$var_radius";
	$hash_out { stroke_thickness }			= "2";
#	$hash_out { stroke_thickness_links }	= "0.1";

	$hash_out { var_increment }				= 1500/10000;
	$hash_out { var_span }					= $var_radius*$var_span;


	$hash_out { var_radius0 }				= $var_radius0;
	$hash_out { var_radius1 }				= $var_radius1;

	$hash_out { var_gap }					= "10";

	$hash_out { var_1Mb_height }			= "40";
	$hash_out { var_100kb_height }			= "20";
	$hash_out { var_10kb_height }			= "10";
	$hash_out { var_1kb_height }			= "5";

	$hash_out { var_factor }				= $var_factor;

#	$hash_out { var_link_width }			= $var_link_width;

	# ----------------------------------

	# Return %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

1;
