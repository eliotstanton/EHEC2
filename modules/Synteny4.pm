package Synteny4;

use strict;
use warnings;
use Circos;

# ---------------------------------------------------------------------------- #

# Subroutine name:	Config
# Description:		This subroutine creates a config file for Circos

# ----------------------------------------------------------

sub Config {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_merged	= @{$hash_in{array_merged}};
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { "file_config" };
	my $file_highlights		= $hash_in { "file_highlights" };
	my $file_ideogram		= $hash_in { "file_ideogram" };
	my $file_image			= $hash_in { "file_image" };
	my $file_karyotype		= $hash_in { "file_karyotype" };
	my $file_ticks			= $hash_in { "file_ticks" };
	my $var_scalar			= scalar @array_merged;
	my $var_increment		= $hash_in { var_increment };
	my $var_radius0			= $hash_in { var_radius0 };
	my $var_radius1			= $hash_in { var_radius1 };
	my $r					= "r";

	# ----------------------------------

	push @array_out, "#$file_out";
	push @array_out, "karyotype = $file_karyotype";
	push @array_out, "<<include $file_ideogram>>";
	push @array_out, "<ideogram>\nshow = no\n</ideogram>";

	push @array_out, "<highlights>";

	for ( my $i = 0; $i < $var_scalar; $i++ ) {

		# Add highlights:
		push @array_out, "\t<highlight>";
		push @array_out, "\t\tfile	= $file_highlights.$i\n\t\tideogram	= no";
		push @array_out, "\t\tr0	= $var_radius0$r\n\t\tr1	= $var_radius1$r";
	 	push @array_out, "\t</highlight>";

		# Add ticks:
#		push @array_out, "\t<highlight>";
#		push @array_out, "\t\tfile	= $file_ticks.$i\n\t\tideogram	= no";
#	 	push @array_out, "\t</highlight>";
	 
	 	$var_radius0	-= $var_increment;
	 	$var_radius1	-= $var_increment;

	}
 
 	push @array_out, "</highlights>";
 
	push @array_out, "<<include $file_image>>";
	push @array_out, "<<include etc/colors_fonts_patterns.conf>>";
	push @array_out, "<<include etc/housekeeping.conf>>";

	# Write @array_out to file:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Highlights
# Description:		This subroutine handles formatting highlights for Synteny.

# ----------------------------------------------------------

sub Highlights {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{ $hash_in { array_merged } };
	my @array_features	= @{ $hash_in { array_features } };
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { file_highlights };
	my $var_thickness	= $hash_in { stroke_thickness };
	my $var_span		= $hash_in { var_span };
	my $var_scalar		= scalar @array_in;

	# ----------------------------------

	# Add highlights for regions:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate down through @array_in formatting data for Circos:
		for ( my $j	= 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			my $var_copy		= $array_in[$i][$j][4];
			my $var_inverted	= $array_in[$i][$j][8];
			my $var_loc0		= $array_in[$i][$j][11];
			my $var_loc1		= $array_in[$i][$j][12];

			# Ignore blank regions:		 
			next if $var_loc0	== 0;

			my $var_border		= $hash_colours { chromosome_border };

			my $var_colour		= Circos::ColoursNew ( "chromosome", $var_copy, $var_scalar );

			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "offset=$var_span," if $var_inverted;
			$var_string			.= "stroke_color=$var_border,";
			$var_string			.= "stroke_thickness=$var_thickness,";
			$var_string			.= "fill_color=$var_colour";

			push @{$array_out[$i]}, $var_string;

		}

	}

	# ----------------------------------

	# Repeat the process for genomic features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		next unless $array_features[$i];

		# Go down through each genomic feature stored in @array_features:
		for ( my $j	= 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			# Define start and end locations of feature block:
			my $var_loc0		= $array_features[$i][$j][0];
			my $var_loc1		= $array_features[$i][$j][1];

			# Define feature type, inverted status and copy number:
			my $var_type		= $array_features[$i][$j][2];
			my $var_inverted	= $array_features[$i][$j][3];
			my $var_copy		= $array_features[$i][$j][4];

			# Define feature colours:
			my $var_border		= $hash_colours { "$var_type\_border" };

			my $var_colour		= Circos::ColoursNew ( $var_type, $var_copy, $var_scalar );

			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "offset=$var_span," if $var_inverted;
			$var_string			.= "stroke_color=$var_border," if $var_border;
			$var_string			.= "stroke_thickness=0,";
			$var_string			.= "fill_color=$var_colour";

			push @{$array_out[$i]}, $var_string;

		}

	}

	# ----------------------------------

	# Write data from @array_out to file:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my @array_temp	= @{$array_out[$i]};

		unshift @array_temp, "# $file_out.$i";

		# Print @array_out to $file_out:
		General::ArrayToFile ( \@array_temp, "$file_out.$i" );

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Karyotype
# Description:		This subroutine produces a formatted Circos karyotype file.

# ----------------------------------------------------------

sub Karyotype {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{$hash_in { array_merged }};
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { "file_karyotype" };
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 1 );

	# ----------------------------------

	my $var_scalar		= scalar @{$array_in[0]} - 1;
	my $var_length		= $array_in[0][$var_scalar][12];

	my $var_string	= "chr - chr0 0 0 $var_length $var_colour\n";

	print "$var_string\n";

	push @array_out, $var_string;

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatFeatures
# Description:		This subroutine create feature blocks with adjusted
#					locations and with feature type data.

# ----------------------------------------------------------

sub FormatFeatures {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_out;

	# ----------------------------------34 PM

	# Iterate through regions for each sequence:
	for ( my $i = 0; $i < scalar @array_merged; $i ++ ) {

#		print "$i\n";

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

#			print "\t$j: @{$array_merged[$i][$j]}\n";

			my $var_loc0		= $array_merged[$i][$j][1];
			my $var_loc1		= $array_merged[$i][$j][2];
			my $var_copy		= $array_merged[$i][$j][4];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc2		= $array_merged[$i][$j][11];
			my $var_loc3		= $array_merged[$i][$j][12];

			my $var_diff	= $var_loc2 - $var_loc0;

			for ( my $k = 0; $k < scalar @{$array_in[$i]}; $k++ ) {

				# Define start and end locations of genomic feature:
				my $var_type	= $array_in[$i][$k][1];
				my $var_loc4	= $array_in[$i][$k][3];
				my $var_loc5	= $array_in[$i][$k][4];
				my $var_orient	= $array_in[$i][$k][5];

				# If genomic feature lies completely within the region:
				if ( $var_loc0 <= $var_loc4 && $var_loc1 >= $var_loc5 ) {

#					print "\t\t$k: @{$array_in[$i][$k]}\n";

					my $var_new0	= $var_loc4 + $var_diff;
					my $var_new1	= $var_loc5 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If region is covered by the genomic feature:
				elsif ( $var_loc0 >= $var_loc4 && $var_loc1 <= $var_loc5 ) {

#					print "\t\t$k: @{$array_in[$i][$k]}\n";

					my $var_new0	= $var_loc2;
					my $var_new1	= $var_loc3;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If genomic feature straddles beginning of region:
				elsif ( $var_loc4 < $var_loc0 && $var_loc5 > $var_loc0 && $var_loc5 < $var_loc1 ) {

#					print "\t\t$k: @{$array_in[$i][$k]}\n";

					my $var_new0	= $var_loc2;
					my $var_new1	= $var_loc5 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If genomic feature straddles the end of the region:
				elsif ( $var_loc4 < $var_loc1 && $var_loc4 > $var_loc0 && $var_loc5 > $var_loc1 ) {

#					print "\t\t$k: @{$array_in[$i][$k]}\n";

					my $var_new0	= $var_loc4 + $var_diff;
					my $var_new1	= $var_loc3;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}


			}

		}

	}

	# ----------------------------------

	# Iterate down through genomic features for each sequence:
#	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#		print "$i\n";

		# Go through each feature and find overlapping regions of the
		# chromosome:
#		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

#			print "\t$j: @{$array_in[$i][$j]}\n";

			# Define start, end, type, and orientation of feature:
#			my $var_loc0	= $array_in[$i][$j][3];
#			my $var_loc1	= $array_in[$i][$j][4];
#			my $var_type	= $array_in[$i][$j][1];
#			my $var_orient	= $array_in[$i][$j][5];

			# Iterate down through regions in @array_merged finding those region
			# blocks that overlap::
#			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				# Define start, end, and inverted status of each regions:
#				my $var_loc2		= $array_merged[$i][$k][1];
#				my $var_loc3		= $array_merged[$i][$k][2];
#				my $var_copy		= $array_merged[$i][$k][4];
#				my $var_inverted	= $array_merged[$i][$k][8];

				# Skip blanks:
#				next if $var_loc2 == 0;

				# Define adjusted locations of region:
#				my $var_loc4	= $array_merged[$i][$k][11];
#				my $var_loc5	= $array_merged[$i][$k][12];

				# Define variables to hold adjusted locations of features:
#				my $var_diff	= $var_loc4-$var_loc2;
#				my $var_new0	= 0;
#				my $var_new1	= 0;

				# If region is entirely covered by the genomic feature:
#				if ( $var_loc2 >= $var_loc0 && $var_loc3 <= $var_loc1 ) {

#					print "\t\t$k: @{$array_merged[$i][$k]}\n";

#					my @array_temp	= ( $var_loc4, $var_loc5, $var_type, $var_inverted, $var_copy, $var_orient );

#					push @{$array_out[$i]}, \@array_temp;

#				}

				# If genomic feature is entirely within a single region:
#				elsif ( $var_loc2 <= $var_loc0 && $var_loc3 >= $var_loc1 ) {

#					print "\t\t$k: @{$array_merged[$i][$k]}\n";

#					$var_new0	= $var_loc0 + $var_diff;
#					$var_new1	= $var_loc1 + $var_diff;

#					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

#					push @{$array_out[$i]}, \@array_temp;

#				}

				# If genomic feature straddles the beginning of a region:
#				if ( $var_loc0 <= $var_loc2 && $var_loc1 > $var_loc2 && $var_loc1 < $var_loc3 ) {

#					print "\t\t$k: @{$array_merged[$i][$k]}\n";
#
#					$var_new0	= $var_loc4;
#					$var_new1	= $var_loc1 + $var_diff;

#					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

#					push @{$array_out[$i]}, \@array_temp;


#				}

				# If genomic feature straddles the end of a region:
#				if ( $var_loc0 > $var_loc2 && $var_loc0 < $var_loc3 && $var_loc1 >= $var_loc3 ) {

#					print "\t\t$k: @{$array_merged[$i][$k]}\n";
#
#					$var_new0	= $var_loc0 + $var_diff;
#					$var_new1	= $var_loc5;

#					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

#					push @{$array_out[$i]}, \@array_temp;


#				}

#			}

#		}

#	}

	# ----------------------------------

	# Store @array_out in %hash_out:
	$hash_out { array_features }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

1;
