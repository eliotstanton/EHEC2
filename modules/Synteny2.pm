package Synteny2;

use strict;
use warnings;
use Algorithm::NeedlemanWunsch;

# ---------------------------------------------------------------------------- #

# File name:		Synteny2.pm
# Date created:		31 October, 2018
# Last modified:	12 January, 2020
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to visualising
#					synteny data.

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
		push @array_out, "\t<highlight>";
		push @array_out, "\t\tfile	= $file_ticks.$i\n\t\tideogram	= no";
	 	push @array_out, "\t</highlight>";
	 
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

# Subroutine name:	Ticks
# Description:		This subroutine takes region data and creates tick marks.

# ----------------------------------------------------------

sub Ticks {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_in		= @{$hash_out{ array_merged }};
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out			= $hash_out { file_ticks };
	my $var_thickness		= $hash_out { stroke_thickness };
	my $var_span0			= $hash_out { var_span };
	my $var_colour			= $hash_colours { ticks };

	my $var_gap				= $hash_out { var_gap };

	my $var_radius0			= $hash_out { var_radius0 };

	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_increment		= $hash_out { var_increment };
	my $var_span			= $hash_out { var_span } * 2;

	# ----------------------------------

	# Iterate through @array_in adding horizontal and vertical marks for ticks:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#		print "$i:\n";

		# Define variables for each size marker:
		my $var_1Mb		= 10**6;
		my $var_100kb	= 10**5;
		my $var_10kb	= 10**4;

		# Add first tick:
		my $var_string		= "chr0 1 1 ";
		$var_string			.= "stroke_color=$var_colour,";
		$var_string			.= "stroke_thickness=$var_thickness,";

		unless ( scalar @array_in == 2 && $i == 1 ) {

			$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
			$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

		}

		else {

			$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
			$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

		}

		push @{$array_out[$i]}, $var_string;

		# Iterate through regions in @array_in making up vertical ticks:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

#			print "\t$j: @{$array_in[$i][$j]}\n";

			# Define original locations:
			my $var_loc0		= $array_in[$i][$j][1];
			my $var_loc1		= $array_in[$i][$j][2];

			# Move on if region is a blank:
			next if $var_loc0 == 0;

			# Define adjusted locations:
			my $var_loc2		= $array_in[$i][$j][11];
			my $var_loc3		= $array_in[$i][$j][12];

			# Define inverted status of region:
			my $var_inverted	= $array_in[$i][$j][8];

			# Define difference between original and adjusted locations:
			my $var_diff	= $var_loc2 - $var_loc0;
		 
			# Handle ticks every 10 kb:
			if ( $var_10kb >= $var_loc0 && $var_10kb < $var_loc1 ) {

				my $var_temp		= $var_10kb + $var_diff;

#				print "\t\t\t$var_10kb - $var_temp\n";

				$var_10kb += 10**4;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_10kb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}

			elsif ( $var_loc0 > $var_10kb ) {

				$var_10kb += 10**4;

				$j--;

				next;

			}

			# Handle ticks every 100kb:
			if ( $var_100kb >= $var_loc0 && $var_100kb < $var_loc1 ) {

				my $var_temp	= $var_100kb + $var_diff;

#				print "\t\t$var_100kb - $var_temp\n";

				$var_100kb += 10**5;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_100kb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}

			# Handle ticks every 1 Mb:
			if ( $var_1Mb >= $var_loc0 && $var_1Mb < $var_loc1 ) {

				my $var_temp	= $var_1Mb + $var_diff;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$var_1Mb += 10**6;

				$j--;

			}

		}

	 	$var_radius0	-= $var_increment;

	}

	# ----------------------------------

	# Iterate through @array_in adding horizontal and vertical marks for ticks:
#	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variables for each size marker:
#		my $var_1Mb		= 10**6;
#		my $var_100kb	= 10**5;
#		my $var_10kb	= 10**4;

		# Add first tick:
#		my $var_string		= "chr0 1 1 ";
#		$var_string			.= "stroke_color=$var_colour,";
#		$var_string			.= "stroke_thickness=$var_thickness,";

#		unless ( scalar @array_in == 2 && $i == 1 ) {

#			$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
#			$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

#		}

#		else {

#			$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
#			$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

#		}

#		push @{$array_out[$i]}, $var_string;

		# Iterate through regions in @array_in making up vertical ticks:
#		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define original locations:
#			my $var_loc0		= $array_in[$i][$j][1];
#			my $var_loc1		= $array_in[$i][$j][2];

			# Move on if region is a blank:
#			next if $var_loc0 == 0;

			# Define adjusted locations:
#			my $var_loc2		= $array_in[$i][$j][11];
#			my $var_loc3		= $array_in[$i][$j][12];

			# Define inverted status of region:
#			my $var_inverted	= $array_in[$i][$j][8];

			# Define difference between original and adjusted locations:
#			my $var_diff	= $var_loc2 - $var_loc0;

			# Handle ticks every 10 kb:
#			if ( $var_10kb >= $var_loc0 && $var_10kb < $var_loc1 ) {

#				my $var_temp	= $var_10kb + $var_diff;

#				$var_10kb += 10**4;

#				my $var_string		= "chr0 $var_temp $var_temp ";

#				$var_string			.= "offset=$var_span0," if $var_inverted;
#				$var_string			.= "stroke_color=$var_colour,";
#				$var_string			.= "stroke_thickness=$var_thickness,";

#				unless ( scalar @array_in == 2 && $i == 1 ) {

#					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

#				}

#				else {

#					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_10kb_height$p-$var_gap$p";

#				}

#				push @{$array_out[$i]}, $var_string;

#				$j--;

#				next;

#			}		 

			# Handle ticks every 100kb:
#			if ( $var_100kb >= $var_loc0 && $var_100kb < $var_loc1 ) {

#				my $var_temp	= $var_100kb + $var_diff;

#				$var_100kb += 10**5;

#				my $var_string		= "chr0 $var_temp $var_temp ";

#				$var_string			.= "offset=$var_span0," if $var_inverted;
#				$var_string			.= "stroke_color=$var_colour,";
#				$var_string			.= "stroke_thickness=$var_thickness,";

#				unless ( scalar @array_in == 2 && $i == 1 ) {

#					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

#				}

#				else {

#					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_100kb_height$p-$var_gap$p";

#				}

#				push @{$array_out[$i]}, $var_string;

#				$j--;

#				next;

#			}		 

			# Handle ticks every 1 Mb:
#			if ( $var_1Mb >= $var_loc0 && $var_1Mb < $var_loc1 ) {

#				my $var_temp	= $var_1Mb + $var_diff;

#				my $var_string		= "chr0 $var_temp $var_temp ";

#				$var_string			.= "offset=$var_span0," if $var_inverted;
#				$var_string			.= "stroke_color=$var_colour,";
#				$var_string			.= "stroke_thickness=$var_thickness,";

#				unless ( scalar @array_in == 2 && $i == 1 ) {

#					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

#				}

#				else {

#					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
#					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

#				}

#				push @{$array_out[$i]}, $var_string;

#				$var_1Mb += 10**6;

#				$j--;

#			}

#		}

#	 	$var_radius0	-= $var_increment;

#	}

	# ----------------------------------

	# Iterate through @array_out and write data to @array_out:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my @array_temp	= @{$array_out[$i]};

		unshift @array_temp, "# $file_out.$i";

		# Print @array_out to $file_out:
		General::ArrayToFile ( \@array_temp, "$file_out.$i" );

	}

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

1;
