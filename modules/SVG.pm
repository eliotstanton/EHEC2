package SVG;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# File name:		SVG.pm
# Date created:		26 August, 2019
# Last modified:	26 August, 2019
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to creating SVG
#					images.

# ---------------------------------------------------------------------------- #
	
# Subroutine name:	Rectangle
# Description:		This subroutine returns SVG formatted rectangles.

# ----------------------------------------------------------

sub Rectangle {

	# Arguments:
	my ( $var_x, $var_y, $var_length, $var_height, $var_colour, $var_stroke )	= @_;

	# ----------------------------------

	my $var_return	= "<rect x=\"$var_x\" y=\"$var_y\" ";
	$var_return		.= "width=\"$var_length\" height=\"$var_height\" ";
	$var_return		.= "style=\"fill:rgb($var_colour);stroke-width:3;stroke:rgb($var_stroke)\" />";

	# ----------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	RegionsSVG
# Description:		This subroutine handles formatting chromosomal regions as
#					SVG formatted rectangles.

# ----------------------------------------------------------

sub RegionsSVG {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{ array_merged }};
	my @array_out;

	# Variables:
	my $var_height		= $hash_out{var_span};
	my $var_factor		= $hash_out{var_factor};
	my $var_y			= 100;
	my $var_increment	= 200;
	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";

	# ----------------------------------

	# Iterate through @array_merged:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][6];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc0		= $array_merged[$i][$j][11];
			my $var_loc1		= $array_merged[$i][$j][12];

			$var_y += $var_height if $var_inverted;

			next if $var_loc0 	== 0;

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			my $var_colour	= Circos::ColoursNew ( "chromosome", $var_copy, $var_scalar );

			my $var_string	= Rectangle ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

			push @{$array_out[$i]}, $var_string;

			$var_y -= $var_height if $var_inverted;

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_regions }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FeaturesSVG
# Description:		This subroutine handles formatting genomic regions as
#					SVG formatted rectangles.

# ----------------------------------------------------------

sub FeaturesSVG { 

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @{$hash_out{ array_features }};
	my @array_out;

	# Variables:
	my $var_height		= $hash_out{var_span};
	my $var_height2		= $hash_out{var_span}/2;
	my $var_factor		= $hash_out{var_factor};
	my $var_y			= 100;
	my $var_increment	= 200;
	my $var_scalar		= scalar @array_features;
	my $var_stroke		= "0,0,0";

	# ----------------------------------

	# Iterate through @array_merged:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		next unless $array_features[$i];

		# Iterate down through features for each sequence:
		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			# Define start/end locations of genomic feature:
			my $var_loc0		= $array_features[$i][$j][0];
			my $var_loc1		= $array_features[$i][$j][1];

			# Define the type of genomic feature, inverted status and copy:
			my $var_type		= $array_features[$i][$j][2];
			my $var_inverted	= $array_features[$i][$j][3];
			my $var_copy		= $array_features[$i][$j][4];
			my $var_orient		= $array_features[$i][$j][5] || "0";

			# Calculate length of genomic feature:
			my $var_length		= $var_loc1 - $var_loc0 + 1;

			# Change $var_height if feature is inverted:
			$var_y += $var_height if $var_inverted;

			# Change height if feature is reverse oriented:
			$var_y += $var_height2 if $var_orient eq "-";

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			# Calculate genomic feature colour:
			my $var_colour	= Circos::ColoursNew ( $var_type, $var_copy, $var_scalar );

			if ( $var_type eq "gene" ) {

				# Create string holding SVG formatted region:
				my $var_string	= Rectangle ( $var_loc0, $var_y, $var_length, $var_height2, $var_colour, $var_stroke );

				# Store $var_string in @array_out:
				push @{$array_out[$i]}, $var_string;

#				print "$var_colour $var_string\n";

			}

			else {

				# Create string holding SVG formatted region:
				my $var_string	= Rectangle ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

				# Store $var_string in @array_out:
				push @{$array_out[$i]}, $var_string;

#				print "$var_string\n";

			}

			# Decrease $var_height if handling an inverted region:
			$var_y -= $var_height if $var_inverted;

			# Change height if feature is reverse oriented:
			$var_y -= $var_height2 if $var_orient eq "-";

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_SVGfeatures }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	TicksSVG
# Description:		This subroutine handles adding tick marks to indicate 
#					lenght of chromosome.

# ----------------------------------------------------------

sub TicksSVG { 

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{ array_merged }};
	my @array_out;

	# Variables:
	my $var_thick		= $hash_out{stroke_thickness}/2;
	my $var_height		= $hash_out{var_span};
	my $var_factor		= $hash_out{var_factor};
	my $var_increment	= 200;

	my $var_gap			= $hash_out{var_gap};
	my $var_y			= 100 - $var_gap;

	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";
	my $var_colour		= "0,0,0";

	$var_scalar			= scalar @{$array_merged[0]};
	my $var_length		= $array_merged[0][$var_scalar-1][2];

	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $var_1kb_height		= $hash_out{var_1kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_1Mb				= 10**6;
	my $var_100kb			= 10**5;
	my $var_10kb			= 10**4;
	my $var_1kb				= 10**3;

	# ----------------------------------

	# Create horizontal lines:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][6];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc0		= $array_merged[$i][$j][11];
			my $var_loc1		= $array_merged[$i][$j][12];

			next if $var_loc0 == 0;

			$var_y += $var_height if $var_inverted;

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			my $var_string	= Rectangle ( $var_loc0, $var_y, $var_length, $var_thick, $var_colour, $var_stroke );

			push @{$array_out[$i]}, $var_string;

			$var_y -= $var_height if $var_inverted;

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	$var_y			= 100 - $var_gap;

	# Create vertical lines:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		my $var_length	= 0;

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_loc1	= $array_merged[$i][$j][2];

			$var_length = $var_loc1 if $var_length < $var_loc1;

		}

		# Iterate through by smallest increment for each individual tick:
		for ( my $j = 0; $j <= $var_length; $j+=$var_1kb ) {

			# Define variables for use outside of loop immediately below:
			my $var_diff		= 0;
			my $var_inverted	= 0;

			# Read through each region and find where to place tick:
			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				# Define locations:
				my $var_loc0	= $array_merged[$i][$k][1];
				my $var_loc1	= $array_merged[$i][$k][2];
			 	$var_inverted	= $array_merged[$i][$k][8];
				my $var_loc2	= $array_merged[$i][$k][11];

				# If start of region contains next hit define difference between
				# local and global location values:
				if ( $var_loc0 < $j && $var_loc1 > $j ) {

					$var_diff	= $var_loc2 - $var_loc0;

					last;

				}

			}

			my $var_x		= ( $j + $var_diff )/$var_factor;

			my $var_tick_height	= $var_1kb_height;

			$var_tick_height	= $var_10kb_height if $j % ( $var_10kb ) == 0;
			$var_tick_height	= $var_100kb_height if $j % ( $var_100kb ) == 0;
			$var_tick_height	= $var_1Mb_height if $j % ( $var_1Mb ) == 0;

			$var_y += $var_height if $var_inverted;

			$var_y			-= $var_tick_height;

			my $var_string	= Rectangle ( $var_x, $var_y, $var_thick,
									 $var_tick_height, $var_colour, $var_stroke );

			$var_y			+= $var_tick_height;

			$var_y -= $var_height if $var_inverted;

#			if ( $var_factor == 1000 ) {

#				next if $var_tick_height == $var_1kb_height;

#				print "$var_string\n";

				push @{$array_out[$i]}, $var_string;

#			}		

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_SVGticks }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	MergeSVG
# Description:		This subroutine handles merging SVG formatted data for
#					printing to file.

# ----------------------------------------------------------

sub MergeSVG {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_regions	= @{$hash_out{array_regions}};
	my @array_features	= @{$hash_out{array_SVGfeatures}};
	my @array_ticks		= @{$hash_out{array_SVGticks}};
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_out;

	# Variables:
	my $var_scalar			= scalar @array_regions;
	my $file_out			= $hash_out{file_svg};
	my $var_image_radius	= $hash_out{ var_image_radius };
	my $var_increment		= $hash_out { var_increment };
	my $var_factor			= $hash_out { var_factor };
	my $var_height			= 400;

	$var_scalar				= scalar @{$array_merged[0]};
	my $var_width			= ( $array_merged[0][$var_scalar-1][2] ) / $var_factor;

	# ----------------------------------

	# Add header information:
	push @array_out, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

	# Add background:
	push @array_out, "<g id=\"bg\">";
	my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
	$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";
	push @array_out, $var_string;
	push @array_out, "</g>";

	# ----------------------------------

	# Add regions from @array_regions:
	for ( my $i = 0; $i < scalar @array_regions; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_regions[$i]}; $j++ ) {

			push @array_out, $array_regions[$i][$j];

		}

	}

	# ----------------------------------

	# Add features from @array_features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# TODO: Add tags:

		next unless $array_features[$i];

		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			push @array_out, $array_features[$i][$j];

		}

	}

	# ----------------------------------

	# Add ticks from @array_ticks:
	for ( my $i = 0; $i < scalar @array_ticks; $i++ ) {

		# TODO: Add tags:

		for ( my $j = 0; $j < scalar @{$array_ticks[$i]}; $j++ ) {

			push @array_out, $array_ticks[$i][$j];

		}

	}

	# ----------------------------------

	push @array_out, "</svg>";

	# ----------------------------------

	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
