package Circos;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# File name:		Circos.pm
# Date created:		27 October, 2018
# Last modified:	31 October, 2018
# Created by:		Eliot Stanton (estanton@wisc.edu)

# Description:		This package contains subroutines specific to creating
#					Circos files.

# ---------------------------------------------------------------------------- #

# Subrotuine name:	Colours
# Description:		This subroutine stores colors used in other subroutines
#					used to create images using Circos.

# ----------------------------------------------------------

sub Colours {

	# Data structures:
	my %hash_out;

	# ----------------------------------

	$hash_out { chromosome_border }	= "grey";
	$hash_out { prophage_border }	= "dred";
	$hash_out { PLE_border }		= "dorange";
	$hash_out { IS_border }			= "dblue";
	$hash_out { rRNA_border }		= "dgreen";
	$hash_out { tRNA_border }		= "lgreen";
	$hash_out { ticks }				= "black";
	$hash_out { rhs_border }		= "purple";


	# ----------------------------------
	# Return %hash_out and end subroutines:
	return \%hash_out;

}

sub ColoursNew {

	# Arguments:
	my ( $var_type, $var_copy, $var_scalar )	= @_;

	# Variables:
	my $var_red		= 0;
	my $var_green	= 0;
	my $var_blue	= 0;
	my $var_return;

	# ----------------------------------

	if ( $var_type eq "chromosome" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 200;
			$var_green	= 200;
			$var_blue	= 200;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 220;
			$var_green	= 220;
			$var_blue	= 220;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "prophage" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 203;
			$var_green	= 24;
			$var_blue	= 29;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 220;
			$var_green	= 90;
			$var_blue	= 70;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "PLE" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 241;
			$var_green	= 105;
			$var_blue	= 19;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 253;
			$var_green	= 160;
			$var_blue	= 70;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "IS" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 3;
			$var_green	= 78;
			$var_blue	= 123;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 5;
			$var_green	= 112;
			$var_blue	= 176;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "rRNA" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 35;
			$var_green	= 132;
			$var_blue	= 67;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 65;
			$var_green	= 171;
			$var_blue	= 93;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "tRNA" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 112;
			$var_green	= 214;
			$var_blue	= 123;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 112;
			$var_green	= 214;
			$var_blue	= 123;

		}

	}

	# ----------------------------------

	elsif ( $var_type eq "rhs" ) {

		if ( $var_copy == $var_scalar ) {

			$var_red	= 180;
			$var_green	= 50;
			$var_blue	= 170;

		}

		elsif ( $var_copy != $var_scalar ) {

			$var_red	= 210;
			$var_green	= 90;
			$var_blue	= 200;

		}

	}

	# ----------------------------------

	# Set red, green, and blue values in $var_return:
	$var_return	= "$var_red,$var_green,$var_blue";

	# ----------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Ideogram
# Description:		This subroutine produced a file with Circos ideogram data.

# ----------------------------------------------------------

sub Ideogram {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { "file_ideogram" };
	my $var_chr_border	= $hash_colours { chromosome_border };
	my $var_stroke		= $hash_in { "stroke_thickness" };
	my $var_radius0		= $hash_in { var_radius0 };
	my $var_radius1		= $hash_in { var_radius1 };
	my $var_span		= ($var_radius0 - $var_radius1)/2;
	my $r				= "r";

	# ----------------------------------

	push @array_out, "#$file_out";

	push @array_out, "<ideogram>";
	push @array_out, "<spacing>\n\tdefault = 0.005r\n\t</spacing>";
	push @array_out, "\tradius	= $var_radius0$r";
	push @array_out, "\tthickness	= $var_span$r\n\tfill	= yes";
	push @array_out, "\tstroke_color = $var_chr_border";
	push @array_out, "\tstroke_thickness = $var_stroke";
	push @array_out, "</ideogram>";

	# Write @array_out to file:
	General::ArrayToFile ( \@array_out, $file_out );
	# ----------------------------------
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Image
# Description:		This subroutine produces a file containing image data for
#					Circos.

# ----------------------------------------------------------

sub Image {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in		= %$hash_in;
	my @array_out;

	# Variables:
	my $file_out	= $hash_in { "file_image" };
	my $file_png	= $hash_in { "file_png" };
	my $var_radius	= $hash_in { var_image_radius };
	my $p			= "p";

	# ----------------------------------
	push @array_out, "# $file_out";
	push @array_out, "<image>";

	push @array_out, "dir                = .";
	push @array_out, "file               = $file_png";
	push @array_out, "svg                = yes";
	push @array_out, "\tradius            = $var_radius$p";
	push @array_out, "angle_offset       = -90";
	push @array_out, "auto_alpha_colors  = yes";
	push @array_out, "auto_alpha_steps   = 5";
	push @array_out, "background         = white";

	push @array_out, "<\/image>";

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );
	# ----------------------------------
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Run
# Description:		This subroutine runs Circos.

# ----------------------------------------------------------

sub Run {

	# Arguments:
	my ( $hash_in )	= @_;
	# Data structures:
	my %hash_in		= %$hash_in;

	# Variables:
	my $file_config	= $hash_in { "file_config" };
	my $file_png	= $hash_in { "file_png" };

	# ----------------------------------

	my $var_return	= system ( "circos -conf $file_config -silent" );

	system ( "xviewer $file_png" ) unless $var_return;

	# ----------------------------------
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
