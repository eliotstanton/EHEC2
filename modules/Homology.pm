package Homology;

use strict;
use warnings;
use JSON::XS;

# ---------------------------------------------------------------------------- #

# File name:		Homology.pm
# Date created:		27 October, 2018
# Last modified:	18 November, 2018
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to calculating
#					homology.

# ---------------------------------------------------------------------------- #

# Subroutine name:	CircosConfig
# Description:		This subroutine creates a config file for Circos

# ----------------------------------------------------------

sub CircosConfig {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { "file_config" };
	my $file_direct_links	= $hash_in { "file_direct_links" };
	my $file_features		= $hash_in { "file_features" };
	my $file_ideogram		= $hash_in { "file_ideogram" };
	my $file_inverted_links	= $hash_in { "file_inverted_links" };
	my $file_image			= $hash_in { "file_image" };
	my $file_karyotype		= $hash_in { "file_karyotype" };
	my $file_ticks			= $hash_in { "file_ticks" };
	my $var_radius0			= $hash_in { var_radius0 };
	my $var_radius1			= $hash_in { var_radius1 };
	my $var_gap				= $hash_in { var_gap };
	my $r					= "r";
	my $p					= "p";


	# ----------------------------------

	push @array_out, "#$file_out";

	push @array_out, "karyotype = $file_karyotype";
	push @array_out, "<<include $file_ideogram>>";
	push @array_out, "<ideogram>\nshow = no\n</ideogram>";

	push @array_out, "<highlights>";
	push @array_out, "\t<highlight>";
	push @array_out, "\t\tfile	= $file_features\nideogram	= no";
	push @array_out, "\t\tr0	= $var_radius0$r\n\t\tr1	= $var_radius1$r";
 	push @array_out, "\t</highlight>";
	push @array_out, "\t<highlight>";
	push @array_out, "\t\tfile	= $file_ticks\n\t\tideogram	= no";
	push @array_out, "\t</highlight>";
	push @array_out, "</highlights>";

	push @array_out, "<links>";
	push @array_out, "\tradius	= $var_radius1$r-$var_gap$p\nribbon	= yes";
 	push @array_out, "<link>\nfile	= $file_direct_links\n</link>" unless $hash_in{d};
 	push @array_out, "<link>\nfile	= $file_inverted_links\n</link>" unless $hash_in{i};
 	push @array_out, "</links>";

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

# Subroutine name:	CircosFeatures
# Description:		This subroutine handles writing highlight files representing
#					genomic features.

# ----------------------------------------------------------

sub CircosFeatures {

	# Arguments:
	my ( $array_in, $array_FASTA, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my @array_FASTA		= @$array_FASTA;
	my %hash_in			= %$hash_in;
	my %hash_colours	= %{Circos::Colours( )};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { file_features };
	my $var_thickness	= $hash_in { stroke_thickness };
	my $var_length		= length $array_FASTA[0][2];
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 2 );
	my $var_stroke		= $hash_colours { "chromosome_border" };

	# ----------------------------------

	# Store file name at beginning of file:
	push @array_out, "#$file_out";

	# Create standing ideogram block(s):
	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

		my $var_length	= length $array_FASTA[$i][2];

		my $var_string	= "chr$i 1 $var_length fill_color=$var_colour";
		$var_string		.= ",stroke_color=$var_stroke";
		$var_string		.= ",stroke_thickness=$var_thickness";

		push @array_out, $var_string;

	}

	# Iterate through features in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_seq		= $array_in[$i][0];
		my $var_type	= $array_in[$i][1];
		my $var_loc0	= $array_in[$i][3];
		my $var_loc1	= $array_in[$i][4];

		my $var_colour	= Circos::ColoursNew ( $var_type, 1, 2 );
		my $var_stroke	= $hash_colours { "$var_type\_border" };

		$var_colour		= "0,0,0" unless $var_colour;
		$var_stroke		= "0,0,0" unless $var_stroke;

		my $var_string	= "chr$var_seq $var_loc0 $var_loc1 ";
		$var_string		.= "fill_color=$var_colour";
		$var_string		.= ",stroke_color=$var_stroke";
		$var_string		.= ",stroke_thickness=0";

#		print "$i: $var_string - $var_type\n";

		push @array_out, $var_string;

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	CircosKaryotype
# Description:		This subroutine produces a formatted Circos karyotype file.

# ----------------------------------------------------------

sub CircosKaryotype {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @$array_in;
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { "file_karyotype" };
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 2 );

	# ----------------------------------

	push @array_out, "#$file_out";

	# Iterate through @array_FASTA writing data to @array_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_header		= $array_in[$i][0];

		my $var_length		= length $array_in[$i][2];

		my $var_string	= "chr - chr$i $i 0 $var_length $var_colour";

		push @array_out, $var_string;

	}

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	CircosLinks
# Description:		This subroutine produces a file with Circos links data.

# ----------------------------------------------------------

sub CircosLinks {

	# Arguments:
	my ( $var_variable, $array_features, $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_features	= @$array_features;
	my @array_in;
	my %hash_colours	= %{Circos::Colours( )};
	my @array_out;

	# Variables:
#	my $var_link_width	= $hash_in { "var_link_width" };
	my $var_thickness	= $hash_in { stroke_thickness_links };
	my $var_min			= $hash_in { var_minimum };
	my $file_in;	
	my $file_out;

	# ----------------------------------

	# Define path for $file_in:
	if ( $var_variable eq "direct" ) {

		$file_in	= $hash_in { file_direct_features };
		$file_out	= $hash_in { file_direct_links };

		$var_variable = 1;

	}

	if ( $var_variable eq "inverted" ) {

		$file_in	= $hash_in { file_inverted_features };
		$file_out	= $hash_in { file_inverted_links };

		$var_variable = (-1);

	}

	@array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	# Iterate through @array_in, merge regions and break if regions don't match:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define sequence numbers:
		my $var_seq0	= $array_in[$i][0];
		my $var_seq1	= $array_in[$i][3];

		# Define the start and end locations of link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Define classID of each repeat:
		# Define class of end locations of link:
		my $var_class0	= $array_in[$i][6];
		my $var_class1	= $array_in[$i][7];		

#		print "$i: @{$array_in[$i]}\n";

		for ( my $j = $i+1; $j < scalar @array_in; $j++ ) {

			# Define sequence numbers:
			my $var_seq2	= $array_in[$j][0];
			my $var_seq3	= $array_in[$j][3];

			# Define the start and end locations of link:
			my $var_loc4	= $array_in[$j][1];
			my $var_loc5	= $array_in[$j][2];
			my $var_loc6	= $array_in[$j][4];
			my $var_loc7	= $array_in[$j][5];

			# Define classID of each repeat:
			# Define class of end locations of link:
			my $var_class2	= $array_in[$j][6];
			my $var_class3	= $array_in[$j][7];		

			last if $var_loc4 > $var_loc0 + 1;

			if ( $var_loc0 + 1 == $var_loc4 && $var_loc1 + 1 == $var_loc5 ) {

				if ( $var_loc2 + $var_variable == $var_loc6 && $var_loc3 + $var_variable == $var_loc7 ) {

					if ( $var_class0 ne $var_class2 || $var_class1 ne $var_class3 ) {

#						print "\t$j: @{$array_in[$j]}\n";

						next;

					}

					# Redefine locations in original element:
					$array_in[$i][2]	= $var_loc5;
					$array_in[$i][5]	= $var_loc7 if $var_variable == 1;
					$array_in[$i][4]	= $var_loc6 if $var_variable == -1;

					# Remove adjacent link from @array_in:
					splice @array_in, $j, 1;

					$j--;

					# Update location variables for next iteration:
					$var_loc0	= $var_loc4;
					$var_loc1	= $var_loc5;
					$var_loc2	= $var_loc6;
					$var_loc3	= $var_loc7;

				}

			}

		}

	}

	# ----------------------------------

	# Remove regions smaller in length than $var_min and resolve discrepancies
	# in classification of repeat regions:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define sequence numbers:
		my $var_seq0	= $array_in[$i][0];
		my $var_seq1	= $array_in[$i][3];

		# Define the start and end locations of link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Define classID of each repeat:
		# Define class of end locations of link:
		my $var_class0	= $array_in[$i][6];
		my $var_class1	= $array_in[$i][7];	

		# Define length of region:
		my $var_length	= $var_loc1 - $var_loc0 + 1;

		# Remove region if it is smaller than $var_min:
		if ( $var_length < $var_min ) {

			splice @array_in, $i , 1;

			$i--;

			next;

		}

		if ( $var_class0 ne $var_class1 ) {	

			if ( $var_class0 eq "chromosome" || $var_class1 eq "chromosome" ) {

				$array_in[$i][6]	= $var_class1 if $var_class0 eq "chromosome";
				$array_in[$i][7]	= $var_class0 if $var_class1 eq "chromosome";

			}

			elsif ( $var_class0 eq "IS" || $var_class1 eq "IS" ) {

				$array_in[$i][6]	= $var_class1 if $var_class0 ne "IS";
				$array_in[$i][7]	= $var_class0 if $var_class1 ne "IS";

			}

			elsif ( $var_class0 eq "prophage" || $var_class1 eq "prophage" ) {

				$array_in[$i][6]	= $var_class1 if $var_class0 eq "PLE";
				$array_in[$i][7]	= $var_class0 if $var_class1 eq "PLE";

			}

		}

	}


	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#		print "$i: @{$array_in[$i]}\n";

		# Define sequence numbers:
		my $var_seq0	= $array_in[$i][0];
		my $var_seq1	= $array_in[$i][3];

		# Define the start and end locations of link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Define class of end locations of link:
		my $var_class0	= $array_in[$i][6];
		my $var_class1	= $array_in[$i][7];

		if ( $var_class0 ne $var_class1 ) {

			print "$i: @{$array_in[$i]}\n";

		}

		# Define length of repeat:
		my $var_length	= $var_loc1 - $var_loc0 - 1;

		# Move on if $var_length is less than $var_link_width:
#		next if $var_length < $var_link_width;

		# Define variables for use in writing link:
		my $var_colour;
		my $var_stroke;

		$var_colour	= Circos::ColoursNew ( $var_class0, 1, 2 );
		$var_stroke	= $hash_colours { "$var_class0\_border" };

		$var_colour		= "0,0,0" unless $var_colour;
		$var_stroke		= "0,0,0" unless $var_stroke;

		# Create string holding data for Circos:
		my $var_string	= "chr$var_seq0 $var_loc0 $var_loc1 ";
		$var_string		.= "chr$var_seq1 $var_loc2 $var_loc3 ";
		$var_string		.= "color=$var_colour";

		# Push $var_string to @array_out for writing to file:
		push @array_out, $var_string;

	}

	# ----------------------------------

	# Write @array_out_links to file:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Ticks
# Description:		This subroutine takes region data and creates tick marks.

# ----------------------------------------------------------

sub CircosTicks {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_in			= @$array_in;
	my %hash_colours		= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $var_length0			= 1;
	my $var_length1			= 0;
	my $file_out			= $hash_out { file_ticks };
	my $var_thickness		= $hash_out { stroke_thickness };
	my $var_colour			= $hash_colours { ticks };
	my $var_gap				= $hash_out { var_gap };
	my $var_radius0			= $hash_out { var_radius0 };
	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_1Mb				= 10**6;
	my $var_100kb			= 10**5;
	my $var_10kb			= 10**4;

	# ----------------------------------

	# Add file name to beginning of file:
	push @array_out, "#$file_out";

	# Add ring highlight(s):
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_length	= length $array_in[$i][2];

		my $var_string	= "chr$i 1 $var_length ";
		$var_string		.= "stroke_color=$var_colour,";
		$var_string		.= "stroke_thickness=$var_thickness,";
		$var_string		.= "r0=$var_radius0$r+$var_gap$p,";
		$var_string		.= "r1=$var_radius0$r+$var_gap$p";

#		push @array_out, $var_string;

	}

	# ----------------------------------

	# Add first tick:
	my $var_string		= "chr0 1 1 ";
	$var_string			.= "stroke_color=$var_colour,";
	$var_string			.= "stroke_thickness=$var_thickness,";
	$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
	$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

	push @array_out, $var_string;

	# Iterate through @array_in adding vertical marks for ticks:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_length2	= length $array_in[$i][2];

		$var_length1 	+= $var_length2;

		my $var_loc		= 1;

		for ( my $j = $var_length0; $j <= $var_length1; $j++ ) {

			if ( $j % $var_10kb == 0 ) {

				my $var_string	= "chr$i $var_loc $var_loc ";
				$var_string		.= "stroke_color=$var_colour,";
				$var_string		.= "stroke_thickness=$var_thickness,";

				if ( $j % $var_1Mb == 0 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

				} 

				elsif ( $j % $var_100kb == 0 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

				} 

				else {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

				} 

				push @array_out, $var_string;

			}

			$var_loc++;

		}

		$var_length0	= $var_length1 + 1;

	}

	# ----------------------------------

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, "$file_out" );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		Links
# Description:		This subroutine creates an array of formatted direct and 
#					inverted links.

# ----------------------------------------------------------

sub Links {

	# Arguments:
	my ( $hash_in, $hash_nmer, $var_nmer )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_nmer		= %$hash_nmer if $hash_nmer;
	my @array_direct;
	my @array_inverted;

	# Variables:
	my $file_direct		= $hash_out { file_direct };
	my $file_inverted	= $hash_out { file_inverted };

	# ----------------------------------

	# End subroutine early if the output files are already present:
	if ( -e $file_direct && $file_inverted ) {

		print "SKIPPING: $file_direct and $file_inverted present\n";

		@array_direct	= @{ General::FileToArray ( $file_direct ) };
		@array_inverted	= @{ General::FileToArray ( $file_inverted ) };

		my $var_scalar	= scalar @array_direct;

		print "\t- $var_scalar direct links\n";

		$var_scalar	= scalar @array_inverted;

		print "\t- $var_scalar inverted links\n";

		# End subroutine:
		return;

	}

	# ----------------------------------

	# Iterate through each sequence stored in %hash_nmer and process locations:
	foreach my $var_seq ( keys %hash_nmer ) {

		my @array_tags	= @{$hash_nmer { $var_seq }};

		# Proceed if more than one location:
		if ( scalar @array_tags > 1 ) {

			# Iterate through locations and generate links:
			for ( my $i = 0; $i < scalar @array_tags - 1; $i++ ) {

				my $var_seq0	= $array_tags[$i][0];
				my $var_loc0	= $array_tags[$i][1];
				my $var_loc1	= $array_tags[$i][1] + $var_nmer;

				for ( my $j = $i+1; $j < scalar @array_tags; $j++ ) {

					my $var_seq1	= $array_tags[$j][0];
					my $var_loc2	= $array_tags[$j][1];
					my $var_loc3	= $array_tags[$j][1] + $var_nmer;

					# Define temporary array:
					my @array_temp;

					# Add lower sequence first:
					if ( $var_seq0 != $var_seq1 ) {

						if ( $var_seq0 < $var_seq1 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					else {

						# Add lower location first:
						if ( $var_loc0 < $var_loc2 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					next if $array_temp[4] <= $array_temp[2];

					# Store link in @array_out:
					push @array_direct, \@array_temp;

				}

			}

		}

		# Generate reverse complement of $var_seq:
		my $var_reverse	= General::ReverseComplement ( $var_seq );

		# If $var_reverse is also present in %hash_in make links:
		if ( $hash_nmer { $var_reverse } ) {

			# Define both arrays of locations:
			my @array_tags1	= @{$hash_nmer{$var_seq}};
			my @array_tags2	= @{$hash_nmer{$var_reverse}};

			# Iterate through arrays making links:
			for ( my $i = 0; $i < scalar @array_tags1; $i++ ) {

				my $var_seq0	= $array_tags1[$i][0];
				my $var_loc0	= $array_tags1[$i][1];
				my $var_loc1	= $array_tags1[$i][1] + $var_nmer;

				for ( my $j = 0; $j < scalar @array_tags2; $j++ ) {

					my $var_seq1	= $array_tags2[$j][0];
					my $var_loc2	= $array_tags2[$j][1];
					my $var_loc3	= $array_tags2[$j][1] + $var_nmer;

					# Define temporary array:
					my @array_temp;

					# Add lower sequence first:
					if ( $var_seq0 != $var_seq1 ) {

						if ( $var_seq0 < $var_seq1 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					else {

						# Add lower location first:
						if ( $var_loc0 < $var_loc2 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					next if $array_temp[4] <= $array_temp[2] && $var_seq0 == $var_seq1;

					# Store link in @array_out:
					push @array_inverted, \@array_temp;

				}

			}

		}

		# Remove $var_seq from %hash_out:
		delete $hash_nmer { $var_seq };

	}

	# ----------------------------------

	# Sort elements in @array_direct by location:
	@array_direct	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_direct;

	# Sort elements in @array_inverted by location:
	@array_inverted	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_inverted;

	# ----------------------------------

	# Convert arrays in @array_temp into strings:
	for ( my $i = 0; $i < scalar @array_direct; $i++ ) {

		my $var_string	= join " ", @{$array_direct[$i]};

		$array_direct[$i]	= $var_string;

	}

	for ( my $i = 0; $i < scalar @array_inverted; $i++ ) {

		my $var_string	= join " ", @{$array_inverted[$i]};

		$array_inverted[$i]	= $var_string;

	}

	# ----------------------------------

	# Write @array_direct to $file_direct:
	General::ArrayToFile ( \@array_direct, $file_direct );

	# Write @array_inverted to $file_inverted:
	General::ArrayToFile ( \@array_inverted, $file_inverted );

	# ----------------------------------

	my $var_scalar	= scalar @array_direct;

	print "\t- $var_scalar direct links\n";

	$var_scalar	= scalar @array_inverted;

	print "\t- $var_scalar inverted links\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #


# Subroutine name:		FeatureData
# Description:			This subroutine attaches feature information to each
#						formatted direct and inverted link.

# --------------------------------------

sub FeatureData {

	# Arguments:
	my ( $hash_in, $var_variable, $array_features )	= @_;

	# Data structures:
	my %hash_out				= %$hash_in;
	my @array_features			= @$array_features;
	my @array_in;

	# Variables:
	my $file_in;
	my $file_out;

	# ----------------------------------

	# Define $file_in based upon $var_variable:
	if ( $var_variable eq "direct" ) {

		$file_in	= $hash_out { file_direct };
		$file_out	= $hash_out { file_direct_features };

	}

	if ( $var_variable eq "inverted" ) {

		$file_in	= $hash_out { file_inverted };
		$file_out	= $hash_out { file_inverted_features };

	}

	# If $file_out already exists:
	if ( -e $file_out ) {

		print "SKIPPING: $file_out present\n";

		# End subroutine:
		return;

	}

	# ----------------------------------

	# Import links from $file_in to @array_in:
	@array_in	= @{ General::FileToArray( $file_in ) };

	# ----------------------------------

	# Iterate through links in @array_in and assign identity for each repeat:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define locations of link:
		my $var_seq0	= $array_in[$i][0];
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_seq1	= $array_in[$i][3];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# If links overlap, remove them:
		if ( $var_loc2 <= $var_loc1 && $var_seq0 == $var_seq1 ) {

			splice @array_in, $i, 1;

			$i--;

			next;

		}

		# Define variables to hold the identity of each link:
		my $var_class0;
		my $var_class1;

		# Iterate through @array_features and classify genomic feature for each
		# link:
		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			# Define class of the genomic feature:
			my $var_class	= $array_features[$j][1];

			# Define locations of genomic feature:
			my $var_loc4	= $array_features[$j][3];
			my $var_loc5	= $array_features[$j][4];

			# Handle repeats on left edge of genomic feature:
			if ( $var_loc0 >= $var_loc4 && $var_loc1 <= $var_loc5 ) {

				$var_class0	= $var_class unless $var_class0;

			}

			# handle repeats on right edge of genomic feature:
			if ( $var_loc2 >= $var_loc4 && $var_loc3 <= $var_loc5 ) {

				$var_class1	= $var_class unless $var_class1;

			}

			# End loop if both $var_class0 and $var_class1 are defined:
			if ( $var_class0 && $var_class1 ) {

				last;

			}

		}

		# Set $var_class0 and $var_class1 to default to chromosome:
		$var_class0			= "chromosome" unless $var_class0;
		$var_class1			= "chromosome" unless $var_class1;

		# Store classes in array for repeat stored in @array_in:
		$array_in[$i][6]	= $var_class0 unless $array_in[$i][6];
		$array_in[$i][7]	= $var_class1 unless $array_in[$i][7];

	}

	# ----------------------------------

	# Save @array_in to $file_inverted_features:
	General::ArrayToFile ( \@array_in, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:		FlattenLinks
# Description:			This subroutine flattens pairs of repeats into repeat
#						regions.

# --------------------------------------

sub FlattenLinks {

	# Arguments:
	my ( $hash_in, $var_variable )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_in;
	my @array_out;

	# Variables:
	my $file_in;
	my $file_out;

	# ----------------------------------

	# Set filepath for $file_in and $file_out:
	if ( $var_variable eq "direct" ) {

		$file_in		= $hash_out { file_direct_merged };
		$file_out		= $hash_out { file_direct_flat };

	}

	elsif ( $var_variable eq "inverted" ) {

		$file_in		= $hash_out { file_inverted_merged };
		$file_out		= $hash_out { file_inverted_flat };

	}

	if ( $var_variable eq "combined" ) {

		$file_out		= $hash_out { file_flat };

	}

	# ----------------------------------

	# If $file_out already exists:
	if ( -e $file_out ) {

		print "SKIPPING: $file_out present\n";

		@array_out		= @{General::FileToArray ( $file_out )};

		my $var_scalar	= scalar @array_out;
		my $var_total	= 0;

		print "\t- $var_scalar flattened $var_variable regions\n";

		for ( my $i = 0; $i < scalar @array_out; $i++ ) {

			my $var_loc0	= $array_out[$i][0];
			my $var_loc1	= $array_out[$i][1];

			my $var_length	= $var_loc1 - $var_loc0 + 1;

			$var_total		+= $var_length;

		}	
		
		print "\t- Total length: $var_total bp\n";

		# Return @array_out and end subroutine:
		return \@array_out;

	}

	# ----------------------------------

	# Import data and split links for inverted or direct:
	if ( $var_variable eq "direct" || $var_variable eq "inverted" ) {

		# Import data from $file_in:
		@array_in	= @{General::FileToArray ( $file_in )};

		# Split the links:
		@array_out	= @{ SplitLinks ( \@array_in ) }; 

	}

	# If combining data import from direct and inverted files:
	else {

		my $file_direct		= $hash_out { file_direct_flat };
		my $file_inverted	= $hash_out { file_inverted_flat };

		my @array_direct	= @{General::FileToArray ( $file_direct )};
		my @array_inverted	= @{General::FileToArray ( $file_inverted )};

		@array_out			= ( @array_direct, @array_inverted );

	}

	# Sort @array_out by leading location:
	@array_out	= sort {

		$a -> [0] <=> $b -> [0] || $b -> [1] <=> $a -> [1]
 
		} @array_out;

	# ----------------------------------

	# Remove duplicates from @array_out:
	@array_out	= @{ RemoveDuplicates ( \@array_out )};

	# Merge overlapping regions
	@array_out	= @{ MergeOverlapping ( \@array_out )};

	# ----------------------------------

	my $var_scalar	= scalar @array_out;
	my $var_total	= 0;

	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my $var_loc0	= $array_out[$i][0];
		my $var_loc1	= $array_out[$i][1];

		my $var_length	= $var_loc1 - $var_loc0 + 1;

		$var_total		+= $var_length;

	}

	print "\t- $var_scalar flattened $var_variable regions\n";
	print "\t- Total length: $var_total bp\n";

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		MergeLinks
# Description:		This subroutine merges adjacent pairs of links.

# ----------------------------------------------------------

sub MergeLinks {

	# Arguments:
	my ( $hash_in, $var_variable, $var_min )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_in;
	my @array_out;

	# Variables:
	my $file_in;
	my $file_out;
	my $var_orientation	= $var_variable;

	# ----------------------------------

	# Set $var_variable equal to one (direct) or minus one (inverted):
	if ( $var_variable eq "direct" ) {

		$var_variable	= "1";

		$file_in		= $hash_out { file_direct };
		$file_out		= $hash_out { file_direct_merged };

	}

	if ( $var_variable eq "inverted" ) {

		$var_variable = "-1" ;

		$file_in		= $hash_out { file_inverted };
		$file_out		= $hash_out { file_inverted_merged };

	}

	# If $file_out already exists:
	if ( -e $file_out ) {

		print "SKIPPING: $file_out present\n";

		@array_out		= @{General::FileToArray ( $file_out )};

		my $var_scalar	= scalar @array_out;

		print "\t- $var_scalar merged $var_orientation links\n";

		# Return @array_out and end subroutine:
		return \@array_out;

	}

	@array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	my $var_scalar2	= scalar @array_in;

	print "$var_scalar2 unmerged $var_orientation links\n";

	# ----------------------------------

	# Sort elements in @array_in by location:
	@array_in	= sort { 

		$a -> [1] <=> $b -> [1] || $a -> [3] <=> $b -> [3] 

	} @array_in;

	# ----------------------------------

	# Iterate down through each link in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define first sequence ID and location of initial link:
		my $var_seq0	= $array_in[$i][0];
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		# Define second sequence ID and location of initial link:
		my $var_seq1	= $array_in[$i][3];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Iterate downward from link at hand finding adjacent links and merging:
		for ( my $j = $i + 1; $j < scalar @array_in; $j++ ) {

			# Define first sequence ID and location of new link:
			my $var_seq2	= $array_in[$j][0];
			my $var_loc4	= $array_in[$j][1];
			my $var_loc5	= $array_in[$j][2];

			# Define second sequence ID and location of new link:
			my $var_seq3	= $array_in[$j][3];
			my $var_loc6	= $array_in[$j][4];
			my $var_loc7	= $array_in[$j][5];

			# End loop if locations are larger than the original:
			last if $var_loc4 > $var_loc0 + 1;

			# End loop is sequence ID varies:
			last if $var_seq0 != $var_seq2;

			# If links are adjacent merge the data:
			if ( $var_loc0 + 1 == $var_loc4 && $var_loc1 + 1 == $var_loc5 ) {

				if ( $var_loc2 + $var_variable == $var_loc6 && $var_loc3 + $var_variable == $var_loc7 ) {

					# Redefine locations in original element:
					$array_in[$i][2]	= $var_loc5;
					$array_in[$i][5]	= $var_loc7 if $var_variable == 1;
					$array_in[$i][4]	= $var_loc6 if $var_variable == -1;

					# Remove adjacent link from @array_in:
					splice @array_in, $j, 1;

					$j--;

					# Update location variables for next iteration:
					$var_loc0	= $var_loc4;
					$var_loc1	= $var_loc5;
					$var_loc2	= $var_loc6;
					$var_loc3	= $var_loc7;

				}

			}

		}

	}

	# ----------------------------------

	# Store links greater than or equal to minimum length in @array_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end location of first link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		# Calculate length of link:
		my $var_length	= $var_loc1 - $var_loc0 + 1;

		# If link is long enough store it in @array_out:
		if ( $var_length >= $var_min ) {

			$array_in[$i][6]	= $var_length;

			push @array_out, $array_in[$i];

		}		

	}

	# ----------------------------------

	# Sort elements in @array_out by size:
	@array_out	= sort { $b -> [6] <=> $a -> [6] } @array_out;

	# Iterate through @array_out and remove links inside other links:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my $var_loc0	= $array_out[$i][1];
		my $var_loc1	= $array_out[$i][2];
		my $var_loc2	= $array_out[$i][4];
		my $var_loc3	= $array_out[$i][5];

		# Remove length information:
		$array_out[$i][6]	= "";

		for ( my $j = $i+1; $j < scalar @array_out; $j++ ) {

			my $var_loc4	= $array_out[$j][1];
			my $var_loc5	= $array_out[$j][2];
			my $var_loc6	= $array_out[$j][4];
			my $var_loc7	= $array_out[$j][5];	

			if ( $var_loc0 <= $var_loc4 && $var_loc1 >= $var_loc5 ) {

				if ( $var_loc2 <= $var_loc6 && $var_loc3 >= $var_loc7 ) {

					splice @array_out, $j, 1;

					$j--;

				}

			}

		}	

	}


	# ----------------------------------

	# Sort elements in @array_out by location and store them back in @array_in:
	@array_in	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_out;

	# Undefine @array_out:
	undef @array_out;

	# ----------------------------------

	# Remove links below minimum size in @array_in and store them in @array_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end of first link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		# Determine length of link:
		my $var_length	= $var_loc1 - $var_loc0 + 1;

		# If link is long enough, store it in @array_out:
		if ( $var_length >= $var_min ) {

			my $var_string	= join " ", @{$array_in[$i]};

			push @array_out, $var_string;

		}

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	my $var_scalar	= scalar @array_out;

	print "\t- $var_scalar merged $var_orientation links\n";

	# ----------------------------------

	# Return @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #


# Subroutine name:	SplitLinks
# Description:		This subroutine splits pairs of links into seperate regions.

# ---------------------------------------------------------------------------- #

sub SplitLinks {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Split links in @array_in into separate regions:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations of both ends of link:
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Create temporary arrays holding locations:
		my @array_temp0	= ( $var_loc0, $var_loc1 );
		my @array_temp1	= ( $var_loc2, $var_loc3 );

		# Store references to temporary arrays in @array_flat:
		push @array_out, \@array_temp0;
		push @array_out, \@array_temp1;

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		Nmer
# Description:		This subroutine stores each possible nmer length nucleotide
#					in a hash.

# ----------------------------------------------------------

sub Nmer {

	# Arguments:
	my ( $hash_in, $array_in, $var_nmer )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{$array_in};
	my %hash_out;

	# Variables:
	my $file_hash		= $hash_in { file_hash };
#	my $file_direct		= $hash_in { file_direct };
#	my $file_inverted	= $hash_in { file_inverted };
	my $var_null		= ("N")x$var_nmer;

	# ----------------------------------

	# Leave subroutine early if files containing direct and inverted repeat
	# locations is already present:
#	if ( -e $file_direct && $file_inverted ) {

#		print "SKIPPING: $file_direct and $file_inverted present\n";

		# End subroutine:
#		return;

#	}

	# ----------------------------------

	# If $file_hash exists load it into memory:
	if ( -e $file_hash ) {

		print "LOADING HASH from $file_hash\n";

		# Open file containing stored hash.
		open ( my $file_read, '<', $file_hash ) or die " ERROR: UNABLE TO OPEN $file_hash!\n";

		my $hash_json = <$file_read>;

		%hash_out = %{decode_json( $hash_json )};

		close $file_read or die " ERROR: UNABLE TO CLOSE $file_hash!\n";

		my $var_scalar	= scalar keys %hash_out;

		print "\t- $var_scalar sequences in hash (nmer=$var_nmer bp)\n";

		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			# Define variable holding FASTA header:
			my $var_header		= $array_in[$i][1];

			# Define variable to hold sequence:
			my $var_sequence	= $array_in[$i][2];

			# Define variable to hold length of sequence:
			my $var_length		= length $var_sequence;

			print "$i: $var_header - $var_length bp\n";

		}

		# End subroutine and return %hash_out:
		return \%hash_out;

	}

	# ----------------------------------

	# Iterate through each FASTA sequence and store nucleotide sequences in 
	# %hash_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variable holding FASTA header:
		my $var_header		= $array_in[$i][1];

		# Define variable to hold sequence:
		my $var_sequence	= $array_in[$i][2];

		# Define variable to hold length of sequence:
		my $var_length		= length $var_sequence;

		# Presize %hash_out to $var_length;
		keys ( %hash_out )	= $var_length;

		print "$i: $var_header - $var_length bp\n";

		my @array_temp;

		# Iterate through sequence character by character creating a hash of
		# sequences of $var_nmer length:
		for ( my $j = 0; $j < $var_length - $var_nmer + 1; $j++ ) {

			my $var_seq		= substr $var_sequence, $j, $var_nmer;

			my @array_temp	= ( $i, $j+1 );

			push ( @{ $hash_out { $var_seq } }, \@array_temp ) unless $var_seq eq $var_null;

		}

	}

	# ----------------------------------

	my $var_scalar	= scalar keys %hash_out;

	print "\t- $var_scalar sequences in hash (nmer=$var_nmer bp)\n";

	# ----------------------------------

	# Remove extraneous entries from the hash:
	while ( my ( $var_seq, $array_tags ) = each %hash_out) {

		# Define array holding locations:
		my @array_tags	= @$array_tags;

		# If only one location is present in @array_tags check to see if a
		# reverse complement is present in the hash:
		if ( scalar @array_tags == 1 ) {

			my $var_reverse	= General::ReverseComplement ( $var_seq );

			# Delete entry if reverse complement isn't present in the hash:
			delete $hash_out { $var_seq } unless $hash_out { $var_reverse };

		}

	}

	# ----------------------------------

	# Store data to file if not already present:
	unless ( -e $file_hash ) {

		print "SAVING\n";

		# Encode %hash_out in JSON format:
		my $hash_json = encode_json( \%hash_out );

		# Open $file_hash to write data:
		open ( my $file_write, '>', $file_hash ) or die " ERROR: UNABLE TO OPEN $file_hash!\n";

		# Print the hash to file.
		print $file_write $hash_json;

		# Close $file_out:
		close $file_write or die " ERROR: UNABLE TO CLOSE $file_hash!\n";

	}

	# ----------------------------------

	# Update the number of sequences in %hash_out:
	$var_scalar	= scalar keys %hash_out;

	# Report the number of sequences in %hash_out:
	print "\t- $var_scalar sequences in hash post-processing\n";

	# ----------------------------------

	# End subroutine and return %hash_out:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	ProcessFeatures
# Description:		This subroutine helps correcting bad coordinates.

# ---------------------------------------------------------------------------- #

sub ProcessFeatures {

	# Arguments:
	my ( $array_in )	= @_;

	# Data-structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Catch and correct genome features in @array_in that have coordinates in 
	# the wrong order:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][3];
		my $var_loc1	= $array_in[$i][4];

		if ( $var_loc1 < $var_loc0 ) {

			$array_in[$i][3]	= $var_loc1;
			$array_in[$i][4]	= $var_loc0;

		}

	}

	# ----------------------------------

	# Add size information to regions in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations of region:
		my $var_loc0		= $array_in[$i][3];
		my $var_loc1		= $array_in[$i][4];

		# Define length of region:
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		$array_in[$i][5]	= $var_length;

	}

	# ----------------------------------

	# Sort regions in @array_in ascending by size:	
	@array_out	= sort { $a -> [5] <=> $b -> [5] } @array_in;

	# ----------------------------------

	# Return @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	RemoveDuplicates
# Description:		This subroutine removes duplicate and enclosed regions.

# ---------------------------------------------------------------------------- #

sub RemoveDuplicates {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Remove duplicate and enclosed regions from @array_out:
	# Remove totally enclosed and duplicate links:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations of region:
		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		# Iterate down and remove any subsequent regions that are identical:
		for ( my $j = $i + 1; $j < scalar @array_in; $j++ ) {

			# Define start and end locations of region:
			my $var_loc2	= $array_in[$j][0];
			my $var_loc3	= $array_in[$j][1];

			if ( $var_loc2 >= $var_loc0 && $var_loc1 >= $var_loc3 ) {

				splice @array_in, $j, 1;

				$j--;

				next;

			}

			elsif ( $var_loc2 > $var_loc1 ) {

				last;

			}

		}

	}

	@array_out	= @array_in;

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	MergeOverlapping
# Description:		This subroutine Merges overlapping regions.

# ---------------------------------------------------------------------------- #

sub MergeOverlapping {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define start and end locations of region:
		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

#		print "$i: $var_loc0 $var_loc1\n";

		# Iterate down through array finding overlapping regions:
		for ( my $j = $i + 1; $j < scalar @array_in; $j++ ) {

			# Define start and end locations of region:
			my $var_loc2	= $array_in[$j][0];
			my $var_loc3	= $array_in[$j][1];

			last if $var_loc2 > $var_loc1 + 1;

			if ( $var_loc2 <= $var_loc1 ) {

				$array_in[$i][1]	= $var_loc3 if $var_loc3 > $var_loc1;

				splice @array_in, $j, 1;

				$i = -1;

				last;

			}

			elsif ( $var_loc1 + 1 == $var_loc2 ) {

#				print "\t$j: $var_loc2 $var_loc3\n";

				$array_in[$i][1]	= $var_loc3 if $var_loc3 > $var_loc1;

				splice @array_in, $j, 1;

				$i = -1;

				last;

			}

		}

	}

	@array_out	= @array_in;

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

1;
