#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib '/etc/perl/modules';
use Check;
use General;
use SVG;

# ---------------------------------------------------------------------------- #

# File name:		GenBankToORFs.pl
# Date created:		04 January, 2020
# Last modified:	20 May, 2020
# Created by:		Eliot Stanton

# Description:		This is a script for converting ORFs in GenBank data into 
#					SVG image data.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_genbank;
my @array_align;
my @array_out;

# Import variables from the command-line provided by the user:
getopts('a:e:g:n:o:s:t:', \%hash_variables);

# Define variables used by this script:
my $file_alignment	= $hash_variables{a};
my $var_stem		= $hash_variables{e} || 0;
my $file_genbank	= $hash_variables{g};
my $var_output		= $hash_variables{o} || ".";
my $var_start		= $hash_variables{s} || 0;
my $var_stop		= $hash_variables{t} || 0;

my $file_out		= "$var_output\/default.$var_stem.svg";
my $var_y			= 100;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nORFs.pl [OPTIONS]

	-a Alignment file from Synteny (required)
	-e Stem number for file name
	-g GenBank file (required)
	-o Output directory (default: .)
	-s Start location
	-t Stop location\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_genbank is missing print $var_help and stop script:
unless ( $file_genbank && $file_alignment ) {

	print "$var_help\n";

	print "\tGenBank file required!\n" unless $hash_variables{g};

	print "\tGenBank file required!\n" unless $hash_variables{a};

	exit;

}

# ---------------------------------------------------------------------------- #

# Store $var_output in %hash_variables:
$hash_variables { var_output }	= $var_output;

# Add more variables to %hash_variables:
%hash_variables			= %{General::Variables ( \%hash_variables )};

my $var_factor			= $hash_variables { var_factor };
my $var_height			= "150";
my $var_thick			= $hash_variables { stroke_thickness }/2;
my $var_gap				= $hash_variables { var_gap };

my $var_1Mb_height		= $hash_variables { var_1Mb_height };
my $var_100kb_height	= $hash_variables { var_100kb_height };
my $var_10kb_height		= $hash_variables { var_10kb_height };
my $var_1kb_height		= $hash_variables { var_1kb_height };

# ---------------------------------------------------------------------------- #

# Check if output directory is present, if it isn't, create it:
Check::Directory ( $var_output, \%hash_variables );

# Import GenBank file path to @array_files:
my @array_files	= $file_genbank;

# Check if GenBank file exists and contains data:
Check::Files ( \@array_files, \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import GenBank data into @array_genbank:
@array_genbank		= @{ General::FileToArray ( $file_genbank ) };

# Import alignment data into @array_alignment:
if ( $file_alignment ) {

	@array_align		= @{ General::FileToArray ( $file_alignment ) };

}

# ---------------------------------------------------------------------------- #

# Reformat GenBank data to hold only coordinates, sense, and annotation
# information:
@array_genbank	= @{ FormatGB ( \@array_genbank, $var_start, $var_stop ) };

# Format GenBank data to conform to alignment data:
@array_genbank	= @{ Align ( \@array_genbank, \@array_align )};

for ( my $i = 0; $i < scalar @array_genbank; $i++ ) {

#	print "$i: @{$array_genbank[$i]}\n";

}

# Translate alignment and coding data into SVG formatting:
@array_out		= @{ Translate ( \@array_align, \@array_genbank )};

# Print data to file:
General::ArrayToFile ( \@array_out, $file_out );

# ---------------------------------------------------------------------------- #

# Subroutine name:	Translate
# Description:		This subroutine translates data into SVG format and adds
#					tick marks.

# ----------------------------------------------------------

sub Translate {

	# Arguments:
	my ( $array_align, $array_genbank )	= @_;

	# Data structures:
	my @array_align		= @$array_align;
	my @array_genbank	= @$array_genbank;
	my @array_out;

	# Variables:
	my $var_colour		= "200,200,200";
	my $var_stroke		= "0,0,0";
	my $var_width		= $array_align[$#array_align][3]/$var_factor;
	my $var_tall		= 500;

	my $var_1Mb			= 10**6;
	my $var_100kb		= 10**5;
	my $var_10kb		= 10**4;
	my $var_1kb			= 10**3;

	# ----------------------------------

	# Add header information:
	push @array_out, "<svg width=\"$var_width\" height=\"$var_tall\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

	# Add background:
	push @array_out, "<g id=\"background\">";
	my $var_string	= "\t<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
	$var_string		.= "height=\"$var_tall\" style=\"fill:rgb(255,255,255)\"/>";
	push @array_out, $var_string;
	push @array_out, "</g>";

	# ----------------------------------

	push @array_out, "<g id=\"aligned regions\">";

	# Format aligned region data:
	for ( my $i = 0; $i < scalar @array_align; $i++ ) {

		my $var_loc0		= $array_align[$i][2];
		my $var_loc1		= $array_align[$i][3];
		my $var_copy		= $array_align[$i][4];
		my $var_invert	 	= $array_align[$i][5];
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		# Move on if region is not conserved:
		next if $var_copy == 1;

		$var_length		/= $var_factor;
		$var_loc0		/= $var_factor;

		$var_y			+= $var_height if $var_invert;

		my $var_string	= SVG::Rectangle ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

		$var_string		= "\t$var_string";

		$var_y			-= $var_height if $var_invert;

		# Add $var_string to @array_out:
		push @array_out, $var_string;

	}

	push @array_out, "</g>";

	# ----------------------------------

	push @array_out, "<g id=\"non-aligned regions\">";

	$var_colour		= "220,220,220";

	# Format non-aligned region data:
	for ( my $i = 0; $i < scalar @array_align; $i++ ) {

		my $var_loc0		= $array_align[$i][2];
		my $var_loc1		= $array_align[$i][3];
		my $var_copy		= $array_align[$i][4];
		my $var_invert	 	= $array_align[$i][5];
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		next unless $var_copy == 1;

		$var_length		/= $var_factor;
		$var_loc0		/= $var_factor;

		$var_y			+= $var_height if $var_invert;

		my $var_string	= SVG::Rectangle ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

		$var_string		= "\t$var_string";

		$var_y			-= $var_height if $var_invert;

		# Add $var_string to @array_out:
		push @array_out, $var_string;

	}

	push @array_out, "</g>";

	# ----------------------------------

	$var_y -= $var_gap;

	$var_colour		= "0,0,0";

	push @array_out, "<g id=\"horizontal tick marks\">";

	# Add horizontal tick marks:
	for ( my $i = 0; $i < scalar @array_align; $i++ ) {

		my $var_loc0		= $array_align[$i][2];
		my $var_loc1		= $array_align[$i][3];
		my $var_invert	 	= $array_align[$i][5];
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		$var_length		/= $var_factor;
		$var_loc0		/= $var_factor;

		$var_y			+= $var_height if $var_invert;

		my $var_string	= SVG::Rectangle ( $var_loc0, $var_y, $var_length, $var_thick, $var_colour, $var_stroke );

		$var_y			-= $var_height if $var_invert;

		$var_string		= "\t$var_string";

		push @array_out, $var_string;

	}

	push @array_out, "</g>";

	# ----------------------------------

	push @array_out, "<g id=\"vertical tick marks\">";

	# Add vertical tick marks:
	for ( my $i = 0; $i < scalar @array_align; $i++ ) {

		my $var_loc0		= $array_align[$i][0];
		my $var_loc1		= $array_align[$i][1];
		my $var_loc2		= $array_align[$i][2];
		my $var_invert	 	= $array_align[$i][5];
		my $var_diff		= $var_loc2 - $var_loc0;

		if ( $var_1kb >= $var_loc0 && $var_1kb <= $var_loc1 ) {

			my $var_x	= $var_1kb + $var_diff;

			$var_x /= $var_factor;

			my $var_tick_height	= $var_1kb_height;

			$var_tick_height	= $var_10kb_height if $var_1kb % ( $var_10kb ) == 0;
			$var_tick_height	= $var_100kb_height if $var_1kb % ( $var_100kb ) == 0;
			$var_tick_height	= $var_1Mb_height if $var_1kb % ( $var_1Mb ) == 0;

			$var_y				-= ( $var_tick_height + $var_thick );
			$var_y				+= $var_height if $var_invert;

			my $var_string	= SVG::Rectangle ( $var_x, $var_y, $var_thick,
									 $var_tick_height, $var_colour, $var_stroke );

			$var_y				+= ( $var_tick_height + $var_thick );
		$var_y					-= $var_height if $var_invert;

			push @array_out, $var_string;

			$var_1kb += 10**3;

			$i--;

		}

		if ( $var_loc0 > $var_1kb ) {

			$var_1kb += 10**3;

			$i--;

		}

	}

	$var_y += $var_gap;

	$var_string		= "\t$var_string";

	push @array_out, "</g>";

	# ----------------------------------

	push @array_out, "<g id=\"coding regions\">";

	# Format coding region data:
	for ( my $i = 0; $i < scalar @array_genbank; $i++ ) {

		$var_colour			= "0,0,0";
		my $var_loc0		= $array_genbank[$i][0];
		my $var_loc1		= $array_genbank[$i][1];
		my $var_sense		= $array_genbank[$i][2];
		my $var_copy		= $array_genbank[$i][3];
		my $var_invert		= $array_genbank[$i][4];
		my $var_tag			= $array_genbank[$i][5];
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		$var_length		/= $var_factor;
		$var_loc0		/= $var_factor;

		$var_y	+= $var_height/2 if $var_sense eq "-";
		$var_y	+= $var_height if $var_invert;

		# Determine colour for unknown/hypothetical genes:
		$var_colour	= "120,120,120" if $var_tag =~ /hypothetical/;
		$var_colour	= "120,120,120" if $var_tag =~ /unknown/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage protein/;
		$var_colour	= "120,120,120" if $var_tag =~ /ORF/;
		$var_colour	= "120,120,120" if $var_tag =~ /Orf/;
		$var_colour	= "120,120,120" if $var_tag =~ /putative/;
		$var_colour	= "120,120,120" if $var_tag =~ /Putative/;
		$var_colour	= "120,120,120" if $var_tag =~ /Uncharacterized/;
		$var_colour	= "120,120,120" if $var_tag =~ /Uncharacterized/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage DNA binding protein Roi/;
		$var_colour	= "120,120,120" if $var_tag =~ /Protein co-occuring with molybdenum cofactor/;
		$var_colour	= "120,120,120" if $var_tag =~ /host killer protein/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage TraR\/YbiI family protein/;
		$var_colour	= "120,120,120" if $var_tag =~ /ea22/;
		$var_colour	= "120,120,120" if $var_tag =~ /ea10/;
		$var_colour	= "120,120,120" if $var_tag =~ /EaA/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage EaE protein/;
		$var_colour	= "120,120,120" if $var_tag =~ /generated by GeneMarkS/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage DNA binding protein/;
		$var_colour	= "120,120,120" if $var_tag =~ /COG1896: Predicted hydrolases of HD/;
		$var_colour	= "120,120,120" if $var_tag =~ /Phage endonuclease/;

		# Determine colour for integration/excision and recombination:
		$var_colour	= "63,155,196" if $var_tag =~ /integrase/;
		$var_colour	= "63,155,196" if $var_tag =~ /excisionase/;
		$var_colour	= "63,155,196" if $var_tag =~ /invertase/;
		$var_colour	= "63,155,196" if $var_tag =~ /Phage exonuclease \(EC 3.1.11.3\)/;
		$var_colour	= "63,155,196" if $var_tag =~ /protein Gam/;
		$var_colour	= "63,155,196" if $var_tag =~ /protein Bet/;
		$var_colour	= "63,155,196" if $var_tag =~ /Phage recombination protein/;
		$var_colour	= "63,155,196" if $var_tag =~ /Phage exonuclease (EC 3.1.11.3)/;
		$var_colour	= "63,155,196" if $var_tag =~ /Phage Cox/;
		$var_colour	= "63,155,196" if $var_tag =~ /Holliday junction resolvase \/ Crossover/;

		# Determine colour for DNA replication genes:
		$var_colour	= "72,130,132" if $var_tag =~ /DNA replication protein O/;
		$var_colour	= "72,130,132" if $var_tag =~ /DNA helicase \(EC 3.6.4.12\)/;
		$var_colour	= "72,130,132" if $var_tag =~ /DNA replication protein P/;
		$var_colour	= "72,130,132" if $var_tag =~ /Phage replication protein GpA, endonuclease/;
		$var_colour	= "72,130,132" if $var_tag =~ /Phage replication protein GpB/;
		$var_colour	= "72,130,132" if $var_tag =~ /DNA primase, phage associated/;
		$var_colour	= "72,130,132" if $var_tag =~ /Replicative helicase RepA/;

		# Determine colour for transcriptional regulators:
		$var_colour	= "137,119,219" if $var_tag =~ /transcriptional regulator/;
		$var_colour	= "137,119,219" if $var_tag =~ /antitermination protein Q/;
		$var_colour	= "137,119,219" if $var_tag =~ /antitermination protein N/;
		$var_colour	= "137,119,219" if $var_tag =~ /antitermination protein N/;
		$var_colour	= "137,119,219" if $var_tag =~ /Phage repressor protein cI/;
		$var_colour	= "137,119,219" if $var_tag =~ /Phage activator protein cII/;
		$var_colour	= "137,119,219" if $var_tag =~ /Phage regulatory protein cIII/;
		$var_colour	= "137,119,219" if $var_tag =~ /Phage antirepressor protein/;
		$var_colour	= "137,119,219" if $var_tag =~ /Phage immunity repressor protein GpC/;
	$var_colour	= "137,119,219" if $var_tag =~ /Transcriptional regulator/;

		# Determine colour for lysis genes:
		$var_colour	= "128,128,10" if $var_tag =~ /endopeptidase Rz/;
		$var_colour	= "128,128,10" if $var_tag =~ /lysozyme R/;
		$var_colour	= "128,128,10" if $var_tag =~ /holin/;
		$var_colour	= "128,128,10" if $var_tag =~ /Bor/;
#		$var_colour	= "128,128,10" if $var_tag =~ /Lom/;
		$var_colour	= "128,128,10" if $var_tag =~ /Phage host-killing protein Kil/;
		$var_colour	= "128,128,10" if $var_tag =~ /Phage lysis regulatory protein/;

		# Determine colour for structural/morphogenesis genes and DNA packaging:
		$var_colour	= "184,107,29" if $var_tag =~ /Phage tail fiber protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage head completion protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /DNA packaging/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage terminase/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage tail/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage head/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage fibritin \(wac\) protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage capsid/;
		$var_colour	= "184,107,29" if $var_tag =~ /major tail protein/;

		$var_colour	= "184,107,29" if $var_tag =~ /Phage baseplate assembly/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage major capsid protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Tail fiber assembly/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage major tail tube protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage portal vertex protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage P2 GpE family protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Head-tail adaptor/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage small tail protein E/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage tape measure protein Mup42/;
		$var_colour	= "184,107,29" if $var_tag =~ /Phage portal protein/;
		$var_colour	= "184,107,29" if $var_tag =~ /Putative terminase small subunit/;

		# Determine colour for virulence/toxin genes:
		$var_colour	= "255,157,35" if $var_tag =~ /Shiga/;
		$var_colour	= "255,157,35" if $var_tag =~ /Lom/;

		# Determine colour for tRNA genes:
		$var_colour	= "128,195,79" if $var_tag =~ /tRNA/;

		# Determine colour for other proteins:
		$var_colour	= "198,178,178" if $var_tag =~ /Phage DNA N-6-adenine methyltransferase/;
		$var_colour	= "198,178,178" if $var_tag =~ /Phage DNA adenine methylase/;
		$var_colour	= "198,178,178" if $var_tag =~ /Adenine DNA methyltransferase/;
		$var_colour	= "198,178,178" if $var_tag =~ /Serine\/threonine protein kinase/;
		$var_colour	= "198,178,178" if $var_tag =~ /Phage exclusion protein ren/;
		$var_colour	= "198,178,178" if $var_tag =~ /HNH homing endonuclease/;
		$var_colour	= "198,178,178" if $var_tag =~ /Phage superinfection exclusion protein B/;
		$var_colour	= "198,178,178" if $var_tag =~ /DNA-damage-inducible protein I/;
		$var_colour	= "198,178,178" if $var_tag =~ /Phage-associated homing endonuclease/;
		$var_colour	= "198,178,178" if $var_tag =~ /Anti-adapter protein IraM/;
		$var_colour	= "198,178,178" if $var_tag =~ /Phage serine\/threonine protein phosphatase NinI/;
		$var_colour	= "198,178,178" if $var_tag =~ /Exonuclease SbcC/;
		$var_colour	= "198,178,178" if $var_tag =~ /Retron-type RNA-directed DNA polymerase/;
		$var_colour	= "198,178,178" if $var_tag =~ /Single-stranded DNA-binding protein/;

		# Determine colour for IS genes:
		$var_colour	= "39,80,147" if $var_tag =~ /Insertion element/;
		$var_colour	= "39,80,147" if $var_tag =~ /insertion sequence/;
		$var_colour	= "39,80,147" if $var_tag =~ /Transposase/;
		$var_colour	= "39,80,147" if $var_tag =~ /IS, phage, Tn; Transposon-related functions/;

		print "$i: @{$array_genbank[$i]}\n" if $var_colour eq "0,0,0";

		my $var_string	= SVG::Rectangle ( $var_loc0, $var_y, $var_length, $var_height/2, $var_colour, $var_stroke );

#		print "\t$var_string\n";

		$var_y	-= $var_height/2 if $var_sense eq "-";
		$var_y	-= $var_height if $var_invert;

		$var_string		= "\t$var_string";

		# Add $var_string to @array_out:
		push @array_out, $var_string;

		$var_string		= "# $var_tag";

		push @array_out, $var_string;

	}

	push @array_out, "</g>";

	# ----------------------------------

	push @array_out, "</svg>";

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Align
# Description:		This subroutine aligns data with the Mauve alignment file
#					if required.

# ----------------------------------------------------------

sub Align {

	# Arguments:
	my ( $array_in, $array_align, $var_number )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_align	= @$array_align;
	my @array_out;

	# ----------------------------------

	# Iterate through @array_in containing coding region coordinates:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];
		my $var_sense	= $array_in[$i][2];
		my $var_tag		= $array_in[$i][5];

#		print "$i: $var_loc0 $var_loc1 $var_sense $var_tag\n";

		my $var_copy	= 1;
		my $var_invert	= 0;

		for ( my $j = 0; $j < scalar @array_align; $j++ ) {

			# Define original and adjusted locations:
			my $var_loc2	= $array_align[$j][0];
			my $var_loc3	= $array_align[$j][1];
			my $var_loc4	= $array_align[$j][2];
			my $var_loc5	= $array_align[$j][3];
			$var_copy		= $array_align[$j][4];
			$var_invert		= $array_align[$j][5];

#			print "\t$j: @{$array_align[$j]}\n";

			# Calculate difference between the original location and the
			# adjusted location:
			my $var_diff	= $var_loc4 - $var_loc2;

			# If coding region sits within a block, adjust the locations of the
			# coding region using $var_diff:
			if ( $var_loc2 <= $var_loc0 && $var_loc3 >= $var_loc1 ) {

#				print "\t$j: @{$array_align[$j]}\n";

				my @array_temp	= ( $var_loc0 + $var_diff, $var_loc1+$var_diff, $var_sense, $var_copy, $var_invert, $var_tag );

#				my @array_temp	= ( $var_loc0, $var_loc1, $var_sense, $var_copy, $var_invert, $var_tag );

				push @array_out, \@array_temp;

#				print "\t@array_temp\n";

				last;

			}

			# If coding region goes past the end of the block:
			if ( $var_loc0 > $var_loc2 && $var_loc0 < $var_loc3 && $var_loc1 > $var_loc3 ) {

#				print "\t$j: @{$array_align[$j]}\n";

				# Adjust the location of the current region:
				$array_in[$i][1]	= $var_loc3;

				# Create a new coding region:
				my @array_temp	= ( $var_loc3 + 1, $var_loc1, $var_sense, $var_copy, $var_invert, $var_tag );

#				print "\t\t$j: @array_temp\n";

				splice @array_in, $i+1, 0, \@array_temp;

				$i--;

				last;

			}

			# If coding region slips over the left end:
			if ( $var_loc0 < $var_loc2 && $var_loc1 > $var_loc2 && $var_loc1 < $var_loc3 ) {

#				print "\t$j: @{$array_align[$j]}\n";

				my @array_temp	= ( $var_loc2, $var_loc1, $var_sense, $var_copy, $var_invert, $var_tag );

#				print "\t\t$j: @array_temp\n";

				$array_in[$i]	= \@array_temp;

				$i--;

				last;

			}

		}

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatGB
# Description:		This subroutine reformats GenBank data to an array holding
#					only coordinates and orientation information.

# ----------------------------------------------------------

sub FormatGB {

	# Arguments:
	my ( $array_in, $var_start, $var_stop )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

#		print "$i: @{$array_in[$i]}\n";

		# If line contains product information:
		if ( $array_in[$i][0] =~ /^\/product="/ ) {

			my @array_temp	= split /["]/, $array_in[$i][0];

			$array_in[$i][0]	= $array_temp[1];

			my $var_string	= join " ", @{$array_in[$i]};

#			print "$i: $var_string\n";

			$array_out[$#array_out][5]	= $var_string;

#			print "$i: @{$array_out[$#array_out]}\n";

		}

		# If line contains coding region locations grab them and store:
		if ( $array_in[$i][0]  ) {

			next unless $array_in[$i][0] =~ /^CDS/ || $array_in[$i][0] =~ /^tRNA/;

			next unless length $array_in[$i][0] == 3 || length $array_in[$i][0] == 4;

#			print "$i: @{$array_in[$i]}\n";

			my @array_temp	= split /[\.,),(]/, $array_in[$i][1];

			my $var_orientation	= "+";

			# Define sense of region:
			if ( $array_temp[0] eq "complement" ) {

				$var_orientation	= "-";

			}

			# Define locations:
			my $var_loc0	= $array_temp[$#array_temp-2];
			my $var_loc1	= $array_temp[$#array_temp];

			@array_temp		= ( $var_loc0, $var_loc1, $var_orientation );

			# Handle start stop flags:
			if ( $var_start && $var_stop ) {

				next unless $var_loc0 < $var_stop;
				next unless $var_loc1 > $var_start;

				push @array_out, \@array_temp;

			}

			# Add temporary array to @array_out:
			else {

				push @array_out, \@array_temp;

			}

		}

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Shift
# Description:		This subroutine moves locations to the beginning as
#					required.

# ----------------------------------------------------------

sub Shift {

	# Arguments:
	my ( $array_in, $var_start, $var_stop )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# Variables:
	my $var_length	= $var_stop - $var_start;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		$var_loc0 -= $var_start;
		$var_loc1 -= $var_start;

		$array_out[$i][0]	= $var_loc0;
		$array_out[$i][1]	= $var_loc1;
		$array_out[$i][2]	= $array_in[$i][2];

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #
