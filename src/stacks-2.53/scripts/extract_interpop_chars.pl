#!/usr/bin/env perl
#
# Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
#
# This file is part of Stacks.
#
# Stacks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stacks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use Bio::SeqIO;

use constant true  => 1;
use constant false => 0;

my $debug        = 0;
my $in_file      = "";
my $barcode_list = "";
my $out_file     = "";
my $pop_limit    = 10;

parse_command_line();

my (@barcodes, %barcode_key, %ids, %characters, %polymorph_level, %site_dist);

load_barcode_list(\@barcodes, \%barcode_key);

my ($id, $barcode);

# 
# Initialize the characters hash.
#
foreach $barcode (@barcodes) {
    $characters{$barcode} = "";
    $polymorph_level{$barcode} = {};
    $polymorph_level{$barcode}->{'poly_loci'} = 0;
    $polymorph_level{$barcode}->{'loci'}  = 0;
    $polymorph_level{$barcode}->{'poly'}  = 0;
    $polymorph_level{$barcode}->{'sites'} = 0;
}

#
# Obtain a list of unique tags to process
#
parse_interspecies_tags($in_file, \%ids);

my $i = 1;
my $num_ids = keys %ids;

foreach $id (keys %ids) {
    print STDERR "Processing tag $i of $num_ids.           \r" if ($i % 1000 == 0);

    process_tag($ids{$id}, \@barcodes, \%characters, \%polymorph_level, \%site_dist);

    $i++;
}

print_results(\%characters, \%barcode_key, \%polymorph_level, \%site_dist);

print STDERR "\n";

sub process_tag {
    my ($tag, $barcodes, $characters, $poly_level, $dist) = @_;

    my (%substack, $bc, $len, $col);

    return if (check_tag_counts($tag) == false);

#     if (scalar(keys %{$tag->{'count'}}) >= 10) {
#         print STDERR "Keeping tag ", $tag->{'id'}, " with tags from ", scalar(keys %{$tag->{'count'}}), " populations.\n";
#     }

    #
    # Make sure enough populations have this marker before recording it.
    #
    return if (scalar(keys %{$tag->{'seqs'}}) < $pop_limit);

    foreach $bc (@{$barcodes}) {
	$substack{$bc} = [];
    }

    foreach $bc (keys %{$tag->{'seqs'}}) {
	my $aref;

	#
	# Store the individual reads as a two-dimensional array.
	#
	@{$aref} = split(//, $tag->{'seqs'}->{$bc});
	$substack{$bc} = $aref;

	#
	# Tally the number of non-homozygous sites and the total 
	# number of sites at each locus.
	#
	$len = scalar(@{$aref});
	$poly_level->{$bc}->{'loci'}++;
	$poly_level->{$bc}->{'sites'} += scalar(@{$aref});
	my $poly = 0;
	foreach $col (0..$len - 1) {
	    if ($aref->[$col] eq "N") {
		$poly_level->{$bc}->{'poly'}++;
		$poly++;
	    }
	}
	if ($poly > 0) {
	    $poly_level->{$bc}->{'poly_loci'}++;
	}
    }

    $len = 0;

    foreach $bc (@{$barcodes}) {
	$len = scalar(@{$substack{$bc}});
	last if ($len > 0);
    }

    foreach $col (0..$len - 1) {
	if (homozygous(\%substack, $col) == false) {
	    record_character(\%substack, $col, $characters, $dist);
	}
    }
}

sub check_tag_counts {
    my ($tag) = @_;
    #
    # Check to make sure that there is only a single tag from
    # each individual sample.
    #
    my ($bc);

    foreach $bc (keys %{$tag->{'count'}}) {
	return false if ($tag->{'count'}->{$bc} > 1);
    }

    return true;
}

sub homozygous {
    my ($substack, $col) = @_;

    my ($bc, %nuc);

    foreach $bc (keys %{$substack}) {
        if (scalar(@{$substack->{$bc}}) > 0) {
            $nuc{$substack->{$bc}->[$col]}++;
        }
    }

    my @keys = sort {$nuc{$b} <=> $nuc{$a}} keys(%nuc);

    return true if (scalar(@keys) == 1);

    return true if (scalar(@keys) == 2 && ($keys[0] eq "N" || $keys[1] eq "N"));

    return false;
}

sub record_character {
    my ($stack, $col, $characters, $dist) = @_;

    my ($bc, $cnt);

    $cnt = 0;
    foreach $bc (keys %{$stack}) {
	if (scalar(@{$stack->{$bc}}) == 0) {
	    $characters->{$bc} .= "N";
	} else {
	    $characters->{$bc} .= $stack->{$bc}->[$col];

	    $cnt++ if ($stack->{$bc}->[$col] ne "N");
	}
    }

    $dist->{$cnt}++;
}

sub print_results {
    my ($characters, $key, $poly_level, $dist) = @_;

    my ($out, $seq, $log_fh, $barcode);

    $out = Bio::SeqIO->new(-file => ">$out_file", -format => "fasta");

    foreach $barcode (sort keys %{$characters}) {
	$seq = Bio::Seq->new('-seq' => $characters->{$barcode},
			     '-display_id' => $key->{$barcode},
			     '-alphabet' => 'dna');
	print STDERR "Writing sequence for population '", $key->{$barcode}, " / ", $barcode, "' with a length ", $seq->length(), "\n";
	$out->write_seq($seq);
    }

    print STDERR "Pop\tTotal Loci\tNon-hom Loci\tTotal Sites\tNon-homozygous sites\n";
    foreach $barcode (sort keys %{$poly_level}) {
	print STDERR 
	    $key->{$barcode}, "\t", 
	    $poly_level->{$barcode}->{'loci'}, "\t",
	    $poly_level->{$barcode}->{'poly_loci'}, "\t",
	    $poly_level->{$barcode}->{'sites'}, "\t",
	    $poly_level->{$barcode}->{'poly'}, "\n";
    }

    print STDERR "Number of Populations\tSites\n";
    foreach $barcode (sort {$b <=> $a} keys %{$dist}) {
	print STDERR 
	    $barcode, "\t", 
	    $dist->{$barcode}, "\n";
    }
}

sub parse_interspecies_tags {
    my ($in_path, $ids) = @_;

    my (@parts, $line, $id, $barcode, $catalog_id);

    open(IN, "<$in_path") or 
	die("Unable to open input file '$in_path'; $!\n");

    while ($line = <IN>) {
	chomp $line;
	@parts = split(/\t/, $line);

        next if ($parts[5] eq "consensus");

	$id = $parts[2];
	($barcode, $catalog_id) = ($parts[7] =~ /^(\d+)\_(\d+)/);
        #($barcode, $catalog_id) = ($parts[4] =~ /^(\d+)\_(\d+)/);

	#print STDERR "Barcode_key: $barcode, $barcode_key{$barcode}\n";
	#next if ($barcode_key{$barcode} eq "CB2");

	if (!defined($ids->{$id})) {
	    $ids->{$id} = {};
            $ids->{$id}->{'id'}       = $parts[2];
	    $ids->{$id}->{'seqs'}     = {};
	    $ids->{$id}->{'count'}    = {};
	    $ids->{$id}->{'batch_id'} = $parts[1];
	}

	$ids->{$id}->{'seqs'}->{$barcode} = $parts[8];
	#$ids->{$id}->{'seqs'}->{$barcode} = $parts[5];
	$ids->{$id}->{'count'}->{$barcode}++;
	#print STDERR "Adding '$parts[8]' with barcode $barcode to $id\n";
    }
}

sub load_barcode_list {
    my ($bl, $bk) = @_;

    open(BC, "<" . $barcode_list) 
	or die("Unable to open barcodes file '$barcode_list': $!\n");

    my ($line, $pop_id, $sample_id);

    while ($line = <BC>) {
	chomp $line;

	next if (length($line) == 0 || substr($line, 0, 1) eq "#");

	($pop_id, $sample_id) = ($line =~ /^(\w+)\t(\d+)$/);

        $bk->{$sample_id} = $pop_id;

	push(@{$bl}, $sample_id);
    }

    close(BC);

    if (scalar(@{$bl}) == 0) {
	print STDERR "Unable to load any barcodes from '$barcode_list'\n";
	usage();
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-f$/) { $in_file      = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $barcode_list = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_file     = shift @ARGV; }
	elsif ($_ =~ /^-p$/) { $pop_limit    = shift @ARGV; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }


}

sub usage {
	print << "EOQ";
extract-interpop-chars.pl -f in_file -b barcodes [-o path] [-p limit] [-d] [-h]
  f: input file containing interspecies tags.
  b: list of barcodes to process (one barcode per population).
  o: output file.
  p: minimum number of populations required before recording a marker.
  h: display this help message.
  d: turn on debug output.

EOQ
    exit(0);
}
