#!/usr/bin/perl
#
# Copyright 2014, Julian Catchen <jcatchen@uoregon.edu>
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
use constant stacks_version => "_VERSION_";

use constant true  => 1;
use constant false => 0;

my $tags_midpt = 11;
my $tags_len   = 13;
my $snps_midpt = 3;
my $snps_len   = 9;

my $debug     = 0;
my $in_path   = "";
my $out_path  = "";
my $gzipped   = false;

parse_command_line();

my (@files, @catalog_files);

build_file_list(\@files, \@catalog_files);

my ($file, $num_files, $i, $key);

$num_files = scalar(@catalog_files);
$i         = 1;
foreach $file (@catalog_files) {
    printf(STDERR "Converting catalog files, file % 2s of % 2s [%s]\n", $i, $num_files, $file);

    convert_tags_file($in_path, $out_path, $file . ".catalog");

    convert_snps_file($in_path, $out_path, $file . ".catalog");

    #
    # Just copy the unchanged *.alleles.tsv files.
    #
    if ($gzipped == true) {
 	`cp $in_path/${file}.catalog.alleles.tsv.gz $out_path/.`;
    } else {
 	`cp $in_path/${file}.catalog.alleles.tsv $out_path/.`;
    }

    $i++;
}

$num_files = scalar(@files);
$i         = 1;
foreach $file (@files) {
    printf(STDERR "Converting sample files, file % 2s of % 2s [%s]\n", $i, $num_files, $file);

    convert_tags_file($in_path, $out_path, $file);

    convert_snps_file($in_path, $out_path, $file);

    convert_matches_file($in_path, $out_path, $file);

    #
    # Just copy the unchanged *.alleles.tsv files.
    #
    if ($gzipped == true) {
 	`cp $in_path/${file}.alleles.tsv.gz $out_path/.`;
    } else {
 	`cp $in_path/${file}.alleles.tsv $out_path/.`;
    }

    $i++;
}

sub convert_matches_file {
    my ($in_path, $out_path, $file) = @_;

    my ($path, $in_fh, $out_fh, $line, @parts);

    $path = $in_path . "/" . $file . ".matches.tsv";
    if ($gzipped) {
	$path .= ".gz";
	open($in_fh, "gunzip -c $path |") or die("Unable to open matches file '$path', $!\n");
    } else {
	open($in_fh, "<$path") or die("Unable to open matches file '$path', $!\n");
    }

    $path = $out_path . "/" . $file . ".matches.tsv";
    open($out_fh, ">$path") or die("Unable to open tags output file '$path', $!\n");

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

	print $out_fh 
	    join("\t", @parts), "\t",
	    "0.0\n"; # Missing column
    }

    close($in_fh);
    close($out_fh);

    `gzip -f $path` if ($gzipped);
}

sub convert_tags_file {
    my ($in_path, $out_path, $file) = @_;

    my ($path, $in_fh, $out_fh, $line, @parts);

    $path = $in_path . "/" . $file . ".tags.tsv";
    if ($gzipped) {
	$path .= ".gz";
	open($in_fh, "gunzip -c $path |") or die("Unable to open tags file '$path', $!\n");
    } else {
	open($in_fh, "<$path") or die("Unable to open tags file '$path', $!\n");
    }

    $path = $out_path . "/" . $file . ".tags.tsv";
    open($out_fh, ">$path") or die("Unable to open tags output file '$path', $!\n");

    while ($line = <$in_fh>) {
        chomp $line;

	if ($parts[6] eq "consensus") {
	    print $out_fh 
		$line, "\t",
		"0.0\n"; # Missing column
        } else {
	    print $out_fh
		$line, "\t",
		"\n", # Missing column
        }
    }

    close($in_fh);
    close($out_fh);

    `gzip -f $path` if ($gzipped);
}

sub convert_snps_file {
    my ($in_path, $out_path, $file) = @_;

    my ($path, $in_fh, $out_fh, $line, @parts);

    $path = $in_path . "/" . $file . ".snps.tsv";
    if ($gzipped) {
	$path .= ".gz";
	open($in_fh, "gunzip -c $path |") or die("Unable to open tags file '$path', $!\n");
    } else {
	open($in_fh, "<$path") or die("Unable to open tags file '$path', $!\n");
    }

    $path = $out_path . "/" . $file . ".snps.tsv";
    open($out_fh, ">$path") or die("Unable to open tags output file '$path', $!\n");

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

	print $out_fh
	    join("\t", @parts[0 .. $snps_midpt]), "\t",
	    "E\t", # Missing column
	    join("\t", @parts[$snps_midpt+1 .. $snps_len-1]), "\n";
    }

    close($in_fh);
    close($out_fh);

    `gzip -f $path` if ($gzipped);
}

sub build_file_list {
    my ($files, $catalog_files) = @_;

    my (@wl, @ls, $line, $prefix);

    @ls = `ls -1 $in_path/*.tags.tsv* 2> /dev/null`;

    if (scalar(@ls) == 0) {
	print STDERR "Unable to locate any input files to process within '$in_path'\n";
	usage();
    }

    foreach $line (@ls) {
	chomp $line;

	if ($line =~ /\.tags\.tsv\.gz$/) {
	    $gzipped  = true;
	    ($prefix) = ($line =~ /$in_path\/(.+)\.tags\.tsv\.gz/);
	} else {
	    ($prefix) = ($line =~ /$in_path\/(.+)\.tags\.tsv/);
	}

        next if ($prefix =~ /catalog/);

	push(@{$files}, $prefix);
    }

    @ls = `ls -1 $in_path/*.catalog.tags.tsv* 2> /dev/null`;

    if (scalar(@ls) == 0) {
	print STDERR "Unable to locate any catalog input files to process within '$in_path'\n";
	usage();
    }

    foreach $line (@ls) {
	chomp $line;

	if ($line =~ /\.catalog\.tags\.tsv\.gz$/) {
	    $gzipped  = true;
	    ($prefix) = ($line =~ /$in_path\/(.+)\.catalog\.tags\.tsv\.gz/);
	} else {
	    ($prefix) = ($line =~ /$in_path\/(.+)\.catalog\.tags\.tsv/);
	}

	push(@{$catalog_files}, $prefix);
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { $in_path  = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $in_path   = substr($in_path, 0, -1)   if (substr($in_path, -1)   eq "/");
    $out_path  = substr($out_path, 0, -1)  if (substr($out_path, -1)  eq "/");

    if ($in_path eq $out_path) {
	print STDERR "Input and output paths cannot be the same.\n";
	usage();
    }
}

sub version {
    print STDERR "convert_stacks.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
convert_stacks.pl -p path -s path -o path [-t type] [-d] [-h]
    p: path to the Stacks output files.
    o: path to output the converted Stacks files.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
