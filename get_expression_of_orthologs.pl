#!/usr/bin/env perl

## author: reubwn Oct 2020

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Sort::Naturally;

my $usage = "
SYNOPSIS:
  Get DE expression (log2FC) data for orthologs inferred from an Orthogroups.txt file.

OPTIONS:
  -i|--orthogroups   [FILE] : Orthogroups.txt file [required]
  -d|--DE_files      [FILE] : DE results file(s) to be parsed and collated
  -p|--padj         [FLOAT] : FDR threshold for defining DE genes [1e-3]
  -c|--logFC        [FLOAT] : log2 fold-change threshold [2]
  -o|--out           [FILE] : prefix for outfiles ['orthologs']
  -h|--help                 : prints this help message
\n";

my ($orthogroups_file, $help);
my (@DE_files, @other_files);
my $method = "DESeq2";
my $padj_threshold = 0.001;
my $logfc_threshold = 2;
my $out_prefix = "orthologs";

GetOptions (
  'i|orthogroups=s'    => \$orthogroups_file,
  'd|DE_files:s{1,}'   => \@DE_files,
  'r|method:s'         => \$method,
  'e|other_files:s{,}' => \@other_files,
  'p|padj:f'           => \$padj_threshold,
  'c|logFC:f'          => \$logfc_threshold,
  'o|out:s'            => \$out_prefix,
  'h|help'             => \$help
);

die $usage if $help;
die $usage unless ($orthogroups_file);

my %features_hash;

## parse log2FC from DE results file(s)
print STDERR "[INFO] Number of DE analysis files: ".scalar(@DE_files)."\n";
foreach my $current_file (@DE_files) {
  print STDERR "[INFO]   $current_file\n";
  open (my $fh, "$current_file") or die $!;
  while (my $line = <$fh>) {
    next if $. == 1; ## header
    chomp $line;
    my @F = split (m/\s+/, $line);

    die "[ERROR] File doesn't look like a DESEq2-style DE results file! (11-cols)\n" if (scalar(@F) != 11);

    ## get the log2FC
    $features_hash{$F[0]} = $F[6];
    ## if the files have been supplied sensibly, each feature should only have 1 log2FC value

  }
  close $fh;
}

## open $out_file
my $out_file = "${out_prefix}.tab";
open (my $OUT, ">$out_file") or die $!;

## parse Orthogroups.txt file
open (my $fh, $orthogroups_file) or die $!;
while (my $line = <$fh>) {
  chomp $line;
  my @F = split (m/\s+/, $line);
  shift @F; ## remove OG name
  next if scalar(@F) < 2; ## skip single-member OGs
  foreach my $t (@F) {
    next unless ($features_hash{$t}); ## skip if there is no expression data for $t
    foreach my $q (@F) {
      next unless ($features_hash{$q}); ## skip if there is no expression data for $q
      ## print in long format, but not for self
      print $OUT join ("\t", $t, $features_hash{$t}, $q, $features_hash{$q}) unless $t eq $q;
    }
  }
}
close $fh;
close $OUT;
