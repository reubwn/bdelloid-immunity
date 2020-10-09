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
    ## if the files have been supplied sensibly, each feature should only have 1 log2FC value
    $features_hash{$F[0]}{log2FC} = $F[6];
    ## get the padj
    $features_hash{$F[0]}{padj} = $F[10];
    ## calculate the negLogPadj
    if ($F[10] == 0) {
      $features_hash{$F[0]}{negLogPadj} = -log(5e-324)/log(10);
    } else {
      $features_hash{$F[0]}{negLogPadj} = -log($F[10])/log(10); ## base-N log of a number is equal to the natural log of that number divided by the natural log of N
    }
  }
  close $fh;
}

## open $out_file
my $out_file = "${out_prefix}.DE_results.tab";
open (my $OUT, ">$out_file") or die $!;
print $OUT "g1\tlog2FC.g1\tpadj.g1\tnegLogPadj.g1\tis_DE.g1\tg2\tlog2FC.g2\tpadj.g2\tnegLogPadj.g2\tis_DE.g2\n"; ## header
my %seen_already;

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
      ## count the number of times a given pair is seen
      my $key = join ("_", nsort($t,$q));
      $seen_already{$key}++;
      ## and skip if the pair has already been seen (e.g. in reverse)
      next if $seen_already{$key} > 1;
      ## print in long format, but not for self!
      unless ($t eq $q) {
        print $OUT join ("\t", $t, $features_hash{$t}{log2FC}, $features_hash{$t}{padj}, $features_hash{$t}{negLogPadj}, ((abs($features_hash{$t}{log2FC}) > $logfc_threshold) && ($features_hash{$t}{padj} < $padj_threshold) ? 1 : 0));
        print $OUT "\t";
        print $OUT join ("\t", $q, $features_hash{$q}{log2FC}, $features_hash{$q}{padj}, $features_hash{$q}{negLogPadj}, ((abs($features_hash{$q}{log2FC}) > $logfc_threshold) && ($features_hash{$q}{padj} < $padj_threshold) ? 1 : 0));
        print $OUT "\n";
      }
    }
  }
}
close $fh;
close $OUT;
