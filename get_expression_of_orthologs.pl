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
