#!/usr/bin/env perl

## author: reubwn Sep 2020

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Sort::Naturally;

use Bio::Seq;
use Bio::SeqIO;

my $usage = "
SYNOPSIS:
  Collate DE results into single file.

USAGE:
  collate_DE_results.pl -t <FASTA> [-p 0.001] [-c 2] file1 [file2...]

OPTIONS:
  -t|--transcripts  [FILE] : transcriptome fasta file [required]
  -d|--DE_files     [FILE] : DE results file(s) to be parsed and collated
  -e|--other_files  [FILE] : other file(s) to be collated
  -p|--padj        [FLOAT] : FDR threshold for defining DE genes [1e-3]
  -c|--logFC       [FLOAT] : log2 fold-change threshold [2]
  -h|--help                : prints this help message
\n";

my ($transcripts_file,$help);
my (@DE_files, @other_files);
my $padj_threshold = 0.001;
my $logfc_threshold = 2;

GetOptions (
  't|transcripts=s'    => \$transcripts_file,
  'd|DE_files:s{1,}'   => \@DE_files,
  'e|other_files:s{,}' => \@other_files,
  'p|padj:f'           => \$padj_threshold,
  'c|logFC:f'          => \$logfc_threshold,
  'h|help'             => \$help
);

die $usage if $help;
die $usage unless ($transcripts_file);

my %features_hash;

## parse feature names from fasta
my $in = Bio::SeqIO -> new ( -file => $transcripts_file, -format => 'fasta' );
while ( my $seq_obj = $in->next_seq() ) {
  $features_hash{$seq_obj->display_id()}{length} = $seq_obj->length();
}
print STDERR "[INFO] Number of sequences in $transcripts_file: ".scalar(keys %features_hash)."\n";

## parse features from DE results file(s)
print STDERR "[INFO] Number of DE analysis files to collate: ".scalar(@DE_files)."\n";
foreach my $current_DE_file (@DE_files) {
  print STDERR "[INFO] $current_DE_file\n";
  open (my $fh, "$current_DE_file") or die $!;
  while (my $line = <$fh>) {
    next if $. == 1; ## header
    chomp $line;
    my @F = split (m/\s+/, $line);
    # if ( exists($features_hash{$F[0]}) ) {
      ## if current gene has DE result
      push ( @{$features_hash{$F[0]}{log2FC}}, $F[6] );
      push ( @{$features_hash{$F[0]}{padj}}, $F[10] );

      ## can't take log of 0 (Inf), so replace with some very small number
      if ($F[10] == 0) {
        push ( @{$features_hash{$F[0]}{negLogPadj}}, -log(5e-324)/log(10) );
      } else {
        push ( @{$features_hash{$F[0]}{negLogPadj}}, -log($F[10])/log(10) ); ## base-N log of a number is equal to the natural log of that number divided by the natural log of N
      }

      ## is DE expressed based on thresholds

    # } else {
    #   ## if not, populate with NAs
    #   push ( @{$features_hash{$F[0]}{log2FC}}, "NA" );
    #   push ( @{$features_hash{$F[0]}{padj}}, "NA" );
    #   push ( @{$features_hash{$F[0]}{negLogPadj}}, "NA" );
    # }
  }
  close $fh;
}

print Dumper (\%features_hash);
