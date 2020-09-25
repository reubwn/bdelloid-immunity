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
  collate_DE_results.pl -t <FASTA> -d file1 [file2...] [-e file3 [file4...]][-p 0.001] [-c 2] [-m col1 [col2...]]

OPTIONS:
  -t|--transcripts   [FILE] : transcriptome fasta file [required]
  -d|--DE_files      [FILE] : DE results file(s) to be parsed and collated
  -e|--other_files   [FILE] : other file(s) to be collated
  -p|--padj         [FLOAT] : FDR threshold for defining DE genes [1e-3]
  -c|--logFC        [FLOAT] : log2 fold-change threshold [2]
  -m|--col_mapping   [FILE] : tab-delim file with replacements for filename <-> column names
  -h|--help                 : prints this help message
\n";

my ($transcripts_file,$col_map_file,$help);
my (@DE_files, @other_files);
my $padj_threshold = 0.001;
my $logfc_threshold = 2;

GetOptions (
  't|transcripts=s'    => \$transcripts_file,
  'd|DE_files:s{1,}'   => \@DE_files,
  'e|other_files:s{,}' => \@other_files,
  'p|padj:f'           => \$padj_threshold,
  'c|logFC:f'          => \$logfc_threshold,
  'm|col_mapping:s{,}' => \$col_map_file,
  'h|help'             => \$help
);

die $usage if $help;
die $usage unless ($transcripts_file);

my %features_hash;

## make %col_mapping
my %col_map;
if ( $col_map_file ) {
  open (my $fh, $col_map_file) or die $!;
  while (<$fh>) {
    chomp;
    my @F = split (m/\s+/, $_);
    $col_map{$F[0]} = $F[1];
  }
  close $fh;
} else {
  foreach (@DE_files, @other_files) {
    my $new_name = `basename $_`;
    $col_map{$_} = $new_name;
  }
}
# my @all_files = (@DE_files, @other_files);
# for my $i (0..$#all_files) {
#   if (scalar(@column_mapping) == scalar(@all_files)) {
#     $col_map{$all_files[$i]} = $column_mapping[$i];
#   } else {
#     $col_map{$all_files[$i]} = $all_files[$i];
#   }
# }

print Dumper (\%col_map);

## parse feature names from fasta
print STDERR "[INFO] Parsing sequences in file: $transcripts_file\n";
my $in = Bio::SeqIO -> new ( -file => $transcripts_file, -format => 'fasta' );
while ( my $seq_obj = $in->next_seq() ) {
  $features_hash{$seq_obj->display_id()}{length} = $seq_obj->length();
}
print STDERR "[INFO] Number of sequences in $transcripts_file: ".scalar(keys %features_hash)."\n";

## parse features from DE results file(s)
print STDERR "[INFO] Number of DE analysis files to collate: ".scalar(@DE_files)."\n";
foreach my $current_file (@DE_files) {
  print STDERR "[INFO] $current_file\n";
  open (my $fh, "$current_file") or die $!;
  while (my $line = <$fh>) {
    next if $. == 1; ## header
    chomp $line;
    my @F = split (m/\s+/, $line);

    die "[ERROR] File doesn't look like a DESEq2-style DE results file! (11-cols)\n" if (scalar(@F) != 11);

    # ## if current gene has DE result
    # push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{log2FC}}, $F[6] );
    # push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{padj}}, $F[10] );
    #
    # ## can't take log of 0 (Inf), so replace with some very small number
    # if ($F[10] == 0) {
    #   push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{negLogPadj}}, -log(5e-324)/log(10) );
    # } else {
    #   push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{negLogPadj}}, -log($F[10])/log(10) ); ## base-N log of a number is equal to the natural log of that number divided by the natural log of N
    # }

    ## if current gene has DE result
    $features_hash{$F[0]}{"$col_map{$current_file}.log2FC"} = $F[6];
    $features_hash{$F[0]}{"$col_map{$current_file}.padj"} = $F[10];

    ## can't take log of 0 (Inf), so replace with some very small number
    if ($F[10] == 0) {
      $features_hash{$F[0]}{"$col_map{$current_file}.negLogPadj"} = -log(5e-324)/log(10);
    } else {
      $features_hash{$F[0]}{"$col_map{$current_file}.negLogPadj"} = -log($F[10])/log(10); ## base-N log of a number is equal to the natural log of that number divided by the natural log of N
    }

    ## is feature DE based on thresholds?
    if ( $F[10] < $padj_threshold ) {
      ## feature is significant
      if ( abs($F[6]) > $logfc_threshold ) {
        ## feature is DE
        push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE}}, "1" );
        ## feature is sig up-regulated
        if ( $F[6] > $logfc_threshold ) {
          push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_up}}, "1" );
        } else {
          push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_up}}, "0" );
        }
        ## feature is sig down-regulated
        if ( $F[6] < -$logfc_threshold ) {
          push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_down}}, "1" );
        } else {
          push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_down}}, "0" );
        }
      }
    } else {
      push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE}}, "0" );
      push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_up}}, "0" );
      push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{is_DE_down}}, "0" );
    }

  }
  close $fh;
}

## parse features from other results file(s)
print STDERR "[INFO] Number of other files to collate: ".scalar(@other_files)."\n";
foreach my $current_file (@other_files) {
  print STDERR "[INFO] $current_file\n";
  open (my $fh, "$current_file") or die $!;
  while (my $line = <$fh>) {
    next if $. == 1; ## header
    chomp $line;
    my @F = split (m/\s+/, $line);

    ## whatever is in the file, we take a 1/0 based on the feature name in col0
    push ( @{$features_hash{$F[0]}{$col_map{$current_file}}{"is_$col_map{$current_file}"}}, "1" );

  }
}

print STDOUT join ("\t", "feature", (nsort values %col_map)) . "\n";
foreach my $feature (nsort keys %features_hash) {
  print STDOUT "$feature\t";
  foreach (nsort values %col_map) {

  }
}


print Dumper (\%features_hash);
