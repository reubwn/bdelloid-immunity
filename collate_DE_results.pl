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
  Can provide a list of genes to '--other_files' that will inherit 1/0 designation in a column with header name taken from '--col_mapping'

USAGE:
  collate_DE_results.pl -t <FASTA> -d file1 [file2...] [-e file3 [file4...]][-p 0.001] [-c 2] [-m col1 [col2...]]

OPTIONS:
  -t|--transcripts   [FILE] : transcriptome fasta file [required]
  -d|--DE_files      [FILE] : DE results file(s) to be parsed and collated
  -r|--method      [STRING] : method used to run DE for format [DESeq2|voom]
  -g|--HGT_files     [FILE] : HGT_locations.txt file
  -e|--other_files   [FILE] : other file(s) to be collated
  -p|--padj         [FLOAT] : FDR threshold for defining DE genes [1e-3]
  -c|--logFC        [FLOAT] : log2 fold-change threshold [2]
  -n|--missing     [STRING] : what to print for genes with no expression data ['NA']
  -m|--col_mapping   [FILE] : tab-delim file with replacements for filename <-> column names
  -o|--out           [FILE] : prefix for outfiles ['collated']
  -h|--help                 : prints this help message
\n";

my ($transcripts_file, $col_map_file, $help);
my (@DE_files, @HGT_files, @other_files);
my $method = "DESeq2";
my $padj_threshold = 0.001;
my $logfc_threshold = 2;
my $missing = "NA";
my $out_prefix = "collated";

GetOptions (
  't|transcripts=s'    => \$transcripts_file,
  'd|DE_files:s{1,}'   => \@DE_files,
  'r|method:s'         => \$method,
  'g|HGT_files:s{,}'   => \@HGT_files,
  'e|other_files:s{,}' => \@other_files,
  'p|padj:f'           => \$padj_threshold,
  'c|logFC:f'          => \$logfc_threshold,
  'n|missing:s'        => \$missing,
  'm|col_mapping:s'    => \$col_map_file,
  'o|out:s'            => \$out_prefix,
  'h|help'             => \$help
);

die $usage if $help;
die $usage unless ($transcripts_file);

my %features_hash;
my $out_file = "${out_prefix}.${method}_P${padj_threshold}_C${logfc_threshold}.DE_results.tab";

## make %col_mapping
my %col_map;
if ( $col_map_file ) {
  open (my $fh, $col_map_file) or die "[ERROR] Can't open mapping file '$col_map_file': $!";
  while (<$fh>) {
    chomp;
    my @F = split (m/\s+/, $_);
    $col_map{$F[0]} = $F[1];
  }
  close $fh;
} else {
  foreach (@DE_files, @HGT_files, @other_files) {
    my $new_name = `basename $_`;
    $col_map{$_} = $new_name;
  }
}
# print Dumper (\%col_map);

## parse feature names from fasta
print STDERR "[INFO] Parsing sequences in file: $transcripts_file\n";
if ( $transcripts_file =~ m/gz$/ ) { ## file is gzipped
  my $in = Bio::SeqIO -> new ( -file => "zcat $transcripts_file |", -format => 'fasta' );
  while ( my $seq_obj = $in->next_seq() ) {
    $features_hash{$seq_obj->display_id()}{length} = $seq_obj->length();
  }
} else {
  my $in = Bio::SeqIO -> new ( -file => $transcripts_file, -format => 'fasta' );
  while ( my $seq_obj = $in->next_seq() ) {
    $features_hash{$seq_obj->display_id()}{length} = $seq_obj->length();
  }
}
print STDERR "[INFO] Number of sequences in $transcripts_file: ".scalar(keys %features_hash)."\n";

## parse features from DE results file(s)
print STDERR "[INFO] Format of DE analysis files is '$method'\n";
print STDERR "[INFO] Number of files to collate: ".scalar(@DE_files)."\n";

foreach my $current_file (@DE_files) {
  print STDERR "[INFO]   $current_file\n";
  open (my $fh, "$current_file") or die $!;
  while (my $line = <$fh>) {
    next if $. == 1; ## header
    chomp $line;
    my @F = split (m/\s+/, $line);

    if ($method =~ m/DESeq2/i) {
      ##
      ## DESeq2 format
      ##

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
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "1";
          ## feature is sig up-regulated
          if ( $F[6] > $logfc_threshold ) {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "1";
          } else {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
          }
          ## feature is sig down-regulated
          if ( $F[6] < -$logfc_threshold ) {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "1";
          } else {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
          }
        } else {
          ## catch features with significant DE but not above FC threshold
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "0";
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
        }
      } else {
        ## catch genes with non-significant DE regardless of FC
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "0";
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
      }

    } elsif ($method =~ m/voom/i) {
      ##
      ## VOOM format
      ##

      ## if current gene has DE result
      $features_hash{$F[0]}{"$col_map{$current_file}.log2FC"} = $F[3];
      $features_hash{$F[0]}{"$col_map{$current_file}.padj"} = $F[6];

      ## can't take log of 0 (Inf), so replace with some very small number
      if ($F[6] == 0) {
        $features_hash{$F[0]}{"$col_map{$current_file}.negLogPadj"} = -log(5e-324)/log(10);
      } else {
        $features_hash{$F[0]}{"$col_map{$current_file}.negLogPadj"} = -log($F[6])/log(10); ## base-N log of a number is equal to the natural log of that number divided by the natural log of N
      }

      ## is feature DE based on thresholds?
      if ( $F[6] < $padj_threshold ) {
        ## feature is significant
        if ( abs($F[3]) > $logfc_threshold ) {
          ## feature is DE
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "1";
          ## feature is sig up-regulated
          if ( $F[3] > $logfc_threshold ) {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "1";
          } else {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
          }
          ## feature is sig down-regulated
          if ( $F[3] < -$logfc_threshold ) {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "1";
          } else {
            $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
          }
        } else {
          ## catch features with significant DE but not above FC threshold
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "0";
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
          $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
        }
      } else {
        ## catch genes with non-significant DE regardless of FC
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE"} = "0";
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_up"} = "0";
        $features_hash{$F[0]}{"$col_map{$current_file}.is_DE_down"} = "0";
      }

    } else {
      die "[ERROR] Method is unspecified! Must be 'DESeq2' or 'voom'\n";
    }
  }
  close $fh;
}

if (scalar(@HGT_files)>0) {
  ## parse features from HGT_locations.txt file(s)
  print STDERR "[INFO] Number of HGT files to collate: ".scalar(@HGT_files)."\n";
  foreach my $current_file (@HGT_files) {
    print STDERR "[INFO]   $current_file\n";

    open (my $fh, "$current_file") or die $!;
    while (my $line = <$fh>) {
      next if $. == 1; ## header
      chomp $line;
      my @F = split (m/\s+/, $line);

      ## based on format of HGT_locations.txt file
      if ($F[7] eq "NA") {
        ## no information (eg no BLAST hit)
        $features_hash{$F[3]}{"is.$col_map{$current_file}"} = "NA";
      } elsif ($F[8] == 2) {
        ## good evidence for HGT
        $features_hash{$F[3]}{"is.$col_map{$current_file}"} = "1";
      } else {
        ## otherwise not-HGT
        $features_hash{$F[3]}{"is.$col_map{$current_file}"} = "0";
      }
    }
    close $fh;
  }
}

if (scalar(@other_files)>0) {
  ## parse features from other results file(s)
  print STDERR "[INFO] Number of other files to collate: ".scalar(@other_files)."\n";
  foreach my $current_file (@other_files) {
    print STDERR "[INFO]   $current_file\n";
    open (my $fh, "$current_file") or die $!;
    while (my $line = <$fh>) {
      chomp $line;
      my @F = split (m/\s+/, $line);
      ## whatever is in the file, we take a 1/0 based on the feature name in col0
      $features_hash{$F[0]}{"is.$col_map{$current_file}"} = "1";
    }
    close $fh;
  }
}

## pull out all the keys
my %keys;
print STDERR "[INFO] Collecting the keys...\n";
foreach (keys %features_hash) {
  foreach (keys %{$features_hash{$_}}) {
    $keys{$_}++;
  }
}
# print STDOUT join ("\t", keys %keys) . "\n";

## print to file
print STDERR "[INFO] Printing results to '$out_file'\n";
open (my $fh, ">$out_file") or die $!;
print $fh join ("\t", "feature", (nsort keys %keys)) . "\n";
foreach my $feature (nsort keys %features_hash) {
  print $fh "$feature";
  my %hash = %{$features_hash{$feature}};
  foreach my $key (nsort keys %keys) {
    if (exists $hash{$key}) {
      print $fh "\t$hash{$key}";
    } else {
      ## pad missing information with $missing (NA by default)
      print $fh "\t$missing";
    }
  }
  print $fh "\n";
}
close $fh;

print STDERR "[INFO] Done! ".`date`;


# print Dumper (\%features_hash);
