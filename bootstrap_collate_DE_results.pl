#!/usr/bin/env perl

## author: reubwn Sep 2020

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Sort::Naturally;

my $usage = "
SYNOPSIS:
  Bootstrap resample DE file.

USAGE:
  bootstrap_collate_DE_results.pl -i <FILE> -n <NUMBER_OF_BOOTS>

OPTIONS:
  -i|--in    [FILE] : input DE file [required]
  -n|--boots  [INT] : number of bootstraps to run [100]
  -h|--help         : help message
\n";

my ($in_file, $help);
my $number_boots = 100;

GetOptions (
  'i|in=s'     => \$in_file,
  'n|boots:i'  => \$number_boots,
  'h|help'     => \$help
);

die $usage if $help;
die $usage unless ($in_file);
