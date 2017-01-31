#!/usr/bin/perl
# extract_row.pl
# given a file, extract a given row

# 07-25-02
# option of saving file or just displaying

use warnings;
use strict;

my $record = "none yet";
my $i;
my $input_file;
my $output_file;
my $row_number;
my $disporsave;
my $n_args;

$n_args=$#ARGV+1;
if ($n_args<2) {
    print "perl extract_row.pl <filename> <row_number> (<disporsave>)\n";
    print " where <disporsave>=1 implies save file (default) and 0 implies just display\n";
    exit;
}

$input_file=$ARGV[0];
$row_number=$ARGV[1];

if ($n_args eq 3) {
    $disporsave=$ARGV[2];
} else {
    $disporsave=1;
}

if ($disporsave eq 1) {
    $output_file=">row${row_number}_${input_file}";
}
$i=0;

open (INPUTFILE,$input_file) || die "couldn't open the file ${input_file}!";

if ($disporsave eq 1) {
    open (OUTPUTFILE,$output_file) || die "couldn't open the file ${output_file}!";
}
while ($i<${row_number}) {
    $record = <INPUTFILE>;
    $i++;
}

close(INPUTFILE);

print "${record}";

if ($disporsave eq 1) {
    print OUTPUTFILE "$record\n";
close(OUTPUTFILE);
}




