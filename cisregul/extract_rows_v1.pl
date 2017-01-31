#!/usr/bin/perl
# extract_rows_v1.pl
# given a file1 extract those lines given by the indices in file2

use warnings;
use strict;

my $record1="none yet";
my $record2="none yet";
my $input_file_name1;
my $input_file_name2;
my $output_file_name;
my $n_args;

$n_args=$#ARGV+1;
if ($n_args<2) {
    print "extract_rows_v1.pl\n";
    print "usage: perl extract_rows_v1.pl datafilename indices\n";
    print "where datafilename contains the data\n";
    print "and   indices contains the indices to extract\n";
    print "good luck!\n";
    print "\ngabriel kreiman, gkreiman\@gnf.org\n";
    exit;
}

$input_file_name1=$ARGV[0];
$input_file_name2=$ARGV[1];
$output_file_name=">${input_file_name1}.extracted";

open (INPUTFILE1, $input_file_name1) || die "couldn't open the file $input_file_name1 !";
open (INPUTFILE2, $input_file_name2) || die "couldn't open the file $input_file_name2 !";
open (OUTFILE,$output_file_name) || die "couldn't open the file for writing!";

my $k=0;
while ($record1 = <INPUTFILE2>) {
    while ($k<$record1) {
	$record2 = <INPUTFILE1>;
        $k++;
    }
    print OUTFILE $record2;
}

close(INPUTFILE1);
close(INPUTFILE2);
close(OUTFILE);











