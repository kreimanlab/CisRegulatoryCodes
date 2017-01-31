#!/usr/bin/perl
#appendfiles.pl
# append file(s) to a given file
# if the original file does not exist, create it and add the contents of the other files to it
#
# 03-08-2004: add possibility of creating the file
#
# gk, last modified 03-08-2004

use warnings;
use strict;

my $pref;
my $suff;
my $n_args;
my $n_files;
my $main_filename;
my $input_filename;
my $record;
my $i;

$n_args=$#ARGV+1;

if ($n_args < 2) {
    print "appendfiles.pl <main_filename> <input_filename1> <input_filename2> ...\n";
    print "contents of input_filename1/2/... are appended to main_filename\n";
    exit;
}
$main_filename=$ARGV[0];
$i=1;

if (-e $main_filename) {
    print "main_filename=$main_filename (appending)\n";
    open (main_file,">>${main_filename}") || die "could not open $main_filename for appending\n";
} else {
    print "main_filename=$main_filename (creating)\n";
    open (main_file,">${main_filename}") || die "could not open $main_filename for writing\n";
}
while ($i<$n_args) {
    $input_filename=$ARGV[$i];
    print "appending i=$i\t$input_filename...\n";
    open (input_file,$input_filename) || die "could not open $input_filename for reading\n";
    while ($record=<input_file>) {
	print main_file "$record";
    }
    close (input_file);
    $i++;
}
close (main_file);
