#!/usr/bin/perl
# generate_background_list_v2.pl
# generate a random background gene list for a given set of transcripts

# new in version 2:
# this is a much simpler version
# here we assume that the background regions are taken from a completely separate list (e.g. from the downstream regions) and therefore we do not need to extract
# a complementary set of genes

# 05-26-2003: allow input sequences to be in fasta format (default now)
# 05-26-2003: use environmental variables for code_dir

use warnings;
use strict;

my $CODE_DIR=$ENV{'CODEDIR'};
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/math_methods.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";

my $cluster_filename;
my $n_cluster=0;
my $n_bckgrnd;
my $log_filename=">log.txt";
my $record;
my @whole_list;
my $whole_list;
my $n_matches;
my $ln;
my @ln_array;
my $background_filename=">info_background.txt";
my $whole_list_filename;
my @ln_bckgrnd;
my @v;
my @splitarray;
my $curr_query;
my @searchstuff;
my $i;
my $n_args;
my $nx=5;
my $seq_filename;
my @seq;
my $curr_seq;
my $background_seq_filename=">seq_background.txt";
my $printseq=0;
my $fastaformat=1;                                        # if fastaformat==1, then the input sequences are in fasta format, the output has same format as the input
my $prog_version=2;                                       # program version

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "usage:\n";
    print "perl generate_background_list_v${prog_version}.pl <n_cluster> <whole_list> (<nx> <seq_filename> <fastaformat>)\n";
    print "\ngenerate a list of random background genes not present in the cluster\n";
    print "the list contains nx as many genes as n_cluster (default nx=5)\n";
    print "if a file with sequences is entered, the sequences are also retrieved and stored into seq_background.txt\n";
    print "e.g. cluster_filename=c:/life/databh/spellman98etal/spellman98etal_clb2cluster.txt\n";
    print "whole_list=c:/life/databh/upstream_sequences/sc/utr5_sc_1000.fasta.upstream_sequences_names.txt\n";
    print "if fastaformat==1, then the input sequences are in fasta format, the output has same format as the input\n";
    print "output is stored in log.txt as well as the background.txt\n";
    gk();
    exit;
}
$n_cluster=$ARGV[0];                     # number of sequences in the cluster
$whole_list_filename=$ARGV[1];           # list of sequences from which the background set is to be selected
if ($n_args>2) {
    $nx=$ARGV[2]; # default = 5
}
if ($n_args>3) {
    $seq_filename=$ARGV[3];
    $printseq=1;
}
if ($n_args>4) {
    $fastaformat=1;
}

@whole_list=read_whole_file_v3($whole_list_filename,1,0);       # list of entries for the background

open (log_file,$log_filename) || die "could not open $log_filename for writing";
$n_bckgrnd=$nx*$n_cluster;
@ln_bckgrnd=randperm($n_bckgrnd+1);

open (background_file,$background_filename) || die "could not open $background_filename for writing";
if ($printseq eq 1) {
    print "reading sequences from $seq_filename...\n";
    if ($fastaformat eq 0) {
	@seq=read_whole_file($seq_filename);
    } else {
	# read sequences in fasta format
	@seq=read_whole_fasta($seq_filename);
    }
    open (background_seq_file,$background_seq_filename) || die "could not open $background_seq_filename for writing";
    print "ready.\n";
}
for ($i=1;$i<=$n_bckgrnd;$i++) {
    $ln=$ln_bckgrnd[$i];
    $record=$whole_list[$ln];
    $record=~s/\>//g;
    print background_file ">$ln\t$record\n";
    if ($printseq eq 1) {
	if ($fastaformat eq 1) {
	    $curr_seq=$seq[${ln}][2];
	    print background_seq_file ">$ln\t$record\n$curr_seq\n";
	} else {
	    $curr_seq=$seq[$ln];
	    print background_seq_file "$curr_seq\n";
	}
    }
}
close (log_file);
close (background_file);
if ($printseq eq 1) {
    close (background_seq_file);
}



