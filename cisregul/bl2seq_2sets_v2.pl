#!/usr/bin/perl
# bl2seq_2sets_v2.pl
# given two sets of sequences (in fasta format)
# blast one against the other
#
# gabriel kreiman
# all rights reserved (c)
# last modified: 05-20-2003
# kreiman@mit.edu

# 05-20-2003: add lib to library directory
# new in v2: use map to identify which sequences to compare
# 02-06-2003: filter repetitive sequences out as an option
# 06-19-2003: retrieve query and subject sequences, with "n" where there is no match
# 10-08-2003: add verbose parameter to control output to screen (to make it faster)

use warnings;
use strict;

my $CODE_DIR=$ENV{'CODEDIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
my $TEMP_DIR=$ENV{'TEMPDIR'};
my $input_filename1;
my $input_file1;
my $input_filename2;
my $input_file2;
my $max_seq_length;
my $i;
my $n_args;
my $record;
my $entry1;
my $entry2;
my $output_filename;
#my $log_filename;
my @splitarray1;
my @splitarray2;
my $n1;
my $n2;
my $seq1_filename=">${TEMP_DIR}/seq1.txt";
my $seq2_filename=">${TEMP_DIR}/seq2.txt";
my $bl2seq_run_command;
my @bl2seq_output;
my $seq_length;
my $j;
my $seq;
my $seq_query;
my $seq_subject;
my $info1;
my $info2;
my @map;
my $score;
my $expect;
my $identities_nt;
my $identities_perc;
my $query_pos;
my @qpos_array;
my @spos_array;
my $subject_pos;
my @score_array;
my $n_hits;
my @expect_array;
my @identities_nt_array;
my @identities_perc_array;
my @seq_query_array;
my @seq_subject_array;
my $n_seq_compared=0;
my $n_seq_hits=0;
my $seq_acgt;
my $l_acgt;
my $l_orig;
my @fasta1;
my @fasta2;
my $n_seqs1;
my $n_seqs2;
my $ind1;
my $ind2;
my $map_filename;
my $n;
my $verbose=0;                                          # set to 1 for screen verbose output

# defaul blast parameters
my $gap_open_cost=5; # -G
my $gap_extend_cost=2; # -E
my $dropoff_value=50; # -X
my $wordsize=11; # -W
my $mismatch_penalty=-2; # -q
my $match_reward=1; # -r
my $expectation_value=10; # -e
my $query_strand=3; # -S 3=both strands
my $expect_threshold=1e-10;
my $filteron=1; # filter repetitive sequences out

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "\nbl2seq_2sets_v2.pl <input_filename1> <input_filename2> (<map_filename> <max_seq_length> <expect_threshold> <filteron> <gap_open_cost> <gap_extend_cost> <dropoff_value> <wordsize> <mismatch_penalty> <match_reward> <expectation_value> <query_strand>) \n";
    print "\nmax_seq_length:\tblast only the last max_seq_length nucleotides (the most 3' nucleotides)\n";
    print "filteron:\t1 to filter repetitive elements out (default=1)\n";
    print "\n";
    print "output:\tbl2seq_results.txt\n\n";
    exit;
}

$input_filename1=$ARGV[0];
$input_filename2=$ARGV[1];
if ($n_args>3) {
    $map_filename=$ARGV[2];
} else {
    $map_filename="none";
}
if ($n_args >= 4) {
    $max_seq_length=$ARGV[3];
} else {
    $max_seq_length=1000000;
}
if ($n_args>= 5) {
    $expect_threshold=$ARGV[4];
}
if ($n_args>= 6) {
    $filteron=$ARGV[5];
}
if ($n_args>=7 ) {
    $gap_open_cost=$ARGV[6];
}
if ($n_args>= 8) {
    $gap_extend_cost=$ARGV[7];
}
if ($n_args>= 9) {
    $dropoff_value=$ARGV[8];
}
if ($n_args>= 10) {
    $wordsize=$ARGV[9];
}
if ($n_args>= 11) {
    $mismatch_penalty=$ARGV[10];
}
if ($n_args>= 12) {
    $match_reward=$ARGV[11];
}
if ($n_args>= 13) {
    $expectation_value=$ARGV[12];
}
if ($n_args>= 14) {
    $query_strand=$ARGV[13];
}

$output_filename=">bl2seq_results.txt";

if ($map_filename ne "none") {
    # read map file name
    open (map_file,$map_filename) || die "could not open $map_filename for reading";
    $i=0;
    while ($record = <map_file>) {
	chomp $record;
	@splitarray1=split /\t/,$record;
	$i++;
	$map[$i][1]=$splitarray1[0];
	$map[$i][2]=$splitarray1[1];
	#print "$i\tmap1=$map[$i][1]\tmap2=$map[$i][2]\n";
    }
    $n=$i;
    close (map_file);
} else {
    $n=fasta_line_count($input_filename1);
    for ($i=1;$i<=$n;$i++) {
	$map[$i][1]=$i;
	$map[$i][2]=$i;
    }
}
print "number of entries to compare = $n\n";

(@fasta1)=read_whole_fasta($input_filename1);
#for ($i=1;$i<=$#fasta1;$i++) {
#    print "$i\t$fasta1[$i][2]\n";
#}
#exit;
(@fasta2)=read_whole_fasta($input_filename2);
open (output_file,$output_filename) || die "could not open the file $output_filename for writing";
$n_seqs1=$#fasta1;
$n_seqs2=$#fasta2;
print "n_seqs1=$n_seqs1\n";
print "n_seqs2=$n_seqs2\n";

if ($filteron eq 1) {
    $bl2seq_run_command="bl2seq -i ${TEMP_DIR}/seq1.txt -j ${TEMP_DIR}/seq2.txt -p blastn -o temp.txt -d 0 -G ${gap_open_cost} -E ${gap_extend_cost} -X ${dropoff_value} -W ${wordsize} -q ${mismatch_penalty} -r ${match_reward} -F T -e ${expectation_value} -S ${query_strand}";
} else {
    $bl2seq_run_command="bl2seq -i ${TEMP_DIR}/seq1.txt -j ${TEMP_DIR}/seq2.txt -p blastn -o temp.txt -d 0 -G -${gap_open_cost} -E ${gap_extend_cost} -X ${dropoff_value} -W ${wordsize} -q ${mismatch_penalty} -r ${match_reward} -F F -e ${expectation_value} -S ${query_strand}";
}
print output_file ">program=\tbl2seq_2sets_v2.pl\n";
print output_file ">\n>input_filename1=\t${input_filename1} ($n_seqs1 sequences)\n";
print output_file ">input_filename2=\t${input_filename2} ($n_seqs2 sequences)\n";
print output_file ">map_filename=$map_filename\n";
print output_file ">max_seq_length=\t${max_seq_length}\n";
print output_file ">\n>blast parameters:\n";
print output_file ">expect_threshold=\t${expect_threshold}\n";
print output_file ">filteron=\t${filteron}\n";
print output_file ">gap_open_cost=\t${gap_open_cost}\n";
print output_file ">gap_extend_cost=\t${gap_extend_cost}\n";
print output_file ">dropoff_value=\t${dropoff_value}\n";
print output_file ">wordsize=\t${wordsize}\n";
print output_file ">mismatch_penalty=\t${mismatch_penalty}\n";
print output_file ">match_reward=\t${match_reward}\n";
print output_file ">expectation_value=\t${expectation_value}\n";
print output_file ">query_strand=\t${query_strand}\n";
print output_file ">bl2seq_run_command=\t${bl2seq_run_command}\n>\n";

#$i=0;
#while ($i<10000000) {
#while ($i<6) {
for ($i=1;$i<=$n_seqs1;$i++) {
    #$i++;
    if ( ($i % 10) eq 0 ) {
	print "processing sequence $i\n";
    }
    $ind1=$map[$i][1];
    $ind2=$map[$i][2];
    print "$i\t$ind1\t$ind2\n";

    if ( ($ind1 > 0) and ($ind2 > 0) ) {
	$entry1="$fasta1[$ind1][1]\n$fasta1[$ind1][2]";
	$entry2="$fasta2[$ind2][1]\n$fasta2[$ind2][2]";
	($seq,$info1)=trimandformat_seq($entry1,$max_seq_length);
	($seq_acgt,$l_acgt,$l_orig)=remove_n_fromseq($seq);
	if ($l_acgt>200) {
	    open (seq_file1,$seq1_filename);
	    print seq_file1 ">${info1}\n${seq}\n";
	    close (seq_file1);
	    
	    ($seq,$info2)=trimandformat_seq($entry2,$max_seq_length);
	    ($seq_acgt,$l_acgt,$l_orig)=remove_n_fromseq($seq);
	    if ($l_acgt>200) {
		open (seq_file2,$seq2_filename);
		print seq_file2 ">${info2}\n${seq}\n";
		print "${info2}\n";
		close (seq_file2);
		
		if ($filteron eq 1) {
		    @bl2seq_output=`bl2seq -i ${TEMP_DIR}/seq1.txt -j ${TEMP_DIR}/seq2.txt -p blastn -d 0 -G ${gap_open_cost} -E ${gap_extend_cost} -X ${dropoff_value} -W ${wordsize} -q ${mismatch_penalty} -r 1 -F T -e ${expectation_value} -S ${query_strand}`;
		} else {
		    @bl2seq_output=`bl2seq -i ${TEMP_DIR}/seq1.txt -j ${TEMP_DIR}/seq2.txt -p blastn -d 0 -G ${gap_open_cost} -E ${gap_extend_cost} -X ${dropoff_value} -W ${wordsize} -q ${mismatch_penalty} -r 1 -F F -e ${expectation_value} -S ${query_strand}`;
		}
		#print output_file "@bl2seq_output\n";
		    
		$seq_query="";
		$seq_subject="";
		($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos)=extract_blast_results_from_array(@bl2seq_output);
		@score_array=split /\n/,$score;
		$n_hits=$#score_array+1;
		print "n_hits=$n_hits\n";
		
		$n_seq_compared++;
		if ($n_hits>0) {
		    $n_seq_hits++;
		    @expect_array=split /\n/,$expect;
		    @identities_nt_array=split /\n/,$identities_nt;
		    @identities_perc_array=split /\n/,$identities_perc;
		    @seq_query_array=split /\n/,$seq_query;
		    @seq_subject_array=split /\n/,$seq_subject;
		    @qpos_array=split /\n/,$query_pos;
		    @spos_array=split /\n/,$subject_pos;
		    print output_file ">$i\t$info1\t$info2\n";
		    for ($j=0;$j<$n_hits;$j++) {
			#print "j=$j\t";
			$score=$score_array[$j];
			$expect=$expect_array[$j];
			#print "expect=$expect\n";
			if ( ($expect>=0) and ($expect<$expect_threshold) ) {
			    $identities_nt=$identities_nt_array[$j];
			    $identities_perc=$identities_perc_array[$j];
			    $seq_query=$seq_query_array[$j];
			    $seq_subject=$seq_subject_array[$j];
			    $query_pos=$qpos_array[$j];
			    $subject_pos=$spos_array[$j];
			    print output_file "$j\t$expect\t$score\t$identities_nt\t$identities_perc\t$query_pos\t$subject_pos\t$seq_query\t$seq_subject\n";
			    if ($verbose eq 1) {
				print "$j\t$expect\t$score\t$identities_nt\t$identities_perc\t$query_pos\t$subject_pos\t$seq_query\t$seq_subject\n";
			    }
			}
		    }
		}
	    } else {
		print "i=$i\tl_acgt=$l_acgt\tmm\n";
	    }
	} else {
	    print "i=$i\tl_acgt=$l_acgt\ths\n";
	    #print "seq=$seq\n";
	    #print "entry1=$entry1\n";
	    #print "ind1=$ind1\n";
	    #print "$fasta1[$ind1][2]\n";
	    #exit;
	}
    } else {
	print "ind1=$ind1\tind2=$ind2\n";
    }
}
close (output_file);
print "processed $i entries\n";
print "compared $n_seq_compared entries\n";
print "got $n_seq_hits pairs with some hit below expectation=${expect_threshold}\n";
print "\ncheck results in bl2seq_results.txt\n";

sub trimandformat_seq{
    # trim and format the sequence
    # usage: ($seq,$info)=trimandformat($entry,$max_seq_length)

    my $entry=$_[0];
    my $max_seq_length=$_[1];
    my @splitarray;
    my $info;
    my $seq_length;
    my $j;
    my $seq;

    ($info,$seq)=split_infoseq($entry);
    #@splitarray=split /\n/,$entry;
    #$info=$splitarray[0];
    #$entry =~ s/\n//g;    # convert to single line
    #$entry =~ s/$info//;  # remove the information
    #$seq_length=length($entry);
    $seq_length=length($seq);
    if ($seq_length>$max_seq_length) {
	$j=$seq_length-$max_seq_length;
	$seq=substr($entry,$j,$max_seq_length);
    } #else {
      #$seq=$entry;
      #}
    $seq=singleline2multipleline($seq,80);
    $info =~s /\>//;

    return ($seq,$info);
}
