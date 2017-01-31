#!/usr/bin/perl
# list_nonredundant_wm_v2.pl
# given a list of weight matrices, create a non-redundant list using compare2wm to compare them
# input matrices are in the following format
# > information
#   1 2 ... n
# a
# c
# g
# t
#
# 10-23-2003: added system_status
# 07-22-2003: use the environment directories
# 07-22-2003: use program specific temporary file names
# 07-22-2003: remove temporary files when we are done
# 
# new in version 2
# 03-04-2004 also allow filtering of the list of motifs (assumes that the information about each motif is already computed)

# motif information line structure:
# n_processed	map_score	motif_length	mean_pos_all	std_pos_all	min_pos_all	max_pos_all	n_seqs	mean_pos_p	std_pos_p	min_pos_p	max_pos_p	n_seqs_p	mean_pos_m	std_pos_m	min_pos_m	max_pos_m	n_seqs_m	consensus	i_tot	i_mean	maxlik_seq	n	max_score	min_score	mean_score	std_score	min_possible_score	max_possible_score	p5	p10	mc	p10	";
#
# motif_length=3
# n_seqs=8
# i_mean=21
# seq_variety=23 

my $CODE_DIR=$ENV{'CODEDIR'};
my $TEMP_DIR=$ENV{'TEMPDIR'};
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";

$temp_filename="nonredwm";
$system_status=0;
$prog_version=2;
# default thresholds
$information_threshold=0.2331;              # discard if the mean information per nt is <= than this value
$seq_variety_threshold=1;                   # discard if the sequence variety is <= than this value
$n_seqs_threshold=5;                        # minimum number of sequences used to build the weight matrix
$min_motif_length_threshold=6;              # minimum motif length
$max_motif_length_threshold=13;             # maximum motif length
$comp_nmatrix_threshold=0.7;                # default matrix similarity threshold

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "list_nonredundant_wm_v${prog_version}.pl\n";
    print "usage\n";
    print "list_nonredundant_wm_v${prog_version}.pl <filename> <comp_nmatrix_threshold> <information_threshold> <seq_variety_threshold> <n_seqs_threshold> <min_moitf_length_threshold> <max_motif_length_threshold>\n";
    print "\nthe file contains > information lines separating the weigh matrices\n";
    print "<filename>\tfile with the weight matrices and information about each motif\n";
    print "<comp_nmatrix_threshold>\tmatrix similarity threshold, 0 <= comp_nmatrix_threshold <= 100 [default=0.70]\n";
    print "<information_threshold>\tminimum information per nt [default=0.2331]\n";
    print "<seq_variety_threshold>\tminimum sequence variety [default=1]\n";
    print "<n_seqs_threshold>\tminimum number of sequences used to build the weight matrix [default=5]\n";
    print "<min_moitf_length_threshold>\tminimum motif length [default=6]\n";
    print "<max_motif_length_threshold>\tmaximum motif length [default=13]\n";
    print "\n";
    print "output = <filename>.<thres>.txt";
    print "log = <filename>.<thres>.log.txt";
    exit;
}

$curr_arg=0;$input_filename=$ARGV[$curr_arg];
$curr_arg++;
if ($n_args>$curr_arg) {
    $comp_nmatrix_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $information_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $seq_variety_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $n_seqs_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $min_motif_length_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $max_motif_length_threshold=$ARGV[$curr_arg];
}
print "read $n_args arguments\n";

$input_filename_noext=get_filename_noext($input_filename);
$output_filename=">${input_filename_noext}.${comp_nmatrix_threshold}.txt";
$log_filename=">${input_filename_noext}.${comp_nmatrix_threshold}.log.txt";
if ($comp_nmatrix_threshold>1) {
    $comp_nmatrix_threshold=$comp_nmatrix_threshold/100;
}

open (output_file,$output_filename) || die "could not open $output_filename";
open (log_file,$log_filename) || die "could not open $log_filename";

$headertext="";
$dateinfo=get_time_string();
$headertext.="%% list_nonredundant_wm_v${prog_version}.pl\n";
$headertext.="%% $dateinfo\n";
$temptext.="perl list_nonredundant_wm_v${prog_version}.pl ";
for ($i=0;$i<$n_args;$i++) {
    $temptext.="$ARGV[$i] ";
}
$headertext.="%% ${temptext}\n";
$headertext.="%% input_filename=${input_filename}\n";
$headertext.="%% comp_nmatrix_threshold=${comp_nmatrix_threshold}\n";
$headertext.="%% information_threshold=${information_threshold}\n";
$headertext.="%% seq_variety_threshold=${seq_variety_threshold}\n";
$headertext.="%% n_seqs_threshold=${n_seqs_threshold}\n";
$headertext.="%% min_motif_length_threshold=${min_motif_length_threshold}\n";
$headertext.="%% max_motif_length_threshold=${max_motif_length_threshold}\n";
print log_file "$headertext\n";

$running_n_motifs=0;
$running_filtered_out=0;
$running_compared_out=0;
$running_incorporated_motifs=0;

$i=0;
$input_file=open_read_file($input_filename);
while ($i<1000000000) {
#while ($i<2) {
    $i++;
    $entry=read_entry_beginning_with($input_file,">");
    if (!$entry) {
	print "end of file\n";
	last;
    }
	
    @split2lines=split /\n/,$entry;
    $n_lines=$#split2lines+1;
    if ($n_lines ne 5) {
	print "ERROR! i=$i\tn_lines=$n_lines\n";
    }
    $info=$split2lines[0];
    $info=~s/\>//g;

    @split_info=split /\t/,$info;
    $motif_length=$split_info[2];
    $n_seqs=$split_info[7];
    $i_mean=$split_info[20];
    $seq_variety=$split_info[22]; 

    $in=0;
    if ($motif_length<$min_motif_length_threshold) {
	print log_file "\ti=$i\tmotif_length=$motif_length<$min_motif_length_threshold\tskipping\n";
    } else {
	if ($motif_length>$max_motif_length_threshold) {
	    print log_file "\ti=$i\tmotif_length=$motif_length>$max_motif_length_threshold\tskipping\n";
	} else {
	    if ($n_seqs<=$n_seqs_threshold) {
		print log_file "\ti=$i\tn_seqs=$n_seqs<=$n_seqs_threshold\tskipping\n";
	    } else {
		if ($i_mean<=$information_threshold) {
		    print log_file "\ti=$i\ti_mean=$i_mean<=$information_threshold\tskipping\n";
		} else {
		    if ($seq_variety<=$seq_variety_threshold) {
			print log_file "\ti=$i\tseq_variety=$seq_variety<=$seq_variety_threshold\tskipping\n";
		    } else {
			$in=1;
		    }
		}
	    }
	}
    }

    if ($in==0) {
	$running_filtered_out++;
    } else {
	$wm="";
	for ($j=1;$j<=4;$j++) {
	    $wm.="$split2lines[$j]\n";
	}
	@splitarray=split /\t/,$split2lines[1];
	$ml=$#splitarray+1;
	
	$running_n_motifs++;
	$sameordiff=0; # evreyone is different unless the contrary is proven
	$j=0;
	while ( ($sameordiff eq 0) and ($j<$running_incorporated_motifs) ) {
	    $j++;
	    $wm1=$arr[$j][1];
	    $ml1=$arr[$j][0];
	    open (temp_file,">${TEMP_DIR}/${temp_filename}1.txt") || die "could not open output file";
	    print temp_file "$wm1";
	    close (temp_file);
	    open (temp_file,">${TEMP_DIR}/${temp_filename}2.txt") || die "could not open output file";
	    print temp_file "$wm";
	    close (temp_file);
	    $exec_code="${CODE_DIR}/cpp/executables/compare2wm.exe ${TEMP_DIR}/${temp_filename}1.txt ${TEMP_DIR}/${temp_filename}2.txt ${ml1} ${ml}\n";
	    #print "$exec_code\n";
	    @exec_results=`$exec_code`;
	    $temp=$exec_results[0];
	    chomp $temp;
	    @splitarray=split /\=/,$temp;
	    $score=$splitarray[1];
	    #print "score=$score\n";
	    
	    if ($score>$comp_nmatrix_threshold) {
		$sameordiff=1;
	    }
	}
	if ($sameordiff eq 0) {
	    #$n_motifs2_new++;
	    $running_incorporated_motifs++;
	    $arr[$running_incorporated_motifs][0]=$ml;
	    $arr[$running_incorporated_motifs][1]=$wm;
	    #print output_file ">$running_incorporated_motifs\t$info\n";
	    print output_file ">$info\n";
	    print output_file "$wm";
	} else {
	    $running_compared_out++;
	    print "i=$i\tj=$j\tscore=$score\tfound similarities...\n";
	    print log_file "i=$i\tj=$j\tscore=$score\tfound similarities...\n";
	    #exit;
	}
    }

    print "$running_n_motifs\t$running_incorporated_motifs\n";
}
close ($input_filename);
unlink ("${TEMP_DIR}/${temp_filename}1.txt");
unlink ("${TEMP_DIR}/${temp_filename}2.txt");

print "\n\n";
print "running_n_motifs=$running_n_motifs\n";
print "running_incorporated_motifs=$running_incorporated_motifs\n";
print "running_filtered_out=${running_filtered_out}\n";
print "running_compared_out=${running_compared_out}\n";
print "system_status=$system_status\n";

print log_file "running_n_motifs=$running_n_motifs\n";
print log_file "running_incorporated_motifs=$running_incorporated_motifs\n";
print log_file "system_status=$system_status\n";

close (output_file);
close (log_file);

