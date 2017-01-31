#!/usr/bin/perl
# process_alignace_output.pl
# comes from: convert_alignace2weightmatrix_v2.pl
# given the output from alignACE, convert to a weight matrix for easy plotting
#
# 03-08-2004: rc is an input to the compute_minfp_threshold subroutine
# 10-23-2003: make temporary files specific to this program and remove them at the end
# 10-23-2003: add system_status output
# 10-21-2003: allow processing of drosophila (dm) files

$program_name="process_alignace_output.pl";

# new in v2
# multiple file output

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
$TEMP_DIR=$ENV{'TEMPDIR'};

require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/math_methods.pl";
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";

$system_status=0;
$fudge=1;
$eval_minfp_threshold="true";
$verbose=1;
$isnew="true";                             # to check whether an entry is true or not
$nt_freq="0.262\t0.238\t0.238\t0.262";
$rc=1;                                     # process both strands
$n_motifs=0;
$n_processed=0;
$n_filtered_out=0;
$n_compared_out=0;
$n_motifs=0;

# thresholds
$mean_information_threshold=0.2331;
$seq_variety_threshold=1;
$n_seqs_threshold=5;
$motif_length_threshold=6;
$in=1;
$comp_nmatrix_threshold=0.75;

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "${program_name}\n";
    print "usage:\n";
    print "${program_name}\t<species> <filename_1> <filename_2> ... <filename_n> 1\n";
    print "${program_name}\t<species> <filename_list> 0\n";
    print "where <filename_list> is a list of files to process\n";
    exit;
}
$species=$ARGV[0];
if ( ($species ne "hs") & ($species ne "mm") & ($species ne "dm") & ($species ne "sc") ) {
    print "what is the name of that species again?\n";
    print "i only know hs or mm\n";
    exit;
}
$isfilelist=$ARGV[$n_args-1];
if ($isfilelist eq 1) {
    $n_files=$n_args-2;
    for ($i=1;$i<=$n_files;$i++) {
	$filelist[$i]=$ARGV[$i+1];
    }
} else {
    $filelist_filename=$ARGV[1];
    @filelist=read_whole_file($filelist_filename);
    $n_files=$#filelist;
}

$log_filename=">process_alignace_output.log.txt";
open (log_file,$log_filename) || die "could not open $log_filename for output";
$output_filename=">process_alignace_output.out.txt";
open (output_file,$output_filename) || die "could not open $output_filename";

for ($curr_file=1;$curr_file<=$n_files;$curr_file++) {

    $input_filename=$filelist[$curr_file];
    $input_filename_noext=get_filename_noext($input_filename);
    $input_file=open_read_file($input_filename);

    if ($verbose) {
	print log_file "processing curr_file=$curr_file\tinput_filename=$input_filename\n";
    }
    $i=0;
    while ($i<1000000000) {
    #while ($i<10) {
	$i++;
	$entry=read_entry_beginning_with($input_file,"Motif");
	if (!$entry) {
	    print log_file "end of file $input_filename\n";
	    last;
	}
	$n_processed++;
	#$n_motifs++;
	#$output_filename_wm=">${input_filename_noext}_m${n_motifs}.wm.txt";
	#$output_filename_pos=">${input_filename_noext}_m${n_motifs}.pos.txt";
	#open (output_file_wm,$output_filename_wm) || die "could not open $output_filename_wm";
	#open (output_file_pos,$output_filename_pos) || die "could not open $output_filename_pos";
	#print "entry=$entry\n";

	if ($verbose) {
	    print log_file "\tprocessing entry=$i\n";
	}
	($pos,$strand,$map_score,$n_seqs,@seq)=extract_info_from_alignace_entry($entry);
	#$firstseq=1;
	print log_file "n_seqs=$n_seqs\n";
	$motif_length=length($seq[0]);
	print log_file "motif_length=$motif_length\n";
	#print "seq[0]=$seq[0]\n";
	#print "seq[1]=$seq[1]\n";
	if ($n_seqs > 0) {
	    (@nmatrix)=compute_frequency_from_aligned_text(@seq);
	    $nmatrix_txt="";
	    for ($j=0;$j<=3;$j++) {
		for ($k=0;$k<$motif_length;$k++) {
		    $nmatrix_txt.="$nmatrix[$j][$k]\t";
		}
		$nmatrix_txt.="\n";
	    }
	    if ($verbose eq 1) {
		print log_file "\nnmatrix\n";
		print log_file "$nmatrix_txt\n";
	    }
	} else {
	    $motif_length=0;
	}

	# compare with previous motifs
	if ( ($n_motifs>0) & ($motif_length>0) ) {
	    ($isnew,$comp_score,$comp_number)=nonredundant($nmatrix_txt,$motif_length,$comp_nmatrix_threshold,$n_motifs,@nmatrix_list);
	} else {
	    if ( ($n_motifs eq 0) & ($motif_length>0) ) {
		$isnew="true";
	    } else {
		$isnew="false";
		$comp_score="nan";
		$com_number="nan";
	    }
	}
	if ($isnew ne "true") {
	    $n_compared_out++;
	    print log_file "found similarities:\tcurr_file=$curr_file\ti=$i\tn_processed=$n_processed\t\tcomp_number=$comp_number\tcomp_score=$comp_score\n";
	} else {
	    @pos_array=split /\n/,$pos;
	    @strand_array=split /\n/,$strand;

	    $n_seqs=$#pos_array+1;
	    ($mean_pos_all,$std_pos_all,$min_pos_all,$max_pos_all)=basic_stats(0,@pos_array);
	    $mean_pos_all=sprintf "%.2f",$mean_pos_all;
	    $std_pos_all=sprintf "%.2f",$std_pos_all;
	    if ($verbose eq 1) {
		print log_file "pos_array[0]=$pos_array[0]\n";
		print log_file "pos_array[1]=$pos_array[1]\n";
		print log_file "pos_array[$n_seqs]=$pos_array[$n_seqs]\n";
		print log_file "mean_pos_all=$mean_pos_all\n";
		print log_file "std_pos_all=$std_pos_all\n";
		print log_file "min_pos_all=$min_pos_all\n";
		print log_file "max_pos_all=$max_pos_all\n";
	    }

	    # split into p and m strands
	    $n_seqs_p=-1;
	    $n_seqs_m=-1;
	    for ($j=0;$j<$n_seqs;$j++) {
		if ($strand_array[$j] eq 1) {
		    $n_seqs_p++;
		    $pos_p[$n_seqs_p]=$pos_array[$j];
		} else {
		    $n_seqs_m++;
		    $pos_m[$n_seqs_m]=$pos_array[$j];
		}
	    }
	    $n_seqs_p++;
	    $n_seqs_m++;
	    ($mean_pos_p,$std_pos_p,$min_pos_p,$max_pos_p)=basic_stats(0,@pos_p);
	    $mean_pos_p=sprintf "%.2f",$mean_pos_p;
	    $std_pos_p=sprintf "%.2f",$std_pos_p;
	    if ($verbose eq 1) {
		print log_file "pos_p[0]=$pos_p[0]\n";
		print log_file "pos_p[1]=$pos_p[1]\n";
		print log_file "pos_p[$n_seqs_p]=$pos_p[$n_seqs_p]\n";
		print log_file "mean_pos_p=$mean_pos_p\n";
		print log_file "std_pos_p=$std_pos_p\n";
		print log_file "min_pos_p=$min_pos_p\n";
		print log_file "max_pos_p=$max_pos_p\n";
	    }
	    ($mean_pos_m,$std_pos_m,$min_pos_m,$max_pos_m)=basic_stats(0,@pos_m);
	    $mean_pos_m=sprintf "%.2f",$mean_pos_m;
	    $std_pos_m=sprintf "%.2f",$std_pos_m;	    

	    (@fudged_matrix)=freqmatrix_addfudge($fudge,$motif_length,@nmatrix);
	    if ($verbose eq 1) {
		print log_file "\nfudged matrix\n";
		for ($j=0;$j<=3;$j++) {
		    for ($k=0;$k<$motif_length;$k++) {
			print log_file "$fudged_matrix[$j][$k]\t";
		    }
		    print log_file "\n";
		}
	    }

	    (@matrix01)=freqmatrix_to_01($motif_length,@fudged_matrix);
	    if ($verbose eq 1) {
		print log_file "\nmatrix01 matrix\n";
		for ($j=0;$j<=3;$j++) {
		    for ($k=0;$k<$motif_length;$k++) {
			$txt=sprintf "%.4f",$matrix01[$j][$k];
			print log_file "$txt\t";
		    }
		    print log_file "\n";
		}
	    }

	    $consensus=wm2consensus($motif_length,@matrix01);

	    (@wm)=matrix01_to_weightmatrix($motif_length,$nt_freq,@matrix01);
	    $wm_text="";
	    for ($j=0;$j<=3;$j++) {
		for ($k=0;$k<$motif_length;$k++) {
		    $txt=sprintf "%.4f",$wm[$j][$k];
		    $wm_text.="$txt\t";
		}
		$wm_text.="\n";
	    }
	    if ($verbose eq 1) {
		print log_file "\nweight matrix\n";
		print log_file "$wm_text\n";
	    }

	    if ($verbose eq 1) {
		print log_file "\n\nsequence scores\n";
	    }
	    for ($k=0;$k<$n_seqs;$k++) {
		$curr_seq=$seq[$k];
		$score=score_seq($curr_seq,$motif_length,@wm);
		$scores[$k]=$score;
		if ($verbose eq 1) {
		    print log_file "$i\t$curr_seq\t$score\n";
		}
	    }
	    ($mean_score,$std_score,$min_score,$max_score)=basic_stats(@scores);

	    if ($verbose eq 1) {
		$mean_score=sprintf "%.4f",$mean_score;
		$std_score=sprintf "%.4f",$std_score;
		$min_score=sprintf "%.4f",$min_score;
		$max_score=sprintf "%.4f",$max_score;
		print log_file "mean_score=$mean_score\n";
		print log_file "std_score=$std_score\n";
		print log_file "min_score=$min_score\n";
		print log_file "max_score=$max_score\n";
	    }

	    ($maxlik_seq,$max_possible_score,$min_possible_score)=maxlik_seq($motif_length,@wm);

	    if ($verbose eq 1) {
		$min_possible_score=sprintf "%.4f",$min_possible_score;
		$max_possible_score=sprintf "%.4f",$max_possible_score;
		print log_file "min_possible_score=$min_possible_score\n";
		print log_file "max_possible_score=$max_possible_score\n";
		print log_file "maxlik_seq=$maxlik_seq\n";
	    }
	    $n=seq_variety($maxlik_seq);

	    if ($verbose eq 1) {
		$score=score_seq($maxlik_seq,$motif_length,@wm);
		print log_file "score of maxlik seq=$score\n";

		print log_file "consensus=$consensus\n";
	    }

	    ($i_tot,$i_mean)=motif_information($nt_freq,$motif_length,@matrix01);
	    if ($verbose eq 1) {
		$i_tot=sprintf "%.4f",$i_tot;
		$i_mean=sprintf "%.4f",$i_mean;
		print log_file "i_tot=$i_tot\n";
		print log_file "i_mean=$i_mean\n";
	    }

	    # filtering
	    $in=1;                             # in by default
	    if ($i_mean<=$mean_information_threshold) {
		$in=0;
		$n_filtered_out++;
	    } else {
		if ($n<=$seq_variety_threshold) {
		    $in=0;
		    $n_filtered_out++;
		} else {
		    if ($n_seqs<=$n_seqs_threshold) {
			$in=0;
			$n_filtered_out++;
		    } else {
			if ($motif_length <= $motif_length_threshold) {
			    $in=0;
			    $n_filtered_out++;
			}
		    }
		}
	    }

	    if ($in eq 1) {
		# compute the 5th percentile and 10th percentiles of scores
		@sorted_scores=sort_matrix(1,$n_seqs,1,0,@scores);
		$i5=ceil(0.05*$n_seqs);
		$i10=ceil(0.10*$n_seqs);
		print log_file "i5=$i5\n";
		print log_file "i10=$i10\n";
		$p5=$sorted_scores[$i5];
		$p5=sprintf "%.4f",$p5;
		$p10=$sorted_scores[$i10];
		$p10=sprintf "%.4f",$p10;
		print log_file "p5=$p5\n";
		print log_file "p10=$p10\n";
		
		if ($eval_minfp_threshold eq "true") {
		    ($max_possible_score2,$min_possible_score2,$mean_score2,$std_score2,$max_score2,$min_score2,@percentiles)=compute_minfp_threshold($motif_length,$species,$n_seqs,$rc,@wm);
		    if ($verbose eq 1) {
			print log_file "min_possible_score2=$min_possible_score2\n";
			print log_file "max_possible_score2=$max_possible_score2\n";
			print log_file "percentiles=@percentiles\n";
		    }
		}
		
		$mc=($p10-$min_possible_score)/($max_possible_score-$min_possible_score);
		$mc=sprintf "%.4f",$mc;
		
		$n_motifs++;
		$nmatrix_list[$n_motifs][0]=$motif_length;
		$nmatrix_list[$n_motifs][1]=$nmatrix_txt;
		print output_file ">$n_processed\t$map_score\t$motif_length\t$mean_pos_all\t$std_pos_all\t$min_pos_all\t$max_pos_all\t$n_seqs\t$mean_pos_p\t$std_pos_p\t$min_pos_p\t$max_pos_p\t$n_seqs_p\t$mean_pos_m\t$std_pos_m\t$min_pos_m\t$max_pos_m\t$n_seqs_m\t$consensus\t$i_tot\t$i_mean\t$maxlik_seq\t$n\t$max_score\t$min_score\t$mean_score\t$std_score\t$min_possible_score\t$max_possible_score\t$p5\t$p10\t$mc\t$p10\t";
		for ($j=0;$j<=4;$j++) {
		    print output_file "$percentiles[$j]\t";
		}
		if ($percentiles[1]>$p10) {
		    print output_file "$percentiles[1]\n";
		} else {
		    print output_file "$p10\n";
		}
		
		for ($j=0;$j<=3;$j++) {
		    for ($k=0;$k<$motif_length;$k++) {
			$txt=sprintf "%.4f",$wm[$j][$k];
			print output_file "$txt\t";
		    }
		    print output_file "\n";
		}
	    }     # close check on in eq 1
	}     # close check on isnew
    }         # close i loop on different entries from the alignace file
    close ($input_file);
    print log_file "n_files=$n_files\n";
    print log_file "n_processed=$n_processed\n";
    print log_file "n_motifs=$n_motifs\n";
    print log_file "n_filtered_out=$n_filtered_out\n";
    print log_file "n_compared_out=$n_compared_out\n";
}             # close curr_file loop on multiple alignace files

print log_file "\nSUMMARY OF RESULTS\n\n";
print log_file "n_files=$n_files\n";
print log_file "n_processed=$n_processed\n";
print log_file "n_motifs=$n_motifs\n";
print log_file "n_filtered_out=$n_filtered_out\n";
print log_file "n_compared_out=$n_compared_out\n";
print log_file "system_status=$system_status\n";
print "n_files=$n_files\n";
print "n_processed=$n_processed\n";
print "n_motifs=$n_motifs\n";
print "n_filtered_out=$n_filtered_out\n";
print "n_compared_out=$n_compared_out\n";
print "system_status=$system_status\n";

close (output_file);
close (log_file);

unlink("${TEMP_DIR}/process_alignace_output.compwm1.txt");
unlink("${TEMP_DIR}/process_alignace_output.compwm2.txt");

#end main

sub extract_info_from_alignace_entry{
    # extract sequence, position, strand from an alignace entry
    # usage:
    # ($pos,$strand,$map_score,$n_seqs,@seq)=extract_info_from_alignace_entry($entry);
    
    my $entry=$_[0];
    my $pos="";
    my $strand="";
    my @seq;
    my $j;
    my @split2lines=split /\n/,$entry;
    my $n_lines=$#split2lines;
    my $curr_line;
    my $map_score="nan";
    my $splitarray;
    my $n_seqs=-1;

    for ($j=0;$j<=$n_lines;$j++) {
	$curr_line=$split2lines[$j];
	if ($curr_line !~ /^[A|C|G|T|N]/) {
	    if ($curr_line =~ /MAP/) {
		@splitarray=split /\:/,$curr_line;
		$map_score=$splitarray[1];
		$map_score=~s/\s+//g;
	    }
	} else {
	    if ($curr_line !~ /^Align/) {
		$n_seqs++;
		@splitarray=split /\t/,$curr_line;
		$seq[$n_seqs]=$splitarray[0];
		$pos.="$splitarray[2]\n";
		$strand.="$splitarray[3]\n";
	    }
	}
    }

    return ($pos,$strand,$map_score,$n_seqs,@seq);
}

sub nonredundant{
    # compare a wm to the list of previous weight matrices and output "true" if it is non-redundant and "false" otherwise
    # usage: 
    # ($isnew,$score,$entry_number)=nonredundant($curr_wm,$curr_ml,$threshold,$n,@wm_list)
    # uses global variables TEMP_DIR and CODE_DIR

    my $curr_wm;                    # current weight matrix
    my $curr_ml;                    # current motif length
    my $threshold;                  # threshold in comparsion
    my @wm_list;
    my $ml;
    my $wm;
    my $i;
    my $n;                           # number of motifs in list
    my $isnew="true";                # everybody is different until the contrary is proven
    my $temp_file1="${TEMP_DIR}/process_alignace_output.compwm1.txt";
    my $temp_file2="${TEMP_DIR}/process_alignace_output.compwm2.txt";
    my @exec_results;
    my @splitarray;
    my $temp;

    ($curr_wm,$curr_ml,$threshold,$n,@wm_list)=@_;

    $i=0;
    while ( ($i<$n) & ($isnew eq "true") ) {
	$i++;
	$ml=$wm_list[$i][0];
	$wm=$wm_list[$i][1];

	open (temp_file,">${temp_file1}") || die "could not open output file ${temp_file1}";
	print temp_file "$wm";
	close (temp_file);
	open (temp_file,">${temp_file2}") || die "could not open output file ${temp_file2}";
	print temp_file "$curr_wm";
	close (temp_file);
	$exec_code="${CODE_DIR}/cpp/executables/compare2wm.exe ${temp_file1} ${temp_file2} ${ml} ${curr_ml}\n";
	#print log_file "wm (ml=$ml)=\n$wm\ncurr_wm (curr_ml=$curr_ml)=\n$curr_wm\nexec_code=$exec_code\n";
	@exec_results=`$exec_code`;
	$temp=$exec_results[0];
	chomp $temp;
	@splitarray=split /\=/,$temp;
	$score=$splitarray[1];
	#print "score=$score\n";
		
	if ($score>$threshold) {
	    $isnew="false";
	}
    }

    return ($isnew,$score,$i);
}
