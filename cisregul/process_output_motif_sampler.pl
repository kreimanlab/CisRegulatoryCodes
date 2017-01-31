#!/usr/bin/perl
# process_output_motif_sampler.pl
# read in the matrix output of MotifSampler.exe
#
# 12-16-2003: make rc an input to this program and also to compute_minfp_threshold

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
$TEMP_DIR=$ENV{'TEMPDIR'};

require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/math_methods.pl";

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "process_output_motif_sampler.pl\n";
    print "usage:\n";
    print "process_output_motif_sampler.pl <matrix_filename> <sites_filename> (<species> <mean_information_threshold> <seq_variety_threshold> <n_seqs_threshold> <comp_nmaitrx_threshold> <rc> <verbose>)\n";
    exit;
}

$matrix_filename=$ARGV[0];
$sites_filename=$ARGV[1];

# default parameters
$mean_information_threshold=0.2331;
$seq_variety_threshold=1;
$n_seqs_threshold=5;
$motif_length_threshold=6;
$comp_nmatrix_threshold=0.75;
$rc=1;
$searchpat_instances='instances\: (\d+)\s+';
$searchpat_seq='site \"(\w+)\"';
$species="mm";
$verbose=1;
$n_motifs=0;
$nt_freq="0.262\t0.238\t0.238\t0.262";
$fudge=1;
$eval_minfp_threshold="true";
$n_filtered_out=0;
$n_compared_out=0;
$system_status=0;

$curr_arg=2;
if ($n_args>$curr_arg) {
    $species=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $mean_information_threshold=$ARGV[$curr_arg];
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
    $comp_nmatrix_threshold=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $rc=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $verbose=$ARGV[$curr_arg];
}

$matrix_file=open_read_file($matrix_filename);
$sites_file=open_read_file($sites_filename);

# open log and output files
open (log_file,">motif_sampler.log.txt") || die "could not open motif_sampler.log.txt for writing";
open (output_file,">motif_sampler.out.txt") || die "could not open motif_sampler.log.txt for writing";

print log_file "process_output_motif_sampler.pl\n";
$t=get_time_string();
print log_file "$t\n";
print log_file "matrix_filename=$matrix_filename\n";
print log_file "sites_filename=$sites_filename\n";
print log_file "mean_information_threshold=$mean_information_threshold\n";
print log_file "seq_variety_threshold=$seq_variety_threshold\n";
print log_file "n_seqs_threshold=$n_seqs_threshold\n";
print log_file "motif_length_threshold=$motif_length_threshold\n";
print log_file "comp_matrix_threshold=$comp_nmatrix_threshold\n";
print log_file "species=$species\n";
print log_file "nt_freq=$nt_freq\n";
print log_file "fudge=$fudge\n";
print log_file "rc=$rc\n";

$i=0;
while ($i<1000000) {
    $i++;
    $matrix_entry=read_entry_beginning_with($matrix_file,"#ID");
    $sites_entry=read_entry_beginning_with($sites_file,"#id");

    if ($verbose) {
	print log_file "processing i=$i\n";
    }

    if (!$matrix_entry) {
	last;
    }
    if (!$sites_entry) {
	last;
    }
    @split_entry=split /\n/,$matrix_entry;
    $n=$#split_entry+1;
    $motif_length=0;

    for ($j=0;$j<$n;$j++) {
	$curr_line=$split_entry[$j];
	if ( ($curr_line =~ /^\#/) | (length($curr_line)<1) ) {
	} else {
	    @split_line=split /\t/,$curr_line;
	    for ($nuc=0;$nuc<=3;$nuc++) {
		$matrix01_orig[$nuc][$motif_length]=$split_line[$nuc];
		#print "j=$j\tnuc=$nuc\tmatrix01_orig[$nuc][$motif_length]=$matrix01_orig[$nuc][$motif_length]\n";
	    }
	    $motif_length++;
	}
    }
    ($consensus,$n_matches)=extract_text_from_array("Consensus",@split_entry);
    ($sampler_score,$n_matches)=extract_text_from_array("Score",@split_entry);
    ($sampler_id,$n_matches)=extract_text_from_array("ID",@split_entry);
    if ($verbose) {
	print log_file "\tconsensus=$consensus\n";
	print log_file "\tsampler_score=$sampler_score\n";
	print log_file "\tsampler_id=$sampler_id\n";
	print log_file "\tmotif_length=$motif_length\n";
    }
    
    # convert matrix to text
    $matrix01_orig_txt="";
    for ($row=0;$row<=3;$row++) {
	for ($column=0;$column<$motif_length;$column++) {
	    $matrix01_orig_txt.="$matrix01_orig[$row][$column]\t";
	}
	$matrix01_orig_txt.="\n";
    }
    if ($verbose eq 1) {
	print log_file "\tmatrix01_orig\n";
	print log_file "$matrix01_orig_txt\n";
    }

    if ($n_motifs>0) {
	# compare with previous motifs
	($isnew,$comp_score,$comp_number)=nonredundant($matrix01_orig_txt,$motif_length,$comp_nmatrix_threshold,$n_motifs,@matrix_list);
    } else {
	$isnew="true";
    }

    if ($isnew ne "true") {
	$n_compared_out++;
	print log_file "\tfound similarities:\tcurrent entry i = $i\tcomp_number=$comp_number\tcomp_score=$comp_score\n";
    } else {
	#@pos_array=split /\n/,$pos;
	#@strand_array=split /\n/,$strand;

	@split_sites_entry=split /\n/,$sites_entry;
	$first_site_line=0;
	$firstline_sites_entry="none";
	while ( ($firstline_sites_entry !~ /id/) & ($first_site_line < $#split_sites_entry) ) {
	    $firstline_sites_entry=$split_sites_entry[$first_site_line];
	    $first_site_line++;
	}
	(@searchstuff)=$firstline_sites_entry=~/$searchpat_instances/;
	if (!@searchstuff) {
	    print log_file "firstline_sites_entry=$firstline_sites_entry\ti could not find the number of sequences\n";
	    exit;
	} else {
	    $n_seqs=$searchstuff[0];
	    if ($verbose) {
		print log_file "\tn_seqs=$n_seqs\n";
	    }

	    # create nmatrix
	    for ($row=0;$row<=3;$row++) {
		for ($column=0;$column<$motif_length;$column++) {
		    $nmatrix[$row][$column]=$matrix01_orig[$row][$column]*$n_seqs;
		}
	    }
	    if ($verbose) {
		$nmatrix_txt="";
		for ($row=0;$row<=3;$row++) {
		    for ($column=0;$column<$motif_length;$column++) {
			$nmatrix_txt.="$nmatrix[$row][$column]\t";
		    }
		    $nmatrix_txt.="\n";
		}
		print log_file "\tnmatrix\n";
		print log_file "$nmatrix_txt\n";
	    }

	    (@fudged_matrix)=freqmatrix_addfudge($fudge,$motif_length,@nmatrix);
	    if ($verbose eq 1) {
		print log_file "\n\tfudged matrix\n";
		for ($j=0;$j<=3;$j++) {
		    for ($k=0;$k<$motif_length;$k++) {
			print log_file "$fudged_matrix[$j][$k]\t";
		    }
		    print log_file "\n";
		}
	    }

	    (@matrix01)=freqmatrix_to_01($motif_length,@fudged_matrix);
	    if ($verbose eq 1) {
		print log_file "\n\tmatrix01 matrix\n";
		for ($j=0;$j<=3;$j++) {
		    for ($k=0;$k<$motif_length;$k++) {
			$txt=sprintf "%.4f",$matrix01[$j][$k];
			print log_file "$txt\t";
		    }
		    print log_file "\n";
		}
	    }

	    $consensus2=wm2consensus($motif_length,@matrix01);
	    if ($verbose) {
		print log_file "\tconsensus2=$consensus2\n";
	    }

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
		print log_file "\n\tweight matrix\n";
		print log_file "$wm_text\n";
	    }

	    if ($verbose eq 1) {
		print log_file "\n\n\ti\tn_motifs\tk\tsequence scores\tpos\tstrand\n";
	    }

	    # get sequences for actual occurrences and score 'em
	    $n_seqs_p=0;
	    $n_seqs_m=0;
	    for ($k=1;$k<=$n_seqs;$k++) {
		$curr_line=$split_sites_entry[$k+$first_site_line-1];
		(@searchstuff)=$curr_line=~/$searchpat_seq/;
		if (!@searchstuff) {
		    print log_file "i=$i\tk=$k\tcould not find sequence in curr_line=$curr_line\n";
		} else {
		    $curr_seq=$searchstuff[0];
		    $score=score_seq($curr_seq,$motif_length,@wm);
		    $scores[$k]=$score;

		    @split_curr_line=split /\s+/,$curr_line;
		    $curr_pos=$split_curr_line[3];
		    $curr_strand=$split_curr_line[6];
		    if ($curr_strand eq "+") {
			$n_seqs_p++;
		    } else {
			$n_seqs_m++;
		    }
		    $pos_array[$k]=$curr_pos;
		    if ($verbose eq 1) {
			print log_file "\t$i\t$n_motifs\t$k\t$curr_seq\t$score\t$curr_pos\t$curr_strand\n";
		    }
		}
	    }
	    ($mean_score,$std_score,$min_score,$max_score)=basic_stats(@scores);
	    ($mean_pos_all,$std_pos_all,$min_pos_all,$max_pos_all)=basic_stats(1,@pos_array);
	    $mean_pos_all=sprintf "%.2f",$mean_pos_all;
	    $std_pos_all=sprintf "%.2f",$std_pos_all;
	    if ($verbose eq 1) {
		#print log_file "pos_array[0]=$pos_array[0]\n";
		print log_file "\tpos_array[1]=$pos_array[1]\n";
		print log_file "\tpos_array[$n_seqs]=$pos_array[$n_seqs]\n";
		print log_file "\tmean_pos_all=$mean_pos_all\n";
		print log_file "\tstd_pos_all=$std_pos_all\n";
		print log_file "\tmin_pos_all=$min_pos_all\n";
		print log_file "\tmax_pos_all=$max_pos_all\n";
	    }

	    if ($verbose eq 1) {
		$mean_score=sprintf "%.4f",$mean_score;
		$std_score=sprintf "%.4f",$std_score;
		$min_score=sprintf "%.4f",$min_score;
		$max_score=sprintf "%.4f",$max_score;
		print log_file "\tmean_score=$mean_score\n";
		print log_file "\tstd_score=$std_score\n";
		print log_file "\tmin_score=$min_score\n";
		print log_file "\tmax_score=$max_score\n";
	    }

	    ($maxlik_seq,$max_possible_score,$min_possible_score)=maxlik_seq($motif_length,@wm);

	    if ($verbose eq 1) {
		$min_possible_score=sprintf "%.4f",$min_possible_score;
		$max_possible_score=sprintf "%.4f",$max_possible_score;
		print log_file "\tmin_possible_score=$min_possible_score\n";
		print log_file "\tmax_possible_score=$max_possible_score\n";
		print log_file "\tmaxlik_seq=$maxlik_seq\n";
	    }
	    $n=seq_variety($maxlik_seq);

	    if ($verbose eq 1) {
		print log_file "\tseq variety=$n\n";
		$score=score_seq($maxlik_seq,$motif_length,@wm);
		print log_file "\tscore of maxlik seq=$score\n";
	    }

	    ($i_tot,$i_mean)=motif_information($nt_freq,$motif_length,@matrix01);
	    if ($verbose eq 1) {
		$i_tot=sprintf "%.4f",$i_tot;
		$i_mean=sprintf "%.4f",$i_mean;
		print log_file "\ti_tot=$i_tot\n";
		print log_file "\ti_mean=$i_mean\n";
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
		$p5=$sorted_scores[$i5];
		$p5=sprintf "%.4f",$p5;
		$p10=$sorted_scores[$i10];
		$p10=sprintf "%.4f",$p10;
		if ($verbose) {
		    print log_file "\ti5=$i5\n";
		    print log_file "\ti10=$i10\n";
		    print log_file "\tp5=$p5\n";
		    print log_file "\tp10=$p10\n";
		}

		if ($eval_minfp_threshold eq "true") {
		    ($max_possible_score2,$min_possible_score2,$mean_score2,$std_score2,$max_score2,$min_score2,@percentiles)=compute_minfp_threshold($motif_length,$species,$n_seqs,$rc,@wm);
		    if ($verbose eq 1) {
			print log_file "\tmin_possible_score2=$min_possible_score2\n";
			print log_file "\tmax_possible_score2=$max_possible_score2\n";
			print log_file "\tpercentiles=@percentiles\n";
		    }
		}
		
		$mc=($p10-$min_possible_score)/($max_possible_score-$min_possible_score);
		$mc=sprintf "%.4f",$mc;
		
		$n_motifs++;
		$matrix_list[$n_motifs][0]=$motif_length;
		$matrix_list[$n_motifs][1]=$matrix01_orig_txt;
		print output_file ">$i\t$sampler_score\t$motif_length\t$mean_pos_all\t$std_pos_all\t$min_pos_all\t$max_pos_all\t$n_seqs\t$sampler_id\t$consensus\tnan\tnan\t$n_seqs_p\tnan\tnan\tnan\tnan\t$n_seqs_m\t$consensus2\t$i_tot\t$i_mean\t$maxlik_seq\t$n\t$max_score\t$min_score\t$mean_score\t$std_score\t$min_possible_score\t$max_possible_score\t$p5\t$p10\t$mc\t$p10\t";
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
	}         # close check on whether the word instance is present
    }             # close check on isnew
}                 # close i loop

print log_file "\nSUMMARY OF RESULTS\n\n";
print log_file "n_processed=$i\n";
print log_file "n_motifs=$n_motifs\n";
print log_file "n_filtered_out=$n_filtered_out\n";
print log_file "n_compared_out=$n_compared_out\n";
print log_file "system_status=$system_status\n";
print "n_processed=$i\n";
print "n_motifs=$n_motifs\n";
print "n_filtered_out=$n_filtered_out\n";
print "n_compared_out=$n_compared_out\n";
print "system_status=$system_status\n";

close (log_file);
close (output_file);

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
    my $temp_file1="${TEMP_DIR}/compwm1.txt";
    my $temp_file2="${TEMP_DIR}/compwm2.txt";
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

