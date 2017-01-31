#!/usr/bin/perl
# cooc_checkresults_v12.pl
# take results from cooc_4m ran on a small random background set and then scan for co-occurrences throughout the genome

# new in version 12
# 02-23-2004 allow to scan as we go
# 02-20-2004 call cooc_single_1234m_v4 (allowing order_constraint=-1)
# 02-18-2004 add ml_threhsold (minimum motif length)
# 01-30-2004 add errors file and warnings file output
# 01-31-2004 call cooc_single_1234m_v3
# 02-02-2004 use max_n_transcripts for 'all' genes
# 02-26-2004 call scan_1m_allseq_v6.exe
# 02-26-2004 added file_label to the names of the temporary files holding the weight matrices
# 02-28-2004 added singlefile,seq_source,read_all (scan_1m_allseq_v6 parameters) as input paramters to this program
# 02-29-2004 make n_transfac_motifs an optional input
# 03-04-2004 call cooc_single_1234m_v5

# new in version 11
# 01-28-2004: call cooc_single_1234_v2 (which incorporates strand information)
# try to incorporate both v9 and v10 (i.e. single species and multiple species)
# make the two species generic (not necessarily hs and mm)
# see full history information below the code
# added sp1 and sp2 as input; by default, sp1="hs" and sp2="mm"
# remove the stuff with the filename_full

########################
# list of subroutines  #
########################

# ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_cl,$logtext)=scan_4motifs_cl($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count);
# ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_all,$logtext)=scan_4motifs_all($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count);
# ($motifscan_filename,$motifstrand_filename)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$motif_index,$motif_scan_threshold_column,$species,$n_transfac_motifs);
# $temp_n=scan_motif($wm,$upstream_filename,$scan_output_filename,$motif_index,$col_index,$ml,$thres,$n_seqs,$get_count,$rc);
# $n_transcripts_occurrences=scan_all($wm,$scan_output_filename,$motif_index,$species,$ml,$thres,$get_count,$singlefile,$seq_source,$read_all);
# ($errorcode,$n_cooc,$n_cooc_tr)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename,$max_ul,$max_el,$max_il,$transcript_map,$sortpos);
# ($p_binom,$p)=computebinom($n_background,$n_cooc_background,$n_cluster,$n_cooc_cluster);
# ($p_binom)=computebinom_dist($n_ortho,$n_close,$p_ortho_distance);
# ($max_n_transcripts)=get_max_n_transcripts($k,$n,$y,$n_all);
# ($n_hits,$n_miss,$n_notfound)=map_orthologous($species1,$ind,$clorall,@mm2hs_map);
# ($n_hits,$n_miss,$n_notfound)=count_genes_withhomologous($species1,$ind,$clorall,@map);
# ($comments,$core)=add_genename($input_filename,@info);
# ($p)=compute_random_overlap($n1,$n2,$n,$nr,$n_iter_random_overlap);

# required libraries
use POSIX qw(ceil floor);
$TEMP_DIR=$ENV{'TEMPDIR'};
$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/math_methods.pl";

# default parameters
$last_modified="02_28_2004";                # date last modified
$get_count=1;                               # 1 to get a line count on the scanned files
$threshold_binom_sp1="0.01";                # report if $p_binom_sp1<$threshold_binom_sp1
$threshold_binom_sp2="0.01";                # report if $p_binom_sp2<$threshold_binom_sp2
$threshold_binom_both="0.01";               # report also if ($p_binom_sp1<$threshold_binom_both) & ($p_binom_sp2<$threshold_binom_both) )
$max_dist=100;                              # default maximum distance (bp)
$min_dist_factor=0.5;                       # default minimum distance factor between motifs (bp)
$verbose=0;                                 # not-verbose by default
$prog_version=12;                           # program version                 
$prog_name="cooc_checkresults_v${prog_version}";
$max_n_transcripts=-1;                      # default maximum number of transcripts (i.e. compute for all transcripts by default)
$cooc_single_threshold_p=0.999;             # only run cooc_single until the binomial probability is cooc_single_threshold_p
$rc=1;                                      # rc=1 to consider the reverse complementary strand when scanning
$info_col_cl_sp1=11;                        #  "              "                       info_cl_sp1 file
$info_col_cl_sp2=6;                         #  "              "                       info_cl_sp2 file
$n_all_withhomologous_map=0;                # number of genes with homologous in the mm2hs_all map (or the corresponding between species gene map)
$n_cl_withhomologous_map=0;                 # number of genes with homologous in the mm2hs_cl map (or the corresponding between species gene map)
$p_ortho_distance=0.1;                      # cut-off for the probability that two genes are separated a distance d apart
$remove_temporary=1;                        # if this is 1, remove temporary files at the end
$order_constraint=1;                        # 1 to impose an order constraint for the motifs 
#$motif_filename="${TEMP_DIR}/cooc_checkresults_wm.txt";     # name for temporary weight matrix file
$motif_scan_threshold_column=38;            # threshold column in the scan output file
$init_t=time();                             # initial time
$upstream_length=7000;                      # now this shows up only in the call to compare_cooc_lists.c
$n_iter_random_overlap=1000;                # number of iterations in call to random_overlap.exe
$sp1="hs";                                  # species 1
$sp2="mm";                                  # species 2
$cooc_c_code_version=5;                     # in call to cooc_single_1234m
$max_upstream_length=5000;                  # maximum upstream length
$max_exon_length=5000;                      # maximum exon length
$max_intron_length=5000;                    # maximum intron length
$exp_id="u";                                # experiment id
$init_i=0;                                  # initial module
$locuslink_only=1;                          # 1 to process entries that contain a valid locus link only
$sortpos=1;                                 # 1 to sort scan output positions
$minian_transcripts=4;                      # minimum number of transcripts where the module must be present
$file_label="1";                            # file label (for temporary output)
$seq_upstream_length=5000;                  # maximum length for the available sequences (as used in the scan program), upstream
$seq_exon_length=5000;                      # ", exon
$seq_intron_length=5000;                    # ", intron
$report_all=1;                              # 1 to report all positions
$n_species=2;                               # default = 2 species, comparison of hs and mm
$n_errors=0;                                # number of errors
$n_warnings=0;                              # number of warnings
$ml_threshold=6;                            # minimum motif length
$singlefile=0;                              # 1 if there is a single file with the sequences for all genes (used in scan_1m_allseq_v6, input parameter)
$seq_source="gbk";                          # source for the sequences for all genes (used in scan_1m_allseq_v6, input parameter)
$read_all=0;                                # 1 to read also the empty lines (used in scan_1m_allseq_v6, input parameter)
$n_transfac_motifs=-1;                      # make this an input, if there is no input and it is < 0, we use the default values

$n_args=$#ARGV+1;

if ($n_args<9) {
    usage();
}

$curr_arg=0;
$input_filename=$ARGV[$curr_arg];                         # file with coocurrence results
$curr_arg++;$wm_filename=$ARGV[$curr_arg];                # file with weight matrices
$curr_arg++;$n_species=$ARGV[$curr_arg];                  # number of species (1 or 2)
$curr_arg++;$order_constraint=$ARGV[$curr_arg];           # 1 if order matters and 0 otherwise
$curr_arg++;$upstream_filename_cl_sp1=$ARGV[$curr_arg];    # file with cluster upstream regions, hs
if ($n_species eq 2) {
    $curr_arg++;$upstream_filename_cl_sp2=$ARGV[$curr_arg];    # file with cluster upstream regions, mm
}
$curr_arg++;$info_filename_cl_sp1=$ARGV[$curr_arg];        # file with info, cluster, hs
if ($n_species eq 2) {
    $curr_arg++;$info_filename_cl_sp2=$ARGV[$curr_arg];
}
$curr_arg++;$transcript_lengths_filename_sp1_cl=$ARGV[$curr_arg]; # file with upstream, exon and intron lengths
if ($n_species eq 2) {
    $curr_arg++;$transcript_lengths_filename_sp2_cl=$ARGV[$curr_arg]; # file with upstream, exon and intron lengths
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $max_upstream_length=$ARGV[$curr_arg];                # maximum upstream length [default = 5000]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $max_exon_length=$ARGV[$curr_arg];                    # maximum exon length [default = 5000]
}
$curr_arg++; 
if ($n_args>$curr_arg) {
    $max_intron_length=$ARGV[$curr_arg];                  # maximum intron length [default = 5000]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $max_dist=$ARGV[$curr_arg];                           # maximum distance between motifs [default = 100]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $min_dist_factor=$ARGV[$curr_arg];                    # minimum distance = min_dist_factor * motif_length [default min_dist_factor = 0.5]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $exp_id=$ARGV[$curr_arg];                             # experiment id (this is used only to give a prefix to the scan files) [default = 'u'
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $init_i=$ARGV[$curr_arg];                             # in case we wish to start at a later position [default = 0]
} 
if ($n_species eq 2) {
    $curr_arg++;
    if ($n_args>$curr_arg) {
	$mm2hs_map_filename_all=$ARGV[$curr_arg];             # map hs to mm, all genes  [default = "default_map"]
    } else {
	$mm2hs_map_filename_all="default_map";
    }
    $curr_arg++;
    if ($n_args>$curr_arg) {
	$mm2hs_map_filename_cl=$ARGV[$curr_arg];              # map hs to mm, cluster [default = "default_map"]
    } else {
	$mm2hs_map_filename_cl="default_map";
    }
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $rc=$ARGV[$curr_arg];                                 # 1 to consider both strands [default = 1]
} 
$curr_arg++;
if ($n_args > $curr_arg) {
    $motif_scan_threshold_column=$ARGV[$curr_arg];                   # threshold column [default = 38]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $info_col_cl_sp1=$ARGV[$curr_arg];                     #  column with locuslink information, set=cluster, species=hs, (from the info_cl_hs file) [default=11]
}
if ($n_species eq 2) {
    $curr_arg++;
    if ($n_args>$curr_arg) {
	$info_col_cl_sp2=$ARGV[$curr_arg];                     #  column with locuslink information, set=cluster, species=hs, (from the info_cl_hs file) [default=6]
    }
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $locuslink_only=$ARGV[$curr_arg];                     # only consider entries with a valid locus link id (in the all sequences scan) if this is 1 [default = 1]
} 
$curr_arg++;                      
if ($n_args>$curr_arg) {
    $sortpos=$ARGV[$curr_arg];                            # sort positions in output of scan if this is 1 [default = 1]
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $n_iter_random_overlap=$ARGV[$curr_arg];              # number of iterations in call to random overlap
} 
$curr_arg++;
if ($n_args>$curr_arg) {
    $minian_transcripts=$ARGV[$curr_arg];                 # minimum number of transcripts
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $file_label=$ARGV[$curr_arg];                         # file label for output files to be able to run multiple copies of the program in parallel
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $sp1=$ARGV[$curr_arg];
}
if ($n_species eq 2) {
    $curr_arg++;
    if ($n_args>$curr_arg) {
	$sp2=$ARGV[$curr_arg];
    }
} else {
    $sp2="none";
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $singlefile=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $seq_source=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $read_all=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $n_transfac_motifs=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $verbose=$ARGV[$curr_arg];                            # verbose output
} 

# before starting all the mess, check that input_filename exists
if (-e $input_filename) {
    # all right
} else {
    print "\n\nERROR!!!\nI could not find $input_filename\nI am going to exit now\n";
    exit;
}

# more parameters and file names
if ( ($sp1 eq "hs") | ($sp1 eq "mm") ) {
    if ($n_transfac_motifs<0) {
	$n_transfac_motifs=144;            # there are n_transfac_motifs transfac motifs, therefore we search in a different location within the scan all subroutine
    }
    $info_filename_all_sp1="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/exon_list_v3.txt";     # by default, assuming sp1="hs";
    $info_col_all_sp1=6;                        # 
    $transcript_lengths_filename_sp1_all="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/lengths.locusid.${sp1}.txt";
    $transcript_map_filename_sp1="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/exon_list_v4.locuslink.${sp1}.txt";
    if ($sp2 ne "none") {
	$info_filename_all_sp2="${DATA_DIR}/db/upstream_sequences/${sp2}/lists/exon_list_v3.txt";
	$info_col_all_sp2=6;
	$transcript_lengths_filename_sp2_all="${DATA_DIR}/db/upstream_sequences/${sp2}/lists/lengths.locusid.${sp2}.txt";
	$transcript_map_filename_sp2="${DATA_DIR}/db/upstream_sequences/${sp2}/lists/exon_list_v4.locuslink.${sp2}.txt";
    }
} else {
    if ($n_species eq 2) {
	print "this seems very strange\nplease double-check and re-run\nexiting now...\n";
	exit;
    }
    if ($sp1 eq "dm") {
	if ($n_transfac_motifs<0) {
	    $n_transfac_motifs=30;
	}
	$info_filename_all_sp1="${DATA_DIR}/db/ensembl/${sp1}/lists/exon_list_v3.txt";     
	$transcript_lengths_filename_sp1_all="${DATA_DIR}/db/ensembl/${sp1}/lists/lengths.locusid.${sp1}.txt";
	$transcript_map_filename_sp1="${DATA_DIR}/db/ensembl/${sp1}/lists/exon_list_v4.locuslink.${sp1}.txt";
	$info_col_all_sp1=5;  
    }
    if ($sp1 eq "sc") {
	if ($n_transfac_motifs<0) {
	    $n_transfac_motifs=29;
	}
	$info_filename_all_sp1="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/exon_list_v3.txt";     
	$transcript_lengths_filename_sp1_all="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/lengths.locusid.${sp1}.txt";
	$transcript_map_filename_sp1="${DATA_DIR}/db/upstream_sequences/${sp1}/lists/exon_list_v4.locuslink.${sp1}.txt";
	$info_col_all_sp1=5;  
    }
}

# create the directory for the scan output if it does not exist
$scan_dir="${TEMP_DIR}/scan_${exp_id}";
if (-e $scan_dir) {
    print "directory $scan_dir already exists\n";
} else {
    $code_exec="mkdir ${scan_dir}";
    $system_status=system($code_exec);
    if ($system_status eq 0) {
	print "$code_exec\tok\n";
    } else {
	print "$code_exec\terror\n";
	exit;
    }
}

#####################
# open output files #
#####################
$log_filename="cooc_checkresults_v${prog_version}.${file_label}.log.txt";
if (-e $log_filename) {
    rename($log_filename,"${log_filename}.old");
    print "renamed $log_filename to ${log_filename}.old\n";
}
open (log_file,">${log_filename}") || die "could not open $log_filename for writing";
$err_filename="cooc_checkresults_v${prog_version}.${file_label}.err.txt";
if (-e $err_filename) {
    rename($err_filename,"${err_filename}.old");
    print "renamed $err_filename to ${err_filename}.old\n";
}
open (err_file,">${err_filename}") || die "could not open $err_filename for writing";
$warnings_filename="cooc_checkresults_v${prog_version}.${file_label}.wrn.txt";
if (-e $warnings_filename) {
    rename($warnings_filename,"${warnings_filename}.old");
    print "renamed $warnings_filename to ${warnings_filename}.old\n";
}
open (wrn_file,">${warnings_filename}") || die "could not open $warnings_filename for writing";
$output_filename=">cooc_checkresults_v${prog_version}.${file_label}.out.txt";
if (-e $output_filename) {
    rename($output_filename,"${output_filename}.old");
    print "renamed $otuput_filename to ${output_filename}.old\n";
}
open (output_file,">${output_filename}") || die "could not open $output_filename for writing";
$rejected_filename="cooc_checkresults_v${prog_version}.${file_label}.rej.txt";
open (rejected_file,">${rejected_filename}") || die "could not open $rejected_filename for writing";
print rejected_file "% record\tn_cooc_all_sp1\tn_cooc_tr_all_sp1\tn_cooc_cl_sp1\tn_cooc_tr_cl_sp1\tn_cooc_all_sp2\tn_cooc_tr_all_sp2\tn_cooc_cl_sp2\tn_cooc_tr_cl_sp2\tp_binom_sp1\tp_binom_sp2\terrorcode_all_sp1\terrorcode_cl_sp1\terrorcode_all_sp2\terrorcode_cl_sp2\n";
# files with positions of co-occurring motifs
$pos_cl_sp1_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}.pos-cl.txt";
open (pos_cl_sp1_file,">${pos_cl_sp1_filename}") || die "could not open $pos_cl_sp1_filename for writing";
print pos_cl_sp1_file "% i1\ti2\ti3\ti4\tp_binom_${sp1}\tn_cooc_tr_cl_${sp1}\tn_cooc_tr_all_${sp1}\terrorcode_cl_${sp1}\terrorcode_all_${sp1}\t1\n";
$pos_all_sp1_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}.pos-all.txt";
open (pos_all_sp1_file,">${pos_all_sp1_filename}") || die "could not open $pos_all_sp1_filename for writing";
print pos_all_sp1_file "% i1\ti2\ti3\ti4\tp_binom_${sp1}\tn_cooc_tr_cl_${sp1}\tn_cooc_tr_all_${sp1}\terrorcode_cl_${sp1}\terrorcode_all_${sp1}\t1\n";

if ($sp2 ne "none") {
    $pos_cl_sp2_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp2}.pos-cl.txt";
    open (pos_cl_sp2_file,">${pos_cl_sp2_filename}") || die "could not open $pos_cl_sp2_filename for writing";
    print pos_cl_sp2_file "% i1\ti2\ti3\ti4\tp_binom_${sp2}\tn_cooc_tr_cl_${sp2}\tn_cooc_tr_all_${sp2}\terrorcode_cl_${sp2}\terrorcode_all_${sp2}\t1\n";
    $pos_all_sp2_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp2}.pos-all.txt";
    open (pos_all_sp2_file,">${pos_all_sp2_filename}") || die "could not open $pos_all_sp2_filename for writing";
    print pos_all_sp2_file "% i1\ti2\ti3\ti4\tp_binom_${sp2}\tn_cooc_tr_cl_${sp2}\tn_cooc_tr_all_${sp2}\terrorcode_cl_${sp2}\terrorcode_all_${sp2}\t1\n";

    if ($n_iter_random_overlap > 0) {
	# files with transcript and position overlap between thw two species
	$both_troverlap_all_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}${sp2}.trov-all.txt";
	open (both_troverlap_all_file,">${both_troverlap_all_filename}") || die "could not open $both_troverlap_all_filename for writing";
	print both_troverlap_all_file "% i1\ti2\ti3\ti4\tn_all_ov\tn_tr_all_ov\tn_all_pos\tn_tr_all_pos\t1\t1\n";
	$both_posoverlap_all_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}${sp2}.posov-all.txt";
	open (both_posoverlap_all_file,">${both_posoverlap_all_filename}") || die "could not open $both_posoverlap_all_filename for writing";
	print both_posoverlap_all_file "% i1\ti2\ti3\ti4\tn_all_ov\tn_tr_all_ov\tn_all_pos\tn_tr_all_pos\t1\t1\n";
	$both_troverlap_cl_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}${sp2}.trov-cl.txt";
	open (both_troverlap_cl_file,">${both_troverlap_cl_filename}") || die "could not open $both_troverlap_cl_filename for writing";
	print both_troverlap_cl_file "% i1\ti2\ti3\ti4\tn_cl_ov\tn_tr_cl_ov\tn_cl_pos\tn_tr_cl_pos\t1\t1\n";
	$both_posoverlap_cl_filename="cooc_checkresults_v${prog_version}.${file_label}.${sp1}${sp2}.posov-cl.txt";
	open (both_posoverlap_cl_file,">${both_posoverlap_cl_filename}") || die "could not open $both_posoverlap_cl_filename for writing";
	print both_posoverlap_cl_file "% i1\ti2\ti3\ti4\tn_cl_ov\tn_tr_cl_ov\tn_cl_pos\tn_tr_cl_pos\t1\t1\n";
    }
}


##############
# thresholds #
##############
($n_modules,$temp_n)=line_count($input_filename);
$n_modules=$n_modules-$temp_n;
if ($n_modules <= 0 ) {
    #$n_modules=1;
    print "WARNING! cooc_checkresults_v12.pl\ninput_filename=$input_filename\nn_modules=${n_modules}\nexiting...\n";
    close (log_file);close (err_file);close (wrn_file);
    close (output_file);close (pos_cl_sp1_file);close (pos_all_sp1_file);close (rejected_file);
    if ($sp2 ne "none") {
	close (pos_all_sp2_file);close (pos_cl_sp2_file);close (both_troverlap_all_file);
	close (both_posoverlap_all_file);close (both_troverlap_cl_file);close (both_posoverlap_cl_file);
    }
} else {
    $threshold_binom_sp1=$threshold_binom_sp1/$n_modules;        # bonferroni-like correction
    $threshold_binom_sp1_c=1.0-$threshold_binom_sp1;
    if ($n_species eq 2) {
	$threshold_binom_sp2=$threshold_binom_sp2/$n_modules;        # bonferroni-like correction
	$threshold_binom_both=sqrt($threshold_binom_sp1);
    }

######################
# header information #
######################
    $dateinfo=get_time_string();
    print log_file "%% perl cooc_checkresults_v${prog_version}.pl ";
    for ($i=0;$i<$n_args;$i++) {
	print log_file "$ARGV[$i] ";
    }
    print log_file "\n";
    $headertext="";
    $headertext.="% cooc_checkresults_v${prog_version}.pl\n";
    $headertext.="% last modified = $last_modified\n";
    $headertext.="% $dateinfo\n";
    $temptext="perl cooc_checkresults_v${prog_version}.pl ";
    for ($i=0;$i<$n_args;$i++) {
	$temptext.="$ARGV[$i] ";
    }
    $headertext.="% ${temptext}\n";
    $headertext.="% input_filename=$input_filename\n";
    $headertext.="% wm_filename=$wm_filename\n";
    $headertext.="% n_species=$n_species\n";
    $headertext.="% order_constraint=$order_constraint\n";
    $headertext.="% upstream_filename_cl_sp1=$upstream_filename_cl_sp1\n";
    $headertext.="% max_upstream_length=$max_upstream_length\n";
    $headertext.="% max_exon_length=$max_exon_length\n";
    $headertext.="% max_intron_length=$max_intron_length\n";
    $headertext.="% max_dist=$max_dist\n";
    $headertext.="% min_dist_factor=$min_dist_factor\n";
    $headertext.="% exp_id=$exp_id\n";
    $headertext.="% init_i=$init_i\n";
    $headertext.="% verbose=${verbose}\n";
    $headertext.="% sp1 = ${sp1}\n";
    $headertext.="% sp2 = ${sp2}\n";
    $headertext.="% n_modules=${n_modules}\n";
    $headertext.="% threshold_binom_sp1 = $threshold_binom_sp1\n";
    $headertext.="% info_filename_all_sp1 = $info_filename_all_sp1\n";
    $headertext.="% info_col_all_sp1 = $info_col_all_sp1\n";
    $headertext.="% info_filename_cl_sp1 = $info_filename_cl_sp1\n";
    $headertext.="% info_col_cl_sp1 = $info_col_cl_sp1\n";
    $headertext.="% transcript_lengths_filename_sp1_cl = $transcript_lengths_filename_sp1_cl\n";
    $headertext.="% transcript_map_filename_sp1 = $transcript_map_filename_sp1\n";
    $headertext.="% transcript_lengths_filename_sp1_all = $transcript_lengths_filename_sp1_all\n";
    if ($n_species eq 2) {
	$headertext.="% upstream_filename_cl_sp2=$upstream_filename_cl_sp2\n";
	$headertext.="% threshold_binom_sp2 = $threshold_binom_sp2\n";
	$headertext.="% threshold_binom_both = $threshold_binom_both\n";
	$headertext.="% info_filename_all_sp2 = $info_filename_all_sp2\n";
	$headertext.="% info_col_all_sp2 = $info_col_all_sp2\n";
	$headertext.="% info_filename_cl_sp1 = $info_filename_cl_sp1\n";
	$headertext.="% info_filename_cl_sp2 = $info_filename_cl_sp2\n";
	$headertext.="% info_col_cl_sp2 = $info_col_cl_sp2\n";
	$headertext.="% transcript_lengths_filename_sp2_cl = $transcript_lengths_filename_sp2_cl\n";
	$headertext.="% transcript_map_filename_sp2 = $transcript_map_filename_sp2\n";
	$headertext.="% transcript_lengths_filename_sp2_all = $transcript_lengths_filename_sp2_all\n";
	$headertext.="% mm2hs_map_all=${mm2hs_map_filename_all}\n";
	$headertext.="% mm2hs_map_cl=${mm2hs_map_filename_cl}\n";
    }
    $headertext.="% locuslink_only = $locuslink_only\n";
    $headertext.="% sortpos = $sortpos\n";
    $headertext.="% rc = $rc\n";
    $headertext.="% motif_scan_threshold_column=$motif_scan_threshold_column\n";
    $headertext.="% minian_transcripts = $minian_transcripts\n";
    $headertext.="% seq_upstream_length = $seq_upstream_length\n";
    $headertext.="% seq_exon_length = $seq_exon_length\n";
    $headertext.="% seq_intron_length = $seq_intron_length\n";
    $headertext.="% ml_threshold = $ml_threshold\n";
    $headertext.="% singlefile = $singlefile\n";
    $headertext.="% seq_source = $seq_source\n";
    $headertext.="% read_all = $read_all\n";
    $headertext.="% scan_dir = $scan_dir\n";
    $headertext.="% report_all = $report_all\n";
    $headertext.="% n_transfac_motifs = $n_transfac_motifs\n";
    $headertext.="%";
    print log_file "$headertext\n";
    print output_file "$headertext\n";
    if ($sp2 ne "none") {
	$headertext="% i\tj\tk\tl\tml_1\tml_2\tml_3\tml_4\tn_cl\tn_tr_cl\tn_bck\tntr_bck\tr_oc\tr_tr\tp_b_oc\tp_b_tr\tec\tn_cooc_all_${sp1}\tn_cooc_tr_all_${sp1}\tn_all_withhomologous_${sp1}\tn_cooc_cl_${sp1}\tn_cooc_tr_cl_${sp1}\tn_cl_withhomologous_${sp1}\tn_cooc_all_${sp2}\tn_cooc_tr_all_${sp2}\tn_all_withhomologous_${sp2}\tn_cooc_cl_${sp2}\tn_cooc_tr_cl_${sp2}\tn_cl_withhomologous_${sp2}\tp_binom_${sp1}\tp_binom_${sp2}\tn_cl_ov\tn_tr_cl_ov\tp_cl_trov\tn_cl_pos\tn_tr_cl_pos\tp_cl_posov\tn_all_ov\tn_tr_all_ov\tp_all_trov\tn_all_pos\tn_tr_all_pos\tp_all_posov\terrorcode_all_${sp1}\terrorcode_cl_${sp1}\terrorcode_all_${sp2}\terrorcode_cl_${sp2}\n";
    } else {
	$headertext="% i\tj\tk\tl\tml_1\tml_2\tml_3\tml_4\tn_cl\tn_tr_cl\tn_bck\tntr_bck\tr_oc\tr_tr\tp_b_oc\tp_b_tr\tec\tn_cooc_all_${sp1}\tn_cooc_tr_all_${sp1}\tn_all_withhomologous_${sp1}\tn_cooc_cl_${sp1}\tn_cooc_tr_cl_${sp1}\tn_cl_withhomologous_${sp1}\tnone\tnone\tnone\tnone\tnone\tnone\tp_binom_${sp1}\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\terrorcode_all_${sp1}\terrorcode_cl_${sp1}\tnone\tnone\n";
    }
    print output_file "$headertext";
    print log_file "$headertext";

#######################################
# get the total number of transcripts #
#######################################
    $n_upstream_all_sp1=get_n_transcripts_all($info_filename_all_sp1,$locuslink_only);
    $n_upstream_cl_sp1=line_count($upstream_filename_cl_sp1);
    print log_file "% n_upstream_all_sp1=$n_upstream_all_sp1\n";
    print log_file "% n_upstream_cl_sp1=$n_upstream_cl_sp1\n";
    if ($sp2 ne "none") {
	$n_upstream_all_sp2=get_n_transcripts_all($info_filename_all_sp2,$locuslink_only);
	$n_upstream_all=max2($n_upstream_all_hs,$n_upstream_all_mm);
	$n_upstream_cl_sp2=line_count($upstream_filename_cl_sp2);
	$n_upstream_cl=max2($n_upstream_cl_sp1,$n_upstream_cl_sp2);
	print log_file "% n_upstream_all_sp2=$n_upstream_all_sp2\n";
	print log_file "% n_upstream_cl_sp2=$n_upstream_cl_sp2\n";
    } else {
	$n_upstream_all=$n_upstream_all_sp1;
	$n_upstream_cl=$n_upstream_cl_sp1;
    }
    print log_file "% n_upstream_all=$n_upstream_all\n";
    print log_file "% n_upstream_cl=$n_upstream_cl\n";

#################################################
# load file with all the gene names information #
#################################################
    @info_all_sp1=read_column_v2($info_filename_all_sp1,$info_col_all_sp1,1,0);    # ignore all header lines if there are any, read $info_col column
    @info_cl_sp1=read_column_v2($info_filename_cl_sp1,$info_col_cl_sp1,1,0);       # ignore all header lines if there are any, read $info_col column
    if ($sp2 ne "none") {
	@info_all_sp2=read_column_v2($info_filename_all_sp2,$info_col_all_sp2,1,0);    # ignore all header lines if there are any, read $info_col column
	@info_cl_sp2=read_column_v2($info_filename_cl_sp2,$info_col_cl_sp2,1,0);       # ignore all header lines if there are any, read $info_col column
    }

######################################
# load file with all weight matrices #
######################################
# e.g. wmatrix_blackshawetal_rpe_500.mm.alignace_10.comp.txt
# >3      133.402 21      280.52  134.43  2.00    475.00  84      277.00  143.82  2.00    471.00  43 284.22   123.70  8.00    475.00  41      gvnnnsnnnnsngvgsngvns
#0.00    0.15    0.27    0.18    0.30    0.00    0.26    0.21    0.18    0.26    0.08    0.26    0.05        0.40    0.00    0.06    0.30    0.02    0.11    0.20    0.01
#0.00    0.36    0.25    0.26    0.20    0.12    0.20    0.27    0.21    0.21    0.27    0.20    0.00        0.18    0.00    0.21    0.32    0.00    0.14    0.19    0.17
#1.00    0.42    0.35    0.42    0.37    0.88    0.43    0.33    0.45    0.36    0.64    0.37    0.95        0.42    0.93    0.73    0.27    0.95    0.75    0.48    0.80
#0.00    0.07    0.13    0.14    0.13    0.00    0.11    0.18    0.15    0.17    0.00    0.17    0.00        0.00    0.07    0.00    0.11    0.02    0.00    0.13    0.02
    open (input_file,$wm_filename) || die "could not open $wm_filename";
    $temp="";
    $i=0;
    while ($record=<input_file>) {
	if ($record =~ /^>/) {
	    $wm[$i]=$temp;                                                    # previous weight matrix
	    $temp="";                                                         # re-initialize temp (which holds the weight matrices
	    $i++;                                                             # current entry number
	    chomp $record;                         
	    @splitarray=split /\t/,$record;       
	    $thres[$i]=$splitarray[$motif_scan_threshold_column]-0.0001;      # threshold, avoid rounding problems by subtracting 0.0001
	    $ml[$i]=$splitarray[2];                                           # motif length
	} else {
	    $temp.="$record";
	}
    }
    $wm[$i]=$temp;
    if ($verbose eq 1) {
	print log_file "read $i weight matrices\n";
	print log_file "wm[1]=$wm[1]\tthres[1]=$thres[1]\tml[1]=$ml[1]\n";
	print log_file "wm[$i]=$wm[$i]\tthres[$i]=$thres[$i]\tml[$i]=$ml[$i]\n";
    }
    close (input_file);    
    $n_motifs=$i;
    print log_file "% n_motifs=$n_motifs\n";
    print log_file "motif:\t";
    for ($i=1;$i<=$n_motifs;$i++) {
	print log_file "$i\t";
    }
    print log_file "\n";
    print log_file "ml:\t";
    for ($i=1;$i<=$n_motifs;$i++) {
	print log_file "$ml[$i]\t";
    }
    print log_file "\n";
    print log_file "thres:\t";
    for ($i=1;$i<=$n_motifs;$i++) {
	print log_file "$thres[$i]\t";
    }
    print log_file "\n";

    if ($sp2 ne "none") {
	# load mm2hs map #
	@mm2hs_map_all=read_matrix_file($mm2hs_map_filename_all,2);       # 09-05-2003 read into a matrix
	$n_mm2hs_map_all=$#mm2hs_map_all;
	if ($mm2hs_map_filename_cl ne "default_map") {
	    @mm2hs_map_cl=read_matrix_file($mm2hs_map_filename_cl,2);         # 09-05-2003 read into a matrix
	    $n_mm2hs_map_cl=$#mm2hs_map_cl;
	} else {
	    $n_mm2hs_map_cl=$n_upstream_cl_sp1;
	    for ($i=1;$i<=$n_mm2hs_map_cl;$i++) {
		$mm2hs_map_cl[$i][1]=$i;
		$mm2hs_map_cl[$i][2]=$i;
	    }
	}
	if ($verbose eq 1) {
	    print log_file "% mm2hs_map_filename_all=$mm2hs_map_filename_all\tn_mm2hs_map_all=$n_mm2hs_map_all\n";
	    print log_file "% mm2hs_map_filename_cl=$mm2hs_map_filename_cl\tn_mm2hs_map_cl=$n_mm2hs_map_cl\n";
	}
	# count the number of genes that have a homologous #
	$n_all_withhomologous_map=0;
	for ($i=1;$i <= $n_mm2hs_map_all ;$i++) {
	    $temp_n=$mm2hs_map_all[$i][1];
	    if ($temp_n>0) {
		$temp_n=$mm2hs_map_all[$i][2];
		if ($temp_n>0) {
		    $n_all_withhomologous_map++;
		} 
	    }
	}	
	$n_cl_withhomologous_map=0;
	for ($i=1;$i <= $n_mm2hs_map_cl ; $i++) {
	    $temp_n=$mm2hs_map_cl[$i][1];
	    if ($temp_n>0) {
		$temp_n=$mm2hs_map_cl[$i][2];
		if ($temp_n>0) {
		    $n_cl_withhomologous_map++;
		} 
	    } 
	}    
	print log_file "% n_mm2hs_map_cl=$n_mm2hs_map_cl\nn_cl_withhomologous=$n_cl_withhomologous_map\n";
	print log_file "% n_mm2hs_map_all=$n_mm2hs_map_all\nn_all_withhomologous_map=$n_all_withhomologous_map\n";
    }

################
# initialize n #
################
    for ($i=0;$i<=$n_motifs;$i++) {
	for ($j=0;$j<=4;$j++) {
	    $n[$i][$j]=-1;               # n contains the number of occurrences for each motif (rows) in the 4 seq files (1=all,sp1; 2=all,sp2; 3=cl,sp1; 4=cl,sp2)
	}
    }

#########################################
# compute the list of max_n_transcripts #
#########################################
    print log_file "% max_n_transcripts\n";
    for ($i=1;$i<=$n_upstream_cl_sp1;$i++) {
	($max_n_transcripts)=get_max_n_transcripts($i,$n_upstream_cl_sp1,$threshold_binom_sp1_c,$n_upstream_all_sp1);
	$max_n_transcripts_array[$i]=$max_n_transcripts;
	print log_file "\t$i\t$max_n_transcripts\n";
    }

###################
# begin the begin #
###################
    open (cooc_file,$input_filename) || die "could not open $input_filename for reading";
    $n_processed=0;
    while ( $record = <cooc_file> ) {
	chomp $record;
	if (is_comment($record) ne "true") {
	    $n_processed++;
	    print "n_processed=$n_processed\trecord=$record\n";

	    if ($n_processed>=$init_i) {
		if ($verbose eq 1) {
		    print log_file "n_processed=$n_processed\trecord=$record\n";
		}

		# default initialization
		$errorcode_cl_sp1=0;$errorcode_all_sp1=0;$n_cooc_all_sp1=-1;$n_cooc_tr_all_sp1=-1;$p_binom_sp1=1;$n_cooc_cl_sp1=-1;$n_cooc_tr_cl_sp1=-1;
		$n_all_withhomologous_sp1=-1;$n_cl_withhomologous_sp1=-1;
		$errorcode_cl_sp2=0;$errorcode_all_sp2=0;$n_cooc_all_sp2=-1;$n_cooc_tr_all_sp2=-1;$p_binom_sp2=1;$n_cooc_cl_sp2=-1;$n_cooc_tr_cl_sp2=-1;
		$n_all_withhomologous_sp2=-1;$n_cl_withhomologous_sp2=-1;

		@splitarray=split /\t/,$record;
		$i1=$splitarray[0];  # index of motif 1
		$i2=$splitarray[1];  # index of motif 2
		$i3=$splitarray[2];  # index of motif 3 (-1 for 2-motif interaction)
		$i4=$splitarray[3];  # index of motif 4 (-1 for 2- and 3-motif interactions)
		$ml1=$splitarray[4]; # motif lengths
		$ml2=$splitarray[5];
		$ml3=$splitarray[6];
		$ml4=$splitarray[7];
		$wm1=$wm[$i1];       # weight matrix for motif i1 
		$t1=$thres[$i1];     # threshold for motif i1
		$wm2=$wm[$i2];       # weight matrix for motif i2
		$t2=$thres[$i2];     # threshold for motif i2
		if ($i3>0) {
		    $wm3=$wm[$i3];
		    $t3=$thres[$i3];
		}
		if ($i4>0) {
		    $wm4=$wm[$i4];
		    $t4=$thres[$i4];
		}
		$n_cooc_cl=$splitarray[8];        # number of co-ocurrences in the cluster set (previously computed in cooc_4m)
		$n_cooc_cl_tr=$splitarray[9];     # number of transcripts with co-occurrences (previously computed in cooc_4m)

		if ($verbose eq 1) {
		    print log_file "\ti1=$i1\ti2=$i2\ti3=$i3\ti4=$i4\tml1=$ml1\tml2=$ml2\tml3=$ml3\tml4=$ml4\n";
		}
		$compute_cooc_cl_sp1=1;
		if ( ($ml1>0) & ($ml1<$ml_threshold) ) {
		    $compute_cooc_cl_sp1=0;print log_file "ml1=$ml1\t$i1 $i2 $i3 $i4 did not pass threhsold = $ml_threshold; skipping...\n";
		}
		if ( ($ml2>0) & ($ml2<$ml_threshold) ) {
		    $compute_cooc_cl_sp1=0;print log_file "ml2=$ml2\t$i1 $i2 $i3 $i4 did not pass threhsold = $ml_threshold; skipping...\n";
		}
		if ( ($ml3>0) & ($ml3<$ml_threshold) ) {
		    $compute_cooc_cl_sp1=0;print log_file "ml3=$ml3\t$i1 $i2 $i3 $i4 did not pass threhsold = $ml_threshold; skipping...\n";
		}
		if ( ($ml4>0) & ($ml4<$ml_threshold) ) {
		    $compute_cooc_cl_sp1=0;print log_file "ml4=$ml4\t$i1 $i2 $i3 $i4 did not pass threhsold = $ml_threshold; skipping...\n";
		}
		if ($compute_cooc_cl_sp1>0) {
		    ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_cl_sp1,$logtext)=scan_4motifs_cl($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$sp1,$get_count);
		    print log_file "$logtext\n";
		    if ($n1<=0) {
			print "warning! i1=$i1\tn1=$n1\n";
			$compute_cooc_cl_sp1=0;
		    }
		    if ( ($i2>0) & ($n2<=0) ) {
			print "warning! i2=$i2\tn2=$n2\n";
			$compute_cooc_cl_sp1=0;
		    }
		    if ( ($i3>0) & ($n3<=0) ) {
			print "warning! i3=$i3\tn3=$n3\n";
			$compute_cooc_cl_sp1=0;
		    }
		    if ( ($i4>0) & ($n4<=0) ) {
			print "warning! i4=$i4\tn4=$n4\n";
			$compute_cooc_cl_sp1=0;
		    }
		    if ($verbose eq 1) {
			print log_file "\t$logtext\n";
		    }
		} else {
		    print wrn_file "$i1\t$i2\t$i3\t$i4\tml1=$ml1\tml2=$ml2\tml3=$ml3\tml4=$ml4\tskipped because motif length is less than threshold\n";
		}
		if ($compute_cooc_cl_sp1>0) {
		    $max_n_transcripts=-1;
		    $n_transcripts=$n_upstream_cl_sp1;
		    $transcript_map="none";
		    $pos_format=-1;
		    ($errorcode_cl_sp1,$n_cooc_cl_sp1,$n_cooc_tr_cl_sp1)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename_sp1_cl,$max_upstream_length,$max_exon_length,$max_intron_length,$transcript_map,$sortpos,$pos_format);
		    ($system_status,$msg)=my_rename("cooc_single_1234m_v${cooc_c_code_version}.log.txt","cooc_cl_${sp1}.txt",0,1);
		    print log_file "\t$msg\n";
		} else {
		    $errorcode_cl_sp1=0;$n_cooc_cl_sp1=0;$n_cooc_tr_cl_sp1=0;
		}    # close check on compute_coo_cl_sp1 > 0
		if ($verbose eq 1) {
		    print log_file "\t\tn_cooc_cl_${sp1}=$n_cooc_cl_sp1\tn_cooc_tr_cl_${sp1}=$n_cooc_tr_cl_sp1\terrorcode_cl_${sp1}=$errorcode_cl_sp1\n";
		}

		if ( ($n_cooc_cl_tr ne $n_cooc_tr_cl_sp1) & ($compute_cooc_cl_sp1 eq 1) ) {
		    $n_errors++;
		    print err_file "after cooc4m\t$record\tn_cooc_cl=${n_cooc_cl}\tn_cooc_cl_tr=${n_cooc_cl_tr}\tn_cooc_cl_sp1=${n_cooc_cl_sp1}\tn_cooc_tr_cl_sp1=${n_cooc_tr_cl_sp1}\n";
		} else {
		    if ($n_cooc_cl ne $n_cooc_cl_sp1) {
			$n_warnings++;
			print wrn_file "after cooc4m\t$record\tn_cooc_cl=${n_cooc_cl}\tn_cooc_cl_tr=${n_cooc_cl_tr}\tn_cooc_cl_sp1=${n_cooc_cl_sp1}\tn_cooc_tr_cl_sp1=${n_cooc_tr_cl_sp1}\n";
		    }
		}

		if ($n_cooc_cl_sp1<=0) {
		    print "n_cooc_cl_sp1=$n_cooc_cl_sp1\n";
		} else {
		    ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_all_sp1,$logtext)=scan_4motifs_all($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$sp1,$get_count,$singlefile,$seq_source,$read_all);
		    if ($verbose eq 1) {
			print log_file "$logtext\n";
		    }
		    if ($compute_cooc_all_sp1>0) {
			#$max_n_transcripts=-1;
			$max_n_transcripts=$max_n_transcripts_array[$n_cooc_tr_cl_sp1];
			if (length($max_n_transcripts)<1) {
			    print "\n\nerror!\tcooc_checkresults_v12.pl\tn_cooc_tr_cl_sp1=${n_cooc_tr_cl_sp1}\tmax_n_transcripts=$max_n_transcripts\n\n";
			    exit;
			}
			print "n_cooc_cl_sp1=$n_cooc_cl_sp1\tn_cooc_tr_cl_sp1=$n_cooc_tr_cl_sp1\tmax_n_transcripts=$max_n_transcripts\n";#exit;
			$n_transcripts=$n_upstream_all_sp1;
			$pos_format=1;
			($errorcode_all_sp1,$n_cooc_all_sp1,$n_cooc_tr_all_sp1)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename_sp1_all,$max_upstream_length,$max_exon_length,$max_intron_length,$transcript_map_filename_sp1,$sortpos,$pos_format);
			($system_status,$msg)=my_rename("cooc_single_1234m_v${cooc_c_code_version}.log.txt","cooc_all_${sp1}.txt",0,1);
			print log_file "\t$msg\n";
			if ($verbose eq 1) {
			    print log_file "\t\tn1=$n1\tn2=$n2\tn3=$n3\tn4=$n4\n";
			    print log_file "\t\tn_cooc_all_${sp1}=$n_cooc_all_sp1\tn_cooc_tr_all_${sp1}=$n_cooc_tr_all_sp1\terrorcode_all_${sp1}=$errorcode_all_sp1\n";
			}
			
			# compute statistics for enrichment
			($p_binom_sp1,$p_sp1)=computebinom($n_upstream_all_sp1,$n_cooc_tr_all_sp1,$n_upstream_cl_sp1,$n_cooc_tr_cl_sp1);
			if ($verbose eq 1) {
			    print log_file "\t\tp_binom_${sp1}=${p_binom_sp1}\n";
			    print log_file "\t\tp_${sp1}=${p_sp1}\n";
			}
		    }       # if compute_cooc_all_sp1>0
		}           # if n_cooc_cl_sp1>0 


		if (${sp2} ne "none") {
		    ########################
		    # scan cl upstream sp2 #
		    ########################
		    ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_cl_sp2,$logtext)=scan_4motifs_cl($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$sp2,$get_count);
		    if ($verbose eq 1) {
			print log_file "\t$logtext\n";
		    }
		    if ($n1<=0) {
			print "warning! i1=$i1\tn1=$n1\tsp2\n";$compute_cooc_cl_sp2=0;
		    }
		    if ( ($i2>0) & ($n2<=0) ) {
			print "warning! i2=$i2\tn2=$n2\tsp2\n";$compute_cooc_cl_sp2=0;
		    }
		    if ( ($i3>0) & ($n3<=0) ) {
			print "warning! i3=$i3\tn3=$n3\tsp2\n";$compute_cooc_cl_sp2=0;
		    }
		    if ( ($i4>0) & ($n4<=0) ) {
			print "warning! i4=$i4\tn4=$n4\tsp2\n";$compute_cooc_cl_sp2=0;
		    }
		    if ($compute_cooc_cl_sp2>0) {
			$max_n_transcripts=-1;
			$n_transcripts=$n_upstream_cl_sp2;
			$transcript_map="none";
			$pos_format=-1;
			($errorcode_cl_sp2,$n_cooc_cl_sp2,$n_cooc_tr_cl_sp2)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename_sp2_cl,$max_upstream_length,$max_exon_length,$max_intron_length,$transcript_map,$sortpos,$pos_format);
			($system_status,$msg)=my_rename("cooc_single_1234m_v${cooc_c_code_version}.log.txt","cooc_cl_${sp2}.txt",0,1);
			print log_file "\t$msg\n";
			if ($verbose eq 1) {
			    print log_file "\tn1=$n1\tn2=$n2\tn3=$n3\tn4=$n4\tn_cooc_cl_${sp2}=${n_cooc_cl_sp2}\tn_cooc_tr_cl_${sp2}=${n_cooc_tr_cl_sp2}\terrorcode_cl_${sp2}=${errorcode_cl_sp2}\tcompute_cooc_cl_${sp2}=${compute_cooc_cl_sp2}\n";
			}
		    }

		    if ($n_cooc_cl_sp2>0) {
			#########################
			# scan all upstream sp2 #
			#########################
			($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_all_sp2,$logtext)=scan_4motifs_all($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$sp2,$get_count,$singlefile,$seq_source,$read_all);
			if ($verbose eq 1) {
			    print log_file "$logtext\n";
			}
			if ($compute_cooc_all_sp2>0) {
			    #$max_n_transcripts=-1;
			    $max_n_transcripts=$max_n_transcripts_array[$n_cooc_tr_cl_sp2];
			    $n_transcripts=$n_upstream_all_sp2;
			    $pos_format=1;
			    ($errorcode_all_sp2,$n_cooc_all_sp2,$n_cooc_tr_all_sp2)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename_sp2_all,$max_upstream_length,$max_exon_length,$max_intron_length,$transcript_map_filename_sp2,$sortpos,$pos_format);
			    ($system_status,$msg)=my_rename("cooc_single_1234m_v${cooc_c_code_version}.log.txt","cooc_all_${sp2}.txt",0,1);
			    print log_file "\t$msg\n";
			    if ($verbose eq 1) {
				print log_file "\t\tn1=$n1\tn2=$n2\tn3=$n3\tn4=$n4\n";
				print log_file "\t\tn_cooc_all_${sp2}=${n_cooc_all_sp2}\tn_cooc_tr_all_${sp2}=${n_cooc_tr_all_sp2}\terrorcode_all_${sp2}=${errorcode_all_sp2}\n";
			    }
			    
			    # compute statistics for enrichment
			    ($p_binom_sp2,$p_sp2)=computebinom($n_upstream_all_sp2,$n_cooc_tr_all_sp2,$n_upstream_cl_sp2,$n_cooc_tr_cl_sp2);
			    if ($verbose eq 1) {
				print log_file "\t\tp_binom_${sp2}=${p_binom_sp2}\n";
				print log_file "\t\tp_${sp2}=${p_sp2}\n";
			    }
			}     # close if ($compute_cooc_all_sp2>0)
		    }         # close if ($compute_cooc_cl_sp2>0)
		}             # close if sp2 ne "none"
		
		$passed_threshold=0;
		if ($sp2 ne "none") {
		    if ($p_binom_sp2<$threshold_binom_sp2) {
			$passed_threshold=1;
		    }
		    if ($p_binom_sp1<$threshold_binom_sp1) {
			$passed_threshold=1;
		    }
		    if ( ($p_binom_sp2 < $threshold_binom_both) and ($p_binom_sp1 <$threshold_binom_both) ) {
			$passed_threshold=1;
		    }
		} else {
		    if ($p_binom_sp1 < $threshold_binom_sp1) {
			$passed_threshold=1;
		    }
		}

		if ($passed_threshold eq 1) {
		    if ($verbose) {
			print log_file "\n\tpassed threshold\n";
		    }

		    $n_cl_ov=-1;$n_tr_cl_ov=-1;$n_cl_pos=-1;$n_tr_cl_pos=-1;$p_cl_trov=-1;$p_cl_posov=-1;$n_cl_withhomologous_sp2=-1;
		    $n_all_ov=-1;$n_tr_all_ov=-1;$n_all_pos=-1;$n_tr_all_pos=-1;$p_all_trov=-1;$p_all_posov=-1;$n_all_withhomologous_sp2=-1;
		    if ($sp2 ne "none") {
			if ( ($n_cooc_tr_cl_sp1>0) & ($n_cooc_tr_cl_sp2>0) ) {
			    # map the sp1 transcripts in the coexpressed set to the orthologous genes in sp2 and re-sorting the file by transcript number
			    if ($verbose eq 1) {
				print log_file "\n\tmapping the $sp1 genes to orthologous $sp2 genes ";
				print log_file "(n_cooc_tr_cl_${sp1}=${n_cooc_tr_cl_sp1}\tn_cooc_tr_cl_${sp2}=${n_cooc_tr_cl_sp2})\n";
			    }
			    if ($sp1=="hs") {
				$ind=1;
			    } else {
				$ind=2;
			    }
			    ($n_cl_withhomologous_sp1,$n_miss,$n_notfound)=map_orthologous($sp1,$ind,"cl",@mm2hs_map_cl);
			    if ($verbose eq 1) {
				print log_file "\t\tn_cl_withhomologous_${sp1}=${n_cl_withhomologous_sp1}\tn_miss=$n_miss\tn_notfound=$n_notfound\n";
			    }

			    # now compare the hs and mm lists (cl)
			    $t0=time();
			    if ($verbose eq 1) {
				print log_file "\t->comparing cooc_cl_${sp1}.mapped.txt to cooc_cl_${sp2}.txt\n";
			    }
			    $code_exec="${CODE_DIR}/cpp/executables/compare_cooc_lists.exe cooc_cl_${sp1}.mapped.txt cooc_cl_${sp2}.txt ${ml1} ${ml2} ${ml3} ${ml4} ${n_upstream_cl} ${upstream_length} 0 ${p_ortho_distance}";
			    if ($verbose eq 1) {
				print log_file "\t\tcode_exec=$code_exec\n";
			    }
			    @output=`${code_exec}`;
			    # get some output statistics
			    ($n_cl_ov,$n_matches)=extract_text_from_array("n_ov",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_ov (cl), n_matches=$n_matches\n";
			    }
			    ($n_tr_cl_ov,$n_matches)=extract_text_from_array("n_tr_ov",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_tr_ov (cl), n_matches=$n_matches\n";
			    }
			    ($n_cl_pos,$n_matches)=extract_text_from_array("n_pos",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_pos (cl), n_matches=$n_matches\n";
			    }
			    ($n_tr_cl_pos,$n_matches)=extract_text_from_array("n_tr_pos",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_tr_pos (cl), n_matches=$n_matches\n";
			    }
			    ($system_status,$n_matches)=extract_text_from_array("system_status",@output);
			    if ($system_status ne 0) {
				print "error!\tcooc_checkresults_v12.pl\ncode_exec=$code_exec\nsystem_status=$system_status\nexiting...\n";
				exit;
			    }
			    
			    $et=time();
			    $et=$et-$t0;
			    if ($verbose eq 1) {
				print log_file "\t\tn_cl_ov=$n_cl_ov\n";
				print log_file "\t\tn_tr_ov=$n_tr_cl_ov\n";
				print log_file "\t\tn_cl_pos=$n_cl_pos\n";
				print log_file "\t\tn_tr_cl_pos=$n_tr_cl_pos\n";
				print log_file "\t\tsystem_status=${system_status}\n";
				print log_file "\t\tcpu time=$et\n";
			    }
			    
			    # copy  the contents of the pos_overlap and tr_overlap files
			    print both_troverlap_cl_file "$i1\t$i2\t$i3\t$i4\t${n_cl_ov}\t${n_tr_cl_ov}\t${n_cl_pos}\t${n_tr_cl_pos}\t1\t1\n";
			    ($comments,$coretext)=add_genename("tr_overlap.txt",@info_cl_sp2);
			    print both_troverlap_cl_file "% $comments\n$coretext"; # note that coretex already has a carriage return at the end		    
			    print both_posoverlap_cl_file "$i1\t$i2\t$i3\t$i4\t${n_cl_ov}\t${n_tr_cl_ov}\t${n_cl_pos}\t${n_tr_cl_pos}\t1\t1\n";
			    ($comments,$coretext)=add_genename("pos_overlap.txt",@info_cl_sp2);
			    print both_posoverlap_cl_file "% $comments\n$coretext";
			    
			    # count the number of mouse genes in cooc_cl_sp2.txt that have an orthologous sp1gene
			    ($n_cl_withhomologous_sp2,$n_miss,$n_notfound)=count_genes_withhomologous("mm",1,"cl",@mm2hs_map_cl);       # this will need to be checked		    
			    # compute statistics for the number of overlapping transcripts
			    $p_cl_trov=compute_random_overlap($n_cl_withhomologous_sp1,$n_cl_withhomologous_sp2,$n_cl_withhomologous_map,$n_tr_cl_ov,$n_iter_random_overlap);
			    
			    ($p_cl_posov)=computebinom_dist($n_tr_cl_ov,$n_tr_cl_pos,$p_ortho_distance);
			    if ($verbose eq 1) {
				print log_file "\tp_cl_trov=$p_cl_trov\n";
				print log_file "\tp_cl_posov=$p_cl_posov\n";
			    }
			} else {         # close check on n_cooc_tr_cl_sp1>0 and n_cooc_tr_cl_sp2>0
			    if ($verbose eq 1) {
				print log_file "\tn_cooc_tr_cl_${sp1}=${n_cooc_tr_cl_sp1}\n\tn_cooc_tr_cl_${sp2}=${n_cooc_tr_cl_sp2}\n\t\tno comparison performed\n";
			    }

			    # print a line anyway in the trov and posov files so that we can more easily use these files later!
			    print both_troverlap_cl_file "$i1\t$i2\t$i3\t$i4\t${n_cl_ov}\t${n_tr_cl_ov}\t${n_cl_pos}\t${n_tr_cl_pos}\n";
			    print both_troverlap_cl_file "% no tr overlap\n"; 

			    print both_posoverlap_cl_file "$i1\t$i2\t$i3\t$i4\t${n_cl_ov}\t${n_tr_cl_ov}\t${n_cl_pos}\t${n_tr_cl_pos}\n";
			    print both_posoverlap_cl_file "% no pos overlap\n";
			}

			# map the sp2 transcripts to the sp1 orthologous and re-sort the file by transcript number
			if ( ($n_cooc_tr_all_sp1 > 0) & ($n_cooc_tr_all_sp2 > 0) ) {
			    if ($verbose eq 1) {
				print log_file "\n\tmapping the sp1 genes to sp2 genes\t(n_cooc_tr_all_${sp1}=${n_cooc_tr_all_sp1}\tn_cooc_tr_all_${sp2}=${n_cooc_tr_all_sp2})\n";
			    }		    
			    if ($sp2 eq "hs") {
				$ind=1;
			    } else {
				$ind=2;
			    }
			    ($n_all_withhomologous_sp1,$n_miss,$n_notfound)=map_orthologous($sp2,$ind,"all",@mm2hs_map_all);
			    if ($verbose eq 1) {
				print log_file "\t\tn_all_withhomologous_sp1=$n_all_withhomologous_sp1\tn_miss=$n_miss\tn_notfound=$n_notfound\n";
			    }

			    # compare the sp1 and sp2 lists
			    $t0=time();
			    if ($verbose eq 1) {
				print log_file "\t->comparing cooc_all_${sp1}.mapped.txt to cooc_all_${sp2}.txt\n";
			    }
			    $code_exec="${CODE_DIR}/cpp/executables/compare_cooc_lists.exe cooc_all_${sp1}.mapped.txt cooc_all_${sp2}.txt ${ml1} ${ml2} ${ml3} ${ml4} ${n_upstream_all} ${upstream_length}";
			    if ($verbose eq 1) {
				print log_file "\t\tcode_exec=$code_exec\n";
			    }
			    @output=`${code_exec}`;
			    # get some output statistics
			    ($n_all_ov,$n_matches)=extract_text_from_array("n_ov",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_ov (all), n_matches=$n_matches\n";
			    }
			    ($n_tr_all_ov,$n_matches)=extract_text_from_array("n_tr_ov",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_tr_ov (all), n_matches=$n_matches\n";
			    }
			    ($n_all_pos,$n_matches)=extract_text_from_array("n_pos",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_pos (all), n_matches=$n_matches\n";
			    }
			    ($n_tr_all_pos,$n_matches)=extract_text_from_array("n_tr_pos",@output);
			    if ($n_matches ne 1) {
				print "\t\twarning searching for n_tr_pos (all), n_matches=$n_matches\n";
			    }
			    ($system_status,$n_matches)=extract_text_from_array("system_status",@output);
			    
			    if ($verbose eq 1) {
				print log_file "\t\tn_all_ov=$n_all_ov\n";
				print log_file "\t\tn_tr_all_ov=$n_tr_all_ov\n";
				print log_file "\t\tn_all_pos=$n_all_pos\n";
				print log_file "\t\tn_tr_all_pos=$n_tr_all_pos\n";
				print log_file "\t\tsystem_status=${system_status}\n";
			    }

			    # copy  the contents of the pos_overlap and tr_overlap files
			    print both_troverlap_all_file "$i1\t$i2\t$i3\t$i4\t${n_all_ov}\t${n_tr_all_ov}\t${n_all_pos}\t${n_tr_all_pos}\t1\t1\n";
			    ($comments,$coretext)=add_genename("tr_overlap.txt",@info_all_sp2);
			    print both_troverlap_all_file "% $comments\n$coretext";
			    print both_posoverlap_all_file "$i1\t$i2\t$i3\t$i4\t${n_all_ov}\t${n_tr_all_ov}\t${n_all_pos}\t${n_tr_all_pos}\t1\t1\n";
			    ($comments,$coretext)=add_genename("pos_overlap.txt",@info_all_sp2);
			    print both_posoverlap_all_file "% $comments\n$coretext";
			    # count the number of sp2 genes in cooc_all_sp2.txt that have an orthologous sp1 gene
			    ($n_all_withhomologous_sp2,$n_miss,$n_notfound)=count_genes_withhomologous("mm",1,"all",@mm2hs_map_all);
			    # compute statistics for the number of overlapping transcripts
			    $p_all_trov=compute_random_overlap($n_all_withhomologous_sp1,$n_all_withhomologous_sp2,$n_all_withhomologous_map,$n_tr_all_ov,$n_iter_random_overlap);
			    ($p_all_posov)=computebinom_dist($n_tr_all_ov,$n_tr_all_pos,$p_ortho_distance);
			    if ($verbose eq 1) {
				print log_file "\t\tp_all_trov=$p_all_trov\n";
				print log_file "\t\tp_all_posov=$p_all_posov\n";
			    }
			} else {        # close check on n_cooc_tr_all_sp1>0 and n_cooc_tr_all_sp2>0
			    
			    # print a line in trov and posov files so that we can more easily process these files
			    print both_troverlap_all_file "$i1\t$i2\t$i3\t$i4\t${n_all_ov}\t${n_tr_all_ov}\t${n_all_pos}\t${n_tr_all_pos}\n";
			    print both_troverlap_all_file "% no tr overlap\n";

			    print both_posoverlap_all_file "$i1\t$i2\t$i3\t$i4\t${n_all_ov}\t${n_tr_all_ov}\t${n_all_pos}\t${n_tr_all_pos}\n";
			    print both_posoverlap_all_file "% no pos overlap\n";
			}
		    }                   # close check on sp2 ne "none"

		    ################################################
		    # if it is an interesting case, then report it #
		    ################################################
		    $record=~s/\t$//;        # remove last tab
		    print output_file "${record}\t${n_cooc_all_sp1}\t${n_cooc_tr_all_sp1}\t${n_all_withhomologous_sp1}\t${n_cooc_cl_sp1}\t${n_cooc_tr_cl_sp1}\t${n_cl_withhomologous_sp1}\t${n_cooc_all_sp2}\t${n_cooc_tr_all_sp2}\t${n_all_withhomologous_sp2}\t${n_cooc_cl_sp2}\t${n_cooc_tr_cl_sp2}\t${n_cl_withhomologous_sp2}\t${p_binom_sp1}\t${p_binom_sp2}\t${n_cl_ov}\t${n_tr_cl_ov}\t${p_cl_trov}\t${n_cl_pos}\t${n_tr_cl_pos}\t${p_cl_posov}\t${n_all_ov}\t${n_tr_all_ov}\t${p_all_trov}\t${n_all_pos}\t${n_tr_all_pos}\t${p_all_posov}\t$errorcode_all_sp1\t$errorcode_cl_sp1\t$errorcode_all_sp2\t$errorcode_cl_sp2\n";
		    print log_file "${record}\t${n_cooc_all_sp1}\t${n_cooc_tr_all_sp1}\t${n_all_withhomologous_sp1}\t${n_cooc_cl_sp1}\t${n_cooc_tr_cl_sp1}\t${n_cl_withhomologous_sp1}\t${n_cooc_all_sp2}\t${n_cooc_tr_all_sp2}\t${n_all_withhomologous_sp2}\t${n_cooc_cl_sp2}\t${n_cooc_tr_cl_sp2}\t${n_cl_withhomologous_sp2}\t${p_binom_sp1}\t${p_binom_sp2}\t${n_cl_ov}\t${n_tr_cl_ov}\t${p_cl_trov}\t${n_cl_pos}\t${n_tr_cl_pos}\t${p_cl_posov}\t${n_all_ov}\t${n_tr_all_ov}\t${p_all_trov}\t${n_all_pos}\t${n_tr_all_pos}\t${p_all_posov}\t$errorcode_all_sp1\t$errorcode_cl_sp1\t$errorcode_all_sp2\t$errorcode_cl_sp2\n";
		    ###############################################################
		    # if it is an interesting case, then store also the positions #
		    ###############################################################		
		    # coexp genes
		    if ($verbose eq 1) {
			print log_file "\treporting the positions in the coexpressed genes\n";
		    }
		    $t0=time();
		    if ($sp2 ne "none") {
			# sp2
			print pos_cl_sp2_file "$i1\t$i2\t$i3\t$i4\t$p_binom_sp2\t$n_cooc_tr_cl_sp2\t$n_cooc_tr_all_sp2\t$errorcode_cl_sp2\t$errorcode_all_sp2\t1\n";
			if (-e "cooc_cl_${sp2}.txt") {
			    ($comments,$coretext)=add_genename("cooc_cl_${sp2}.txt",@info_cl_sp2);
			    print pos_cl_sp2_file "% $comments\n$coretext";
			}
		    }
		    print pos_cl_sp1_file "$i1\t$i2\t$i3\t$i4\t$p_binom_sp1\t$n_cooc_tr_cl_sp1\t$n_cooc_tr_all_sp1\t$errorcode_cl_sp1\t$errorcode_all_sp1\t1\n";
		    if (-e "cooc_cl_${sp1}.txt") {
			($comments,$coretext)=add_genename("cooc_cl_${sp1}.txt",@info_cl_sp1);
			print pos_cl_sp1_file "% $comments\n$coretext";
		    }
		    $et=time();
		    $et=$et-$t0;
		    if ($verbose eq 1) {
			print log_file "\t\tcpu time=$et (for reporting positions in set=cluster)\n";
		    }
		    # all genes
		    $t0=time();
		    if ($verbose eq 1) {
			print log_file "\treporting the positions in all genes\n";
		    }
		    if ($sp2 ne "none") {
			print pos_all_sp2_file "$i1\t$i2\t$i3\t$i4\t$p_binom_sp2\t$n_cooc_tr_cl_sp2\t$n_cooc_tr_all_sp2\t$errorcode_cl_sp2\t$errorcode_all_sp2\t1\n";
			if (-e "cooc_all_${sp2}.txt") {
			    ($comments,$coretext)=add_genename("cooc_all_${sp2}.txt",@info_all_sp2);
			    print pos_all_sp2_file "% $comments\n$coretext";
			}
		    }
		    print pos_all_sp1_file "$i1\t$i2\t$i3\t$i4\t$p_binom_sp1\t$n_cooc_tr_cl_sp1\t$n_cooc_tr_all_sp1\t$errorcode_cl_sp1\t$errorcode_all_sp1\t1\n";
		    if (-e "cooc_all_${sp1}.txt") {
			($comments,$coretext)=add_genename("cooc_all_${sp1}.txt",@info_all_sp1);
			print pos_all_sp1_file "% $comments\n$coretext";
		    }
		    $et=time();
		    $et=$et-$t0;
		    if ($verbose eq 1) {
			print log_file "\t\tcpu time=$et (for reporting positions in set=all)\n";
		    }
		} else {
		    # did not pass the p_binom thresholds, print out to the rejected list
		    print rejected_file "$record\t$n_cooc_all_sp1\t$n_cooc_tr_all_sp1\t$n_cooc_cl_sp1\t$n_cooc_tr_cl_sp1\t$n_cooc_all_sp2\t$n_cooc_tr_all_sp2\t$n_cooc_cl_sp2\t$n_cooc_tr_cl_sp2\t$p_binom_sp1\t$p_binom_sp2\t$errorcode_all_sp1\t$errorcode_cl_sp1\t$errorcode_all_sp2\t$errorcode_cl_sp2\n";
		}

		# remove temporary files
		if ($remove_temporary eq 1) {
		    $unlinked_files=unlink("compare_cooc_lists.log.txt","pos_overlap.txt","cooc_all_${sp1}.tempmap.txt","cooc_all_${sp1}.mapped.txt","cooc_cl_${sp1}.tempmap.txt","cooc_cl_${sp1}.mapped.txt","tr_overlap.txt","cooc_cl_${sp1}.txt","cooc_all_${sp2}.txt","cooc_all_${sp1}.txt","cooc_cl_${sp2}.txt");
		    print "unlinked $unlinked_files files\n";
		    if (-e "cooc_cl_${sp1}.txt") {
			print "i could not delete cooc_cl_${sp1}.txt\n";
			exit;
		    }
		}
	    } else {
		print "skipping n_processed=$n_processed (init_i=$init_i)\n";
	    }     # close check on n_processed>=init_i
	    $partial_t=time();
	    $partial_t-=$init_t;
	    if ($verbose eq 1) {
		print log_file "\n\tpartial_t=$partial_t\n";
	    }    
	}  else {
	    # this is a comment line, just print the information afte replacing ">" by "%"
	    $record=~s/\>/\%/;
	    print output_file "$record\n";
	}         # close check on whether record is a comment file 
    }             # close the loop on entries from cooc_file

    $finit_t=time();
    $finit_t=$finit_t-$init_t;
    if ($verbose eq 1) {
	print log_file "total processing time = $finit_t secs\n";
    }
    print "total processing time = $finit_t\nnow closing the files and saying bye-bye.\n";

    close (cooc_file);
    close (log_file);
    close (err_file);
    close (wrn_file);
    close (output_file);
    close (pos_cl_sp1_file);
    close (pos_all_sp1_file);
    close (rejected_file);
    if ($sp2 ne "none") {
	close (pos_all_sp2_file);
	close (pos_cl_sp2_file);
	close (both_troverlap_all_file);
	close (both_posoverlap_all_file);
	close (both_troverlap_cl_file);
	close (both_posoverlap_cl_file);
    }

    print "system_status=0\n";
}

# end main

sub scan_motif{
    # usage:
    # $n_transcripts_occurrences=scan_motif($wm,$upstream_filename,$scan_output_filename,$motif_index,$species,$ml,$thres,$get_count,$rc);
    #
    # scan the occurrences of motif given by weight matrix $wm in the $upstream_filename file and writes the output to $scan_output_filename
    # note that in the current version (v6) we directly use the output from the c program
    # note also that here we require the number of sequences used to compute the motif model
    # if get_count is 1 then we do a word count and return the number of lines, otherwise, return -1 if there was an error and 0 if everything was ok
    #
    # scan_motif_v8/9 <background_sequence> <weight_matrix> <motif_length> <threshold> <max_upstream_length> <rc> <sortpos> <verbose>
    # 09-04-2003: changed name of temporary weight matrix file to cooc_checkresults_wm.txt
    # 09-05-2003: add computation time information
    # 12-28-2003: clean up, store also strand information, call scan_motif_v9, remove temporary files
    # 02-02-2004: use scan_motif_version=10
    # 02-16-2004: use scan_motif_version=11
    # 02-26-2004: add file_label in the name of the temporary file with the weight matrix
    # 03-08-2004: add exp_id to the name of the temporary file
    #
    # uses global variables: @n, $sp1, $sp2, $file_label log_file, verbose, $exp_id

    my $code_exec;                    # code to be executed
    my $col_index;                    # column to store the number of occurrences in the global variable @n
    my $et;                           # end time
    my $get_count=$_[7];              # if get_count is 1 then we do a word count and return the number of lines 
    my $ml=$_[5];                     # motif length
    my $motif_filename="${TEMP_DIR}/cooc_checkresults_wm.${exp_id}.${file_label}.txt";
    my $motif_index=$_[3];            # motif index
    my $n_comments;                   # used to count the number of lines in the scan output
    my $n_transcripts_occurrences=-1; # number of transcripts with motif occurrences
    my $rc=$_[8];                     # 1 to consider both strands, 0 otherwise
    my $scan_motif_version=11;        # version of the c code for scan_motif
    my $scan_output_filename=$_[2];   # output filename
    my $scanlog_filename;             # scan log filename
    my $sortpos=1;                    # sort positions in call to scan code
    my $species=$_[4];                # species
    my $strand_filename;              # scan strand filename
    my $system_status;                # system status
    my $t0=time();                    # initial time
    my $thres=$_[6];                  # threshold 
    my $upstream_filename=$_[1];      # sequence filename
    my $upstream_length=50000;        # use the actual sequence length
    my $verbose_output=-1;            # use this verbose mode in call to scan file
    my $wm=$_[0];                     # weight matrix (text, \t and \n delimited)
    my ($i,$j,$k);                    # loop indices
    my ($system_status,$msg);         # warning, used in call to my_rename
    my @scan_output;                  # output of the scan program

    if ($verbose_output>=0) {
	$scanlog_filename=$scan_output_filename;
	$scanlog_filename=~s/txt/log.txt/;
    }
    $strand_filename=$scan_output_filename;
    $strand_filename=~s/txt/strand.txt/;

    # column in the @n variable where the number of transcripts with occurrences is stored
    if ($species eq $sp1) {
	$col_index=3;
    } else {
	$col_index=4;
    }

    if ($verbose eq 1) {
	print log_file "\tscan_motif\n\t\tupstream_filename=$upstream_filename\n\t\tscan_output_filename=$scan_output_filename\n\t\tmotif_index=$motif_index\n\t\tcol_index=$col_index\n";
    }

    if (-e ${scan_output_filename}) {
	# cool, no need to scan
	if ($get_count eq 1) {
	    $n_transcripts_occurrences=$n[$motif_index][$col_index];
	    if ($n_transcripts_occurrences<0) {
		($n_transcripts_occurrences,$n_comments)=line_count($scan_output_filename);
		$n_transcripts_occurrences=$n_transcripts_occurrences-$n_comments;
		$n[$motif_index][$col_index]=$n_transcripts_occurrences;
	    } 
	}
    } else {
	# we do need to scan
	open (temp_file,">${motif_filename}") || die "could not open $motif_filename for writing";
	print temp_file "$wm\n";
	close (temp_file);
	
	$code_exec="${CODE_DIR}/cpp/executables/scan_motif_v${scan_motif_version}.exe ${upstream_filename} ${motif_filename} ${ml} ${thres} ${upstream_length} ${rc} ${sortpos} ${verbose_output}";
	if ($verbose eq 1) {
	    print log_file "\t\tcode_exec=$code_exec\n";
	}
	# output is stored in scan_motif_v${prog_version}.out.txt and scan_motif_v${prog_version}.log.txt
	@scan_output=`$code_exec`;
	($system_status,$n_matches)=extract_text_from_array("system_status",@scan_output);
	if ($system_status ne 0) {
	    print "WARNING!!!\n\tcode_exec=$code_exec\n\tsystem_status=$system_status\n";
	    $n_errors++;
	    print err_file "scan_motif\t${scan_output_filename}\t$code_exec\tsystem_status=$system_status\n";
	    $n_transcripts_occurrences=-1;
	}	
	($system_status,$msg)=my_rename("scan_motif_v${scan_motif_version}.out.txt",$scan_output_filename,0,1);
	if ($verbose_output eq 1) {
	    print log_file "$msg\n";
	}
	($system_status,$msg)=my_rename("scan_motif_v${scan_motif_version}.strand.txt",$strand_filename,0,1);
	if ($verbose_output eq 1) {
	    print log_file "$msg\n";
	}
	if ($verbose_output eq 1) {
	    ($system_status,$msg)=my_rename("scan_motif_v${scan_motif_version}.log.txt",$scanlog_filename,0,1);
	    if ($verbose_output eq 1) {
		print log_file "$msg\n";
	    }
	}

	if ($get_count eq 1) {
	    ($n_transcripts_occurrences,$n_matches)=extract_text_from_array("n_transcripts_hits",@scan_output);
	    $n[$motif_index][$col_index]=$n_transcripts_occurrences;
	}
    }

    if ($verbose eq 1) {
	print log_file "\t\tcount=$n_transcripts_occurrences transcripts\n";
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    if ($verbose_output<=0) {            	# remove temporary files
	$n_unlinked=unlink($motif_filename,"scan_motif_v${scan_motif_version}.log.txt");
    }

    return $n_transcripts_occurrences;
}

sub scan_all{
    # call scan_1m_allseq_v6.exe
    # usage:
    # $n_transcripts_occurrences=scan_all($wm,$scan_output_filename,$motif_index,$species,$ml,$thres,$get_count,$singlefile,$seq_source,$read_all);
    #
    # scan the occurrences of motif given by weight matrix $wm in all the sequences for a given species
    # note that in the current version (v6) we directly use the output from the c program
    # if get_count is 1 then we do a word count and return the number of lines, otherwise, return -1 if there was an error and 0 if everything was ok
    #
    # scan_1m_allseq_v6 usage:
    # scan_1m_allseq_v6.exe <motif_filename> <species> <motif_threshold> <motif_length> ( <max_upstream_length> <max_exon_length> <max_intron_length> <locuslink_only> <rc> <sortpos> <min_dist_factor> <singlefile> <seq_source> <read_all> <verbose> )
    #
    # global variables: @n, $CODE_DIR, $max_upstream_length, $max_exon_length, $max_intron_length, $locuslink_only, $rc, $sortpos, $min_dist_factor, $verbose, $file_label, log_file, $sp1
    #
    # 09-04-2003: changed name of temporary weight matrix file to cooc_checkresults_wm.txt
    # 12-28-2003: call scan_1m_allseq_v4.exe 
    # 02-26-2004: call scan_1m_allseq_v6.exe
    # 02-26-2004: add file_label to the name of the temporary file with weight matrix
    # 03-08-2004: add exp_id to the name of the temporary file with weight matrix

    my $code_exec;                  # code to execute
    my $col_index;                  # in the current version, we fill in the column depending on the species
    my $ml;                         # motif length
    my $motif_index;                # this is used only to check and/or fill in the global variable @n
    my $scan_output_filename;       # where the scan results will be stored 
    my $system_status;              # system status (0 = ok)
    my $thres;                      # threshold
    my $wm;                         # weight matrix 
    my $get_count;                  # if get_count is 1 then we do a word count and return the number of lines 
    my $scanlog_filename;           # log filename
    my $strand_filename;            # strand filename
    my $n_transcripts_occurrences=-1;
    my $motif_filename="${TEMP_DIR}/cooc_checkresults_wm.${exp_id}.${file_label}.txt";
    my $species;                    # hs/mm/dm/sc
    my $verbose_output=-1;          # >0 for verbose output
    my @scan_output;                # used for the output of the scan program
    my $t0=time();                  # used to compute the elapsed time
    my $et;                         # used to compute the elapsed time
    my $msg;                        # used for my_rename and run_system
    my $c_code_version=6;           # version of the c code to call
    my $seq_source;                 # e.g. ensembl, gbk
    my $read_all;                   # 1 to read all lines including those with 0 length 
    my $singlefile;                 # 1 if all the sequences are in a single file (e.g. yeast) and 0 otherwise (e.g. hs and mm)

    ($wm,$scan_output_filename,$motif_index,$species,$ml,$thres,$get_count,$singlefile,$seq_source,$read_all)=@_;
    if ($verbose_output>=0) {
	$scanlog_filename=$scan_output_filename;
	$scanlog_filename=~s/txt/log.txt/;
    }
    $strand_filename=$scan_output_filename;
    $strand_filename=~s/txt/strand.txt/;

    # column in the @n variable where the number of transcripts with occurrences is stored
    if ($species eq ${sp1}) {
	$col_index=1;
    } else {
	$col_index=2;
    }

    if ($verbose eq 1) {
	print log_file "\tscan_motif\n\t\tscan_output_filename=$scan_output_filename\n\t\tmotif_index=$motif_index\n\t\tcol_index=$col_index\n\t\tml=$ml\n\t\tspecies=$species\n";
    }

    if (-e ${scan_output_filename}) {
	# cool, no need to scan
	if ($get_count eq 1) {
	    $n_transcripts_occurrences=$n[$motif_index][$col_index];
	    if ($n_transcripts_occurrences<0) {
		($n_transcripts_occurrences,$n_comments)=line_count($scan_output_filename);
		$n_transcripts_occurrences=$n_transcripts_occurrences-$n_comments;
		$n[$motif_index][$col_index]=$n_transcripts_occurrences;
	    } 
	}
    } else {
	# we do need to scan
	if ($verbose eq 1) {
	    print log_file "\t\twe need to scan. $scan_output_filename does not exist\n";
	}
	print "\t\twe need to scan. $scan_output_filename does not exist\n";

	# print out weight matrix
	open (temp_file,">${motif_filename}") || die "could not open $motif_filename for writing";
	print temp_file "$wm\n";
	close (temp_file);
	
	if ($species eq "dm") {
	    $seq_source="ensembl";
	    $read_all=1;
	}
	if ($species eq "sc") {
	    $singlefile=1;
	    $read_all=1;
	}

	# scan_1m_allseq_v6.exe <motif_filename> <species> <motif_threshold> <motif_length> ( <max_upstream_length> <max_exon_length> <max_intron_length> <locuslink_only> <rc> <sortpos> <min_dist_factor> <singlefile> <seq_source> <read_all> <verbose> )
	$code_exec="${CODE_DIR}/cpp/executables/scan_1m_allseq_v${c_code_version}.exe ${motif_filename} ${species} ${thres} ${ml} ${max_upstream_length} ${max_exon_length} ${max_intron_length} ${locuslink_only} ${rc} ${sortpos} ${min_dist_factor} ${singlefile} ${seq_source} ${read_all} ${verbose_output}";
	#$code_exec="${CODE_DIR}/cpp/executables/scan_1m_allseq.exe ${motif_filename} ${species} ${thres} ${ml} $max_upstream_length $max_exon_length $max_intron_length $locuslink_only $rc $sortpos $min_dist_factor $verbose_output";

	if ($verbose eq 1) {
	    print log_file "\t\tcode_exec=$code_exec\n";
	}
	print "\t\tcode_exec=$code_exec\n";
	# output is stored in scan_1m_allseq_v?.out.txt, scan_1m_allseq_v?.strand.txt (and, if verbose_output>=0) in scan_1m_allseq_v?.log.txt)
	@scan_output=`$code_exec`;
	($system_status,$n_matches)=extract_text_from_array("system_status",@scan_output);
	if ($system_status ne 0) {
	    print "WARNING!!!\n\tcode_exec=$code_exec\n\tsystem_status=$system_status\n";
	    $n_errors++;
	    print err_file "scan_all\t${scan_output_filename}\t$code_exec\tsystem_status=$system_status\n";
	    $n_transcripts_occurrences=-1;
	} else {
	    ($n_transcripts_occurrences,$n_matches)=extract_text_from_array("n_transcripts_hits",@scan_output);
	    ($system_status,$msg)=my_rename("scan_1m_allseq_v${c_code_version}.out.txt",$scan_output_filename,0,1);
	    if ($verbose_output>=0) {
		print log_file "\t\t$msg\n";
	    }
	    ($system_status,$msg)=my_rename("scan_1m_allseq_v${c_code_version}.strand.txt",$strand_filename,0,1);
	    if ($verbose_output>=0) {
		print log_file "\t\t$msg\n";
	    }
	    if ($verbose_output>=1) {
		($system_status,$msg)=my_rename("scan_1m_allseq_v${c_code_version}.log.txt",$scanlog_filename,0,1);
		if ($verbose_output>=0) {
		    print log_file "\t\t$msg\n";
		}
	    }
	    $n[$motif_index][$col_index]=$n_transcripts_occurrences;
	}
    }

    if ($verbose eq 1) {
	print log_file "\t\tcount=$n_transcripts_occurrences\n";
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    if ($verbose_output<=0) {           	# remove temporary files
	$n_unlinked=unlink($motif_filename,"scan_1m_allseq_v${scan_motif_version}.log.txt");
    }

    return $n_transcripts_occurrences;
}

sub cooc4m{
    # usage:
    # ($errorcode,$n_cooc,$n_cooc_tr)=cooc4m($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename,$max_ul,$max_el,$max_il,$transcript_map,$sortpos,$pos_format);
    # call cooc_single_2/3/4m_v3
    # call in verbose=0 mode (we will use the log.txt output from the cooc_single c code)
    # cooc_single_4m_v2 <file1> <file2> <file3> <file4> <ml_1> <ml_2> <ml_3> <ml_4> <n_transcripts> (<max_dist> <min_dist_factor> <verbose> <report_all> <max_n_transcripts>)
    # cooc_single_4m_v3 <filename1> <filename2> <filename3> <filename4> <ml_1> <ml_2> <ml_3> <ml_4> <n_transcripts> (<max_dist> <min_dist_factor> <verbose> <report_all> <order_constraint> <max_n_transcripts> <neg2pos> <sortpos>)
    # cooc_single_1234_v1.exe
    # cooc_single_1234m_v1,exe <max_n_motifs> <filename1> ( <filename2> <filename3> <filename4> ) <ml_1> ( <ml_2> <ml_3> <ml_4> ) <n_transcripts> <transcript_lengths_filename> ( <max_ul> <max_el> <max_il> <max_dist> <min_dist_factor> <report_all> <order_constraint> <max_n_transcripts> <neg2pos> <sortpos> <transcript_map_filename> <verbose>)
    # cooc_single_1234m_v3 <max_n_motifs> <filename1_pos> <filename1_strand> ( <filename2_pos> <filename2_strand> <filename3_pos> <filename3_strand> <filename4_pos> <filename4_strand>) <ml_1> ( <ml_2> <ml_3> <ml_4> ) <n_transcripts> <transcript_lengths_filename> ( <max_ul> <max_el> <max_il> <max_dist> <min_dist_factor> <report_all> <order_constraint> <max_n_transcripts> <sortpos> <transcript_map_filename> <pos_format> <filteron> <verbose>)
    #
    # 01-31-2004 added pos_format as input (+1 for positions with respect to TSS but with upstream > 0 and -1 for positions with respect to 3' end)
    # 01-31-2004 call cooc_single_1234_v3.exe (measures positions with respect to TSS)
    # 01-28-2004 call cooc_single_1234_v2.exe (incorporates strand information)
    # 11-10-2003 call cooc_single_1234_v1.exe
    # 11-10-2003 added report_all, transcript_lengths_filename, max_ul, max_el, max_il, transcript_map, sortpos as input
    # 08-27-2003 always add the max_n_transcripts variable, if it < 0, then the cooc_single program will ignore it

    my ($ml1,$ml2,$ml3,$ml4);
    my ($motifscan_filename1,$motifscan_filename2,$motifscan_filename3,$motifscan_filename4,$motifstrand_filename1,$motifstrand_filename2,$motifstrand_filename3,$motifstrand_filename4);
    my $max_n_transcripts;
    my $n_transcripts;
    my $order_constraint;
    my $errorcode;
    my $n_cooc;
    my $n_cooc_tr;
    my $code_exec;
    my $search_pattern;
    my @cooc_output;
    my $i;
    my $temp;
    my $temp1;
    my $temp2;
    my @splitarray;
    my $cooc_c_version=5;                    # version of the c code to call
    my $t0=time();
    my $et;
    my $system_status;
    my $report_all;
    my $transcript_lengths_filename;
    my $max_ul;
    my $max_el;
    my $max_il;
    my $transcript_map_filename;
    my $sortpos;
    my $function_verbose=1;                  # determines whether this subroutine is verbose or not
    my $cooc_verbose=1;                      # determines whether the call to cooc_single_1234m_v1.exe is verbose or not
    my $filteron=1;

    ($ml1,$ml2,$ml3,$ml4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$max_n_transcripts,$n_transcripts,$order_constraint,$report_all,$transcript_lengths_filename,$max_ul,$max_el,$max_il,$transcript_map_filename,$sortpos)=@_;

    if ($function_verbose eq 1) {
	print log_file "\tcooc4m\nmotifscan_filename1=$motifscan_filename1\tmotifstrand_filename1=$motifstrand_filename1\nmotifscan_filename2=$motifscan_filename2\tmotifstrand_filename2=$motifstrand_filename2\nmotifscan_filename3=$motifscan_filename3\tmotifstrand_filename3=${motifstrand_filename3}\nmotifscan_filename4=$motifscan_filename4\tmotifstrand_filename4=$motifstrand_filename4\n";
	print log_file "\t\tml1=$ml1\tml2=$ml2\tml3=$ml3\tml4=$ml4\n";
	print log_file "\t\tmax_n_transcripts=$max_n_transcripts\n";
	print log_file "\t\tn_transcripts=$n_transcripts\n";
	print log_file "\t\torder_constraint=$order_constraint\n";
	print log_file "\t\treport_all=$report_all\ttranscript_lengths_filename=$transcript_lengths_filename\n";
	print log_file "\t\tmax_ul=$max_ul\tmax_el=$max_el\tmax_il=$max_il\n";
	print log_file "\t\ttranscript_map_filename=$transcript_map_filename\n";
	print log_file "\t\tsortpos=$sortpos\n";
    }


    # cooc_single_1234m_v1,exe <max_n_motifs> <filename1> ( <filename2> <filename3> <filename4> ) <ml_1> ( <ml_2> <ml_3> <ml_4> ) <n_transcripts> <transcript_lengths_filename> ( <max_ul> <max_el> <max_il> <max_dist> <min_dist_factor> <report_all> <order_constraint> <max_n_transcripts> <neg2pos> <sortpos> <transcript_map_filename> <verbose>)
    # cooc_single_1234m_v2 <max_n_motifs> <filename1_pos> <filename1_strand> ( <filename2_pos> <filename2_strand> <filename3_pos> <filename3_strand> <filename4_pos> <filename4_strand>) <ml_1> ( <ml_2> <ml_3> <ml_4> ) <n_transcripts> <transcript_lengths_filename> ( <max_ul> <max_el> <max_il> <max_dist> <min_dist_factor> <report_all> <order_constraint> <max_n_transcripts> <neg2pos> <sortpos> <transcript_map_filename> <verbose>)
    
    if ( $ml3 < 0 ) {	    # 2 motifs
	$code_exec="${CODE_DIR}/cpp/executables/cooc_single_1234m_v${cooc_c_version}.exe 2 ${motifscan_filename1} ${motifstrand_filename1} ${motifscan_filename2} ${motifstrand_filename2} ${ml1} ${ml2} ${n_transcripts} ${transcript_lengths_filename} ${max_ul} ${max_el} ${max_il} ${max_dist} ${min_dist_factor} ${report_all} ${order_constraint} ${max_n_transcripts} ${sortpos} ${transcript_map_filename} ${pos_format} ${filteron} ${cooc_verbose}";
    } else {
	if ($ml4<0) {	    # 3 motifs 
	    $code_exec="${CODE_DIR}/cpp/executables/cooc_single_1234m_v${cooc_c_version}.exe 3 ${motifscan_filename1} ${motifstrand_filename1} ${motifscan_filename2} ${motifstrand_filename2} ${motifscan_filename3} ${motifstrand_filename3} ${ml1} ${ml2} ${ml3} ${n_transcripts} ${transcript_lengths_filename} ${max_ul} ${max_el} ${max_il} ${max_dist} ${min_dist_factor} ${report_all} ${order_constraint} ${max_n_transcripts} ${sortpos} ${transcript_map_filename} ${pos_format} ${filteron} ${cooc_verbose}";
	} else {	    # 4 motifs 
	    $code_exec="${CODE_DIR}/cpp/executables/cooc_single_1234m_v${cooc_c_version}.exe 4 ${motifscan_filename1} ${motifstrand_filename1} ${motifscan_filename2} ${motifstrand_filename2} ${motifscan_filename3} ${motifstrand_filename3} ${motifscan_filename4} ${motifstrand_filename4} ${ml1} ${ml2} ${ml3} ${ml4} ${n_transcripts} ${transcript_lengths_filename} ${max_ul} ${max_el} ${max_il} ${max_dist} ${min_dist_factor} ${report_all} ${order_constraint} ${max_n_transcripts} ${sortpos} ${transcript_map_filename} ${pos_format} ${filteron} ${cooc_verbose}";
	}
    }
    if ($function_verbose eq 1) {
	print log_file "\t\tcode_exec=\t$code_exec\n";
    }
    @cooc_output=`$code_exec`;
    
    $errorcode=-1;
    $n_cooc=-1;
    $n_cooc_tr=-1;
    
    ($n_cooc_tr,$n_matches)=extract_text_from_array("n_tr_cooc",@cooc_output);
    if ($n_matches ne 1) {
	print "\t\tWARNING! n_matches=$n_matches in search for n_tr_cooc\n";
    }
    ($n_cooc,$n_matches)=extract_text_from_array("n_cooc",@cooc_output);
    if ($n_matches ne 1) {
	print "\t\tWARNING! n_matches=$n_matches in search for n_cooc\n";
    }
    ($errorcode,$n_matches)=extract_text_from_array("errorcode",@cooc_output);
    if ($n_matches ne 1) {
	print "\t\tWARNING! n_matches=$n_matches in search for errorcode\n";
    }
    ($system_status,$n_matches)=extract_text_from_array("system_status",@cooc_output);
    if ($n_matches ne 1) {
	print "\t\tWARNING! n_matches=$n_matches in search for system_status\n";
    }
    if ($system_status ne 0) {
	print "\n\nERROR!!!\n\ncode_exec=\n$code_exec\nsystem_status=$system_status\n";
	exit;
    }
    if ($verbose eq 1) {
	print log_file "\t\tn_cooc=$n_cooc\tn_cooc_tr=$n_cooc_tr\terrorcode=$errorcode\n";
    }

    chomp $n_cooc;
    chomp $n_cooc_tr;
    chomp $errorcode;

    if ($verbose eq 1) {
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    return ($errorcode,$n_cooc,$n_cooc_tr);
}

sub computebinom{
    # usage:
    # ($p_binom,$p)=computebinom($n_background,$n_cooc_background,$n_cluster,$n_cooc_cluster);
    # computes the probability of obtaining $n_cooc_cluster hits out of $n_cluster trials assuming a probability given by $n_cooc_background out of $n_background
    my $n_background=$_[0];
    my $n_cooc_background=$_[1];
    my $n_cluster=$_[2];
    my $n_cooc_cluster=$_[3];
    my $p;
    my $temp;
    my $code_exec;
    my @binocdf_output;
    my @splitarray;
    my $p_binom=-1;

    if ($verbose eq 1) {
	print log_file "\tcomputebinom $n_background $n_cooc_background $n_cluster $n_cooc_cluster\n";
    }
    if ($n_cooc_background>-1) {
	if ($n_cooc_background>=10) {
	    $p=(${n_cooc_background}-${n_cooc_cluster})/${n_background};     # because we assume here that the cluster transcripts are included in the background set
	} else {
	    $p=$n_cooc_background/$n_background;
	}
	if ( ($p<=0) & ($n_cooc_cluster>0) ) {		
	    print log_file "error! n_cooc_background=${n_cooc_background} and n_cooc_cluster=${n_cooc_cluster}\n";
	    print log_file "it should be n_cooc_cluster<=n_cooc_background\n";
	    $p=1.0/${n_background};
	}
	if ($verbose eq 1) {
	    print log_file "\t\tn_cooc_background=$n_cooc_background\n";
	    print log_file "\t\tp=$p\n";
	}
	
	#$temp=$n_cooc_cluster-1;  # no need to do this, this is now taken care within the binocdf.exe program
	#$code_exec="${CODE_DIR}/cpp/code/misc/binocdf ${n_cluster} ${temp} ${p}";
	$code_exec="${CODE_DIR}/cpp/executables/binocdf.exe ${n_cluster} ${n_cooc_cluster} ${p}";
	if ($verbose eq 1) {
	    print log_file "\t\tcode_exec=$code_exec\n";
	}
	@binocdf_output=`$code_exec`;
	#if ($verbose eq 1) {
	#    print log_file "binocdf_output=@binocdf_output";
	#}
	$temp=$binocdf_output[4];
	chomp $temp;
	@splitarray=split /\=/,$temp;
	$p_binom=$splitarray[2];
	chomp $p_binom;
	if ($p_binom<0) {
	    $p_binom=0.00000000000000001;
	}
    }

    return ($p_binom,$p);
}

sub computebinom_dist{
    # given a set of $n_ortho orthologous transcripts, we are interested here in estimating the probability that $n_close of them occur close to the corresponding 
    # orthologous genes. for one single pair, the probability of them being close is $p_ortho_distance
    # usage:
    # ($p_binom)=computebinom($n_ortho,$n_close,$p_ortho_distance);

    my $n_ortho=$_[0];                 # total number of orthologous transcripts
    my $n_close=$_[1];                 # number of mm transcripts that are close to the hs orthologous ones      
    my $p_ortho_distance=$_[2];        # probability cut-off that a given single mm gene is close to its orthologous gene
    my $n_cooc_cluster=$_[3];
    my $temp;
    my $code_exec;
    my @binocdf_output;
    my @splitarray;
    my $p_binom=-1;
    my $t0=time();
    my $et;

    # take care of special cases
    if ($n_close eq 0) {
	$p_binom=1.0;
    }
    if ($n_ortho eq 0) {
	$p_binom=1.0;
    }
    if ($p_ortho_distance eq 0) {
	$p_binom=0.0;
    }
    if ($p_ortho_distance eq 1.0) {
	$p_binom=1;
    }

    if ($verbose eq 1) {
	print log_file "\tcomputebinom_dist\n";
    }
    
    if ($p_binom < 0) {
	# none of the special cases apply
	$code_exec="${CODE_DIR}/cpp/executables/binocdf.exe ${n_ortho} ${n_close} ${p_ortho_distance}";
	if ($verbose eq 1) {
	    print log_file "\t\t$code_exec\n";
	}
	@binocdf_output=`$code_exec`;
	#if ($verbose eq 1) {
	#    print log_file "binocdf_output=@binocdf_output";
	#}
	$temp=$binocdf_output[4];
	chomp $temp;
	@splitarray=split /\=/,$temp;
	$p_binom=$splitarray[2];
	chomp $p_binom;
	if ($p_binom<0) {
	    $p_binom=0.00000000000000001;
	}
    }

    if ($verbose eq 1) {
	print log_file "\t\tp_binom=$p_binom\n";
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    return ($p_binom);
}

sub get_max_n_transcripts{
    # get the max_n_transcripts that one needs to run for cooc_single
    # usage: ($max_n_transcripts)=get_max_n_transcripts($k,$n,$y,$n_all);
    # where 
    #    $n is the total number of transcripts in the cluster
    #    $k is the number of transcripts containing the module within the cluster 
    #    $y is the threshold probability
    #    $n_all is the total number of transcripts in the genome
    # returns the maximum number of transcripts in the genome that can contain the module such that the binomial probability does not exceed $y
    # 12-29-2003: clean up the code, use CODE_DIR global variable

    my $k=$_[0];
    my $n=$_[1];
    my $y=$_[2];
    my $n_all=$_[3];
    my $r=0.01;
    my $code_exec;
    my $max_n_transcripts;
    my @code_output;
    my $i=0;
    my $temp;
    my @splitarray;
    my $p;

    # invbinocdf2 <n> <k> <y> (<r>)

    $code_exec="${CODE_DIR}/${code_exec}/cpp/executables/invbinocdf2.exe ${n} ${k} ${y} ${r}";
    @code_output=`$code_exec`;

    $max_n_transcripts=$n_all;
    while ($i<=$#code_output) {
	$temp=$code_output[$i];
	chomp $temp;
	if ($temp =~ /^p/) {
	    @splitarray=split /\=/,$temp;
	    $p=$splitarray[1];
	    $i=$#code_output+1;
	    if ($p<=0) {
		$p=$r;
	    }
	    $max_n_transcripts=ceil($p * $n_all);
	}
	$i++;
    }

    return $max_n_transcripts;
}

sub map_orthologous{
    # map the orthologous genes in the upstream regions with cl/all genes from one species to the other 
    # usage: ($n_hits,$n_miss,$n_notfound)=map_orthologous($species1,$ind,$clorall,@mm2hs_map);
    # where mm2hs map file contains in each entry the following information:  index_of_mm_gene  index_of_hs_gene
    #       species1 indicates the species to be mapped
    #       ind      indicates whether the species to be mapped corresponds to the first (ind=1) or second (ind=2) column in the map
    #
    # 09-04-2003: we are not calling the remove_redundancies program, therefore, we do not use here the 'nored' files
    # 09-05-2003: we now map the other way around. we convert the hs gene numbers to mm gene numbers. 
    # 09-05-2003: we use a matrix as map
    # 09-05-2003: we change a few file names
    # 12-29-2003: back to making it generic for cl/all
    # 12-29-2003: use species1 and ind to allow mapping in any direction

    my $clorall;                                         # clorall = 'cl' or 'all'
    my $species1;                                        # first species
    my $ind;                                             # whether species1 corresponds to the first (ind=1) or second (ind=2) entry in the map file
    my $ind2;                                            # index for the other species
    my @map;                                             # map in the format index1 index2 where 1 and 2 correspond to the two species (ind says which is which)
    my $record;
    my $sp2_gene;                                        # gene in species 2
    my $sp1_gene;                                        # gene in species 1
    my $n_hits=0;
    my $n_miss=0;
    my $n_notfound=0;
    my $deja_vu="";
    my $t0=time();
    my $et;
    my $previous_gene_number=-1;
    my $renamed_files;
    my $original_filename;                                # original file with cooc_all_hs information
    my $tempmap_filename;                                 # temporary file mapping the hs gene numbers to mm gene numbers
    my $sorted_filename;                                  # mapped file with the results from cooc_all_hs but with orthologous mm gene numbers instead of hs 
    my ($system_status,$msg);                             # used in my_rename

    ($species1,$ind,$clorall,@map)=@_;
    print "species1=$species1\n";
    print "ind=$ind\n";
    print "clorall=$clorall\n";
    $ind2=3-$ind;
    $original_filename="cooc_${clorall}_${species1}.txt";
    $tempmap_filename="cooc_${clorall}_${species1}.tempmap.txt";
    $sorted_filename="cooc_${clorall}_${species1}.mapped.txt";

    if ($verbose eq 1) {
	print log_file "\t\tmapping ${species1} gene numbers in ${original_filename} to orthologous gene numbers...\n";
    }

    open (cooc_sp1_tempmap_file,">${tempmap_filename}") || die "could not open $tempmap_filename for output";
    open (cooc_sp1_orig_file,$original_filename) || die "could not open $original_filename for reading";

    while ($record=<cooc_sp1_orig_file>) {
	if ($record =~ /^\%/) {
	    print cooc_sp1_tempmap_file "$record";
	} else {
	    chomp $record;
	    $sp1_gene=get_column($record,0);          # gene number species 1
	    $sp2_gene=$map[$sp1_gene][$ind2];         # orthologous gene number

	    if ($previous_gene_number eq $sp1_gene) {
	    } else {
		$previous_gene_number=$sp1_gene;
		if ($sp2_gene > -1) {
		    $n_hits++;
		} else {
		    $n_miss++;
		}
	    }

	    $record=~s/^${sp1_gene}/${sp2_gene}/;
	    print cooc_sp1_tempmap_file "${record}\n";
	}
    }
    close (cooc_sp1_orig_file);
    close (cooc_sp1_tempmap_file);

    # sort file
    if ($verbose eq 1) {
	print log_file "\t\tsorting $tempmap_filename\n";
    }
    sort_file_num($tempmap_filename,1);
    ($system_status,$msg)=my_rename("cooc_${clorall}_${species1}.tempmap.sortedbycol1.txt",$sorted_filename,0,1);

    if ($verbose eq 1) {
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    return ($n_hits,$n_miss,$n_notfound);
}

# sub map_orthologous_cl{
#     # map the orthologous genes in the upstream regions of coexpressed genes from mouse to human 
#     # usage: ($n_hits,$n_miss,$n_notfound)=map_orthologous_cl(@mm2hs_map_cl);
#     # where the mm2hs_cl map file contains in each entry the following information:
#     #   index_of_mm_gene  index_of_hs_gene (for orthologous pairs within the coexpressed genes
#     # 09-04-2003: we are not calling the remove_redundancies program, therefore, we do not use here the 'nored' files
#     # 09-05-2003: we now map the other way around. we convert the hs gene numbers to mm gene numbers. 
#     # 09-05-2003: we use a matrix as map
#     # 09-05-2003: we change a few file names
#     # 09-29-2003: corrected bug so that we now close cooc_cl_hs_orig_file (and not cooc_cl_mm_orig_file!)

#     my @map=@_;
#     my $record;
#     my $mouse_gene;
#     my $human_gene;
#     my $n_notfound=0; # not found in map
#     my $n_miss=0;     # found in map but orthologous is -1
#     my $n_hits=0;     # found in map and orthologous not -1
#     #my $deja_vu="none";
#     my $previous_gene=0;
#     my $renamed_files;
#     my $original_filename="cooc_cl_hs.txt";             # original file with cooc_cl_hs information
#     my $tempmap_filename="cooc_cl_hs.tempmap.txt";      # temporary file mapping the hs gene numbers to mm gene numbers
#     my $tempmap_filename2=">cooc_cl_hs.tempmap.txt";    # temporary file mapping the hs gene numbers to mm gene numbers
#     my $sorted_filename="cooc_cl_hs.mapped.txt";        # mapped file with the results from cooc_cl_hs but with orthologous mm gene numbers instead of hs 

#     my $t0=time();
#     my $et;

#     if ($verbose eq 1) {
# 	print log_file "\t\tmapping hs gene numbers in cooc_cl_hs.txt to orthologous mm gene numbers...\n";
#     }

#     open (cooc_cl_hs_tempmap_file,$tempmap_filename2) || die "could not open $tempmap_filename2 for output";
#     open (cooc_cl_hs_orig_file,$original_filename) || die "could not open $original_filename for reading";

#     while ($record=<cooc_cl_hs_orig_file>) {
# 	if ($record =~ /^\%/) {
# 	    print cooc_cl_hs_tempmap_file "$record";
# 	} else {
# 	    chomp $record;
# 	    $human_gene=get_column($record,0);
# 	    $mouse_gene=$map[$human_gene][2];     # orthologous mm gene number

# 	    if ($previous_gene_number eq $human_gene) {
# 	    } else {
# 		$previous_gene_number=$human_gene;
# 		if ($mouse_gene > -1) {
# 		    $n_hits++;
# 		} else {
# 		    $n_miss++;
# 		}
# 	    }

# 	    $record=~s/^${human_gene}/${mouse_gene}/;
# 	    print cooc_cl_hs_tempmap_file "${record}\n";
# 	}
#     }

#     close (cooc_cl_hs_tempmap_file);
#     close (cooc_cl_hs_orig_file);
#     #printf("ready\n");

#     # sort file
#     if ($verbose eq 1) {
# 	print log_file "\t\tsorting cooc_cl_hs.tempmap.txt by orthologous mm transcript...\n";
#     }
#     #sort_file_num("cooc_cl_mm.nored.txt",1);
#     sort_file_num($tempmap_filename,1);
#     #rename("cooc_cl_mm.nored.txt","cooc_cl_mm.nored.bak2");
#     #rename($cooc_filename,$cooc_filename_bak2);
#     #rename("cooc_cl_mm.nored.sortedbycol1.txt","cooc_cl_mm.nored.txt");
#     #rename("cooc_cl_mm.sortedbycol1.txt",$cooc_filename);
#     $renamed_files=rename("cooc_cl_hs.tempmap.sortedbycol1.txt",$sorted_filename);
#     if ($renamed_files ne 1) {
# 	print "ERROR!!! renaming cooc_cl_hs.tempmap.sortedbycol1.txt to $sorted_filename, renamed_files = $renamed_files\n";
#     }

#     if ($verbose eq 1) {
# 	$et=time();
# 	$et=$et-$t0;
# 	print log_file "\t\tcpu time=$et\n";
#     }

#     return ($n_hits,$n_miss,$n_notfound);
# }

sub count_genes_withhomologous{
    # count the number of genes with homologous genes
    # usage: ($n_hits,$n_miss,$n_notfound)=count_genes_withhomologous($species1,$ind,$clorall,@map);
    # where map is the map from sp1 to sp2 genes
    #       species1 indicates the species to be mapped
    #       ind      indicates whether the species to be mapped corresponds to the first (ind=1) or second (ind=2) column in the map
    # 
    # 09-08-2003: convert to count_mm_withhomologous, the count of hs genes with homologous is done within the map_orthologous subroutines
    # 09-04-2003: there is no call to remove redundancies program, therefore, remove the 'nored' from the file names
    # 09-05-2003: map is now a matrix with two columns, one for each species, entries are assumed to be in order and exhaustive in the first column (hs)
    # 12-29-2003: allow to move in any direction

    my @map;
    my $record;
    my $sp1_gene;
    my $sp2_gene;
    my $n_hits=0;             # human gene found in orthology map and there is an orthologous gene
    my $n_miss=0;             # human gene found in orthology map but there is no orthologous 
    my $n_notfound=0;         # human gene not found in the orthology map
    my $clorall;              # 'cl' or 'all'
    my $t0=time();
    my $et;
    my $previous_gene=-1;
    my $count_filename;       # file with entries to count
    my $species1;
    my $ind;
    my $ind2;

    ($species1,$ind,$clorall,@map)=@_;
    $ind2=3-$ind;

    $count_filename="cooc_${clorall}_${species1}.txt";
    #if ($allorcl eq 1) {
    #$count_mm_filename="cooc_all_mm.txt";
    #} else {
    #$count_mm_filename="cooc_cl_mm.txt";
    #}

    open (count_file,$count_filename) || die "could not open $count_filename for reading";
    if ($verbose eq 1) {
	print log_file "\tcounting ${species1} genes with homologous genes, clorall=$clorall, count_filename=$count_filename\n";
    }

    while ($record=<count_file>) {
	if ($record !~ /^\%/) {	    # if it is not a comment line
	    chomp $record;
	    $sp1_gene=get_column($record,0);
	    if ($sp1_gene ne $previous_gene) {
		$previous_gene=$sp1_gene;
		$sp2_gene=$map[$sp1_gene][$ind2];
		if ($sp2_gene > -1) {
		    $n_hits++;
		} else {
		    $n_miss++;
		}
	    }
	}
    }    
    close (count_file);

    if ($verbose eq 1) {
	print log_file "\t\tn_hits=$n_hits\n";
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    return ($n_hits,$n_miss,$n_notfound);
}

sub add_genename{
    # given an output file with the following format
    # gene-number position(s)
    # return a text variable with the following format
    # gene-name gene-number position(s)
    # (plus the corresponding comments)
    # usage:
    # ($comments,$core)=add_genename($input_filename,$info_filename,$info_col);
    # 
    # 09-08-2003: pre-load all the info arrays and use this as input to this file (so that we do not need to re-load the info every time...)

    my $input_filename;
    my $comments="";
    my $core="";
    my @splitarray;
    my @info;
    my $n_cols;
    my $i;
    my $curr_line;
    my $record;
    my $curr_index;
    my $ll;
    my $t0=time();
    my $et;
    my $cleanll;

    ($input_filename,@info)=@_;

    if ($verbose eq 1) {
	#print log_file "\tadd_genename for input=$input_filename, info=$info_filename and col=$info_col\n";
	print log_file "\tadd_genename for input=$input_filename\n";
    }
    
    open (results_file,$input_filename) || die "could not open $input_filename for reading";
    while ($record=<results_file>) {
	chomp $record;
	if (is_header($record)) {
	    $record=~s/\%//;
	    $comments.="$record\t";
	} else {
	    if (length($record)>1) {
		@splitarray=split /\t/,$record;
		$curr_index=$splitarray[0];
		$ll=$info[$curr_index];
		$cleanll=ll_clean($ll);
		$curr_line="${cleanll}\t$curr_index";
		for ($i=1;$i<=$#splitarray;$i++) {
		    $curr_line.="\t$splitarray[$i]";
		}
		$core.="$curr_line\n";
	    }
	}
    }
    close (results_file);

    if ($verbose eq 1) {
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }

    return ($comments,$core);
}

sub compute_random_overlap{
    # call random_overlap.exe
    # usage: ($p)=compute_random_overlap($n1,$n2,$n,$nr,$n_iter);
    # usage of c code: random_overlap <n1> <n2> <n> <n_iter> <nr>
    # given a set of n balls numbered 1 ... n, we extract n1 without replacement and randomly. we write down the numbers and put them back in the bag. then we extract n2 in the same fasion. what is the probability that we will get nr or more of intersection between the two cases? analyze this using n_iter of simulations, we return the total distribution of number of intersecting balls and the probability

    my $n1=$_[0];
    my $n2=$_[1];
    my $n=$_[2];
    my $nr=$_[3];
    my $n_iter=$_[4];
    my $code_exec;
    my $p;
    my @code_output;
    my $t0=time();
    my $et;
    my $call_code_exec=1;        # do not call random_overlap.exe if it is a trivial case by assigning call_code_exec=0

    if ($verbose eq 1) {
	print log_file "\tcompute_random_overlap\tn1=$n1\tn2=$n2\tn=$n\tnr=$nr\n";
    }

    # take care of special cases, like n1 / n2 / nr null
    if ($nr eq 0) {
	$p=1;             # for sure we will get 0 or more overlapping genes
	$call_code_exec=0;
    } else {
	if ($n2 eq 0) {
	    print log_file "\t\tWARNING!!! call to compute_random_overlap, this is strange, n2=$n2 and nr=$nr. please check. assigning p=0\n";
	    $p=0;
	    $call_code_exec=0;
	    $n_errors++;
	    print err_file "compute_random_overlap\t$code_exec\tn2=$n2\n";
	}
	if ($n1 eq 0) {
	    print log_file "\t\tWARNING!!! call to compute_random_overlap, this is strange, n2=$n2 and nr=$nr. please check. assigning p=0\n";
	    $p=0;
	    $call_code_exec=0;
	    $n_errors++;
	    print err_file "compute_random_overlap\t$code_exec\tn1=$n1\n";
	}
    }
    print "\t\tcall_code_exec=$call_code_exec (compute_random_overlap)\n";
    if ($verbose eq 1) {
	print log_file "\t\tcall_code_exec=$call_code_exec\n";
    }
    if ($call_code_exec eq 1) {
	$code_exec="${CODE_DIR}/cpp/executables/random_overlap.exe ${n1} ${n2} ${n} ${n_iter} ${nr}";
	@code_output=`$code_exec`;
	
	if ($verbose eq 1) {
	    print log_file "\t\tcode_exec=$code_exec\n";
	}
	
	$p=$code_output[0];
	chomp $p;
	if ($verbose eq 1) {
	    print log_file "\t\tp=$p\n";
	}
    }

    if ($verbose eq 1) {
	$et=time();
	$et=$et-$t0;
	print log_file "\t\tcpu time=$et\n";
    }
    return ($p);
}

sub get_n_transcripts_all{
    # get the total number of transcripts
    # usage:
    # $n_upstream_all=get_n_transcripts_all($list_filename,$ll_only);

    my $list_filename=$_[0];
    my $ll_only=$_[1];
    my $record;
    my $n=0;

    open (list_file,$list_filename) || die "could not open $list_filename for reading";
    while ($record=<list_file>) {
	if (is_comment($record) eq "true") {
	} else {
	    if ($ll_only eq 1) {
		if ($record =~ /LocusID/) {
		    $n++;
		}
	    } else {
		$n++;
	    }
	}
    }
    close (list_file);
    
    return $n;
}

sub ll_clean{
    # we want the locus link ID to be a number, therefore we clean the text here
    # if it is not a proper locus ID (e.g. interim ID, LOC entries, etc.) we convert to a negative number
    # usage: $cleanll=ll_clean($ll);

    my $ll=$_[0];
    my $cleanll;
    
    $cleanll=$ll;
    if ($cleanll =~ /ll/) {
	$cleanll=~s/ll\.//;
    }
    if ($cleanll =~ /LL/) {
	$cleanll=~s/LL\.//;
    }
    if ($cleanll =~ /LOC/) {
	$cleanll=~s/LOC//;
	$cleanll=-$cleanll;
    }
    if ($cleanll =~ /LocusID\:/) {
	$cleanll=~s/LocusID\://;
    }
    if ($cleanll =~ /InterimID\:/) {
	$cleanll=~s/InterimID\://;
	$cleanll=-$cleanll;
    }

    return $cleanll;
}

sub create_scan_all_filename{
    # generate the motifscan_filename by taking into account the motif index and using either the transfac filename or a temporary scan of the current experiment id
    # usage:
    # ($motifscan_filename,$motifstrand_filename)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$motif_index,$motif_scan_threshold_column,$species,$n_transfac_motifs);
    # 01-28-2004: added strand filename as output 
    
    my $exp_id;                         # experiment id (used in the file name unless it is a TRANSFAC motif
    my $max_upstream_length;            # maximum upstream length
    my $motif_index;                    # motif index
    my $motif_scan_threshold_column;    # motif scan threshold column
    my $species;                        # species (hs/mm)
    my $n_transfac_motifs;              # number of TRANSFAC motifs
    my $motifscan_filename;             # output filename
    my $motifstrand_filename;           # output strand filename
    my $seq_upstream_length;
    my $seq_exon_length;
    my $seq_intron_length;

    ($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$motif_index,$motif_scan_threshold_column,$species,$n_transfac_motifs)=@_;
    
    if ($motif_index > $n_transfac_motifs) {    
	$motifscan_filename="${TEMP_DIR}/scan_${exp_id}/ms_${exp_id}_${species}_all_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${motif_index}_${motif_scan_threshold_column}.txt";
	$motifstrand_filename="${TEMP_DIR}/scan_${exp_id}/ms_${exp_id}_${species}_all_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${motif_index}_${motif_scan_threshold_column}.strand.txt";
    } else {
	$motifscan_filename="${DATA_DIR}/db/transfac/scan/ms_transfac_${species}_all_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${motif_index}_${motif_scan_threshold_column}.txt";
	$motifstrand_filename="${DATA_DIR}/db/transfac/scan/ms_transfac_${species}_all_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${motif_index}_${motif_scan_threshold_column}.strand.txt";
    }
    return ($motifscan_filename,$motifstrand_filename);
}

sub scan_4motifs_cl
{
    # scan_4motifs_cl
    # usage:
    # ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_cl,$logtext)=scan_4motifs_cl($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count);
    # uses global variables: scan_dir, exp_id, seq_upstream_length, seq_exon_length, seq_intron_length, motif_scan_thresholds_column, 
    #                        upstream_filename_cl_sp1,upstream_filename_cl_sp2,sp1,sp2,rc
    # 01-28-2004: add motifstrand_filename as output
    # last modified: 01-28-2004

    my $species;
    my $get_count;
    my $logtext="";                  # output variable
    my ($motifscan_filename1,$motifscan_filename2,$motifscan_filename3,$motifscan_filename4);
    my ($motifstrand_filename1,$motifstrand_filename2,$motifstrand_filename3,$motifstrand_filename4);
    my ($n1,$n2,$n3,$n4);
    my ($i1,$i2,$i3,$i4);
    my ($wm1,$wm2,$wm3,$wm4);
    my ($ml1,$ml2,$ml3,$ml4);
    my ($t1,$t2,$t3,$t4);
    my $compute_cooc_cl=1;           # output variable
    my $upstream_filename;
    
    ($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count)=@_;
    $n1=-1;$n2=-1;$n3=-1;$n4=-1;

    if ($species eq $sp1) {
	$upstream_filename=$upstream_filename_cl_sp1;
    } else {
	$upstream_filename=$upstream_filename_cl_sp2;
    }

    $logtext.="\t-> processing scan, set=cluster, species=${species}\n";	    
    $motifscan_filename1="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i1}_${motif_scan_threshold_column}.txt";
    $motifstrand_filename1="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i1}_${motif_scan_threshold_column}.strand.txt";
    $n1=scan_motif($wm1,$upstream_filename,$motifscan_filename1,$i1,$species,${ml1},${t1},$get_count,$rc);
    if ($n1<=0) {
	$compute_cooc_cl_hs=0;
    } else {
	$motifscan_filename2="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i2}_${motif_scan_threshold_column}.txt";
	$motifstrand_filename2="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i2}_${motif_scan_threshold_column}.strand.txt";
	$n2=scan_motif($wm2,$upstream_filename,$motifscan_filename2,$i2,$species,${ml2},${t2},$get_count,$rc);
	if ($n2<=0) {
	    $compute_cooc_cl=0;
	} else {
	    if ($i3>0) {
		$motifscan_filename3="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i3}_${motif_scan_threshold_column}.txt";
		$motifstrand_filename3="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i3}_${motif_scan_threshold_column}.strand.txt";
		$n3=scan_motif($wm3,$upstream_filename,$motifscan_filename3,$i3,$species,${ml3},${t3},$get_count,$rc);	    
		if ($n3<=0) {
		    $compute_cooc_cl=0;
		} else {
		    if ($i4>0) {
			$motifscan_filename4="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i4}_${motif_scan_threshold_column}.txt";
			$motifstrand_filename4="${scan_dir}/ms_${exp_id}_${species}_cl_${seq_upstream_length}_${seq_exon_length}_${seq_intron_length}_${i4}_${motif_scan_threshold_column}.strand.txt";
			$n4=scan_motif($wm4,$upstream_filename,$motifscan_filename4,$i4,$species,${ml4},${t4},$get_count,$rc);
			if ($n4<=0) {
			    $compute_cooc_cl=0;
			}
		    }
		}
	    }
	}
    }

    #print "motifstrand_filename1=$motifstrand_filename1\tn1=$n1\n";
    #print "motifstrand_filename2=$motifstrand_filename2\tn2=$n2\n";
    #print "motifstrand_filename3=$motifstrand_filename3\tn3=$n3\n";
    #print "motifstrand_filename4=$motifstrand_filename4\tn4=$n4\n";

    $logtext.="n1=$n1\tn2=$n2\tn3=$n3\tn4=$n4\n";
    $logtext.="compute_cooc_cl=$compute_cooc_cl\n";
    return ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_cl,$logtext);
}

sub scan_4motifs_all
{
    # scan_4motifs_all
    # usage:
    # ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_all,$logtext)=scan_4motifs_all($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count,$singlefile,$seq_source,$read_all);
    # uses global variables: scan_dir, exp_id, seq_upstream_length, seq_exon_length, seq_intron_length, motif_scan_thresholds_column, 
    #                        upstream_filename_cl_sp1, rc, n_transfac_motifs
    # 02-28-2004: added singlefile,seq_source and read_all (new parameteres required in scan_1m_allseq_v6.exe

    my $species;
    my $get_count;
    my $logtext="";                  # output variable
    my ($motifscan_filename1,$motifscan_filename2,$motifscan_filename3,$motifscan_filename4);
    my ($motifstrand_filename1,$motifstrand_filename2,$motifstrand_filename3,$motifstrand_filename4);
    my ($n1,$n2,$n3,$n4);
    my ($i1,$i2,$i3,$i4);
    my ($wm1,$wm2,$wm3,$wm4);
    my ($ml1,$ml2,$ml3,$ml4);
    my ($t1,$t2,$t3,$t4);
    my $compute_cooc_all=1;           # output variable
    my $singlefile;
    my $seq_source;
    my $read_all;
    
    ($i1,$i2,$i3,$i4,$wm1,$wm2,$wm3,$wm4,$t1,$t2,$t3,$t4,$ml1,$ml2,$ml3,$ml4,$species,$get_count,$singlefile,$seq_source,$read_all)=@_;
    $n1=-1;$n2=-1;$n3=-1;$n4=-1;

    $logtext.="\t-> processing scan, set=all, species=${species}\n";	    
    ($motifscan_filename1,$motifstrand_filename1)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$i1,$motif_scan_threshold_column,$species,$n_transfac_motifs);
    $n1=scan_all($wm1,$motifscan_filename1,$i1,$species,$ml1,$t1,$get_count,$singlefile,$seq_source,$read_all);
    if ($n1<=0) {
	$compute_cooc_all=0;
    } else {
	($motifscan_filename2,$motifstrand_filename2)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$i2,$motif_scan_threshold_column,$species,$n_transfac_motifs);
	$n2=scan_all($wm2,$motifscan_filename2,$i2,$species,$ml2,$t2,$get_count,$singlefile,$seq_source,$read_all);
	if ($n2<=0) {
	    $compute_cooc_all=0;
	} else {
	    if ($i3>0) {
		($motifscan_filename3,$motifstrand_filename3)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$i3,$motif_scan_threshold_column,$species,$n_transfac_motifs);
		$n3=scan_all($wm3,$motifscan_filename3,$i3,$species,$ml3,$t3,$get_count,$singlefile,$seq_source,$read_all);
		if ($n3<=0) {
		    $compute_cooc_all=0;
		} else {
		    if ($i4>0) {
			($motifscan_filename4,$motifstrand_filename4)=create_scan_all_filename($exp_id,$seq_upstream_length,$seq_exon_length,$seq_intron_length,$i4,$motif_scan_threshold_column,$species,$n_transfac_motifs);
			$n4=scan_all($wm4,$motifscan_filename4,$i4,$species,$ml4,$t4,$get_count,$singlefile,$seq_source,$read_all);
			if ($n4<=0) {
			    $compute_cooc_all=0;
			}
		    }      # close check on i4>0
		}    # close check on n3<=0
	    }        # close check on i3>0
	}            # close check on n2<=0
    }                # close check on n1<=0

    $logtext.="n1=$n1\tn2=$n2\tn3=$n3\tn4=$n4\n";
    $logtext.="compute_cooc_all=$compute_cooc_all\n";
    return ($n1,$n2,$n3,$n4,$motifscan_filename1,$motifstrand_filename1,$motifscan_filename2,$motifstrand_filename2,$motifscan_filename3,$motifstrand_filename3,$motifscan_filename4,$motifstrand_filename4,$compute_cooc_all,$logtext);
}

sub usage{
    print "cooc_checkresults_v${prog_version}.pl\n";
    print "last modified = $last_modified\n";
    print "usage:\n";
    print "cooc_checkresults_v${prog_version}.pl <input_filename> <wm_filename> <n_species> <order_constraint> <upstream_filename_cl_sp1> (<upstream_filename_cl_sp2>) <info_filename_cl_sp1> (<info_filename_cl_sp2>) <transcript_lengths_filename_sp1_cl> (<transcript_lengths_filename_sp2_cl>) <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> (<mm2hs_map_all> <mm2hs_map_cl>) <rc> <threshold_column> <info_col_cl_sp1> (<info_col_cl_sp2>) <locuslink_only> <sortpos> <n_iter_random_overlap> <minian_transcripts> <file_label> <sp1> <sp2> <singlefile> <seq_source> <read_all> <n_transfac_motifs> <verbose>)\n";
    print "\noutput=\n\t(i) cooc_checkresults_v${prog_version}.log.txt\n\t(ii) cooc_checkresults_v${prog_version}.out.txt\n\t(iii) cooc_checkresults_v${prog_version}.mm.pos-cl.txt cooc_checkresults_v${prog_version}.hs.pos-cl.txt cooc_checkresults_v${prog_version}.mm.pos-all.txt cooc_checkresults_v${prog_version}.hs.pos-all.txt\n\t(iv) cooc_checkresults_v${prog_version}.mm.pos-cl-full.txt cooc_checkresults_v${prog_version}.hs.pos-cl-full.txt cooc_checkresults_v${prog_version}.mm.pos-all-full.txt cooc_checkresults_v${prog_version}.hs.pos-all-full.txt\n\t(v) cooc_checkresults_v${prog_version}.rej.txt\n\t(vi) cooc_checkresults_v${prog_version}.hsmm.trov-all.txt cooc_checkresults_v${prog_version}.hsmm.trov-cl.txt cooc_checkresults_v${prog_version}.hsmm.trov-cl.txt cooc_checkresults_v${prog_version}.hsmm.posov-cl.txt\n\t(vii) cooc_checkresults_v${prog_version}.err.txt cooc_checkresults_v${prog_version}.wrn.txt\n";
    print "init_i = initial i (default = 0)\n";
    print "verbose = verbose log output (default=0)\n";
    print "exp_id = experiment id (for nomenclature for the scan files only)\n";
    print "rc = 1 to consider the reverse complementary strand when scanning [default=0]\n";
    print "threshold_column = column within the wm_filename that contains the threshold for the motif [default=38]\n";
    print "info_col_cl_hs = column with locuslink information, set=cluster, species=hs, (from the info_cl_hs file) [default=11]\n";
    print "info_col_cl_mm = column with locuslink information, set=cluster, species=mm, (from the info_cl_mm file) [default=6]\n";
    print "locuslink_only = only consider entries with a valid locus link id if this is 1 (in the all sequences scan) [default=1]\n";
    print "sortpos = sort positions in call to scan programs if this is 1 [default=1]\n";
    print "n_iter_random_overlap = number of iteratinos in call to random_overlap.exe [default=1000]\n";
    print "sp1 = species 1 [default=hs]\n";
    print "sp2 = species 2 [default=mm]\n";
    print "singlefile = 1 if all the sequences are in a single file (e.g. yeast) and 0 otherwise (e.g. hs and mm) [default=0]\n";
    print "seq_source = ensembl, gbk [default=gbk]\n";
    print "read_all = to read all lines including those with 0 length [default=0]\n"; 
    print "n_transfac_motifs = for motifs <= n_transfac_motifs, search for scan all data in the transfac directory; if this is -1 then use the default values [default=-1]\n";
    gk();
    exit;
}

# new in version 9
# 11-10-2003: use temp_dir/scan_${exp_id} to store the scan results (rather than just temp_dir)
# 11-10-2003: changed the nomenclature of scan files to specify all 3 lengths and use seq_upstream/exon/intron_length variable
# 11-10-2003: add a file label to output files to be able to run multiple copies in parallel
# 11-10-2003: call cooc_single_1234_v1.exe
# 11-10-2003: use curr_arg to make parameter input simpler
# 11-10-2003: minian_transcripts is optional input now

# new in version 8
# 10-19-2003: added motif_scan_threshold_column to the file names for the scan cluster output
# 10-07-2003: separate the motifs between the transfac motifs and the new ones to speed things up
# 09-29-2003: added the number of renamed_files to check that renaming is working!!!
# 09-08-2003: solved plenty of bugs
# 09-08-2003: make number of iterations in compute_random_overlap an input (n_iter_random_overlap)
# 09-05-2003: use mm2hs_map_all and mm2hs_map_cl as a matrix
# 09-04-2003: only process the "all" sequences case if there was some overlap in hs or mm
# 09-04-2003: major change: do not call remove redundancies program, this is now done within the cooc_single programs
# 08-27-2003: add sortpos to scan_motif sub
# 08-17-2003: call scan_1m_allseq (this version uses separate files for each chromosome and considers the upstream, exon and intron regions)
# 08-17-2003: do not need names for upstream files for all transcripts

# 07-24-2003: call cooc_single_?m_v3
# 07-24-2003: use thresholds from the same motifs file
# 07-24-2003: add order_constraint
# 07-24-2003: call scan_motif_v7

# 05-17-2003: add count orthologous genes
# 05-16-2003: add print gene names
# 05-15-2003: add output for cooc_all_hs and cooc_all_mm (positions of motifs in all genes)
# 05-15-2003: call scan_motif_v5
# 05-15-2003: store also the log file from the scan(s)

# 03-25-2003 run invbinocdf and run the cooc_single program up to a maximum number of transcripts only
# 03-20-2003 changed the c:/temp directory to a more general temporary directory
# 03-20-2003 changed the c:/life directory to a more general code directory

# new in v6 (major revision):
# call cooc_single_2/3/4m_v2 (using different scan input format, not using the wc on the scan file, processing transript by transcript)
# call compare_cooc_list (to compare the hs and mm results
# call scan_motif_v4 (scan fasta formatted file)

# new in v5:
# save the positions for the cluster for both mice and human

# new in v3: major changes
# now we scan both in mice and in human, both in cluster as well as in background

# new in v2: attempt to store temporarily the scans in case they need to be reused, this should make things work faster...
# 03-06-2003: subtract the co-occurrences in the cluster transcripts before computing the binomial probabilities
# 03-07-2003: make it work for 4 motifs as well
# 03-07-2003: change the readout columns based on the new format
# 03-24-2003: changed the temporary weight matrix to be stored in  the current directory from which the code is run, this avoids overlaps in the temp directory when running multiple copies of this code in parallel

