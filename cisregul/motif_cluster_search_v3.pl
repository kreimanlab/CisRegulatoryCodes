#!/usr/bin/perl
# motif_cluster_search_v3.pl
# run together all the programs of the motif search

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
$TEMP_DIR=$ENV{'TEMPDIR'};

# new in version 3
# 02-20-2004: allow adding additional motifs
# 02-20-2004: allow scanning all sequences first
# 02-29-2004: added n_transfac_motifs as input (if < 0, then use the default values in cooc_checkresults_v12.pl)
# 03-03-2004: call cooc_checkresults_v9 (faster version)
# 03-04-2004: call list_nonredundant_wm_v2.pl 
# 03-04-2004: added more parameters to filter motifs in list_nonredundant_wm_v2.pl
# 03-11-2004: implement the scanning of all sequences before calling cooc_checkresults
# 03-17-2004: create the consensus file and the matlab file

# list of executable functions called:
# motif search: (1) alignACE, (2) MotifSampler.linux.exe, (3) meme
# (4) cooc_4m_v${cooc_c_code_version}.exe 
# (5) scan_allseq_v${scan_allseq_version}.exe 
# (6) compare_cooc_outputs_v${c_code_version}.exe 

require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/motif_cluster_search_lib.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";

$n_args=$#ARGV+1;
$program_name="motif_cluster_search";
$program_version=3;
$last_modified="03_13_2004";
if ($n_args<1) {
    print "${program_name}_${program_version}.pl\n";
    print "usage:\n";
    print "${program_name}_${program_version}.pl <parameters_filename>\n";
    gk();
    print "last modified: ${last_modified}\n";
    exit;
}

# default parameters
$search_type="ll";                    # assume that we start with a list of locus link identifiers (this should be changed if the input is different)
$max_ul_length=2500;
$max_exon1_length=2500;
$max_intron1_length=2500;
$print_notfound=1;
$motif_comparison_threshold=70;       # threhold spearman correlation to determine that two motifs are very similar
$process_alignace=1;
$process_motifsampler=1;
$process_meme=0;
$rc=1;                                # 1 to consider both strands when scanning
$sortpos=1;                           # 1 to sort positions in call to scan function
$scan_threshold_column=38;            # column where the threshold is; used in call to scan function
$process_bl2seq=2;                    # to process bl2seq versions of the sequences
$include_noorth=1;                    # include genes with no orthologous sequence
$additional_motifs="none";            # this can be used to add a list of motifs that the investigator suspects to be involved
$n_transfac_motifs=-1;                # use default values for n_transfac_motifs (unless this is an input to this program)

# default parameters that do not show up in the parameters file
$nt_bckg1=20;                        # for the initial exploration step, size of the background set / size of the cluster set

$parameters_filename=$ARGV[0];
($log_text,$input_filename,$exp_name,$exp_id,$species1,$species2,$search_type,$max_ul_length,$max_exon1_length,$max_intron1_length,$motif_comparison_threshold,$rc,$sortpos,$scan_threshold_column,$process_bl2seq,$p,$min_dist_factor,$init_max_n_motifs,$finit_max_n_motifs,$init_max_dist,$finit_max_dist,$init_order_constraint,$finit_order_constraint,$init_max_seq_length,$finit_max_seq_length,$min_transcripts,$locuslink_only,$n_iter_random_overlap,$source,$include_noorth,$process_motifsampler,$additional_motifs,$n_transfac_motifs,$information_threshold,$seq_variety_threshold,$n_seqs_threshold,$min_motif_length_threshold,$max_motif_length_threshold)=read_parameters_from_file($parameters_filename);           

# open log file
($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
$Month++;
$log_filename=">${program_name}.${program_version}.${Year}.${Month}.${Day}.${Hour}.${Minute}.log.txt";
open ($log_file,$log_filename) || die "could not open $log_filename for output\n";
$t=get_time_string();

print $log_file "% ${program_name}_v${program_version}.pl\n";
print $log_file "% $t\n";
print $log_file "% parameters_filename=$parameters_filename\n";
print $log_file "$log_text\n";

############################################
# create directories if they do not exist  #
############################################
$gen_dir="${DATA_DIR}/databh/${exp_name}";
$upstream_dir="${gen_dir}/upstream_sequences";
$analysis_dir="${gen_dir}/analysis";
$scan_dir="${analysis_dir}/scan";
$cooc_dir="${analysis_dir}/cooc";
$check_dir="${cooc_dir}/check";
$old_dir="${gen_dir}/old";
$scanall_dir="${TEMP_DIR}/scan_${exp_id}";
print $log_file "\nchecking/creating directories...\n";
if (-e $upstream_dir) {
    print $log_file "\t$upstream_dir already exists...\n";
} else {
    $code_exec="mkdir ${upstream_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $analysis_dir) {
    print $log_file "\t$analysis_dir already exists...\n";
} else {
    $code_exec="mkdir ${analysis_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $scan_dir) {
    print $log_file "\t$scan_dir already exists...\n";
} else {
    $code_exec="mkdir ${scan_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $cooc_dir) {
    print $log_file "\t$cooc_dir already exists...\n";
} else {
    $code_exec="mkdir ${cooc_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $check_dir) {
    print $log_file "\t$check_dir already exists...\n";
} else {
    $code_exec="mkdir ${check_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $old_dir) {
    print $log_file "\t$old_dir already exists...\n";
} else {
    $code_exec="mkdir ${old_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
if (-e $scanall_dir) {
    print $log_file "\t$scanall_dir already exists...\n";
} else {
    $code_exec="mkdir ${scanall_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}

##################################################################################
# ll2nm2ug.pl <search_filename> <search_col> <species> <searchtype>              #
# to add further information about the entries                                   #
#   output is stored in query.txt (hits) and log.txt (miss+multiple matches)     # 
#   searchtype is ll / nm / ug                                                   #
# if the ${exp_id}.list.txt file already exists, then this step is skipped                 #
##################################################################################
$query_output_filename="${gen_dir}/${exp_id}.list.txt";
$log_output_filename="${gen_dir}/ll2nm2ug.log.txt";
if (-e $query_output_filename) {
    print $log_file "it seems that we already called ll2nm2ug.pl [$query_output_filename]; skipping...\n";
} else {
    print $log_file "$query_output_filename does not exist. calling ll2nm2ug.pl\n";
    $search_column=1;           # assume that the identifiers to search for are in the first column (this should be changed if the input is different)
    $code_exec="perl ${DATA_DIR}/db/locus_link/code/ll2nm2ug.pl ${input_filename} ${search_column} ${species1} ${search_type}";
    print $log_file "calling ll2nm2ug.pl\n$code_exec\n";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\tsystem_status=$syst_stat\n";
    if ($syst_stat ne 0) {
	print "error!\n";
	exit;
    }
    ($syst_stat,$warning_msg)=my_rename("query.txt",$query_output_filename,0,1);         # non-verbose mode, exit on error
    print $log_file "$warning_msg\n"; 
    ($syst_stat,$warning_msg)=my_rename("log.txt",$log_output_filename,0,1);             # non-verbose mode, exit on error
    print $log_file "$warning_msg\n";
}
($n,$n_comments)=line_count($query_output_filename);
$n_transcripts_cl=$n-$n_comments;
if ($n_transcripts_cl<=0) {
    print "error!\nmotif_cluster_search_v1.pl\nll2nm2ug.pl\nn_transcripts_cl=$n_transcripts_cl\n";
    exit;
}
print $log_file "n_transcripts_cl=${n_transcripts_cl}\n";

# get the sequence source depending on the species
if ($source eq "none") {
    if ($species1 eq "hs") {
	$source="ncbi";
    } 
    if ($species1 eq "mm") {
	$source="ncbi";
    }
    if ($species1 eq "dm") {
	$source="ensembl";
    }
    if ($species1 eq "sc") {
	#$source="ncbi";
	$source="sgd";
    }
}
if ($source eq "none") {
    print "\n\nERROR! I do not recognize species1=$species. Please fix...\n";
    exit;
}
$source_original=$source;
if ($source =~ /bl2seq/) {
    $source="${source}2";
}

##################################
# retrieving sequences           #
# call retrieve_upstream_v3.pl   #
# retrieve_upstream_v3.pl <ids_filename> <species> (<max_upstream> <max_exon> <max_intron> <print_notfound> <separator>)
##################################
if ($source_original =~ /bl2seq/) {
    $source="ncbi";                      # because to retrieve the sequences we want the ncbi files to retrieve the sequences
}
$log_upstream_sp1_filename="${upstream_dir}/${exp_id}.${species1}.${source}.log.txt";
$upstream_filename_sp1="${upstream_dir}/${exp_id}.${species1}.${source}.fasta";
$lengths_upstream_sp1_filename="${upstream_dir}/${exp_id}.${species1}.${source}.lengths.txt";
$upstream_filename_motifsearch_sp1=$upstream_filename_sp1;
if ($source_original =~ /bl2seq/) {
    $upstream_filename_scan_sp1="${upstream_dir}/${exp_id}.${species1}.bl2seq2.fasta";
    $lengths_upstream_scan_sp1_filename="${upstream_dir}/${exp_id}.${species1}.bl2seq2.lengths.txt";
} else {
    $upstream_filename_scan_sp1=$upstream_filename_sp1;
    $lengths_upstream_scan_sp1_filename=$lengths_upstream_sp1_filename;
}
#if (-e $log_upstream_sp1_filename) {
if (-e $upstream_filename_sp1) {
    print $log_file "it seems that we already retrieved the sequences for ${species1} [$upstream_filename_sp1], skipping...\n";
} else {
    print $log_file "$upstream_filename_sp1 not found\n";
    print $log_file "retrieving upstream sequences for species1=${species1}...\n";

    # get the list of entries for species 1
    $list_sp1_filename="${gen_dir}/${exp_id}.list-ll-${species1}.txt";
    if ( ($species1 eq "hs") | ($species1 eq "mm") ) {
	#$code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 5";   # extract the locus link id (hs) to get the upstream sequences
	$code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 2";   # extract the locus link id (hs) to get the upstream sequences
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	($syst_stat,$warning_msg)=my_rename("${gen_dir}/col2_${exp_id}.list.txt",$list_sp1_filename,0,1);
	# } else {
	#if ($species1 eq "mm") {
	#    $code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 2";   # extract the gene identifier to get the upstream sequences
	#    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	#    ($syst_stat,$warning_msg)=my_rename("${gen_dir}/col2_${exp_id}.list.txt",$list_sp1_filename,0,1);
    } else {
	$code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 1";   # extract the gene identifier to get the upstream sequences
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("${gen_dir}/col1_${exp_id}.list.txt",$list_sp1_filename,0,1);
	print $log_file "\t$warning_msg\n";
    }
    
    # get the upstream sequences
    if ($species1 eq "dm") {
	$code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/merge_upexin.pl ${list_sp1_filename} ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_upstream.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_exon1.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_intron1.info.txt ${DATA_DIR}/db/ensembl/dm/dm.ensembl.upstream.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.exon1.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.intron1.fasta ${max_ul_length} ${max_exon1_length} ${max_intron1_length} 1 0 none";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("merge_upexin.out.txt",$upstream_filename_sp1,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("log.txt",$log_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("merge_upexin.lengths.txt",$lengths_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
    }
    if ($species1 eq "sc") {
	$code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/merge_upexin.pl ${list_sp1_filename} ${DATA_DIR}/db/upstream_sequences/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream_sequences/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream_sequences/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${max_ul_length} ${max_exon1_length} ${max_intron1_length} 1 0 none";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("merge_upexin.out.txt",$upstream_filename_sp1,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("log.txt",$log_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("merge_upexin.lengths.txt",$lengths_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
    }
    if ( ($species1 eq "hs") | ($species1 eq "mm") ) {
	$code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/retrieve_upstream_v3.pl ${list_sp1_filename} ${species1} ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${print_notfound}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.out.txt",$upstream_filename_sp1,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.log.txt",$log_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
	($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.lengths.txt",$lengths_upstream_sp1_filename,0,1);
	print $log_file "\t$warning_msg";
    }
}         # close check on whether we already retrieved the upstream sequences for species 1
$info_sp1_filename="${upstream_dir}/${exp_id}.${species1}.${source}.list.txt";
if (-e $info_sp1_filename) {
} else {
    # create a list of sequences
    $code_exec="grep \"\>\" ${upstream_filename_sp1} \> ${info_sp1_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    $code_exec="perl -pi.bak -e \'s/\>//' ${info_sp1_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    $n_unlinked=unlink("${info_sp1_filename}.bak");
}
if ( ($species1 eq "hs")  | ($species1 eq "mm") ) {
    $info_sp1_filename_matlab="${upstream_dir}/${exp_id}.${species1}.${source}.listr.txt";           # reformat the list so that it can be read in matlab
    if (-e $info_sp1_filename_matlab) {
    } else {
	$code_exec="cp ${info_sp1_filename} ${info_sp1_filename_matlab}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\\n/endline/\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\\s\+//g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/endline/\\n/g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\|/\\t/g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\"//g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/GI\:\\d+\,LocusID\://g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\,MIM\:\\d+//g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\,MGI\:\\d+//g\' ${info_sp1_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$n_unlinked=unlink("${info_sp1_filename_matlab}.bak");
    }
}

#####################
# process species 2 #
#####################
if ($species2 ne "none") {
    $log_upstream_sp2_filename="${upstream_dir}/${exp_id}.${species2}.${source}.log.txt";
    $upstream_filename_sp2="${upstream_dir}/${exp_id}.${species2}.${source}.fasta";
    $lengths_upstream_sp2_filename="${upstream_dir}/${exp_id}.${species2}.${source}.lengths.txt";
    $upstream_filename_motifsearch_sp2=$upstream_filename_sp2;
    if ($source_original =~ /bl2seq/) {
	$upstream_filename_scan_sp2="${upstream_dir}/${exp_id}.${species2}.bl2seq2.fasta";
	$lengths_upstream_scan_sp2_filename="${upstream_dir}/${exp_id}.${species2}.bl2seq2.lengths.txt";
    } else {
	$upstream_filename_scan_sp2=$upstream_filename_sp2;
	$lengths_upstream_scan_sp2_filename=$lengths_upstream_sp2_filename;
    }
    if (-e $upstream_filename_sp2) {
	print $log_file "it seems that we already retrieved the sequences for ${species2} [$upstream_filename_sp2], skipping...\n";
    } else {
	print $log_file "$upstream_filename_sp2 not found\n";
	print $log_file "retrieving upstream sequences for species2=${species2}...\n";
	# get the list of entries to retrive the upstream sequences
	$list_sp2_filename="${gen_dir}/${exp_id}.list-ll-${species2}.txt";
	if ( ($species2 eq "hs") | ($species2 eq "mm") ) {
	    $code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 5";   # extract the locus link id (hs) to get the upstream sequences
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    ($syst_stat,$warning_msg)=my_rename("${gen_dir}/col5_${exp_id}.list.txt",$list_sp2_filename,0,1);
	    #} else {
	    #if ($species2 eq "mm") {
	    #$code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 2";   # extract the gene identifier to get the upstream sequences
	    #$code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 5";   # extract the gene identifier to get the upstream sequences
	    #($syst_stat,$warning_msg)=run_system($code_exec,1);
	    #($syst_stat,$warning_msg)=my_rename("${gen_dir}/col2_${exp_id}.list.txt",$list_sp2_filename,0,1);
	} else {
	    $code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${gen_dir}/${exp_id}.list.txt 1";   # extract the gene identifier to get the upstream sequences
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    ($syst_stat,$warning_msg)=my_rename("${gen_dir}/col1_${exp_id}.list.txt",$list_sp2_filename,0,1);
	}
	# get the upstream sequences
	if ($species2 eq "dm") {
	    $code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/merge_upexin.pl ${list_sp2_filename} ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_upstream.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_exon1.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_intron1.info.txt ${DATA_DIR}/db/ensembl/dm/dm.ensembl.upstream.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.exon1.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.intron1.fasta ${max_ul_length} ${max_exon1_length} ${max_intron1_length} 1 0 none";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	    ($syst_stat,$warning_msg)=my_rename("merge_upexin.out.txt",$upstream_filename_sp2,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("log.txt",$log_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("merge_upexin.lengths.txt",$lengths_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	}
	if ($species2 eq "sc") {
	    $code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/merge_upexin.pl ${list_sp2_filename} ${DATA_DIR}/db/upstream_sequences/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream_sequences/sc/lists/exon_list.txt ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${DATA_DIR}/db/upstream_sequences/sc/sc.sgd.upstream.fasta ${max_ul_length} ${max_exon1_length} ${max_intron1_length} 1 0 none";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	    ($syst_stat,$warning_msg)=my_rename("merge_upexin.out.txt",$upstream_filename_sp2,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("log.txt",$log_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("merge_upexin.lengths.txt",$lengths_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	}
	if ( ($species2 eq "hs") | ($species2 eq "mm") ) {
	    $code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/retrieve_upstream_v3.pl ${list_sp2_filename} ${species2} ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${print_notfound}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	    ($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.out.txt",$upstream_filename_sp2,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.log.txt",$log_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	    ($syst_stat,$warning_msg)=my_rename("retrieve_upstream_v3.lengths.txt",$lengths_upstream_sp2_filename,0,1);
	    print $log_file "\t$warning_msg";
	}
    }         # close check on whether we already retrieved the upstream sequences for species 2
    $info_sp2_filename="${upstream_dir}/${exp_id}.${species2}.${source}.list.txt";
    if (-e $info_sp2_filename) {
    } else {
	# create a list of sequences
	$code_exec="grep \"\>\" ${upstream_filename_sp2} \> ${info_sp2_filename}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	$code_exec="perl -pi.bak -e \'s/\>//' ${info_sp1_filename}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	$n_unlinked=unlink("${info_sp2_filename}.bak");
    }
    $info_sp2_filename_matlab="${upstream_dir}/${exp_id}.${species2}.${source}.listr.txt";           # reformat the list so that it can be read in matlab
    if (-e $info_sp2_filename_matlab) {
    } else {
	$code_exec="cp ${info_sp2_filename} ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\\n/endline/\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\\s\+//g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/endline/\\n/g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\|/\\t/g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\"//g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/GI\:\\d+\,LocusID\://g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\,MIM\:\\d+//g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\,MGI\:\\d+//g\' ${info_sp2_filename_matlab}";($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	$n_unlinked=unlink("${info_sp2_filename_matlab}.bak");
    }
}             # close check on whether species2 exists

###################################
# comparing sp1 and sp2 sequences #
###################################
# 02-17-2004: call preporth4alignace_v3.pl
# 02-26-2004: added species to the preporth4alignace log file 
# 02-26-2004: run preporth4alignace for the second species
if ($species2 ne "none") {
    $preporth4alignace_version=3;
    $log_bl2seq_filename_sp1="${upstream_dir}/preporth4alignace.${exp_id}.${species1}.log.txt";
    $log_bl2seq_filename_sp2="${upstream_dir}/preporth4alignace.${exp_id}.${species2}.log.txt";
    if ( (-e $log_bl2seq_filename_sp1) & (-e $log_bl2seq_filename_sp2) ) {
	print $log_file "We already performed the comparison of the ${species1} and ${species2} sequences [$log_bl2seq_filename_sp1]; skipping...\n";
    } else {
	# sequence comparison
	print $log_file "I could not find the file $log_bl2seq_filename_sp1\n";
	$bl2seq_output_filename="${upstream_dir}/${exp_id}.${species1}-bl2seq-${species2}.txt";
	if (-e $bl2seq_output_filename) {
	    print $log_file "\tWe already compared the sequences in the two speices [$bl2seq_output_filename]\n";
	} else {
	    print $log_file "\tI could not find $bl2seq_output_filename\n";
	    print_title("comparing ${species1} and ${species2} sequences...",$log_file); 
	    $map_filename="none";$max_seq_length=50000;$expect_threshold=1;
	    $filteron=1;$gap_open_cost=1;$gap_extend_cost=2;$dropoff_value=50;$word_size=11;$mismatch_penalty=-2;$match_reward=1;$expectation_value=10;$query_strand=3;
	    print $log_file "sp1-sp2 comparison parameters:\n";
	    print $log_file "map_filename=$map_filename\n";print $log_file "expect_threshold=$expect_threshold\n";print $log_file "filteron=$filteron\n";
	    print $log_file "gap_open_cost=$gap_open_cost\n";print $log_file "gap_extend_cost=$gap_extend_cost\n";print $log_file "dropoff_value=$dropoff_value\n";
	    print $log_file "word_size=$word_size\n";print $log_file "mismatch_penalty=$mismatch_penalty\n";print $log_file "match_reward=$match_reward\n";
	    print $log_file "expectation_value=$expectation_value\n";print $log_file "query_strand=$query_strand\n";

	    $code_exec="perl ${CODE_DIR}/perl/bl2seq_2sets_v2.pl ${upstream_filename_sp1} ${upstream_filename_sp2} ${map_filename} ${max_seq_length} ${expect_threshold} ${filteron} ${gap_open_cost} ${gap_extend_cost} ${dropoff_value} ${word_size} ${mismatch_penalty} ${match_reward} ${expectation_value} ${query_strand}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\t$warning_msg\n";
	    $bl2seq_output_filename="${upstream_dir}/${exp_id}.${species1}-bl2seq-${species2}.txt";
	    ($syst_stat,$warning_msg)=my_rename("bl2seq_results.txt",$bl2seq_output_filename,0,1);
	    print $log_file "\t$warning_msg";
	}

	#####################################################################################################################
	# preporth4alignace_v3.pl <sequence_filename> <bl2seq_results_filename> (<species> <map_filename> <include_noorth>) #
	#####################################################################################################################
	if (-e $log_bl2seq_filename_sp1) {
	    print log_file "\t$log_bl2seq_filename_sp1 already exists\n";
	} else {
	    # preporth4alignace, species 1
	    $code_exec="perl ${CODE_DIR}/perl/preporth4alignace_v${preporth4alignace_version}.pl ${upstream_filename_sp1} ${bl2seq_output_filename} 1 nomap ${include_noorth}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\t$warning_msg\n";
	    $orth1_filename_sp1="${upstream_dir}/${exp_id}.${species1}.bl2seq1.fasta";      # with ns wherever there is no match
	    $orth2_filename_sp1="${upstream_dir}/${exp_id}.${species1}.bl2seq2.fasta";      # keeps all nucleotides within long conserved sequence stretches
	    $orth3_filename_sp1="${upstream_dir}/${exp_id}.${species1}.bl2seq3.fasta";      # removes all ns (useful to speed up motif search)
	    ($syst_stat,$warning_msg)=my_rename("orth1.txt",$orth1_filename_sp1);
	    ($syst_stat,$warning_msg)=my_rename("orth2.txt",$orth2_filename_sp1);
	    ($syst_stat,$warning_msg)=my_rename("orth3.txt",$orth3_filename_sp1);
	    ($syst_stat,$warning_msg)=my_rename("preporth4alignace.log.txt",$log_bl2seq_filename_sp1);
	    $upstream_filename_motifsearch_sp1=$orth3_filename_sp1;
	} 
	if (-e $lengths_upstream_scan_sp1_filename) {
	    print log_file "\t$lengths_upstream_scan_sp1_filename already exists\n";
	} else {
	    if ($source_original =~ /bl2seq/) {
		$upstream_filename_scan_sp1=$orth2_filename_sp1;
		# get the lengths for the orth2_filename_sp1 sequences
		$code_exec="grep index ${log_bl2seq_filename_sp1} \> index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg";
		$code_exec="perl -pi.bak -e \'s/\>//' index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg";
		$code_exec="perl ${CODE_DIR}/perl/extract_column.pl index_in.txt 2";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg";
		($n,$n_comments)=line_count("col2_index_in.txt");
		print $log_file "\n\tcol2_index_in\t$n lines\t$n_comments comment lines\n";
		$code_exec="perl ${CODE_DIR}/perl/extract_rows_v1.pl ${lengths_upstream_sp1_filename} col2_index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		($syst_stat,$warning_msg)=my_rename("${lengths_upstream_sp1_filename}.extracted",$lengths_upstream_scan_sp1_filename,0,1);print $log_file "$warning_msg";
		$n_unlinked=unlink("index_in.txt","index_in.txt.bak","col2_index_in.txt");
	    }
	}

	if (-e $log_bl2seq_filename_sp2) {
	    print log_file "\t$log_bl2seq_filename_sp2 already exists\n";
	} else {
	    print $log_file "\n\tcalling preporth4alignace_v${preporth4alignace_version}.pl species=${species2}\n";
	    # preporth4alignace, species 2
	    $code_exec="perl ${CODE_DIR}/perl/preporth4alignace_v${preporth4alignace_version}.pl ${upstream_filename_sp2} ${bl2seq_output_filename} 2 nomap ${include_noorth}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\t$warning_msg\n";
	    $orth1_filename_sp2="${upstream_dir}/${exp_id}.${species2}.bl2seq1.fasta";      # with ns wherever there is no match
	    $orth2_filename_sp2="${upstream_dir}/${exp_id}.${species2}.bl2seq2.fasta";      # keeps all nucleotides within long conserved sequence stretches
	    $orth3_filename_sp2="${upstream_dir}/${exp_id}.${species2}.bl2seq3.fasta";      # removes all ns (useful to speed up motif search)
	    ($syst_stat,$warning_msg)=my_rename("orth1.txt",$orth1_filename_sp2);
	    ($syst_stat,$warning_msg)=my_rename("orth2.txt",$orth2_filename_sp2);
	    ($syst_stat,$warning_msg)=my_rename("orth3.txt",$orth3_filename_sp2);
	    ($syst_stat,$warning_msg)=my_rename("preporth4alignace.log.txt",$log_bl2seq_filename_sp2);
	    $upstream_filename_motifsearch_sp2=$orth3_filename_sp2;
	} 
	if (-e $lengths_upstream_scan_sp2_filename) {
	    print log_file "\t$lengths_upstream_scan_sp2_filename already exists\n";
	} else {
	    if ($source_original =~ /bl2seq/) {
		$upstream_filename_scan_sp2=$orth2_filename_sp2;
		# get the lengths for the orth2_filename_sp1 sequences
		$code_exec="grep index ${log_bl2seq_filename_sp2} \> index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		$code_exec="perl -pi.bak -e \'s/\>//' index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		$code_exec="perl ${CODE_DIR}/perl/extract_column.pl index_in.txt 2";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		($n,$n_comments)=line_count("col2_index_in.txt");
		print $log_file "col2_index_in\t$n lines\t$n_comments comment lines\n";
		$code_exec="perl ${CODE_DIR}/perl/extract_rows_v1.pl ${lengths_upstream_sp2_filename} col2_index_in.txt";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		($syst_stat,$warning_msg)=my_rename("${lengths_upstream_sp2_filename}.extracted",$lengths_upstream_scan_sp2_filename,0,1);
		print $log_file "$warning_msg";
		$n_unlinked=unlink("index_in.txt","index_in.txt.bak","col2_index_in.txt");
	    }
	}
    }      # close check on whether $log_bl2seq_filename_sp2 exists (i.e. whether we already run preporth4alignace)
    if ($process_bl2seq>0) {
	$info_sp1_filename_orig="${upstream_dir}/${exp_id}.${species1}.bl2seq2.list.txt";
	$orth2_filename_sp1="${upstream_dir}/${exp_id}.${species1}.bl2seq2.fasta";
	if (-e $info_sp1_filename_orig) {
	    print $log_file "$info_sp1_filename already exists\n";
	} else {        # create a list of sequences
	    $code_exec="grep \"\>\" ${orth2_filename_sp1} \> ${info_sp1_filename_orig}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    $code_exec="perl -pi.bak -e \'s/\>//' ${info_sp1_filename_orig}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	}
	($n,$n_comments)=line_count(${info_sp1_filename_orig});
	$n_transcripts_cl_orig=$n_transcripts_cl;
	$n_transcripts_cl=$n-$n_comments;
	print $log_file "$info_sp1_filename_orig: n_transcripts = $n_transcripts_cl (n_transcripts_cl_orig=$n_transcripts_cl_orig)\n";
    }
}          # close check on nspecies eq 2
if ($source_original =~ /bl2seq/) {
    $source="bl2seq2";
    print $log_file "source_original=$source_original\n";
    print $log_file "source=$source\n";
} 

############################
# motif search algorithms  #
############################
if ( ($species1 eq "hs") | ($species1 eq "mm") ) {
    $species12="${species1}${species2}";
} else {
    $species12=$species1;
}
if ($process_alignace eq 1) {
    $alignace_dir="${analysis_dir}/alignace";
    if (-e $alignace_dir) {
	print $log_file "\tdirectory ${alignace_dir} already exists...\n";
    } else {
	$code_exec="mkdir $alignace_dir";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
    }
    $alignace_subdir="${analysis_dir}/alignace/original_files";
    if (-e $alignace_subdir) {
	print $log_file "\tdirectory ${alignace_subdir} already exists...\n";
    } else {
	$code_exec="mkdir $alignace_subdir";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
    }

    $alignace_list_filename="${analysis_dir}/alignace/original_files/list.alignace.${exp_id}.${species12}.txt";
    if (-e $alignace_list_filename) {
	print $log_file "it seems that we already run alignACE to search for new motifs [$alignace_list_filename]; skipping...\n";
    } else {
	print $log_file "$alignace_list_filename does not exist\nprocessing alignACE...\n";
	print_title("processing alignACE",$log_file); 
	#$code_exec="cp ${upstream_dir}/*.fasta ${DATA_DIR}/util/alignace";
	#($syst_stat,$warning_msg)=run_system($code_exec,1);
	#print $log_file "\t$warning_msg\n";

	if ($species1 eq "sc") {
	    $gcback=0.38;
	} else {
	    $gcback=0.45;
	}
	print $log_file "gcback=$gcback\n";

	$code_exec="${DATA_DIR}/util/alignace/general_alignace_v5.sh ${exp_id} 10 ${species1} ${source} ${gcback} ${upstream_dir}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	if ($species2 ne "none") {
	    $code_exec="${DATA_DIR}/util/alignace/general_alignace_v5.sh ${exp_id} 10 ${species2} ${source} ${gcback} ${upstream_dir}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    
	    $code_exec="${DATA_DIR}/util/alignace/general_alignace_v5.sh ${exp_id} 10 ${species1} bl2seq3 ${gcback} ${upstream_dir}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    $code_exec="${DATA_DIR}/util/alignace/general_alignace_v5.sh ${exp_id} 10 ${species2} bl2seq3 ${gcback} ${upstream_dir}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	}

	$code_exec="mv ${upstream_dir}/${exp_id}*alignace*.txt ${alignace_subdir}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
    }

    #############################################################
    # process motif search output files and compute thresholds  # 
    #############################################################
    $alignace_motifs_filename="${alignace_dir}/${exp_id}.${species12}.alignace.txt";
    if (-e $alignace_motifs_filename) {
	print $log_file "we already processed the output of alignace [$alignace_motifs_filename]; skipping...\n";
    } else {
	print $log_file "$alignace_motifs_filename does not exist...\n";
	print_title("processing alignace output...",$log_file); 
	#$alignace_output_list="${alignace_subdir}/list.alignace.hsmm.txt";      alignace_list_filename above
	$code_exec="ls -1 ${alignace_subdir}/${exp_id}* > ${alignace_list_filename}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	$code_exec="perl ${DATA_DIR}/util/alignace/process_alignace_output.pl ${species1} ${alignace_list_filename} 0";
	@code_output=`${code_exec}`;
	($syst_stat,$n_matches)=extract_text_from_array("system_status",@code_output);
	if ($syst_stat ne 0) {
	    print $log_file "ERROR!\ncode_exec=$code_exec\nsystem_status=$syst_stat\nexiting...\n";
	    exit;
	}
	($n_files,$n_matches)=extract_text_from_array("n_files",@code_output);
	($n_motifs_processed,$n_matches)=extract_text_from_array("n_processed",@code_output);
	($n_motifs_output,$n_matches)=extract_text_from_array("n_motifs",@code_output);
	($n_filtered_out,$n_matches)=extract_text_from_array("n_filtered_out",@code_output);
	($n_compared_out,$n_matches)=extract_text_from_array("n_compared_out",@code_output);

	print $log_file "\tn_files=$n_files\n";
	print $log_file "\tn_motifs_procesed=$n_motifs_processed\n";
	print $log_file "\tn_motifs_output=$n_motifs_output\n";
	print $log_file "\tn_filtered_out=$n_filtered_out\n";
	print $log_file "\tn_compared_out=$n_compared_out\n";

	($syst_stat,$warning_msg)=my_rename("process_alignace_output.out.txt",$alignace_motifs_filename);
	print $log_file "\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("process_alignace_output.log.txt","${alignace_dir}/${exp_id}.${species12}.alignace.log.txt");
	print $log_file "\t$warning_msg\n";
    }
}
################
# motifsampler #
################
if ($process_motifsampler eq 1) {
    $motifsampler_dir="${analysis_dir}/motifsampler";
    if (-e $motifsampler_dir) {
	print $log_file "\tdirectory ${motifsampler_dir} already exists...\n";
    } else {
	$code_exec="mkdir $motifsampler_dir";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
    }
    $motifsampler_subdir="${motifsampler_dir}/original_files";
    if (-e $motifsampler_subdir) {
	print $log_file "\tdirectory ${motifsampler_subdir} already exists...\n";
    } else {
	$code_exec="mkdir ${motifsampler_subdir}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
    }
    $motifsampler_motifs_filename="${motifsampler_dir}/${exp_id}.${species12}.motifsampler.txt";
    if (-e $motifsampler_motifs_filename) {
	print $log_file "i already processed motifsampler [$motifsampler_motifs_filename]; skipping...\n";
    } else {
	print $log_file "$motifsampler_motifs_filename does not exist...\n";
	$motif_sampler_output1="${motifsampler_subdir}/${exp_id}.${species1}.matrix.${source}.txt";
	$motif_sampler_output2="${motifsampler_subdir}/${exp_id}.${species1}.output.${source}.txt";
	if (-e $motif_sampler_output1) {
	    print $log_file "it seems that we already run motifsampler to search for new motifs [$motif_sampler_output1]; skipping...\n";
	} else {
	    print $log_file "$motif_sampler_output1 does not exist\nprocessing motifsampler...\n";
	    print_title("processing motif_sampler",$log_file); 
	    
	    if ($species1 eq "sc") {
		$background_file="need_to_fill_in";
	    } else {
		if ($species1 eq "dm") {
		    $background_file="dm_order3.txt";
		} else {
		    if ($species1 eq "hs") {
			$background_file="hs_order3.12195339.813.txt";
		    } else {
			if ($species1 eq "mm") {
			    $background_file="mm_order3.8692673.1066.txt";
			}
		    }
		}
	    }
	    $background_file="${DATA_DIR}/util/motif_sampler/${background_file}";
	    if ($source =~ /bl2seq/) { 
		$seqs_filename="${upstream_dir}/${exp_id}.${species1}.bl2seq3.fasta";                   # re-define in case it was not defined before
	    } else {
		$seqs_filename=$upstream_filename_motifsearch_sp1;
	    }
	    print $log_file "background_file=$background_file\n";
	    $code_exec="${DATA_DIR}/util/motif_sampler/MotifSampler.linux.exe -f ${seqs_filename} -b ${background_file} -n 10 -o ${motif_sampler_output2} -m ${motif_sampler_output1}";
	    print $log_file "$code_exec\n";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\t$warning_msg\n";
	}

	#############################################################
	# process motif search output files and compute thresholds  # 
	#############################################################
	print $log_file "$motifsampler_motifs_filename does not exist...\n";
	print_title("processing motif_sampler output...",$log_file); 
	# in /data/kreiman/util/motif_sampler
	# process_output_motif_sampler.pl <matrix_filename> <sites_filename> (<species> <mean_information_threshold> <seq_variety_threshold> <n_seqs_threshold> <comp_nmaitrx_threshold> <rc> <verbose>)
	$code_exec="perl ${DATA_DIR}/util/motif_sampler/process_output_motif_sampler.pl ${motif_sampler_output1} ${motif_sampler_output2} ${species1}";
	@code_output=`${code_exec}`;
	($syst_stat,$n_matches)=extract_text_from_array("system_status",@code_output);
	if ($syst_stat ne 0) {
	    print $log_file "ERROR!\ncode_exec=$code_exec\nsystem_status=$syst_stat\nexiting...\n";
	    exit;
	}
	($n_files,$n_matches)=extract_text_from_array("n_files",@code_output);
	($n_motifs_processed,$n_matches)=extract_text_from_array("n_processed",@code_output);
	($n_motifs_output,$n_matches)=extract_text_from_array("n_motifs",@code_output);
	($n_filtered_out,$n_matches)=extract_text_from_array("n_filtered_out",@code_output);
	($n_compared_out,$n_matches)=extract_text_from_array("n_compared_out",@code_output);
	print $log_file "\tn_files=$n_files\n";
	print $log_file "\tn_motifs_procesed=$n_motifs_processed\n";
	print $log_file "\tn_motifs_output=$n_motifs_output\n";
	print $log_file "\tn_filtered_out=$n_filtered_out\n";
	print $log_file "\tn_compared_out=$n_compared_out\n";

	($syst_stat,$warning_msg)=my_rename("motif_sampler.out.txt",$motifsampler_motifs_filename);
	print $log_file "\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("motif_sampler.log.txt","${motifsampler_dir}/${exp_id}.${species12}.motifsampler.log.txt");
	print $log_file "\t$warning_msg\n";
    }
}

##########################################################
# merge with transfac database if n_transfac_motifs >= 0 #
##########################################################
$motifs_filename="${analysis_dir}/${exp_id}.${species12}.txt";
if (-e $motifs_filename) {
    print $log_file "already merged motifs [$motifs_filename]; skipping...\n";
} else {
    if ($n_transfac_motifs>=0) {
	$code_exec="cp ${DATA_DIR}/db/transfac/transfac_db6.0_${species12}_wm.txt $motifs_filename";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";
    } 
    $code_exec="perl ${CODE_DIR}/perl/appendfiles.pl ${motifs_filename}";
    if ($process_alignace eq 1) {
	$code_exec.=" ${alignace_motifs_filename}";
    } 
    if ($process_motifsampler eq 1) {
	$code_exec.=" ${motifsampler_motifs_filename}";
    }
    if ($additional_motifs ne "none") {
	$code_exec.=" ${additional_motifs}";
    }
    if ( ($process_alignace eq 1) | ($process_motifsampler eq 1) ) {
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";
    }
}

###############################
# generate non-redundant list #
###############################
# 03-04-2004: call list_nonredundant_wm_v2.pl
$motifs_nonredundant_filename="${analysis_dir}/${exp_id}.${species12}.${motif_comparison_threshold}.txt";
if (-e $motifs_nonredundant_filename) {
    print $log_file "already computed non-redundant list of motifs [$motifs_nonredundant_filename]; skipping...\n";
} else {
    print $log_file "$motifs_nonredundant_filename does not exist\n";
    print_title("computing non-redundant list of motifs...",$log_file); 
    #$code_exec="perl ${CODE_DIR}/perl/list_nonredundant_wm.pl ${motifs_filename} ${motif_comparison_threshold}";
    $code_exec="perl ${CODE_DIR}/perl/list_nonredundant_wm_v2.pl ${motifs_filename} ${motif_comparison_threshold} ${information_threshold} ${seq_variety_threshold} ${n_seqs_threshold} ${min_motif_length_threshold} ${max_motif_length_threshold}";
    print $log_file "$code_exec\n";
    print "computing non-redundant list...\n";
    @code_output=`$code_exec`;
    ($syst_stat,$n_matches)=extract_text_from_array("system_status",@code_output);
    if ($syst_stat ne 0) {
	print $log_file "ERROR!\ncode_exec=$code_exec\nsystem_status=$syst_stat\nexiting...\n";
	exit;
    }
    ($processed_motifs,$n_matches)=extract_text_from_array("running_n_motifs",@code_output);
    ($n_motifs,$n_matches)=extract_text_from_array("running_incorporated_motifs",@code_output);
    ($running_filtered_out,$n_matches)=extract_text_from_array("running_filtered_out",@code_output);
    ($running_compared_out,$n_matches)=extract_text_from_array("running_compared_out",@code_output);
    print $log_file "\tprocessed motifs=$processed_motifs\n";
    print $log_file "\trunning_filtered_out=$running_filtered_out\n";
    print $log_file "\trunning_compared_out=$running_compared_out\n";
    print $log_file "\tn_motifs=$n_motifs\n";
}
# get the number of motifs
($n,$n_comments)=line_count($motifs_nonredundant_filename);
$n_motifs=$n_comments;      # just count the number of ">"
$n5=$n/5;
if ( ($n5 ne $n_motifs) | ($n_motifs <= 0) ) {
    print "error!\nmotif_cluster_search_v1.pl\ngenerate non-redundant list\nn_motifs=$n_motifs\tn5=$n5\tn=$n\n";
    exit;
}
print $log_file "number of motifs=$n_motifs\n";
$consensus_filename="${analysis_dir}/${exp_id}.${species12}.${motif_comparison_threshold}.consensus.txt";
if (-e $consensus_filename) {
} else {
    $code_exec="grep \"\>\" ${motifs_nonredundant_filename} \> ${consensus_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
    $code_exec="perl -pi.bak -e \'s/\>//' ${consensus_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
}
$matlab_filename="${analysis_dir}/${exp_id}.${species12}.${motif_comparison_threshold}.matlab.txt";
if (-e $matlab_filename) {
} else {
    $code_exec="cp ${motifs_nonredundant_filename} ${matlab_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
    $code_exec="perl -pi.bak -e \'s/\>/\%/' ${matlab_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
}

###################
# background list #
################### 
if ($source_original =~ /bl2seq/) {
    $source="ncbi";
}
$info_rb1_sp1_filename="${upstream_dir}/${exp_id}_rb1.${species1}.${source}.list.txt";
$upstream_bckg1_sp1_filename="${upstream_dir}/${exp_id}_rb1.${species1}.${source}.fasta";
$lengths_upstream_bckg1_sp1_filename="${upstream_dir}/${exp_id}_rb1.${species1}.${source}.lengths.txt";
if ($source_original =~ /bl2seq/) {
    $upstream_bckg1_sp1_filename_scan="${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.fasta";
    $lengths_upstream_bckg1_scan_sp1_filename="${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.lengths.txt";
    $info_rb1_scan_sp1_filename="${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.list.txt";
} else {
    $upstream_bckg1_sp1_filename_scan=$upstream_bckg1_sp1_filename;
    $lengths_upstream_bckg1_scan_sp1_filename=$lengths_upstream_bckg1_sp1_filename;
    $info_rb1_scan_sp1_filename=$info_rb1_sp1_filename;
}
if (-e $info_rb1_sp1_filename) {
    print $log_file "already created the background set [$info_rb1_sp1_filename]; skipping...\n";
} else {
    print $log_file "$info_rb1_sp1_filename does not exist...\n";
    print_title("preparing background list...",$log_file);
    $extract_column=5;                     # column to extract from the exon_list_v4.locuslink file
    if ( ($species1 eq "hs") | (${species1} eq "mm") | ($species1 eq "sc") ) {
	$code_exec="perl ${CODE_DIR}/perl/generate_background_list.pl ${query_output_filename} ${DATA_DIR}/db/upstream_sequences/${species1}/lists/exon_list_v4.locuslink.${species1}.txt $nt_bckg1 $extract_column";
    } 
    if ($species1 eq "dm") {
	$code_exec="perl ${CODE_DIR}/perl/generate_background_list.pl ${query_output_filename} ${DATA_DIR}/db/ensembl/${species1}/lists/exon_list_v4.locuslink.${species1}.txt $nt_bckg1 $extract_column";
    }
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "$warning_msg\n";
    ($n_transcripts_bck,$n_comments)=line_count("info_background.txt");
    print $log_file "(preliminary) n_transcripts_bck=$n_transcripts_bck\n";
    ($syst_stat,$warning_msg)=my_rename("info_background.txt","${upstream_dir}/info_background.txt",0,1);print $log_file "$warning_msg\n";
    $code_exec="perl -pi.bak -e \'s/\>//' ${upstream_dir}/info_background.txt";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "$warning_msg\n";
    # retrieve sequences
    $code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${upstream_dir}/info_background.txt 2";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    $genelist_filename="${upstream_dir}/col2_info_background.txt";                   # list of genes whose sequences we want to extract
    $core_output_filename="${upstream_dir}/${exp_id}_rb1.${species1}.${source}";     # core for the output file names
    $print_notfound=0;                                                               # do not print out if we cannot retrieve the sequences 
    retrieve_sequences($genelist_filename,$species1,$max_ul_length,$max_exon1_length,$max_intron1_length,$core_output_filename,$print_notfound);
    # create a list of background sequences
    $code_exec="grep \"\>\" ${upstream_bckg1_sp1_filename} \> ${info_rb1_sp1_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    $code_exec="perl -pi.bak -e \'s/\>//' ${info_rb1_sp1_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    
    # move the previous background information files to the old directory
    ($syst_stat,$warning_msg)=my_rename("${upstream_dir}/info_background.txt","${old_dir}/info_background.txt",0,1);print $log_file "$warning_msg\n";
    ($syst_stat,$warning_msg)=my_rename("${upstream_dir}/col2_info_background.txt","${old_dir}/col2_info_background.txt",0,1);print $log_file "$warning_msg\n";   
    # remove the .bak files
    $n_unlinked=unlink("${upstream_dir}/info_background.bak","${upstream_dir}/${exp_id}_rb1.${species1}.list.txt.bak");
    print $log_file "unlinked $n_unlinked files\n";
}

if ($source_original =~ /bl2seq/) {
    if ( (-e "${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.list.txt") & (-e "${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.lengths.txt") ) {
	print $log_file "i already compared the sequences for the background [${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.list.txt]\n";
    } else {
	print $log_file "could not find ${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.list.txt nor ${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.lengths.txt\n";
	if (-e "${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt") {
	    print $log_file "already created ${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt skipping...\n";
	} else {
	    print $log_file "i could not find the file ${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt\n";
	    if (-e "${upstream_dir}/${exp_id}_rb1.${species1}.homol-list.txt") {
	    } else {
		open (temp_file_in,$info_rb1_sp1_filename) || die "could not open $info_rb1_sp1_filename for reading";
		$temp_filename="motif_cluster_search_v3.temp.txt";
		open (temp_file_out,">${temp_filename}") || die "could not open temp.txt for output";
		$searchpat='LocusID:(\d+)"';
		$i=0;
		while ($record=<temp_file_in>) {
		    chomp $record;
		    (@searchstuff)=$record=~/$searchpat/; 
		    if (!@searchstuff) {
			print "\nERROR!\nmotif_cluster_search_v3.pl\nretrieving locus id entries from info_rb1_sp1_filename=$info_rb1_sp1_filename\nobtained line without LocusID\n$record\n";
			exit;
		    }
		    print temp_file_out "$searchstuff[0]\n";
		    $i++;
		}
		close (temp_file_in);
		close(temp_file_out);
		print $log_file "\tprocessed $i entries from $info_rb1_sp1_filename\n";
		$code_exec="perl ${DATA_DIR}/db/locus_link/code/ll2nm2ug.pl ${temp_filename} 1 ${species1} ll";
		($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
		($syst_stat,$warning_msg)=my_rename("query.txt","${upstream_dir}/${exp_id}_rb1.${species1}.homol-list.txt",0,1);
		($syst_stat,$warning_msg)=my_rename("log.txt","${upstream_dir}/${exp_id}_rb1.${species1}.homol-list-log.txt",0,1);
	    }
	    $code_exec="perl ${CODE_DIR}/perl/extract_column.pl ${upstream_dir}/${exp_id}_rb1.${species1}.homol-list.txt 5";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "\t$warning_msg\n";
	}
	($n,$n_comments)=line_count("${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt");
	if ($n == 0) {
	    print "error! motif_cluster_search_v3.pl\n${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt is empty\nplease check and re-run\n";
	    exit;
	}
	# retrieve upstream sequences for rb1
	if (-e "${upstream_dir}/${exp_id}_rb1.${species1}-${species2}.ncbi.fasta") {
	    print $log_file "\ti already retrieved the upstream sequences for species2 [${upstream_dir}/${exp_id}_rb1.${species1}-${species2}.ncbi.fasta]\n";
	} else {
	    print $log_file "\t${upstream_dir}/${exp_id}_rb1.${species1}-${species2}.ncbi.fasta does not exist; retrieving upstream sequences...\n";

	    $genelist_filename="${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt";                   # list of genes from sp2 to extract sequences
	    $core_output_filename="${upstream_dir}/${exp_id}_rb1.${species1}-${species2}.${source}";              # core for the output file names
	    $print_notfound=1;                                                                                    # print out even if we cannot retrieve the sequences 
	    retrieve_sequences($genelist_filename,$species2,$max_ul_length,$max_exon1_length,$max_intron1_length,$core_output_filename,$print_notfound);
	}
	# sequences comparison
	if (-e "${upstream_dir}/${exp_id}_rb1.${species1}-bl2seq-${species2}.txt") {
	    print $log_file "\ti already compared the rb1 sequences for the two species [${upstream_dir}/${exp_id}_rb1.${species1}-bl2seq-${species2}.txt]\n";
	} else {
	    print $log_file "\ti could not find ${upstream_dir}/${exp_id}_rb1.${species1}-bl2seq-${species2}.txt; comparing the rb1 sequences for $species1 and $species2\n";
	    # comparison of the species1 and species2 sequences
	    $map_filename="none";$max_seq_length=50000;$expect_threshold=1;$filteron=1;
	    $gap_open_cost=1;$gap_extend_cost=2;$dropoff_value=50;$word_size=11;$mismatch_penalty=-2;$match_reward=1;$expectation_value=10;$query_strand=3;
	    print $log_file "sp1-sp2 comparison parameters:\n";
	    print $log_file "map_filename=$map_filename\n";
	    print $log_file "expect_threshold=$expect_threshold\n";
	    print $log_file "filteron=$filteron\n";
	    print $log_file "gap_open_cost=$gap_open_cost\n";
	    print $log_file "gap_extend_cost=$gap_extend_cost\n";
	    print $log_file "dropoff_value=$dropoff_value\n";
	    print $log_file "word_size=$word_size\n";
	    print $log_file "mismatch_penalty=$mismatch_penalty\n";
	    print $log_file "match_reward=$match_reward\n";
	    print $log_file "expectation_value=$expectation_value\n";
	    print $log_file "query_strand=$query_strand\n";
	    $code_exec="perl ${CODE_DIR}/perl/bl2seq_2sets_v2.pl ${upstream_bckg1_sp1_filename} ${upstream_dir}/${exp_id}_rb1.${species1}-${species2}.ncbi.fasta ${map_filename} ${max_seq_length} ${expect_threshold} ${filteron} ${gap_open_cost} ${gap_extend_cost} ${dropoff_value} ${word_size} ${mismatch_penalty} ${match_reward} ${expectation_value} ${query_strand}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "$warning_msg\n";
	    ($syst_stat,$warning_msg)=my_rename("bl2seq_results.txt","${upstream_dir}/${exp_id}_rb1.${species1}-bl2seq-${species2}.txt",0,1);
	    print $log_file "$warning_msg\n";
	}
	# prep4alignace
	if (-e "${upstream_dir}/preporth4alignace.${exp_id}_rb1.bl2seq.${species1}.log.txt") {
	    print $log_file "i already prepared the sequences [${upstream_dir}/preporth4alignace.${exp_id}_rb1.bl2seq.${species1}.log.txt]\n";
	} else {
	    print $log_file "i could not find ${upstream_dir}/preporth4alignace.${exp_id}_rb1.bl2seq.${species1}.log.txt; calling preporth4alignace";
	    # prepare orthologous sequences for motif search or for scanning
	    $preporth4alignace_version=3;
	    $code_exec="perl ${CODE_DIR}/perl/preporth4alignace_v${preporth4alignace_version}.pl ${upstream_dir}/${exp_id}_rb1.${species1}.ncbi.fasta ${upstream_dir}/${exp_id}_rb1.${species1}-bl2seq-${species2}.txt 1 nomap ${include_noorth}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "$code_exec\t$warning_msg\n";
	    ($syst_stat,$warning_msg)=my_rename("orth2.txt",$upstream_bckg1_sp1_filename_scan,0,1);          # here we do not need orth1.txt and orth3.txt
	    ($syst_stat,$warning_msg)=my_rename("preporth4alignace.log.txt","${upstream_dir}/preporth4alignace.${exp_id}_rb1.bl2seq.${species1}.log.txt",0,1);
	}
	# get the sequence lengths
	$code_exec="grep index ${upstream_dir}/preporth4alignace.${exp_id}_rb1.bl2seq.${species1}.log.txt > index_in.txt";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\>//\' index_in.txt";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";
	$code_exec="perl ${CODE_DIR}/perl/extract_column.pl index_in.txt 2";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";
	$code_exec="perl ${CODE_DIR}/perl/extract_rows_v1.pl $lengths_upstream_bckg1_sp1_filename col2_index_in.txt";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";
	($syst_stat,$warning_msg)=my_rename("${lengths_upstream_bckg1_sp1_filename}.extracted",${lengths_upstream_bckg1_scan_sp1_filename},0,1);
	# get the list of entries in the bl2seq2 file
	$code_exec="grep \"\>\" $upstream_bckg1_sp1_filename_scan \> $info_rb1_scan_sp1_filename";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";
	$code_exec="perl -pi.bak -e \'s/\>//\' ${info_rb1_scan_sp1_filename}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);print $log_file "$code_exec\t$warning_msg\n";	
	# to unlink: orth1.txt orth3.txt index_in.txt index_in.txt.bak upstream_sequences/col5_ $temp_filename
	$n_unlinked=unlink("orth1.txt","orth3.txt","index_in.txt","index_in.txt.bak","col2_index_in.txt","${upstream_dir}/col5_${exp_id}_rb1.${species1}.homol-list.txt","${upstream_dir}/${exp_id}_rb1.${species1}.bl2seq2.list.txt.bak",${temp_filename});
	print log_file "unlinked $n_unlinked files [target = 7 files]\n";
    }        # close check on whether the lengths and list files exist
}            # close check on whether source contains bl2seq

($n,$n_comments)=line_count(${info_rb1_sp1_filename});
$n_transcripts_bck=$n-$n_comments;
print $log_file "$info_rb1_sp1_filename: n_transcripts = $n_transcripts_bck\n";
($n,$n_comments)=line_count(${info_rb1_scan_sp1_filename});
$n_transcripts_bck=$n-$n_comments;
if ($n_transcripts_bck<=0) {
    print "error!\nmotif_cluster_search_v1.pl\nbackground list\ninfo_rb1_sp1_filename=$info_rb1_sp1_filename\nn=$n\tn_comments=$n_comments\nn_transcripts_bck=$n_transcripts_bck\n";
    exit;
}
print $log_file "$info_rb1_scan_sp1_filename: n_transcripts = $n_transcripts_bck\n";
print $log_file "n_transcripts_bck=${n_transcripts_bck}\n";

####################
# scan sequences   #
####################
# usage:
# perl scan_upstream_sequences_v2.pl <matrix_filename> <upstream_filename> (<threshold_column> <ul> <rc> <sortpos> <ulelil_lengths_filename>)
if ( ($species2 ne "none") & ($process_bl2seq>0) ) {
    $toprocess="bl2seq${process_bl2seq}";
} else {
    $toprocess=$source;
}
$scan_output_filename_sp1="${scan_dir}/${exp_id}.${species1}.scan.${toprocess}.${scan_threshold_column}.txt";
$scan_output_filename_strand_sp1="${scan_dir}/${exp_id}.${species1}.scan.${toprocess}.${scan_threshold_column}.strand.txt";
if ( (-e $scan_output_filename_sp1) & (-e $scan_output_filename_strand_sp1) ) {
    print $log_file "already scanned sequences\n\tpositions=$scan_output_filename_sp1\n\tstrand=$scan-output_filename_strand_sp1\nskipping...\n";
} else { 
    print $log_file "i could not find ${scan_output_filename_sp1} and/or ${scan_output_filename_strand_sp1}\n"; 
    print_title("scan_sequences",$log_file);
    $max_upstream_length=1000000;       # set to a very high value to consider the whole sequence
    $code_exec="perl ${CODE_DIR}/perl/scan_upstream_sequences_v2.pl ${motifs_nonredundant_filename} ${upstream_filename_scan_sp1} ${scan_threshold_column} ${max_upstream_length} ${rc} ${sortpos} ${lengths_upstream_scan_sp1_filename}";
    print $log_file "code_exec=$code_exec\n";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "$warning_msg\n";

    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.out.txt",${scan_output_filename_sp1});
    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.log.txt","${scan_dir}/${exp_id}.${species1}.scan.${toprocess}.${scan_threshold_column}.log.txt");
    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.strand.txt",${scan_output_filename_strand_sp1});
}
$scan_output_filename_sp1_rb1="${scan_dir}/${exp_id}_rb1.${species1}.scan.${toprocess}.${scan_threshold_column}.txt";
$scan_output_filename_strand_sp1_rb1="${scan_dir}/${exp_id}_rb1.${species1}.scan.${toprocess}.${scan_threshold_column}.strand.txt";
if ( (-e $scan_output_filename_sp1_rb1) & (-e $scan_output_filename_strand_sp1_rb1) ) {
    print $log_file "already scanned background sequences\n\tpositions=${scan_output_filename_sp1_rb1}\n\tstrands=${scan_output_filename_strand_rb1}\nskipping...\n";
} else { 
    print $log_file "i could not find the file ${scan_output_filename_sp1_rb1} and/or ${scan_output_filename_strand_sp1_rb1}\n";
    print_title("scan sequences background",$log_file);
    $max_upstream_length=1000000;       # set to a very high value to consider the whole sequence
    $code_exec="perl ${CODE_DIR}/perl/scan_upstream_sequences_v2.pl ${motifs_nonredundant_filename} ${upstream_bckg1_sp1_filename_scan} ${scan_threshold_column} ${max_upstream_length} ${rc} ${sortpos} ${lengths_upstream_bckg1_scan_sp1_filename}";
    print $log_file "code_exec=$code_exec\n";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "$warning_msg\n";

    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.out.txt",${scan_output_filename_sp1_rb1});
    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.log.txt","${scan_dir}/${exp_id}_rb1.${species1}.scan.${toprocess}.${scan_threshold_column}.log.txt");
    ($syst_stat,$warning_msg)=my_rename("scan_upstream_sequences.strand.txt",${scan_output_filename_strand_sp1_rb1});
}

##################
# co-occurrences #
##################
print_title("preliminary search for modules",$log_file);
$cooc_c_code_version=9;                                      # changed on 03-02-2004 to version 8
$report_pos=0;
$max_dist_array[1]=25;$max_dist_array[2]=50;$max_dist_array[3]=100;$max_dist_array[4]=200;
$max_seq_length_array[1]=1000;$max_seq_length_array[2]=2000;$max_seq_length_array[3]=5000;
$exp_id_bck=$exp_id;
$timestamp="${Month}${Day}.${Hour}${Minute}";
$include_exon_intron=1;
$cooc_text="cooc-${toprocess}";      # for file naming purposes
$verbose_output=0;
$pos_format=1;

# 01-30-2004: added filteron
# cooc_4m_v6 <filename_pos_cl> <filename_pos_bck> <filename_strands_cl> <filename_strands_bck> <n_motifs> <n_transcripts_cl> <n_transcripts_bck> (<max_dist> <minian_tr> <minian_p_tr> <min_dist_factor> <max_ul> <max_el> <max_il> <init_i,j,k> <finit_i,j,k> <file_label> <report_pos> <order_constraint> <neg2pos> <sortpos> <max_n_motifs> <transcript_lengths_filename_cl> <transcript_lengths_filename_bck> <filteron> <verbose>)
# cooc_4m_v7 <filename_pos_cl> <filename_pos_bck> <filename_strands_cl> <filename_strands_bck> <n_motifs> <n_transcripts_cl> <n_transcripts_bck> (<max_dist> <minian_tr> <minian_p_tr> <min_dist_factor> <max_ul> <max_el> <max_il> <init_i,j,k> <finit_i,j,k> <file_label> <report_pos> <order_constraint> <sortpos> <max_n_motifs> <transcript_lengths_filename_cl> <transcript_lengths_filename_bck> <filteron> <pos_format> <verbose>)
 
for ($index_max_seq_length=$init_max_seq_length;$index_max_seq_length<=$finit_max_seq_length;$index_max_seq_length++) {
    $max_seq_length=$max_seq_length_array[$index_max_seq_length];
    $filteron=0;
    if ($max_seq_length<$max_ul_length) {
	$filteron=1;
    }
    if ($include_exon_intron eq 1) {
	$max_upstream_length=$max_seq_length;
	$max_exon_length=$max_seq_length;
	$max_intron_length=$max_seq_length;
    } else {
	$max_upstream_length=$max_seq_length;
	$max_exon_length=0;
	$max_intron_length=0;
	$filteron=1;
    }
    for ($max_n_motifs=$init_max_n_motifs;$max_n_motifs<=$finit_max_n_motifs;$max_n_motifs++) {
	for ($index_max_dist=$init_max_dist;$index_max_dist<=$finit_max_dist;$index_max_dist++) {
	    $max_dist=$max_dist_array[$index_max_dist];	    
	    for ($order_constraint=$init_order_constraint;$order_constraint<=$finit_order_constraint;$order_constraint++) {

		$filename_new="${cooc_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${cooc_text}.out.txt";
		if (-e $filename_new) {
		    print $log_file "$filename_new already exists; skipping...\n";
		} else {
		    print $log_file "$filename_new does not exist, calling cooc_4m_v${cooc_c_code_version}\n";
		    #$code_exec="${CODE_DIR}/cpp/executables/cooc_4m_v${cooc_c_code_version}.exe ${scan_output_filename_sp1} ${scan_output_filename_sp1_rb1} ${scan_output_filename_strand_sp1} ${scan_output_filename_strand_sp1_rb1} ${n_motifs} ${n_transcripts_cl} ${n_transcripts_bck} ${max_dist} ${min_transcripts} ${p} ${min_dist_factor} ${max_upstream_length} ${max_exon_length} ${max_intron_length} 1 1 1 1 ${n_motifs} ${n_motifs} ${n_motifs} ${n_motifs} ${timestamp} ${report_pos} ${order_constraint} ${neg2pos} ${sortpos} ${max_n_motifs} ${lengths_upstream_scan_sp1_filename} ${lengths_upstream_bckg1_scan_sp1_filename} ${filteron} ${verbose_output}";
		    $code_exec="${CODE_DIR}/cpp/executables/cooc_4m_v${cooc_c_code_version}.exe ${scan_output_filename_sp1} ${scan_output_filename_sp1_rb1} ${scan_output_filename_strand_sp1} ${scan_output_filename_strand_sp1_rb1} ${n_motifs} ${n_transcripts_cl} ${n_transcripts_bck} ${max_dist} ${min_transcripts} ${p} ${min_dist_factor} ${max_upstream_length} ${max_exon_length} ${max_intron_length} 1 1 1 1 ${n_motifs} ${n_motifs} ${n_motifs} ${n_motifs} ${timestamp} ${report_pos} ${order_constraint} ${sortpos} ${max_n_motifs} ${lengths_upstream_scan_sp1_filename} ${lengths_upstream_bckg1_scan_sp1_filename} ${filteron} ${pos_format} ${verbose_output}";
		    print $log_file "code_exec=$code_exec\n";
		    ($syst_stat,$warning_msg)=run_system($code_exec,1);
		    print $log_file "$warning_msg\n";
		    
		    $filename_orig="c.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.cooc.${timestamp}.log.txt";
		    $filename_new="${cooc_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${cooc_text}.log.txt";
		    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);print log_file "$warning_msg\n";
		    print $log_file "\t$warninig_msg\n";

		    $filename_orig="c.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.cooc.${timestamp}.out.txt";
		    $filename_new="${cooc_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${cooc_text}.out.txt";
		    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);print log_file "$warning_msg\n";
		    print $log_file "\t$warninig_msg\n";		    
		    ($n,$n_comments)=line_count($filename_new);
		    $n=$n-$n_comments;
		    print $log_file "\t\t$n modules\n";

		    if ($report_pos eq 1) {
			$filename_orig="c.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.cooc.${timestamp}.pos.txt";
			$filename_new="${cooc_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${cooc_text}.pos.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);print log_file "$warning_msg\n";
		    } else {
			# remove the pos file
			$filename_orig="c.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.cooc.${timestamp}.pos.txt";
			$n_unlinked=unlink($filename_orig);
			print $log_file "\tunlinked $n_unlinked pos files\n";
		    }
		} # close check on whether current parameter combination was already ran
	    }     # close order constraint loop
	}         # close index_max_dist loop
    }             # close max_n_motifs loop
}                 # close max_seq_length loop

#########################################
# whole genome scan for selected motifs #
#########################################
$list_cooc_motifs_dir="${cooc_dir}/list_cooc_motifs";
if (-e $list_cooc_motifs_dir) {
    print $log_file "$list_cooc_motifs_dir already exists\n";
} else {
    $code_exec="mkdir ${list_cooc_motifs_dir}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
$list_cooc_motifs_filename="${list_cooc_motifs_dir}/list_cooc_motifs_${exp_id}.txt";
if (-e $list_cooc_motifs_filename) {
    print $log_file "already prepared the list of all motifs in cooc [${list_cooc_motifs_filename}]\n";
} else {
    $code_exec="perl ${CODE_DIR}/perl/list_cooc_motifs.pl ${exp_name} ${exp_id} ${init_order_constraint} ${min_transcripts} ${source}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
    ($syst_stat,$warning_msg)=my_rename("list_cooc_motifs.txt",$list_cooc_motifs_filename,0,1);     # non-verbose mode exit on error
    print $log_file "$warning_msg\n";    
}
if ( ($species12 eq "hsmm") | ($species12 eq "mmhs") ) {
    $mm2hs_map_all="${DATA_DIR}/db/upstream_sequences/hs/hs2mm_map.txt";
    #$mm2hs_map_cl="${gen_dir}/hs2mm_map.txt";             # this needs to be created NEED TO WORK ON THIS
    $mm2hs_map_cl="default_map";
    $info_col_cl_sp1=5;$info_col_cl_sp2=2;$nspecies=2;$seq_source="gbk";$read_all=0;$singlefile=0;$max_transcripts_allocate=20000;
} else {
    if ($species1 eq "sc") {
	$info_col_cl_sp1=3;$seq_source="gbk";$read_all=1;$singlefile=1;$max_transcripts_allocate=10000;
    } 
    if ($species1 eq "dm") {
	$info_col_cl_sp1=2;$seq_source="ensembl";$read_all=1;$singlefile=1;$max_transcripts_allocate=15000;
    }
    $nspecies=1;
}
($n,$n_comments)=line_count($list_cooc_motifs_filename);
$n_list_cooc=$n-$n_comments;

# scan_allseq_v3.exe <motif_filename> <species> <motif_threshold_column> <file_prefix> (<max_upstream_length> <max_exon_length> <max_intron_length> <locuslink_only> <rc> <sortpos> <min_dist_factor> <singlefile> <seq_source> <read_all> <max_n_transcripts> <init_motif> <finit_motif> <verbose> <interesting_filename>)
$scan_allseq_version=3;
# perform a check on the number of scan files to see whether we need to call scan_allseq or not
$list_scanall_filename1="list_scanall_${exp_id}_${species1}_${timestamp}.txt";
$code_exec="ls ${scanall_dir}/ms_${exp_id}_${species1}_all* > ${list_scanall_filename1}";
($syst_stat,$warning_msg)=run_system($code_exec,0);
print $log_file "\t$warning_msg\n";
($n,$n_comments)=line_count(${list_scanall_filename1});
$n=$n-$n_comments;
if ($n>=$n_list_cooc) {
    print $log_file "number of entries in ${list_scanall_filename1} = $n\tnumber of entries in $list_cooc_motifs_filename = $n_list_cooc\tno need to scan\n";
} else {
    print $log_file "number of entries in ${list_scanall_filename1} = $n\tnumber of entries in $list_cooc_motifs_filename = $n_list_cooc\twe need to scan\n";
    $code_exec="${CODE_DIR}/cpp/executables/scan_allseq_v${scan_allseq_version}.exe ${motifs_nonredundant_filename} ${species1} $scan_threshold_column ms_${exp_id} ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${locuslink_only} ${rc} ${sortpos} ${min_dist_factor} ${singlefile} ${seq_source} ${read_all} ${max_transcripts_allocate} 1 ${n_motifs} 0 ${list_cooc_motifs_filename}";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
    if ($nspecies eq 2) {
	$list_scanall_filename2="list_scanall_${exp_id}_${species2}_${timestamp}.txt";
	$code_exec="ls ${scanall_dir}/ms_${exp_id}_${species2}_all* > ${list_scanall_filename2}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";
	($n,$n_comments)=line_count(${list_scanall_filename1});
	$n=$n-$n_comments;
	if ($n>=$n_list_cooc) {
	    print $log_file "number of entries in ${list_scanall_filename2} = $n\tnumber of entries in $list_cooc_motifs_filename = $n_list_cooc\tno need to scan\n";
	} else {
	    $code_exec="${CODE_DIR}/cpp/executables/scan_allseq_v${scan_allseq_version}.exe $species2 $scan_threshold_column ms_${exp_id} ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${locusllink_only} ${rc} ${sortpos} ${min_dist_factor} ${singlefile} ${seq_source} ${read_all} 20000 1 ${n_motifs} ${list_cooc_motifs_filename}";
	    ($syst_stat,$warning_msg)=run_system($code_exec,1);
	    print $log_file "\t$warning_msg\n";
	}
	$n_unlinked=unlink($list_scanall_filename2);
	if ($n_unlinked eq 1) { 
	    print $log_file "\t\tremoved $list_scanall_filename2\n";
	}
    }
    # move files to the temp directory
    $code_exec="mv -f ms_${exp_id}*.txt $scanall_dir";
    ($syst_stat,$warning_msg)=run_system($code_exec,1);
    print $log_file "\t$warning_msg\n";
}
$n_unlinked=unlink($list_scanall_filename1);
if ($n_unlinked eq 1) {
    print $log_file "removed $list_scanall_filename1\n";
}

######################
# whole genome check #
######################
print_title("whole genome search",$log_file);
# cooc_checkresults_v12.pl <input_filename> <wm_filename> <n_species> <order_constraint> <upstream_filename_cl_sp1> (<upstream_filename_cl_sp2>) <info_filename_cl_sp1> (<info_filename_cl_sp2>) <transcript_lengths_filename_sp1_cl> (<transcript_lengths_filename_sp2_cl>) <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> (<mm2hs_map_all> <mm2hs_map_cl>) <rc> <threshold_column> <info_col_cl_sp1> (<info_col_cl_sp2>) <locuslink_only> <sortpos> <n_iter_random_overlap> <minian_transcripts> <file_label> <sp1> <sp2> <singlefile> <seq_source> <read_all> <n_transfac_motifs> <verbose>
# 01-28-2004: call cooc_checkresults_v11.pl
# 01-30-2004: call cooc_checkresults_v12.pl (with errors and warnings file)
# 02-28-2004 add singlefile, seq_source and read_all as parameters to cooc_checkresults_v12.pl
# 02-29-2004 add <n_transfac_motifs> as parameter to cooc_checkresults_v12.pl
$cooc_checkresults_version=12;
$init_i=1;
$check_text="check-${toprocess}";
$verbose_output=0;

for ($index_max_seq_length=$init_max_seq_length;$index_max_seq_length<=$finit_max_seq_length;$index_max_seq_length++) {
    $max_seq_length=$max_seq_length_array[$index_max_seq_length];
    if ($include_exon_intron eq 1) {
	$max_upstream_length=$max_seq_length;
	$max_exon_length=$max_seq_length;
	$max_intron_length=$max_seq_length;
    } else {
	$max_upstream_length=$max_seq_length;
	$max_exon_length=0;
	$max_intron_length=0;
    }	
    print $log_file "index_max_seq_length=$index_max_seq_length\tmax_seq_length=$max_seq_length\n";
    for ($max_n_motifs=$init_max_n_motifs;$max_n_motifs<=$finit_max_n_motifs;$max_n_motifs++) {
	print $log_file "\tmax_n_motifs=$max_n_motifs\n";
	for ($index_max_dist=$init_max_dist;$index_max_dist<=$finit_max_dist;$index_max_dist++) {
	    $max_dist=$max_dist_array[$index_max_dist];
	    print $log_file "\t\tindex_max_dist=$index_max_dist\tmax_dist=$max_dist\n";
	    for ($order_constraint=$init_order_constraint;$order_constraint<=$finit_order_constraint;$order_constraint++) {
		print $log_file "\t\t\torder_constraint=$order_constraint\n";
		
		$input_filename="${cooc_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.cooc-${toprocess}.out.txt";
		print $log_file "\t\t\tinput_filename=$input_filename\n";

		if (-e $input_filename) {
		    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.out.txt";
		    if (-e $filename_new) {
			print $log_file "$filename_new already exists; skipping...\n";
		    } else {
			if ($nspecies eq 2) {
			    $code_exec="perl ${CODE_DIR}/perl/cooc_checkresults_v${cooc_checkresults_version}.pl ${input_filename} ${motifs_nonredundant_filename} ${nspecies} ${order_constraint} ${upstream_filename_scan_sp1} ${upstream_filename_scan_sp2} ${query_output_filename} ${query_output_filename} ${lengths_upstream_scan_sp1_filename} ${lengths_upstream_scan_sp2_filename} ${max_upstream_length} ${max_exon_length} ${max_intron_length} ${max_dist} ${min_dist_factor} ${exp_id} ${init_i} ${mm2hs_map_all} ${mm2hs_map_cl} ${rc} ${scan_threshold_column} ${info_col_cl_sp1} ${info_col_cl_sp2} ${locuslink_only} ${sortpos} ${n_iter_random_overlap} ${min_transcripts} ${timestamp} ${species1} ${species2} ${singlefile} ${seq_source} ${read_all} ${n_transfac_motifs} ${verbose}";
			} else {			    
			    #print "scan_threshold_column=$scan_threshold_column\n";
			    #print "info_col_cl_sp1=$info_col_cl_sp1\n";
			    #print "locuslink_only=$locuslink_only\n";
			    $code_exec="perl ${CODE_DIR}/perl/cooc_checkresults_v${cooc_checkresults_version}.pl ${input_filename} ${motifs_nonredundant_filename} ${nspecies} ${order_constraint} ${upstream_filename_scan_sp1} ${query_output_filename} ${lengths_upstream_scan_sp1_filename} ${max_upstream_length} ${max_exon_length} ${max_intron_length} ${max_dist} ${min_dist_factor} ${exp_id} ${init_i} ${rc} ${scan_threshold_column} ${info_col_cl_sp1} ${locuslink_only} ${sortpos} ${n_iter_random_overlap} ${min_transcripts} ${timestamp} ${species1} ${species2} ${verbose}";
			}			
			
			print $log_file "code_exec=$code_exec\n";		    
			($syst_stat,$warning_msg)=run_system($code_exec,1);
			print $log_file "$warning_msg\n";

			$filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.log.txt";
			$filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.log.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);     # non-verbose mode exit on error
			print $log_file "$warning_msg\n";
			
			$filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.err.txt";
			$filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.err.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);     # non-verbose mode exit on error
			print $log_file "$warning_msg\n";

			$filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.wrn.txt";
			$filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.wrn.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);     # non-verbose mode exit on error
			print $log_file "$warning_msg\n";

			$filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.out.txt";
			$filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.out.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			print $log_file "$warning_msg\n";
			
			$filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.rej.txt";
			$filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.rej.txt";
			($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			print $log_file "$warning_msg\n";
			
			if ($nspecies eq 2) {
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.mm.pos-cl.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.mm-pos-cl.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.hs.pos-cl.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.hs-pos-cl.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.mm.pos-all.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.mm-pos-all.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.hs.pos-all.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.hs-pos-all.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species12}.trov-all.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.trov-all.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species12}.posov-all.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.posov-all.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";

			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species12}.trov-cl.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.trov-cl.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";

			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species12}.posov-cl.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.posov-cl.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			} else {         # else to check on the number of species
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species1}.pos-cl.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.pos-cl.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			    
			    $filename_orig="cooc_checkresults_v${cooc_checkresults_version}.${timestamp}.${species1}.pos-all.txt";
			    $filename_new="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.pos-all.txt";
			    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
			    print $log_file "$warning_msg\n";
			}     # close check on the number of species
		    }         # close check on whether we had already ran cooc_checkresults for this set of parameters
		} else {
		    print $log_file "\t\t\t$input_filename does not exist. skipping...\n";
		}             # close check on existance of input_filename
	    }                 # close loop on order constraint
	}                     # close loop on max_dist
    }                         # close loop on max_n_motifs
}                             # close loop on index_max_seq_length


#########################
# compare check outputs #
#########################
print_title("comparing check outputs...",$log_file);
# compare_cooc_outputs_v2 <n_files> <max_n_motifs> <n_motifs> <filename_1> <filename_2> ...  <filename_n>
$c_code_version=2;
for ($max_n_motifs=$init_max_n_motifs;$max_n_motifs<=$finit_max_n_motifs;$max_n_motifs++) {
    $n_files=0;                   # number of files to compare
    $allfiles="";                 # list of files to compare 
    $run_comparison=1;            # if any file is missing, set this to 0 and do not run the comparison
    for ($index_max_dist=$init_max_dist;$index_max_dist<=$finit_max_dist;$index_max_dist++) {
	$max_dist=$max_dist_array[$index_max_dist];
	for ($order_constraint=$init_order_constraint;$order_constraint<=$finit_order_constraint;$order_constraint++) {
	    $curr_filename="${check_dir}/${exp_id}.${max_dist}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.${check_text}.out.txt";
	    if (-e $curr_filename) {
		$n_files++;
		$allfiles.="${curr_filename} ";
	    } else {
		$run_comparison=0;
		print $log_file "$curr_filename is missing; skipping comparison...\n";
	    }
	}
    }
    if ($run_comparison) {
	$code_exec="${CODE_DIR}/cpp/executables/compare_cooc_outputs_v${c_code_version}.exe ${n_files} ${max_n_motifs} ${n_motifs} ${allfiles}";
	($syst_stat,$warning_msg)=run_system($code_exec,1);
	print $log_file "\t$warning_msg\n";

	$filename_orig="compare_cooc_outputs.log.txt";
	$filename_new="${check_dir}/compare_cooc_outputs.${exp_id}.${max_n_motifs}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.log.txt";
	($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
	print $log_file "$warning_msg\n";

	$filename_orig="compare_cooc_outputs.mul.txt";
	if (-e $filename_orig) {
	    $filename_new="${check_dir}/compare_cooc_outputs.${exp_id}.${max_n_motifs}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.mul.txt";
	    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
	    print $log_file "$warning_msg\n";
	}
	$filename_orig="compare_cooc_outputs.uni.txt";
	if (-e $filename_orig) {
	    $filename_new="${check_dir}/compare_cooc_outputs.${exp_id}.${max_n_motifs}.${scan_threshold_column}.${max_upstream_length}.${max_exon_length}.${max_intron_length}.uni.txt";
	    ($syst_stat,$warning_msg)=my_rename($filename_orig,$filename_new,0,1);
	    print $log_file "$warning_msg\n";
	}
    }
}

print_title("COLORIN COLORADO, ESTE CUENTO SE HA ACABADO",$log_file);

close ($log_file);

###########
# history #
###########
# OLD cooc_checkresults_v11.pl <input_filename> <wm_filename> <n_species> <order_constraint> <upstream_filename_cl_sp1> (<upstream_filename_cl_sp2>) <info_filename_cl_sp1> (<info_filename_cl_sp2>) <transcript_lengths_filename_sp1_cl> (<transcript_lengths_filename_sp2_cl>) <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> (<mm2hs_map_all> <mm2hs_map_cl>) <rc> <threshold_column> <info_col_cl_sp1> (<info_col_cl_sp2>) <locuslink_only> <sortpos> <n_iter_random_overlap> <minian_transcripts> <file_label> <sp1> <sp2> <verbose>)
# OLD cooc_checkresults_v9.pl <input_filename> <wm_filename> <order_constraint> <upstream_filename_cl_hs> <upstream_filename_cl_mm> <info_filename_cl_hs> <info_filename_cl_mm> <transcript_lengths_filename_hs_cl> <transcript_lengths_filename_mm_cl> <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> <mm2hs_map_all> <mm2hs_map_cl> <rc> <threshold_column> <info_col_cl_hs> <info_col_cl_mm> <locuslink_only> <sortpos> <n_iter_random_overlap> <file_label> <verbose>)
# OLD cooc_checkresults_v10.pl <input_filename> <wm_filename> <order_constraint> <upstream_filename_cl> <info_filename_cl> <transcript_lengths_filename_cl> <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> <rc> <threshold_column> <info_col_cl> <locuslink_only> <sortpos> <minian_transcripts> <file_label> <verbose>)
# OLD cooc_checkresults_v12.pl <input_filename> <wm_filename> <n_species> <order_constraint> <upstream_filename_cl_sp1> (<upstream_filename_cl_sp2>) <info_filename_cl_sp1> (<info_filename_cl_sp2>) <transcript_lengths_filename_sp1_cl> (<transcript_lengths_filename_sp2_cl>) <max_upstream_length> <max_exon_length> <max_intron_length> <max_dist> <min_dist_factor> (<exp_id> <init_i> (<mm2hs_map_all> <mm2hs_map_cl>) <rc> <threshold_column> <info_col_cl_sp1> (<info_col_cl_sp2>) <locuslink_only> <sortpos> <n_iter_random_overlap> <minian_transcripts> <file_label> <sp1> <sp2> <verbose>)

# new in version 2
# 02-20-2004: call cooc_4m_v8.exe (allowing order_constraint=-1)
# 01-31-2004: call cooc_4m_v7.exe (converting positions to TSS format)
# 01-28-2004: call cooc_checkresults_v12.pl (a single program for 1 or 2 species)
# 02-14-2004: added motifsampler motif search algorithm

# 01-12-2004: changed cooc_c_code_version to 6
# 12-08-2003: added the co-occurrences part
# 12-05-2003: continued to work from the background list part
# 10-23-2003: renamed temporary files so that they are specific to this program; also added unlink at the end for the temporary files
