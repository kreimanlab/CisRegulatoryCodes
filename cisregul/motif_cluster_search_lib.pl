#!/usr/bin/perl
# motif_cluster_search_lib.pl
# methods used in motif_cluster_search.pl

# ($system_status,$warning_msg)=run_system($code_exec,$exit_on_error);
# print_title($text,$$log_file);

sub retrieve_sequences
{
    # retrieve_sequences
    # call retrieve_upstream_v3.pl or merge_upexin.pl to retrieve sequences
    # usage:
    # retrieve_sequences($genelist_filename,$sp,$max_ul_length,$max_exon1_length,$max_intron1_length,$core_output_filename,$print_notfound);

    my $genelist_filename=$_[0];               # file with a single column with the gene identifiers 
    my $sp=$_[1];                              # species
    my $max_ul_length=$_[2];                   # maximum upstream length
    my $max_exon1_length=$_[3];                # maximum exon1 length
    my $max_intron1_length=$_[4];              # maximum intron1 length
    my $core_output_filename=$_[5];            # core output filename
    my $print_notfound=$_[6];
    my $load2memory=0;
    my $separator="none";
    my $code_exec;
    my $system_status;

    if ($sp eq "dm") {
	$code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/merge_upexin.pl ${genelist_filename} ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_upstream.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_exon1.info.txt ${DATA_DIR}/db/ensembl/dm/lists/ensembl_dm_intron1.info.txt ${DATA_DIR}/db/ensembl/dm/dm.ensembl.upstream.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.exon1.fasta ${DATA_DIR}/db/ensembl/dm/dm.ensembl.intron1.fasta ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${print_notfound} ${load2memory} ${separator}";
	($system_status,$warning_msg)=run_system($code_exec,1);       # exit on error
	if ($system_status ne 0) {
	    print "\n\nmotif_cluster_search_lib.pl\n";
	    print "retrieve_sequences\n";
	    print "code_exec=$code_exec\n";
	    print "system_status=$system_status\n\n\n";
	    exit;
	}
	print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	($system_status,$warning_msg)=my_rename("merge_upexin.out.txt","${core_output_filename}.fasta",0,1);
	print $log_file "\t$warning_msg";
	($system_status,$warning_msg)=my_rename("log.txt","${core_output_filename}.log.txt",0,1);
	print $log_file "\t$warning_msg";
	($system_status,$warning_msg)=my_rename("merge_upexin.lengths.txt","${core_output_filename}.lengths.txt",0,1);
	print $log_file "\t$warning_msg";
    } else {
	$code_exec="perl ${DATA_DIR}/db/upstream_sequences/code/retrieve_upstream_v3.pl ${genelist_filename} ${sp} ${max_ul_length} ${max_exon1_length} ${max_intron1_length} ${print_notfound}";
	($system_status,$warning_msg)=run_system($code_exec,1);
	if ($system_status ne 0) {
	    print "\n\nmotif_cluster_search_lib.pl\n";
	    print "retrieve_sequences\n";
	    print "code_exec=$code_exec\n";
	    print "system_status=$system_status\n\n\n";
	    exit;
	}
	print $log_file "\tcode_exec=$code_exec\t$warning_msg\n";
	($system_status,$warning_msg)=my_rename("retrieve_upstream_v3.out.txt","${core_output_filename}.fasta",0,1);
	print $log_file "\t$warning_msg";
	($system_status,$warning_msg)=my_rename("retrieve_upstream_v3.log.txt","${core_output_filename}.log.txt",0,1);
	print $log_file "\t$warning_msg";
	($system_status,$warning_msg)=my_rename("retrieve_upstream_v3.lengths.txt","${core_output_filename}.lengths.txt",0,1);
	print $log_file "\t$warning_msg";
    }
}

sub read_parameters_from_file
{
    # read the parameters from a file
    # usage:
    # @params=read_parameters_from_file($parameters_filename);
    # parameters have the following format
    # parameter_name=parameter_value

    my $parameters_filename=$_[0];
    my $i;
    my @splitarray;
    my $input_filename="none";
    my $experiment_name="none";
    my $search_type="none";
    my $species1="none";
    my $species2="none";
    my $max_ul_length="none";
    my $max_exon1_length="none";
    my $max_intron1_length="none";
    my $motif_comparison_threshold="none";
    my $log_text="";
    my $rc=1;
    my $sortpos=1;
    my $scan_threshold_column=38;
    my $process_bl2seq=2;
    my $p=0.01;
    my $min_dist_factor=0.5;
    my $init_max_n_motifs=2;
    my $finit_max_n_motifs=4;
    my $init_max_dist=1;
    my $finit_max_dist=4;
    my $init_order_constraint=0;
    my $finit_order_constraint=1;
    my $init_max_seq_length=1;
    my $finit_max_seq_length=3;
    my $min_transcripts=4;
    my $locuslink_only=1;
    my $n_iter_random_overlap=10;
    my $include_noorth=1;
    my $source="ncbi";
    my $process_motifsampler=0;
    my $additional_motifs="none";
    my $n_transfac_motifs=-1;
    my $information_threshold=0.2331;              # discard if the mean information per nt is <= than this value
    my $seq_variety_threshold=1;                   # discard if the sequence variety is <= than this value
    my $n_seqs_threshold=5;                        # minimum number of sequences used to build the weight matrix
    my $min_motif_length_threshold=6;              # minimum motif length
    my $max_motif_length_threshold=13;             # maximum motif length

    open (input_file,$parameters_filename) || die "could not open $parameters_filename for reading";
    while ($record=<input_file>) {
	if ($record=~/input_filename/) {
	    $input_filename=get_param_from_record($record);
	}
	if ($record=~/experiment_name/) {
	    $experiment_name=get_param_from_record($record);
	}
	if ($record=~/experiment_id/) {
	    $experiment_id=get_param_from_record($record);
	}
	if ($record=~/species1/) {
	    $species1=get_param_from_record($record);
	}
	if ($record=~/species2/) {
	    $species2=get_param_from_record($record);
	}
	if ($record=~/search_type/) {
	    $search_type=get_param_from_record($record);         # ll, nm, ug
	}
	if ($record=~/max_ul_length/) {
	    $max_ul_length=get_param_from_record($record);
	}
	if ($record=~/max_exon1_length/) {
	    $max_exon1_length=get_param_from_record($record);
	}
	if ($record=~/max_intron1_length/) {
	    $max_intron1_length=get_param_from_record($record);
	}	
	if ($record=~/motif_comparison_threshold/) {
	    $motif_comparison_threshold=get_param_from_record($record);
	}	
	if ($record=~/^rc/) {
	    $rc=get_param_from_record($record);
	}	
	if ($record=~/sortpos/) {
	    $sortpos=get_param_from_record($record);
	}	
	if ($record=~/scan_threshold_column/) {
	    $scan_threshold_column=get_param_from_record($record);
	}	
	if ($record=~/process_bl2seq/) {
	    $process_bl2seq=get_param_from_record($record);
	}	
	if ($record=~/^p\=/) {
	    $p=get_param_from_record($record);
	}
	if ($record=~/min_dist_factor/) {
	    $min_dist_factor=get_param_from_record($record);
	}
	if ($record=~/^init_max_n_motifs/) {
	    $init_max_n_motifs=get_param_from_record($record);
	}
	if ($record=~/finit_max_n_motifs/) {
	    $finit_max_n_motifs=get_param_from_record($record);
	}
	if ($record=~/^init_max_dist/) {
	    $init_max_dist=get_param_from_record($record);
	}
	if ($record=~/finit_max_dist/) {
	    $finit_max_dist=get_param_from_record($record);
	}
	if ($record=~/^init_order_constraint/) {
	    $init_order_constraint=get_param_from_record($record);
	}
	if ($record=~/finit_order_constraint/) {
	    $finit_order_constraint=get_param_from_record($record);
	}
	if ($record=~/^init_max_seq_length/) {
	    $init_max_seq_length=get_param_from_record($record);
	}
	if ($record=~/finit_max_seq_length/) {
	    $finit_max_seq_length=get_param_from_record($record);
	}
	if ($record=~/min_transcripts/) {
	    $min_transcripts=get_param_from_record($record);
	}
	if ($record=~/locuslink_only/) {
	    $locuslink_only=get_param_from_record($record);
	}
	if ($record=~/n_iter_random_overlap/) {
	    $n_iter_random_overlap=get_param_from_record($record);
	}
	if ($record=~/include_noorth/) {
	    $include_noorth=get_param_from_record($record);
	}
	if ($record=~/source/) {
	    $source=get_param_from_record($record);
	}
	if ($record=~/process_motifsampler/) {
	    $process_motifsampler=get_param_from_record($record);
	}
	if ($record=~/additional_motifs/) {
	    $additional_motifs=get_param_from_record($record);
	}
	if ($record=~/n_transfac_motifs/) {
	    $n_transfac_motifs=get_param_from_record($record);
	}
	if ($record=~/information_threshold/) {
	    $information_threshold=get_param_from_record($record);
	}
        if ($record=~/seq_variety_threshold/) {
	    $seq_variety_threshold=get_param_from_record($record);
	}
	if ($record=~/n_seqs_threshold/) {
	    $n_seqs_threshold=get_param_from_record($record);
	}
	if ($record=~/min_motif_length_threshold/) {
	    $min_motif_length_threshold=get_param_from_record($record);
	}
	if ($record=~/max_motif_length_threshold/) {
	    $max_motif_length_threshold=get_param_from_record($record);
	}
    }

    $log_text.="% input_filename=$input_filename\n";
    $log_text.="% experiment_name=$experiment_name\n";
    $log_text.="% experiment_id=$experiment_id\n";
    $log_text.="% species1=$species1\n";
    $log_text.="% species2=$species2\n";
    $log_text.="% search_type=$search_type\n";
    $log_text.="% max_ul_length=$max_ul_length\n";
    $log_text.="% max_exon1_length=$max_exon1_length\n";
    $log_text.="% max_intron1_length=$max_intron1_length\n";
    $log_text.="% motif_comparison_threshold=$motif_comparison_threshold\n";
    $log_text.="% sortpos=$sortpos\n";
    $log_text.="% rc=$rc\n";
    $log_text.="% scan_threshold_column=$scan_threshold_column\n";
    $log_text.="% process_bl2seq=$process_bl2seq\n";
    $log_text.="% p=$p\n";
    $log_text.="% min_dist_factor=$min_dist_factor\n";
    $log_text.="% init_max_n_motifs=$init_max_n_motifs\n";
    $log_text.="% finit_max_n_motifs=$finit_max_n_motifs\n";
    $log_text.="% init_max_dist=$init_max_dist\n";
    $log_text.="% finit_max_dist=$finit_max_dist\n";
    $log_text.="% init_order_constraint=$init_order_constraint\n";
    $log_text.="% finit_order_constraint=$finit_order_constraint\n";
    $log_text.="% init_max_seq_length=$init_max_seq_length\n";
    $log_text.="% finit_max_seq_length=$finit_max_seq_length\n";
    $log_text.="% min_transcripts=$min_transcripts\n";
    $log_text.="% locuslink_only=$locuslink_only\n";
    $log_text.="% n_iter_random_overlap=$n_iter_random_overlap\n";
    $log_text.="% source=${source}\n";
    $log_text.="% include_noorth=$include_noorth\n";
    $log_text.="% additional_motifs=${additional_motifs}\n";
    $log_text.="% process_motifsampler=$process_motifsampler\n";
    $log_text.="% n_transfac_motifs=${n_transfac_motifs}\n";
    $log_text.="% information_threshold=${information_threshold}\n";
    $log_text.="% seq_variety_treshold=${seq_variety_threshold}\n";
    $log_text.="% n_seqs_threshold=${n_seqs_threshold}\n";
    $log_text.="% min_motif_length_threshold=${min_motif_length_threshold}\n";
    $log_text.="% max_motif_length_threshold=${max_motif_length_threshold}\n";

    return ($log_text,$input_filename,$experiment_name,$experiment_id,$species1,$species2,$search_type,$max_ul_length,$max_exon1_length,$max_intron1_length,$motif_comparison_threshold,$rc,$sortpos,$scan_threshold_column,$process_bl2seq,$p,$min_dist_factor,$init_max_n_motifs,$finit_max_n_motifs,$init_max_dist,$finit_max_dist,$init_order_constraint,$finit_order_constraint,$init_max_seq_length,$finit_max_seq_length,$min_transcripts,$locuslink_only,$n_iter_random_overlap,$source,$include_noorth,$process_motifsampler,$additional_motifs,$n_transfac_motifs,$information_threshold,$seq_variety_threshold,$n_seqs_threshold,$min_motif_length_threshold,$max_motif_length_threshold);
}

sub get_param_from_record{
    my $param;
    my $record=$_[0];
    my @splitarray;
    
    chomp $record;
    @splitarray=split /\=/,$record;
    $param=$splitarray[1];

    return $param;
}


sub run_system
{
    # run a system command
    # usage:
    # ($system_status,$warning_msg)=run_system($code_exec,$exit_on_error);
    my $code_exec=$_[0];       # code to run
    my $exit_on_error=$_[1];   # if exit_on_error=1, then exit if system status non-zero
    my $system_status=0;
    my $warning_msg="ok";

    $system_status=system($code_exec);
    if ($system_status ne 0) {
	if ($exit_on_error eq 1) {
	    print "error!\ncode_exec=$code_exec\nsystem_status=$system_status\n";
	    exit;
	} else {
	    $warning_msg="error!\ncode_exec=$code_exec\nsystem_status=$system_status\n";
	}
    } else {
	$warning_msg="${code_exec}\tok\n";
	print "$code_exec\n";
    }

    return ($system_status,$warning_msg);
}

sub print_title
{
    # print_title
    # usage:
    # print_title($text,$$log_file);
    
    my $i;
    my $n;
    my $text=$_[0];
    my $log_file=$_[1];
    
    $n=length($text);
    for ($i=1;$i<=$n;$i++) {
	print $log_file "*";
    }
    print $log_file "\n";
    print $log_file "$text\n";
    for ($i=1;$i<=$n;$i++) {
	print $log_file "*";
    }
    print $log_file "\n";

}    
1;         # it took me about 1 hour to realize i had to do this...
