rule create_seq_dict:
    input:
        ref = config["reference_w_ctrl"]
    output:
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    log:
        seq_dict_log = config['resultsdir'] + "/logs/0_prepare_files/SeqDictCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard CreateSequenceDictionary \
            -R {input.ref} \
            -O {output.seq_dict} \
        2> {log.seq_dict_log}
        """

rule create_all_target_interval_list:
    input:
        all_targets_bed = config["targetregions_all"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_all_target_interval_list
    log:
        target_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/AllTargetIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.all_targets_bed} \
            -O {output_file_all_target_interval_list} \
            -SD {input.seq_dict} \
        2> {log.target_interval_creation_log}
        """

rule create_meth_target_interval_list:
    input:
        meth_targets_bed = config["targetregions_meth"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_meth_target_interval_list
    log:
        target_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/MethTargetIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.meth_targets_bed} \
            -O {output_file_meth_target_interval_list} \
            -SD {input.seq_dict} \
        2> {log.target_interval_creation_log}
        """

rule create_snp_target_interval_list:
    input:
        snp_targets_bed = config["targetregions_snp"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_snp_target_interval_list
    log:
        target_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/SNPTargetIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.snp_targets_bed} \
            -O {output_file_snp_target_interval_list} \
            -SD {input.seq_dict} \
        2> {log.target_interval_creation_log}
        """

rule create_all_bait_interval_list:
    input:
        all_bait_bed = config["bait_regions_all"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_all_bait_interval_list
    log:
        bait_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/AllBaitIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.all_bait_bed} \
            -O {output_file_all_bait_interval_list} \
            -SD {input.seq_dict} \
        2> {log.bait_interval_creation_log}
        """

rule create_meth_bait_interval_list:
    input:
        meth_bait_bed = config["bait_regions_meth"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_meth_bait_interval_list
    log:
        bait_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/MethBaitIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.meth_bait_bed} \
            -O {output_file_meth_bait_interval_list} \
            -SD {input.seq_dict} \
        2> {log.bait_interval_creation_log}
        """

rule create_snp_bait_interval_list:
    input:
        snp_bait_bed = config["bait_regions_snp"],
        seq_dict = config["sequence_dictionary_output_dir"] + "/" + config["genome_version"] + ".dict"
    output:
        output_file_snp_bait_interval_list
    log:
        bait_interval_creation_log = config['resultsdir'] + "/logs/0_prepare_files/SNPBaitIntervalListCreation.log"
    conda:
        "envs/twist_target.yaml"
    shell:
        """
        picard BedToIntervalList \
            -I {input.snp_bait_bed} \
            -O {output_file_snp_bait_interval_list} \
            -SD {input.seq_dict} \
        2> {log.bait_interval_creation_log}
        """