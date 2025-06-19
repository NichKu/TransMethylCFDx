rule index_ref_with_ctrl:
    input:
        ref_w_ctrl=config["reference_w_ctrl"],
    output:
        config["reference_w_ctrl"] + ".bwameth.c2t",
        config["reference_w_ctrl"] + ".bwameth.c2t.amb",
        config["reference_w_ctrl"] + ".bwameth.c2t.ann",
        config["reference_w_ctrl"] + ".bwameth.c2t.bwt",
        config["reference_w_ctrl"] + ".bwameth.c2t.pac",
        config["reference_w_ctrl"] + ".bwameth.c2t.sa",
        config["reference_w_ctrl"] + ".fai",
    log:
        bwameth_index_log=config["resultsdir"] + "/logs/0_prepare_files/BwaIndex.log",
        faidx_log=config["resultsdir"] + "/logs/0_prepare_files/SamtoolsFaidxIndex.log",
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        bwameth.py index {input.ref_w_ctrl} 2> {log.bwameth_index_log}
        samtools faidx {input.ref_w_ctrl} 2> {log.faidx_log}
        """


rule index_ref_without_ctrl:
    input:
        ref_wo_ctrl=config["reference_wo_ctrl"],
    output:
        config["reference_wo_ctrl"] + ".bwameth.c2t",
        config["reference_wo_ctrl"] + ".bwameth.c2t.amb",
        config["reference_wo_ctrl"] + ".bwameth.c2t.ann",
        config["reference_wo_ctrl"] + ".bwameth.c2t.bwt",
        config["reference_wo_ctrl"] + ".bwameth.c2t.pac",
        config["reference_wo_ctrl"] + ".bwameth.c2t.sa",
        config["reference_wo_ctrl"] + ".fai",
    log:
        bwameth_index_log=config["resultsdir"] + "/logs/0_prepare_files/BwaIndex.log",
        faidx_log=config["resultsdir"] + "/logs/0_prepare_files/SamtoolsFaidxIndex.log",
    conda:
        "../envs/twist_target.yaml"
    shell:
        """
        bwameth.py index {input.ref_wo_ctrl} 2> {log.bwameth_index_log}
        samtools faidx {input.ref_wo_ctrl} 2> {log.faidx_log}
        """
