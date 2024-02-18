rule pileup:
    input:
        bam_cons_sort = config['resultsdir'] + "/results/5_dedup_consensus/{sample}.cons.sort.bam"
    output:
        pileup = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}_cons.pileup.tsv"
    params:
        snp_bed = config["targetregions_snp"],
        reference = config['reference_wo_ctrl']
    log:
        pileup_cons_log = config['resultsdir'] + "/logs/8_dd-cfDNA/{sample}.pileup_cons.log"
    conda:
         "../envs/twist_target.yaml"
    shell:
        """
        samtools mpileup \
            -l {params.snp_bed} \
            --fasta-ref {params.reference} \
            -o /dev/stdout \
            {input.bam_cons_sort} \
        2> {log.pileup_cons_log} | \
        python3 /data/users/nkueng/bioinformatic_analysis_twist/workflow/scripts/parse_pileup_file.py \
            | tee {output.pileup} 
        """

rule pileup_dedup:
    input:
        bam_dedup = config['resultsdir'] + "/results/5_dedup_consensus/{sample}_filt_dedup_umitools.bam"
    output:
        pileup = config['resultsdir'] + "/results/8_dd-cfDNA/{sample}_dedup.pileup.tsv"
    params:
        snp_bed = config["targetregions_snp"],
        reference = config['reference_wo_ctrl']
    log:
        pileup_dedup_log = config['resultsdir'] + "/logs/8_dd-cfDNA/{sample}.pileup_dedup.log"
    conda:
         "../envs/twist_target.yaml"
    shell:
        """
        samtools mpileup \
            -l {params.snp_bed} \
            --fasta-ref {params.reference} \
            -o /dev/stdout \
            --min-MQ 20 \
            {input.bam_dedup} \
        2> {log.pileup_dedup_log} | \
        python3 /data/users/nkueng/bioinformatic_analysis_twist/workflow/scripts/parse_pileup_file.py \
            | tee {output.pileup} 
        """