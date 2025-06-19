rule calculate_fraction:
    input:
        expand(
            config["resultsdir"] + "/results/8_dd-cfDNA/{sample}_cons.pileup.tsv",
            sample=SAMPLES,
        ),
    output:
        expand(
            config["resultsdir"] + "/results/8_dd-cfDNA/{sample}_dd-cfDNAfraction.png",
            sample=SAMPLES,
        ),
        summary=config["resultsdir"]
        + "/results/8_dd-cfDNA/"
        + config["project_run_name"],
    shell:
        """
        
        """
