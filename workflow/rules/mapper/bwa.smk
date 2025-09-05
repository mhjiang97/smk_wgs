rule bwa_index:
    conda:
        "../../envs/tiddit.yaml"
    input:
        fasta=config["fasta"],
    output:
        index=files_bwa_index,
    log:
        "logs/bwa_index.log",
    shell:
        """
        bwa index {input.fasta} > {log} 2>&1
        """
