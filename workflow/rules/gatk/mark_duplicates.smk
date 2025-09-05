rule mark_duplicates:
    conda:
        "../../envs/gatk.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        bam=protected(f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.bam"),
        bai=temp(f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.bai"),
        bai_renamed=protected(f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.bam.bai"),
        metrics=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.metrics.txt",
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"{MAPPER}/{wildcards.sample}",
    log:
        "logs/{sample}/mark_duplicates.log",
    shell:
        """
        {{ gatk MarkDuplicates \\
            --java-options \"-Xmx{resources.mem_mb}M -XX:-UsePerfData\"  \\
            --INPUT {input.bam} \\
            --OUTPUT {output.bam} \\
            --REFERENCE_SEQUENCE {input.fasta} \\
            --METRICS_FILE {output.metrics} \\
            --REMOVE_DUPLICATES false \\
            --CREATE_INDEX true \\
            --ASSUME_SORT_ORDER coordinate \\
            --TMP_DIR {resources.tmpdir}

        cp {output.bai} {output.bai_renamed}
        touch {output.bai_renamed}; }} \\
        1> {log} 2>&1
        """
