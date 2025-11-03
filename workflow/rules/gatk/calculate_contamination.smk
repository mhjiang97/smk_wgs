rule get_pileup_summaries:
    conda:
        "../../envs/gatk.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        bed=config["intervals"],
        fasta=config["fasta"],
        resource_germline=config["resource_germline"],
    output:
        table="mutect2/{sample}/{sample}.pileups.table",
    params:
        args=get_extra_arguments("get_pileup_summaries"),
    log:
        "logs/{sample}/get_pileup_summaries.log",
    threads: 1
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"mutect2/{wildcards.sample}",
    shell:
        """
        gatk GetPileupSummaries \\
            {params.args} \\
            --java-options "-Xmx{resources.mem_mb}M -XX:-UsePerfData" \\
            --input {input.bam} \\
            --variant {input.resource_germline} \\
            --reference {input.fasta} \\
            --intervals {input.bed} \\
            --output {output.table} \\
            --tmp-dir {resources.tmpdir} \\
            > {log} 2>&1
        """


rule calculate_contamination:
    conda:
        "../../envs/gatk.yaml"
    input:
        table="mutect2/{sample}/{sample}.pileups.table",
    output:
        table="mutect2/{sample}/{sample}.contamination.table",
        segmentation="mutect2/{sample}/{sample}.segmentation.table",
    params:
        args=get_extra_arguments("calculate_contamination"),
    log:
        "logs/{sample}/calculate_contamination.log",
    threads: 1
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"mutect2/{wildcards.sample}",
    shell:
        """
        gatk CalculateContamination \\
            {params.args} \\
            --java-options "-Xmx{resources.mem_mb}M -XX:-UsePerfData" \\
            --input {input.table} \\
            --output {output.table} \\
            --tumor-segmentation {output.segmentation} \\
            --tmp-dir {resources.tmpdir} \\
            > {log} 2>&1
        """
