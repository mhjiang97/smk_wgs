rule mutect2:
    conda:
        "../../../envs/gatk.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        bed=config["intervals"],
        fasta=config["fasta"],
        resource_germline=config["resource_germline"],
        pon=config["pon"],
    output:
        vcf="mutect2/{sample}/{sample}.raw.vcf",
        f1r2="mutect2/{sample}/{sample}.f1r2.tar.gz",
    params:
        min_coverage=config["min_coverage"],
        min_quality_base=config["min_quality_base"],
        args=get_extra_arguments("mutect2"),
    log:
        "logs/{sample}/mutect2.log",
    threads: 1
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"mutect2/{wildcards.sample}",
    shell:
        """
        gatk Mutect2 \\
            {params.args} \\
            --java-options "-Xmx{resources.mem_mb}M -XX:-UsePerfData" \\
            --native-pair-hmm-threads {threads} \\
            -R {input.fasta} \\
            -I {input.bam} \\
            -O {output.vcf} \\
            -L {input.bed} \\
            --f1r2-tar-gz {output.f1r2} \\
            --callable-depth {params.min_coverage} \\
            --min-base-quality-score {params.min_quality_base} \\
            --f1r2-min-bq 20 \\
            --germline-resource {input.resource_germline} \\
            --panel-of-normals {input.pon} \\
            --tmp-dir {resources.tmpdir} \\
            > {log} 2>&1
        """


rule filter_mutect_calls:
    conda:
        "../../../envs/gatk.yaml"
    input:
        vcf="mutect2/{sample}/{sample}.raw.vcf",
        bed=config["intervals"],
        fasta=config["fasta"],
        table_contamination="mutect2/{sample}/{sample}.contamination.table",
        table_segmentation="mutect2/{sample}/{sample}.segmentation.table",
        table_artifact="mutect2/{sample}/{sample}.artifactprior.tar.gz",
    output:
        vcf=temp("mutect2/{sample}/{sample}.filtered.vcf"),
    params:
        args=get_extra_arguments("filter_mutect_calls"),
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"mutect2/{wildcards.sample}",
    log:
        "logs/{sample}/filter_mutect_calls.log",
    shell:
        """
        gatk FilterMutectCalls \\
            {params.args} \\
            --java-options "-Xmx{resources.mem_mb}M -XX:-UsePerfData" \\
            --reference {input.fasta} \\
            --variant {input.vcf} \\
            --output {output.vcf} \\
            --intervals {input.bed} \\
            --orientation-bias-artifact-priors {input.table_artifact} \\
            --contamination-table {input.table_contamination} \\
            --tumor-segmentation {input.table_segmentation} \\
            --tmp-dir {resources.tmpdir} \\
            > {log} 2>&1
        """


rule format_mutect2:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        vcf="mutect2/{sample}/{sample}.filtered.vcf",
    output:
        vcf=protected("mutect2/{sample}/{sample}.vcf"),
    log:
        "logs/{sample}/format_mutect2.log",
    shell:
        """
        bcftools annotate \\
            --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
            {input.vcf} \\
            1> {output.vcf} 2> {log}
        """
