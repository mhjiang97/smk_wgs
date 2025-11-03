rule tiddit:
    conda:
        "../../../envs/tiddit.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        fasta=config["fasta"],
        index=ancient(files_bwa_index),
    output:
        vcf=temp("tiddit/{sample}/tmp.vcf"),
        tab=protected("tiddit/{sample}/{sample}.ploidies.tab"),
        tmp=directory("tiddit/{sample}/{sample}_tiddit"),
    params:
        prefix="tiddit/{sample}/{sample}",
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
    threads: 1
    log:
        "logs/{sample}/tiddit.log",
    shell:
        """
        tiddit \\
            --sv \\
            --threads {threads} \\
            --ref {input.fasta} \\
            --bam {input.bam} \\
            -n 2 \\
            -p {params.min_reads} \\
            -r {params.min_reads} \\
            -q {params.min_quality_mapping} \\
            -o {params.prefix} \\
            1> {log} 2>&1
        """


rule format_tiddit:
    input:
        vcf="tiddit/{sample}/tmp.vcf",
    output:
        vcf=protected("tiddit/{sample}/{sample}.vcf"),
    log:
        "logs/{sample}/format_tiddit.log",
    shell:
        """
        bcftools annotate \\
            --set-id '%ID\\-%CHROM\\_%POS' \\
            {input.vcf} \\
            1> {output.vcf} \\
            2> {log}
        """
