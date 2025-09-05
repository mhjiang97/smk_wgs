rule wham:
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        fasta=config["fasta"],
    output:
        vcf=protected("wham/{sample}/tmp.vcf"),
    params:
        min_quality_mapping=config["min_quality_mapping"],
        chroms=",".join([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]),
    threads: 1
    log:
        "logs/{sample}/wham.log",
    shell:
        """
        {{ export EXCLUDE={params.chroms}
        whamg \\
            -x {threads} \\
            -f {input.bam} \\
            -a {input.fasta} \\
            -m {params.min_quality_mapping} \\
            -c {params.chroms}; }} \\
        1> {output.vcf} 2> {log}
        """


rule format_wham:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="wham/{sample}/tmp.vcf",
        fai=fai_fasta,
    output:
        vcf="wham/{sample}/{sample}.vcf",
    threads: 1
    log:
        "logs/{sample}/format_wham.log",
    shell:
        """
        {{ bcftools reheader \\
            --threads {threads} \\
            --fai {input.fai} \\
            {input.vcf} \\
            | bcftools annotate \\
                --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
                - \\
                > {output.vcf}; }} \\
        1> {log} 2>&1
        """
