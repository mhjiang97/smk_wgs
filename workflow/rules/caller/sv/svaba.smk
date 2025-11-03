rule svaba:
    conda:
        "../../../envs/svaba.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        fasta=config["fasta"],
        dbsnp=config["dbsnp"],
    output:
        vcf=protected("svaba/{sample}/{sample}.svaba.sv.vcf"),
    params:
        dir="svaba/{sample}",
        min_reads=config["min_reads"],
    threads: 1
    log:
        "logs/{sample}/svaba.log",
    shell:
        """
        {{ bam=$(realpath {input.bam})
        fasta=$(realpath {input.fasta})
        dbsnp=$(realpath {input.dbsnp})

        cd {params.dir} || exit 1

        svaba run \\
            -p {threads} \\
            -t ${{bam}} \\
            -G ${{fasta}} \\
            -L {params.min_reads} \\
            -D ${{dbsnp}} \\
            -a {wildcards.sample}; }} \\
        1> {log} 2>&1
        """


rule format_svaba:
    conda:
        "../../../envs/r.yaml"
    input:
        vcf="svaba/{sample}/{sample}.svaba.sv.vcf",
    output:
        vcf="svaba/{sample}/{sample}.vcf",
        rdata="svaba/{sample}/{sample}.RData",
    params:
        genome=config["genome"],
        libs_r=config.get("libs_r", None),
    log:
        "logs/{sample}/format_svaba.log",
    script:
        "../../scripts/add_svtype_svlen.R"
