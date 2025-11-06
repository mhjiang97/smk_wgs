rule convert_vcf_to_annovar:
    input:
        vcf="{caller}/{sample}/{sample}.vcf",
    output:
        annovar=temp("{caller}/{sample}/av.{sample}.avinput"),
    params:
        prefix="{caller}/{sample}/av",
    log:
        "logs/{sample}/convert_vcf_to_annovar.{caller}.log",
    shell:
        """
        convert2annovar.pl \\
            -format vcf4 \\
            -outfile {params.prefix} \\
            -allsample \\
            -includeinfo \\
            {input.vcf} \\
            1> {log} 2>&1
        """


rule annovar:
    input:
        unpack(get_annovar_inputs),
    output:
        tmp=temp(f"{{caller}}/{{sample}}/{{sample}}.{GENOME2}_multianno.txt"),
        annovar=protected("{caller}/{sample}/{sample}.annovar.tsv"),
    params:
        **get_annovar_arguments(),
        prefix="{caller}/{sample}/{sample}",
        cache=config["cache_annovar"],
        genome=GENOME2,
    threads: 1
    resources:
        tmpdir=lambda wildcards: f"{wildcards.caller}/{wildcards.sample}",
    log:
        "logs/{sample}/annovar.{caller}.log",
    shell:
        """
        {{ table_annovar.pl \\
            {input.annovar} \\
            {params.cache} \\
            --thread {threads} \\
            --buildver {params.genome} \\
            --outfile {params.prefix} \\
            --protocol {params.protocol} \\
            --operation {params.operation} \\
            --tempdir {resources.tmpdir} \\
            --nastring NA \\
            --polish \\
            --remove

        cp {output.tmp} {output.annovar}; }} \\
        1> {log} 2>&1
        """
