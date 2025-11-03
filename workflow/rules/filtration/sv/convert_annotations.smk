rule convert_sv_vep:
    conda:
        "../../../envs/vcf2maf.yaml"
    input:
        vcf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vep.vcf",
        fasta=config["fasta"],
    output:
        maf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vep.maf",
    params:
        version=config["version_vep"],
        genome=config["genome"],
        cache=config["cache_vep"],
        species=config["species"],
    log:
        "logs/{sample}/convert_sv_vep.{caller_sv}.{type_sv}.log",
    shell:
        """
        vcf2maf.pl \\
            --input-vcf {input.vcf} \\
            --output-maf {output.maf} \\
            --ncbi-build {params.genome} \\
            --cache-version {params.version} \\
            --ref-fasta {input.fasta} \\
            --vcf-tumor-id {wildcards.sample} \\
            --tumor-id {wildcards.sample} \\
            --vep-data {params.cache} \\
            --species {params.species} \\
            --inhibit-vep \\
            1> {log} 2>&1
        """


rule convert_sv_snpeff:
    conda:
        "../../../envs/snpeff.yaml"
    input:
        vcf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.snpeff.vcf",
    output:
        tsv="{caller_sv}/{sample}/merged/{sample}.{type_sv}.snpeff.tsv",
    params:
        fields_common=FIELDS_COMMON,
        fields_fmt=lambda wildcards: CALLER2FMTS[wildcards.caller_sv],
    log:
        "logs/{sample}/convert_sv_snpeff.{caller_sv}.{type_sv}.log",
    shell:
        """
        {{ fields_fmt=$(printf "GEN[*].%s " {params.fields_fmt})

        SnpSift extractFields \\
            -s ";" -e "." \\
            {input.vcf} {params.fields_common} ${{fields_fmt}} \\
            | sed '1s/GEN\\[\\*\\]\\.//g ; 1s/ANN\\[\\*\\]\\.//g ; 1s/\\[\\*\\]//g' \\
                > {output.tsv}; }} \\
        > {log} 2>&1
        """
