rule convert_mutation_vep:
    conda:
        "../../../envs/vcf2maf.yaml"
    input:
        vcf="{caller_mutation}/{sample}/{sample}.{mutation}.vep.vcf",
        fasta=config["fasta"],
    output:
        maf="{caller_mutation}/{sample}/{sample}.{mutation}.vep.maf",
    params:
        version=config["version_vep"],
        genome=config["genome"],
        cache=config["cache_vep"],
        species=config["species"],
    log:
        "logs/{sample}/convert_mutation_vep.{caller_mutation}.{mutation}.log",
    shell:
        """
        vcf2maf.pl \\
            --input-vcf {input.vcf} \\
            --output-maf {output.maf} \\
            --ref-fasta {input.fasta} \\
            --vcf-tumor-id {wildcards.sample} \\
            --tumor-id {wildcards.sample} \\
            --ncbi-build {params.genome} \\
            --cache-version {params.version} \\
            --vep-data {params.cache} \\
            --species {params.species} \\
            --inhibit-vep \\
            1> {log} 2>&1
        """


rule convert_mutation_snpeff:
    conda:
        "../../../envs/snpeff.yaml"
    input:
        vcf="{caller_mutation}/{sample}/{sample}.{mutation}.snpeff.vcf",
    output:
        tsv="{caller_mutation}/{sample}/{sample}.{mutation}.snpeff.tsv",
    params:
        fields_common=FIELDS_COMMON,
        fields_fmt=get_convert_snpeff_arguments,
    log:
        "logs/{sample}/convert_mutation_snpeff.{caller_mutation}.{mutation}.log",
    shell:
        """
        {{ SnpSift extractFields \\
            -s ";" -e "." \\
            {input.vcf} {params.fields_common} {params.fields_fmt} \\
            | sed '1s/GEN\\[\\*\\]\\.//g ; 1s/ANN\\[\\*\\]\\.//g ; 1s/\\[\\*\\]//g' \\
            > {output.tsv}; }} \\
        > {log} 2>&1
        """
