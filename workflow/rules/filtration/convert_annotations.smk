rule convert_vep:
    conda:
        "../../envs/vcf2maf.yaml"
    input:
        vcf="{caller}/{sample}/merged/{sample}.{type_sv}.vep.vcf",
        fasta=config["fasta"],
    output:
        maf="{caller}/{sample}/merged/{sample}.{type_sv}.vep.maf",
    params:
        version=config["version_vep"],
        genome=config["genome"],
        cache=config["cache_vep"],
        species=config["species"],
    log:
        "logs/{sample}/convert_vep.{caller}.{type_sv}.log",
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


rule convert_snpeff:
    conda:
        "../../envs/snpeff.yaml"
    input:
        vcf="{caller}/{sample}/merged/{sample}.{type_sv}.snpeff.vcf",
    output:
        tsv="{caller}/{sample}/merged/{sample}.{type_sv}.snpeff.tsv",
    params:
        **FIELDS_SNPEFF,
        fields_gt=lambda wildcards: CALLER2FIELD[wildcards.caller],
    log:
        "logs/{sample}/convert_snpeff.{caller}.{type_sv}.log",
    shell:
        """
        {{ fields_standard="{params.STANDARD}"
        fields_ann=$(printf "ANN[*].%s " {params.ANN})
        fields_lof=$(printf "LOF[*].%s " {params.LOF})
        fields_nmd=$(printf "NMD[*].%s " {params.NMD})
        fields_fmt=$(printf "GEN[*].%s " {params.fields_gt})

        SnpSift extractFields \\
            -s "," \\
            -e "." \\
            {input.vcf} \\
            ${{fields_standard}} ${{fields_ann}} ${{fields_lof}} ${{fields_nmd}} ${{fields_fmt}} \\
            | sed '1s/GEN\\[\\*\\]\\.//g ; 1s/ANN\\[\\*\\]\\.//g ; 1s/\\[\\*\\]//g' \\
                > {output.tsv}; }} \\
        > {log} 2>&1
        """
