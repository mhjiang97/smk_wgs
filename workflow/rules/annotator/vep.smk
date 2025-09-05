rule vep:
    conda:
        "../../envs/vep.yaml"
    input:
        vcf="{caller}/{sample}/{sample}.vcf",
        fasta=config["fasta"],
        dir_cache=path_cache_vep,
    output:
        vcf=protected("{caller}/{sample}/{sample}.vep.vcf"),
        html="{caller}/{sample}/{sample}.vep.html",
    params:
        cache=config["cache_vep"],
        version=config["version_vep"],
        genome=config["genome"],
        species=config["species"],
        max_size=config["max_size_vep"],
    threads: 1
    log:
        "logs/{sample}/vep.{caller}.log",
    shell:
        """
        vep \\
            -i {input.vcf} -o {output.vcf} --stats_file {output.html} \\
            --species {params.species} --assembly {params.genome} \\
            --cache_version {params.version} --fasta {input.fasta} \\
            --dir_cache {params.cache} --fork {threads} \\
            --force_overwrite --cache --vcf --everything --filter_common \\
            --per_gene --total_length --offline --format vcf --dont_skip \\
            --max_sv_size {params.max_size} \\
            1> {log} 2>&1
        """
