rule snpeff:
    conda:
        "../../envs/snpeff.yaml"
    input:
        vcf="{caller}/{sample}/{sample}.vcf",
        fasta=config["fasta"],
        dir_cache=path_cache_snpeff,
    output:
        vcf=protected("{caller}/{sample}/{sample}.snpeff.vcf"),
        html="{caller}/{sample}/{sample}.snpeff.html",
    params:
        cache=config["cache_snpeff"],
        version=config["version_snpeff"],
        genome=config["genome"],
    resources:
        mem_mb=1,
    log:
        "logs/{sample}/snpeff.{caller}.log",
    shell:
        """
        snpEff -Xmx{resources.mem_mb}M \\
            -nodownload -v -lof -canon \\
            -dataDir {params.cache} -s {output.html} \\
            {params.genome}.{params.version} {input.vcf} \\
            1> {output.vcf} \\
            2> {log}
        """
