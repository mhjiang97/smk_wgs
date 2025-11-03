rule separate_types:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        vcf="{caller_sv}/{sample}/{sample}.duphold.filtered.vcf",
    output:
        vcf="{caller_sv}/{sample}/{sample}.{type_sv}.vcf",
    log:
        "logs/{sample}/separate_types.{caller_sv}.{type_sv}.log",
    shell:
        """
        bcftools filter \\
            -i "INFO/SVTYPE ~ '{wildcards.type_sv}'" \\
            {input.vcf} \\
            1> {output.vcf} 2> {log}
        """
