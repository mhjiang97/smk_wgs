rule extract_mutations_annotations:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        ids_snv="{caller_mutation}/{sample}/{sample}.snvs.ids",
        ids_indel="{caller_mutation}/{sample}/{sample}.indels.ids",
        vcf_anno="{caller_mutation}/{sample}/{sample}.{annotator}.vcf",
    output:
        snvs_anno="{caller_mutation}/{sample}/{sample}.snvs.{annotator}.vcf",
        indels_anno="{caller_mutation}/{sample}/{sample}.indels.{annotator}.vcf",
    params:
        min_reads=config["min_reads"],
        min_coverage=config["min_coverage"],
    log:
        "logs/{sample}/extract_mutations_annotations.{caller_mutation}.{annotator}.log",
    shell:
        """
        {{ bcftools filter -i 'ID=@{input.ids_snv}' {input.vcf_anno} > {output.snvs_anno}
        bcftools filter -i 'ID=@{input.ids_indel}' {input.vcf_anno} > {output.indels_anno}; }} \\
        1> {log} 2>&1
        """
