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


rule extract_annovar_annotations:
    conda:
        "../../../envs/python.yaml"
    input:
        anno="{caller}/{sample}/{sample}.annovar.tsv",
        tmp="{caller}/{sample}/av.{sample}.avinput",
        vcf_snv="{caller}/{sample}/{sample}.snvs.vcf",
        vcf_indel="{caller}/{sample}/{sample}.indels.vcf",
    output:
        anno_snv="{caller}/{sample}/{sample}.snvs.annovar.tsv",
        anno_indel="{caller}/{sample}/{sample}.indels.annovar.tsv",
    log:
        "logs/{sample}/extract_annovar_annotations.{caller}.log",
    script:
        "../../../scripts/extract_annovar_annotations.py"
