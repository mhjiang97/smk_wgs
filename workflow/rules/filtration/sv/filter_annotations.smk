rule filter_annotations:
    conda:
        "../../../envs/r.yaml"
    input:
        annotsv="{caller_sv}/{sample}/{sample}.annotsv.tsv",
        maf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vep.maf",
        tsv="{caller_sv}/{sample}/merged/{sample}.{type_sv}.snpeff.tsv",
        tab="survivor/{sample}/{sample}.{type_sv}.tab",
        vep="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vep.vcf",
        snpeff="{caller_sv}/{sample}/merged/{sample}.{type_sv}.snpeff.vcf",
        vcf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vcf",
    output:
        vep="{caller_sv}/{sample}/merged/filtered/{sample}.{type_sv}.vep.vcf",
        snpeff="{caller_sv}/{sample}/merged/filtered/{sample}.{type_sv}.snpeff.vcf",
        vcf="{caller_sv}/{sample}/merged/filtered/{sample}.{type_sv}.vcf",
        ids=temp("{caller_sv}/{sample}/merged/filtered/tmp.{type_sv}"),
        rdata="{caller_sv}/{sample}/merged/filtered/{sample}.{type_sv}.RData",
        table="{caller_sv}/{sample}/merged/filtered/{sample}.{type_sv}.tsv",
    params:
        libs_r=config.get("libs_r", None),
        callers=sorted(CALLERS_SV),
        terms_relative=config["terms_relative"],
    threads: 1
    log:
        "logs/{sample}/filter_annotations.{caller_sv}.{type_sv}.log",
    script:
        "../../../scripts/filter_annotations.R"
