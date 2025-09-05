rule filter_annotations:
    conda:
        "../../envs/r.yaml"
    input:
        annotsv="{caller}/{sample}/{sample}.annotsv.tsv",
        maf="{caller}/{sample}/merged/{sample}.{type_sv}.vep.maf",
        tsv="{caller}/{sample}/merged/{sample}.{type_sv}.snpeff.tsv",
        tab="survivor/{sample}/{sample}.{type_sv}.tab",
        vep="{caller}/{sample}/merged/{sample}.{type_sv}.vep.vcf",
        snpeff="{caller}/{sample}/merged/{sample}.{type_sv}.snpeff.vcf",
        vcf="{caller}/{sample}/merged/{sample}.{type_sv}.vcf",
    output:
        vep="{caller}/{sample}/merged/filtered/{sample}.{type_sv}.vep.vcf",
        snpeff="{caller}/{sample}/merged/filtered/{sample}.{type_sv}.snpeff.vcf",
        vcf="{caller}/{sample}/merged/filtered/{sample}.{type_sv}.vcf",
        ids=temp("{caller}/{sample}/merged/filtered/tmp.{type_sv}"),
        rdata="{caller}/{sample}/merged/filtered/{sample}.{type_sv}.RData",
        table="{caller}/{sample}/merged/filtered/{sample}.{type_sv}.tsv",
    params:
        libs_r=config.get("libs_r", None),
        callers=sorted(CALLERS),
        terms_relative=config["terms_relative"],
    threads: 1
    log:
        "logs/{sample}/filter_annotations.{caller}.{type_sv}.log",
    script:
        "../../scripts/filter_annotations.R"
