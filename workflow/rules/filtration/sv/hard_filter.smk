rule hard_filter_svs:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        vcf="{caller_sv}/{sample}/{sample}.duphold.vcf",
    output:
        vcf="{caller_sv}/{sample}/{sample}.duphold.filtered.vcf",
    params:
        min_size=config["min_size"],
        min_reads=config["min_reads"],
        min_coverage=config["min_coverage"],
        min_dhffc=config["min_dhffc"],
        max_dhbfc=config["max_dhbfc"],
    log:
        "logs/{sample}/hard_filter_svs.{caller_sv}.log",
    script:
        "../../../scripts/hard_filter.sh"
