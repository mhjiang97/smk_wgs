rule survivor_final:
    priority: 10
    conda:
        "../../../envs/survivor.yaml"
    input:
        vcfs=expand(
            "{caller}/{{sample}}/merged/filtered/{{sample}}.{{type_sv}}.vcf",
            caller=sorted(CALLERS_SV),
        ),
    output:
        list="survivor/{sample}/final/vcfs.{type_sv}.list",
        vcf_tmp=temp("survivor/{sample}/final/{sample}.{type_sv}.tmp.vcf"),
        tab=temp("survivor/{sample}/final/{sample}.{type_sv}.rename.tab"),
        vcf="survivor/{sample}/final/{sample}.{type_sv}.merged.vcf",
    params:
        min_size=config["min_size"],
        callers=sorted(CALLERS_SV),
        distance_sv=lambda wildcards: config["distance_sv"][wildcards.type_sv],
        n_callers=lambda wildcards: max(config["n_callers"][wildcards.type_sv] - 1, 1),
        consider_type=lambda wildcards: format_survivor_parameters(
            config["consider_type"][wildcards.type_sv]
        ),
        consider_strand=lambda wildcards: format_survivor_parameters(
            config["consider_strand"][wildcards.type_sv]
        ),
        estimate_distance=lambda wildcards: format_survivor_parameters(
            config["estimate_distance"][wildcards.type_sv]
        ),
    log:
        "logs/{sample}/survivor_final.{type_sv}.log",
    script:
        "../../../scripts/survivor.sh"
