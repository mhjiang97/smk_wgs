rule extract_sv_ids:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        survivor="survivor/{sample}/{sample}.{type_sv}.merged.vcf",
    output:
        tab="survivor/{sample}/{sample}.{type_sv}.tab",
    log:
        "logs/{sample}/extract_sv_ids.{type_sv}.log",
    shell:
        """
        {{ bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%ID]' {input.survivor} \\
            | awk -F'\\t' -v OFS='\\t' '{{
                printf "%s\\t%s\\t%s\\t%s", $1, $2, $3, $4
                for (i = 5; i <= NF; i++) {{split($i, arr, "="); printf "\\t%s", arr[2]}}
                printf "\\n"
            }}' \\
            > {output.tab}; }} \\
        1> {log} 2>&1
        """


rule extract_sv_annotations:
    conda:
        "../../../envs/bcftools.yaml"
    input:
        vep="{caller_sv}/{sample}/{sample}.vep.vcf",
        snpeff="{caller_sv}/{sample}/{sample}.snpeff.vcf",
        vcf="{caller_sv}/{sample}/{sample}.{type_sv}.vcf",
        tab="survivor/{sample}/{sample}.{type_sv}.tab",
    output:
        ids=touch(temp("survivor/{sample}/{caller_sv}.{type_sv}.id")),
        vep="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vep.vcf",
        snpeff="{caller_sv}/{sample}/merged/{sample}.{type_sv}.snpeff.vcf",
        vcf="{caller_sv}/{sample}/merged/{sample}.{type_sv}.vcf",
    params:
        callers=sorted(CALLERS_SV),
    log:
        "logs/{sample}/extract_sv_annotations.{caller_sv}.{type_sv}.log",
    shell:
        """
        {{ IFS=" " read -r -a callers_unsorted <<< "{params.callers}"
        list_string=$(printf "%s\\n" "${{callers_unsorted[@]}}")
        mapfile -t callers < <(echo "${{list_string}}" | sort)

        caller_column=$(python -c "callers = '${{callers[*]}}'.split(); map_caller_column = {{callers[i]: i+5 for i in range(len(callers))}}; print(map_caller_column['{wildcards.caller}'])")

        cut -f ${{caller_column}} {input.tab} \\
            | {{ grep -v 'NaN' || true; }} \\
            > {output.ids}

        bcftools view -i "ID=@{output.ids}" {input.vep} > {output.vep}
        bcftools view -i "ID=@{output.ids}" {input.snpeff} > {output.snpeff}
        bcftools view -i "ID=@{output.ids}" {input.vcf} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
