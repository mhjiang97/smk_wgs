rule annotsv:
    shadow:
        "minimal"
    conda:
        "../../envs/annotsv.yaml"
    input:
        **get_annotsv_cache_outputs(),
        vcf="{caller_sv}/{sample}/{sample}.vcf",
    output:
        tsv=touch(protected("{caller_sv}/{sample}/{sample}.annotsv.tsv")),
        tsv_unannotated=touch(
            protected("{caller_sv}/{sample}/{sample}.annotsv.unannotated.tsv")
        ),
    params:
        dir="{caller_sv}/{sample}",
        dir_cache=config["cache_annotsv"],
        genome=format_genome(config["genome"]),
        max_size=config["max_size_annotsv"],
        size_chunk=config["size_chunk"],
    log:
        "logs/{sample}/annotsv.{caller_sv}.log",
    shell:
        """
        {{ if [ {wildcards.caller_sv} == "wham" ]; then
            input={input.vcf}

            bcftools filter -i "INFO/SVTYPE != 'BND' & ABS(INFO/SVLEN) > {params.max_size}" {input.vcf} > ${{input%.*}}.long.vcf
            bcftools filter -e "INFO/SVTYPE != 'BND' & ABS(INFO/SVLEN) > {params.max_size}" {input.vcf} > ${{input%.*}}.short.vcf

            bcftools view -h ${{input%.*}}.long.vcf > {params.dir}/header.vcf

            vars_per_file={params.size_chunk}

            bcftools view -H ${{input%.*}}.long.vcf | \\
                split -l ${{vars_per_file}} - {params.dir}/chunk_

            chunk_num=1
            for chunk_file in {params.dir}/chunk_*; do
                cat {params.dir}/header.vcf ${{chunk_file}} > ${{input%.*}}.long.${{chunk_num}}.vcf
                ((chunk_num++))
            done

            n_files=$(find ${{input%.*}}.long.*.vcf | wc -l)

            for i in $(seq 1 ${{n_files}}); do

                input_long=${{input%.*}}.long.${{i}}.vcf
                output_long=${{input%.*}}.long.${{i}}.annotsv.tsv

                AnnotSV \\
                    -genomeBuild {params.genome} \\
                    -annotationsDir {params.dir_cache} \\
                    -SVinputFile ${{input_long}} \\
                    -outputFile ${{output_long}} \\
                    -SVminSize 1 \\
                    -overwrite 1 || \\
                    {{ echo "AnnotSV failed for chunk ${{i}}, skipping..."
                    touch ${{output_long}}
                    continue; }}
            done

            AnnotSV \\
                -genomeBuild {params.genome} \\
                -annotationsDir {params.dir_cache} \\
                -SVinputFile ${{input%.*}}.short.vcf \\
                -outputFile ${{input%.*}}.short.annotsv.tsv \\
                -SVminSize 1 \\
                -overwrite 1

            mv ${{input%.*}}.short.annotsv.tsv {output.tsv}
            for f in ${{input%.*}}.long.*.annotsv.tsv; do
                tail -n+2 ${{f}} >> {output.tsv}
            done

            if {{ ls ${{input%.*}}.long.*.annotsv.unannotated.tsv 1> /dev/null 2>&1; }} \\
                || {{ ls ${{input%.*}}.short.annotsv.unannotated.tsv 1> /dev/null 2>&1; }} ; then
                find {params.dir} \\
                    -name "{wildcards.sample}.long.*.annotsv.unannotated.tsv" \\
                    -o \\
                    -name "{wildcards.sample}.short.annotsv.unannotated.tsv" \\
                    2>/dev/null \\
                    | xargs cat \\
                    > {output.tsv_unannotated}
            fi

        else
            AnnotSV \\
                -genomeBuild {params.genome} \\
                -annotationsDir {params.dir_cache} \\
                -SVinputFile {input.vcf} \\
                -outputFile {output.tsv} \\
                -SVminSize 1 \\
                -overwrite 1
        fi; }} \\
        1> {log} 2>&1
        """
