rule fastqc:
    conda:
        "../../envs/fastqc.yaml"
    input:
        unpack(get_fastq_files),
    output:
        dir=branch(
            TO_CLEAN_FQ,
            then=directory("fastqc/fastp/{sample}"),
            otherwise=directory("fastqc/{sample}"),
        ),
    threads: 1
    log:
        "logs/{sample}/fastqc.log",
    shell:
        """
        {{ [ -d {output.dir} ] || mkdir -p {output.dir}
        fastqc -o {output.dir} -t {threads} {input.fq_1} {input.fq_2}; }} \\
        1> {log} 2>&1
        """


rule multiqc:
    conda:
        "../../envs/multiqc.yaml"
    input:
        fastqcs=branch(
            TO_CLEAN_FQ,
            then=expand("fastqc/fastp/{sample}", sample=SAMPLES),
            otherwise=expand("fastqc/{sample}", sample=SAMPLES),
        ),
    output:
        dir=branch(
            TO_CLEAN_FQ,
            then=directory("multiqc/fastp"),
            otherwise=directory("multiqc"),
        ),
        html=branch(
            TO_CLEAN_FQ,
            then="multiqc/fastp/multiqc_report.html",
            otherwise="multiqc/multiqc_report.html",
        ),
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc -o {output.dir} --force {input.fastqcs} 1> {log} 2>&1
        """
