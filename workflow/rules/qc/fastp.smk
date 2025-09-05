rule fastp:
    conda:
        "../../envs/fastp.yaml"
    input:
        fq_1=f"{DIR_DATA}/{{sample}}{SUFFIX_READ_1}",
        fq_2=f"{DIR_DATA}/{{sample}}{SUFFIX_READ_2}",
    output:
        fq_1=f"fastp/{{sample}}/{{sample}}{SUFFIX_READ_1}",
        fq_2=f"fastp/{{sample}}/{{sample}}{SUFFIX_READ_2}",
        fq_1_unpaired="fastp/{sample}/{sample}_R1.unpaired.fq.gz",
        fq_2_unpaired="fastp/{sample}/{sample}_R2.unpaired.fq.gz",
        fq_failed="fastp/{sample}/{sample}.failed.fq.gz",
        json="fastp/{sample}/{sample}.json",
        html="fastp/{sample}/{sample}.html",
    threads: 1
    log:
        "logs/{sample}/fastp.log",
    shell:
        """
        fastp \\
            --thread {threads} \\
            --detect_adapter_for_pe \\
            --in1 {input.fq_1} --in2 {input.fq_2} \\
            --out1 {output.fq_1} --out2 {output.fq_2} \\
            --unpaired1 {output.fq_1_unpaired} --unpaired2 {output.fq_2_unpaired} \\
            --failed_out {output.fq_failed} \\
            --json {output.json} \\
            --html {output.html} \\
            1> {log} 2>&1
        """
