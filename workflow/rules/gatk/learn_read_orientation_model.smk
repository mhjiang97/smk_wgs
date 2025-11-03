rule learn_read_orientation_model:
    conda:
        "../../envs/gatk.yaml"
    input:
        f1r2="mutect2/{sample}/{sample}.f1r2.tar.gz",
    output:
        table="mutect2/{sample}/{sample}.artifactprior.tar.gz",
    params:
        args=get_extra_arguments("learn_read_orientation_model"),
    log:
        "logs/{sample}/learn_read_orientation_model.log",
    threads: 1
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"mutect2/{wildcards.sample}",
    shell:
        """
        gatk LearnReadOrientationModel \\
            {params.args} \\
            --java-options "-Xmx{resources.mem_mb}M -XX:-UsePerfData" \\
            --input {input.f1r2} \\
            --output {output.table} \\
            --tmp-dir {resources.tmpdir} \\
            > {log} 2>&1
        """
