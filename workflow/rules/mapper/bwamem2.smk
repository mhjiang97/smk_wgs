rule bwamem2_index:
    conda:
        "../../envs/bwamem2.yaml"
    input:
        fasta=config["fasta"],
    output:
        index=files_bwamem2_index,
    params:
        prefix=config["index_bwamem2"],
    threads: 1
    log:
        "logs/create_bwamem2_index.log",
    shell:
        """
        bwa-mem2 index -p {params.prefix} {input.fasta} 1> {log} 2>&1
        """


rule bwamem2:
    conda:
        "../../envs/bwamem2.yaml"
    input:
        unpack(get_fastq_files),
        index=ancient(files_bwamem2_index),
    output:
        bam=protected("bwamem2/{sample}/{sample}.sorted.bam"),
        bai=protected("bwamem2/{sample}/{sample}.sorted.bam.csi"),
    params:
        prefix=config["index_bwamem2"],
        read_group='"@RG\\tID:{sample}\\tSM:{sample}\\tLB:DNA\\tPL:ILLUMINA"',
        penalty_mismatch=config["penalty_mismatch"],
    threads: 1
    log:
        "logs/{sample}/bwamem2.log",
    shell:
        """
        {{ bwa-mem2 mem -t {threads} -Y -B {params.penalty_mismatch} -R {params.read_group} {params.prefix} {input.fq_1} {input.fq_2} \\
            | samtools sort -@ {threads} --write-index -o {output.bam} -; }} \\
        1> {log} 2>&1
        """
