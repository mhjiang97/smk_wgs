rule manta_config:
    input:
        fasta=config["fasta"],
    output:
        config="manta/config.ini",
    params:
        min_reads=config["min_reads"],
    log:
        "logs/manta_config.log",
    shell:
        """
        {{ echo "[manta]"
        echo "referenceFasta = {input.fasta}"
        echo "minCandidateVariantSize = 8"
        echo "rnaMinCandidateVariantSize = 1000"
        echo "minEdgeObservations = 3"
        echo "graphNodeMaxEdgeCount = 10"
        echo "minCandidateSpanningCount = {params.min_reads}"
        echo "minScoredVariantSize = 50"
        echo "minDiploidVariantScore = 10"
        echo "minPassDiploidVariantScore = 20"
        echo "minPassDiploidGTScore = 15"
        echo "minSomaticScore = 10"
        echo "minPassSomaticScore = 30"
        echo "enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 1"
        echo "enableRemoteReadRetrievalForInsertionsInCancerCallingModes = 0"
        echo "useOverlapPairEvidence = 0"
        echo "enableEvidenceSignalFilter = 1"; }} \\
            1> {output.config} \\
            2> {log}
        """


rule manta:
    conda:
        "../../envs/manta.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.md.recal.bam",
        fasta=config["fasta"],
        config="manta/config.ini",
    output:
        script="manta/{sample}/runWorkflow.py",
        workspace=directory("manta/{sample}/workspace"),
        vcf=protected("manta/{sample}/results/variants/tumorSV.vcf.gz"),
        vcf_renamed=protected("manta/{sample}/{sample}.vcf"),
    threads: 1
    log:
        "logs/{sample}/manta.log",
    shell:
        """
        {{ configManta.py \\
            --config {input.config} \\
            --referenceFasta {input.fasta} \\
            --tumorBam {input.bam} \\
            --runDir $(dirname {output.script})

        {output.script} -j {threads}

        gunzip -c {output.vcf} > {output.vcf_renamed}; }} \\
        1> {log} 2>&1
        """
