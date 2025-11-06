# *--------------------------------------------------------------------------* #
# * Configuration                                                            * #
# *--------------------------------------------------------------------------* #
include: "utils.smk"


configfile: "config/config.yaml"


pepfile: "config/pep/config.yaml"


if config["dir_run"] and config["dir_run"] is not None:

    workdir: config["dir_run"]


# *--------------------------------------------------------------------------* #
# * Constants                                                                * #
# *--------------------------------------------------------------------------* #
CALLERS_SV = ["gridss", "manta", "svaba", "tiddit", "wham"]
CALLERS_SV.sort()
CALLERS_MUTATION = ["mutect2"]
ANNOTATORS = ["vep", "snpeff", "annovar", "annotsv"]
MUTATIONS = ["snvs", "indels"]
CALLER2FMTS = {
    "gridss": ["GT", "AF"],
    "manta": ["PR", "SR"],
    "svaba": ["GT", "LO", "DR", "SR", "PL", "LR", "GQ", "DP", "AD"],
    "tiddit": ["GT", "CN", "COV", "DV", "RV", "LQ", "RR", "DR"],
    "wham": ["GT", "DP", "SP"],
    "mutect2": ["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"],
}
FIELDS_COMMON = (
    "CHROM POS ID REF ALT QUAL "
    "ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID "
    "ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C "
    "ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS "
    "ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS "
    "LOF[*].GENE LOF[*].GENEID LOF[*].NUMTR LOF[*].PERC "
    "NMD[*].GENE NMD[*].GENEID NMD[*].NUMTR NMD[*].PERC"
)
PROTOCOLS_UCSC = ["cytoBand"]
TYPES_SV = ["DEL", "INS", "DUP", "INV", "BND"]
MAPPER = config["mapper"]
SUFFIX_READ_1, SUFFIX_READ_2 = config["suffixes_fastq"]
SPECIES = config["species"]
GENOME2 = convert_genome(config["genome"])
DF_SAMPLE = pep.sample_table
SAMPLES = DF_SAMPLE["sample_name"]
DIR_DATA = config["dir_data"]
TO_CLEAN_FQ = config["clean_fq"]
TO_RUN_FASTQC = config["run_fastqc"]
TO_RUN_MULTIQC = config["run_multiqc"]
TO_CALL_MUTATIONS = config["mutation"]


# *--------------------------------------------------------------------------* #
# * Wildcard constraints                                                     * #
# *--------------------------------------------------------------------------* #
wildcard_constraints:
    sample=r"|".join(SAMPLES),
    caller_sv=r"|".join(CALLERS_SV),
    caller_mutation=r"|".join(CALLERS_MUTATION),
    caller=r"|".join(CALLERS_SV + CALLERS_MUTATION),
    type_sv=r"|".join(TYPES_SV),
    annotator=r"|".join(ANNOTATORS),
    mutation=r"|".join(MUTATIONS),


# *--------------------------------------------------------------------------* #
# * Files and directories required by a few rules                            * #
# *--------------------------------------------------------------------------* #
files_bwamem2_index = multiext(
    config["index_bwamem2"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
)
files_bwa_index = multiext(config["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa")

fai_fasta = f"{config['fasta']}.fai"

path_cache_snpeff = (
    f"{config['cache_snpeff']}/{config['genome']}.{config['version_snpeff']}"
)
path_cache_vep = (
    f"{config['cache_vep']}/{config['species']}/{config['version_vep']}_{config['genome']}"
)


# *--------------------------------------------------------------------------* #
# * Additional validation for config parameters                              * #
# *--------------------------------------------------------------------------* #
perform_validations_with_rich(
    config,
    workflow.source_path("../envs/vep.yaml"),
    workflow.source_path("../envs/annotsv.yaml"),
    [
        "fasta",
        "polymorphism_known",
        "dbsnp",
        "bed_exclude",
        "jar_gridss",
    ],
)
