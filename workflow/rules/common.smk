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
CALLERS = ["gridss", "manta", "svaba", "tiddit", "wham"]
CALLERS.sort()

CALLER2FIELD = {
    "gridss": ["GT", "AF"],
    "manta": ["PR", "SR"],
    "svaba": ["GT", "LO", "DR", "SR", "PL", "LR", "GQ", "DP", "AD"],
    "tiddit": ["GT", "CN", "COV", "DV", "RV", "LQ", "RR", "DR"],
    "wham": ["GT", "DP", "SP"],
}

FIELDS_SNPEFF = {
    "ANN": [
        "ALLELE",
        "EFFECT",
        "IMPACT",
        "GENE",
        "GENEID",
        "FEATURE",
        "FEATUREID",
        "BIOTYPE",
        "RANK",
        "HGVS_C",
        "HGVS_P",
        "CDNA_POS",
        "CDNA_LEN",
        "CDS_POS",
        "CDS_LEN",
        "AA_POS",
        "AA_LEN",
        "DISTANCE",
        "ERRORS",
    ],
    "LOF": ["GENE", "GENEID", "NUMTR", "PERC"],
    "NMD": ["GENE", "GENEID", "NUMTR", "PERC"],
    "STANDARD": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
}

TYPES_SV = ["DEL", "INS", "DUP", "INV", "BND"]

MAPPER = config["mapper"]

ANNOTATORS = ["vep", "snpeff", "annotsv"]

SUFFIX_READ_1, SUFFIX_READ_2 = config["suffixes_fastq"]

SPECIES = config["species"]

DF_SAMPLE = pep.sample_table
SAMPLES = DF_SAMPLE["sample_name"]

DIR_DATA = config["dir_data"]

TO_CLEAN_FQ = config["clean_fq"]
TO_RUN_FASTQC = config["run_fastqc"]
TO_RUN_MULTIQC = config["run_multiqc"]


# *--------------------------------------------------------------------------* #
# * Wildcard constraints                                                     * #
# *--------------------------------------------------------------------------* #
wildcard_constraints:
    sample=r"|".join(SAMPLES),
    caller=r"|".join(CALLERS),
    type_sv=r"|".join(TYPES_SV),


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
