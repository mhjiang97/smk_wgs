libs_r <- snakemake@params[["libs_r"]]
if (!is.null(libs_r)) {
  .libPaths(c(libs_r, .libPaths()))
}
snakemake@source("utils.R")
load_pkg(
  c(
    "vroom", "dplyr", "tibble", "glue", "tidyr", "purrr", "GenomicRanges",
    "stringr", "BiocParallel", "parallel", "logger", "progress"
  )
)
logger::log_appender(logger::appender_file(snakemake@log[[1]]))

logger::log_info("Reading Snakemake input, output, and parameters")
file_annotsv <- snakemake@input[["annotsv"]]
file_maf <- snakemake@input[["maf"]]
file_tsv <- snakemake@input[["tsv"]]
file_tab <- snakemake@input[["tab"]]
vep <- snakemake@input[["vep"]]
snpeff <- snakemake@input[["snpeff"]]
vcf <- snakemake@input[["vcf"]]
vep_final <- snakemake@output[["vep"]]
snpeff_final <- snakemake@output[["snpeff"]]
vcf_final <- snakemake@output[["vcf"]]
path_ids <- snakemake@output[["ids"]]
table <- snakemake@output[["table"]]
callers <- sort(snakemake@params[["callers"]])
terms_relative <- snakemake@params[["terms_relative"]]
threads <- snakemake@threads
s <- snakemake@wildcards[["sample"]]
t <- snakemake@wildcards[["caller"]]
v <- snakemake@wildcards[["type_sv"]]

logger::log_info("Reading AnnotSV annotations")
annotsv <- file_annotsv |>
  my_vroom(
    n_col = 120, quote = "\'", delim = "\t", col_names = TRUE,
    num_threads = threads
  )

logger::log_info("Reading SV ID information")
id_info <- file_tab |>
  my_vroom(col_names = FALSE, num_threads = threads) |>
  setNames(c("chr", "pos", "ref", "alt", callers))

logger::log_info("Adding overlap flags to AnnotSV annotations")
if (nrow(annotsv) != 0) {
  tsv_annotsv <- annotsv |>
    dplyr::filter(ID %in% na.omit(id_info[[t]])) |>
    add_overlap_flag_annotsv(sv_type = v, workers = threads)
} else {
  tsv_annotsv <- annotsv
}

logger::log_info("Reading VEP annotations")
maf_vep <- file_maf |>
  my_vroom(
    na_append = ".", delim = "\t", comment = "#", num_threads = threads
  ) |>
  tidyr::drop_na(Start_Position)

logger::log_info("Reading SnpEff annotations")
tsv_snpeff <- file_tsv |>
  my_vroom(delim = "\t", num_threads = threads)

logger::log_info(glue("Applying filters on {s} | {t} | {v} ..."))
ids <- filter_svid(
  caller = t, sv_type = v, annotsv_result = tsv_annotsv, vep_result = maf_vep,
  snpeff_result = tsv_snpeff, terms_relative = terms_relative, stringent = FALSE
)

list_annotation <- list(
  VEP = maf_vep, SnpEff = tsv_snpeff, AnnotSV = tsv_annotsv
)

annotation_modified <- lapply(seq_along(list_annotation), function(i) {
  x <- list_annotation[[i]]
  if ("vcf_id" %in% colnames(x)) {
    x <- x |>
      dplyr::mutate(ID = vcf_id)
  }
  if (length(ids) != 0 && nrow(x) != 0) {
    tmp <- x[x$ID %in% ids, ]
    colnames(tmp) <- glue::glue("{names(list_annotation)[i]}.{colnames(tmp)}")
    tmp
  } else {
    NULL
  }
}) |>
  dplyr::bind_rows()

logger::log_info("Writing the filtered annotations")
annotation_modified |>
  vroom::vroom_write(table, na = "", num_threads = threads)

logger::log_info("Filtering VCF files")
ids |>
  tibble::as_tibble() |>
  vroom::vroom_write(path_ids, col_names = FALSE, num_threads = threads)
system(glue::glue('bcftools filter -i "ID=@{path_ids}" {vep} > {vep_final}'))
system(glue::glue('bcftools filter -i "ID=@{path_ids}" {snpeff} > {snpeff_final}'))
system(glue::glue('bcftools filter -i "ID=@{path_ids}" {vcf} > {vcf_final}'))

logger::log_info("Saving the whole R session")
save.image(snakemake@output[["rdata"]])
