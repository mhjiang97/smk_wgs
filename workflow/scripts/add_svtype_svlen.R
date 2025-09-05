libs_r <- snakemake@params[["libs_r"]]
if (!is.null(libs_r)) {
  .libPaths(c(libs_r, .libPaths()))
}
snakemake@source("utils.R")
load_pkg(
  c(
    "VariantAnnotation", "StructuralVariantAnnotation", "stringr", "dplyr",
    "tibble"
  )
)
logger::log_appender(logger::appender_file(snakemake@log[[1]]))

logger::log_info("Reading Snakemake input, output, and parameters")
vcf_path <- snakemake@input[["vcf"]]
vcf_path_formatted <- snakemake@output[["vcf"]]
rdata <- snakemake@output[["rdata"]]
genome <- snakemake@params[["genome"]]
sample <- snakemake@wildcards[["sample"]]

if (genome == "GRCh37") genome <- "hg19"

caller <- vcf_path |>
  dirname() |>
  dirname()

logger::log_info("Reading the VCF file")
vcf <- readVcf(vcf_path, genome)

df_header <- vcf |>
  header() |>
  info() |>
  as.data.frame() |>
  dplyr::bind_rows(
    data.frame(
      row.names = c("SVLEN"),
      Number = c("1"),
      Type = c("Integer"),
      Description = c("SVLEN")
    )
  ) |>
  as("DataFrame") |>
  unique()
info(header(vcf)) <- df_header

logger::log_info("Extracting breakpoint ranges and SV types")
grs <- breakpointRanges(vcf, inferMissingBreakends = FALSE)

logger::log_info("Determining SV types and lengths")
types_sv <- get_sv_type(grs)
df <- tibble::tibble(
  id = grs$sourceId,
  sv_type = types_sv,
  length = grs$svLen
)

if (caller %in% c("gridss")) {
  df <- df |>
    mutate(
      length = case_when(
        sv_type == "INV" ~ abs(length),
        sv_type == "INS" ~ grs$insLen[match(id, grs$sourceId)],
        TRUE ~ length
      )
    )
}

if (caller %in% c("svaba")) {
  header(vcf)@samples <- sample
  rownames(vcf@colData) <- sample
}

info(vcf[df$id])$SVTYPE <- df$sv_type
info(vcf)$SVLEN <- NA_integer_
info(vcf[df$id])$SVLEN <- df$length

logger::log_info("Writing the formatted VCF file")
writeVcf(vcf, vcf_path_formatted)

logger::log_info("Saving the environment to the RData file")
save.image(rdata)
