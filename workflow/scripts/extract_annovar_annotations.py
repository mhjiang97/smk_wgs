"""Extract ANNOVAR annotations."""

# /// script
# dependencies = [
#   "polars",
# ]
# ///

# %%
import logging

import polars as pl


# %%
def configure_logger(console_output=None, log_file=None, log_level="INFO"):
    logger = logging.getLogger()

    log_levels = {
        "CRITICAL": logging.CRITICAL,
        "ERROR": logging.ERROR,
        "WARNING": logging.WARNING,
        "INFO": logging.INFO,
        "DEBUG": logging.DEBUG,
        "NOTSET": logging.NOTSET,
    }

    log_level_constant = log_levels.get(log_level.upper(), logging.INFO)
    logger.setLevel(log_level_constant)

    formatter = logging.Formatter(
        "%(levelname)s [%(asctime)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    if logger.hasHandlers():
        logger.handlers.clear()

    if console_output:
        ch = logging.StreamHandler()
        ch.setLevel(log_level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(log_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def filter_annotations(
    IDs: list[str], annovar_input: pl.LazyFrame, annovar_output: pl.LazyFrame
):
    df_fil = annovar_input.filter(pl.col("id").is_in(IDs))

    cols = list(OCOL2NCOL.values())
    cols.remove("id")
    res = annovar_output.join(df_fil, on=cols, how="right").select(
        [
            "Chr",
            "Start",
            "End",
            "id",
            pl.all().exclude(["Chr", "Start", "End", "id"]),
        ]
    )

    return res


# %%
path_anno = snakemake.input[0]  # type: ignore # noqa: F821
path_tmp = snakemake.input[1]  # type: ignore # noqa: F821
path_snv = snakemake.input[2]  # type: ignore # noqa: F821
path_indel = snakemake.input[3]  # type: ignore # noqa: F821
path_out_snv = snakemake.output[0]  # type: ignore # noqa: F821
path_out_indel = snakemake.output[1]  # type: ignore # noqa: F821
path_log = snakemake.log[0]  # type: ignore # noqa: F821
sample = snakemake.wildcards["sample"]  # type: ignore # noqa: F821

OCOL2NCOL = {
    "column_1": "Chr",
    "column_2": "Start",
    "column_3": "End",
    "column_4": "Ref",
    "column_5": "Alt",
    "column_8": "id",
}
COLS = [
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "qual",
    "filter",
    "info",
    "format",
    sample,
]


# %%
def main():
    configure_logger(console_output=False, log_file=path_log, log_level="INFO")

    df_anno = pl.scan_csv(
        path_anno,
        separator="\t",
        null_values=["NA"],
        infer_schema=False,
    )
    df_tmp = (
        pl.scan_csv(path_tmp, separator="\t", has_header=False, infer_schema=False)
        .rename(OCOL2NCOL)
        .select(OCOL2NCOL.values())
    )
    df_snv = pl.scan_csv(
        path_snv, separator="\t", comment_prefix="#", new_columns=COLS
    )
    df_indel = pl.scan_csv(
        path_indel, separator="\t", comment_prefix="#", new_columns=COLS
    )

    ids_snv = df_snv.select("id").collect().to_series().to_list()
    logging.info(f"Filtering ANNOVAR annotations for {len(ids_snv)} SNVs...")
    df_anno_snv = filter_annotations(ids_snv, df_tmp, df_anno).collect()
    logging.info(f"Writing SNV annotations to {path_out_snv}...")
    df_anno_snv.write_csv(
        path_out_snv, separator="\t", null_value="NA", include_header=True
    )

    ids_indel = df_indel.select("id").collect().to_series().to_list()
    logging.info(f"Filtering ANNOVAR annotations for {len(ids_indel)} indels...")
    df_anno_indel = filter_annotations(ids_indel, df_tmp, df_anno).collect()
    logging.info(f"Writing indel annotations to {path_out_indel}...")
    df_anno_indel.write_csv(
        path_out_indel, separator="\t", null_value="NA", include_header=True
    )


# %%
if __name__ == "__main__":
    main()

# %%
