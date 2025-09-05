load_pkg <- function(pkgs) {
  for (pkg in pkgs) {
    if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))) {
      install.packages(pkg)
      if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
        suppressMessages(require(pkg, character.only = TRUE))
      }
    }
  }
}

get_sv_type <- function(grs, ins_threshold = 0.7) {
  pgr <- partner(grs)

  types_sv <- case_when(
    as.logical(seqnames(grs) != seqnames(pgr)) ~ "BND",

    as.logical(strand(grs) == strand(pgr)) ~ "INV",

    grs$insLen >= abs(grs$svLen) * ins_threshold ~ "INS",

    as.logical(xor(start(grs) < start(pgr), strand(grs) == "-")) ~ "DEL",

    TRUE ~ "DUP"
  )

  types_sv
}

my_vroom <- function(
  file, na_append = NULL, guess_lines = 10, n_col = NULL, ...
) {
  na_origin <- c("", "NA", "NaN")
  myna <- c(na_origin, na_append)

  if (is.null(n_col)) {
    n_col <- vroom::vroom(
      file, n_max = guess_lines, col_types = vroom::cols(), na = myna,
      guess_max = Inf, ...
    ) |>
      ncol()
  }

  dat <- vroom::vroom(
    file, col_types = paste0(rep_len("c", n_col), collapse = ""), na = myna, ...
  )

  dat
}

nothing <- function(x) {
  x
}

create_grs <- function(string) {
  if (is.na(string)[1]) return(as("1:1-1", "GRanges"))
  as(string, "GRanges")
}

calculate_overlap_flag <- function(coord, sv) {
  gr <- str_replace_all(coord, "chr", "") |>
    str_split(";") |>
    unlist() |>
    create_grs()

  pcts_gr <- width(pintersect(gr, sv)) / width(gr)
  pcts_sv <- width(pintersect(gr, sv)) / width(sv)
  any(pcts_gr >= 0.9 & pcts_sv >= 0.8)
}

add_overlap_flag_annotsv <- function(tsv, sv_type, workers = detectCores()) {
  columns <- c(
    "B_gain_coord", "B_loss_coord", "po_B_loss_allG_coord",
    "B_ins_coord", "po_B_gain_allG_coord", "B_inv_coord"
  )
  flag_columns <- c(
    "flag_gain", "flag_loss", "flag_po_loss",
    "flag_ins", "flag_po_gain", "flag_inv"
  )

  for (i in seq_along(flag_columns)) {
    tsv[[flag_columns[i]]] <- rep(FALSE, nrow(tsv))
  }

  if (sv_type == "INS") {
    gr_sv <- GRanges(
      seqnames = as.character(tsv$SV_chrom),
      ranges = IRanges(
        as.numeric(tsv$SV_start),
        as.numeric(tsv$SV_start) + as.numeric(tsv$SV_length) - 1
      )
    )
  } else {
    gr_sv <- GRanges(
      seqnames = as.character(tsv$SV_chrom),
      ranges = IRanges(as.numeric(tsv$SV_start), as.numeric(tsv$SV_end))
    )
  }

  pb <- progress::progress_bar$new(
    format = "Processing overlap flags [:bar] :percent in :elapsed ETA: :eta",
    total = length(columns),
    clear = FALSE
  )

  register(MulticoreParam(workers = workers))

  for (j in seq_along(columns)) {
    pb$tick(tokens = list(column = flag_columns[j]))

    tsv[[flag_columns[j]]] <- bplapply(
      seq_len(nrow(tsv)), function(i) {
        calculate_overlap_flag(tsv[[columns[j]]][i], gr_sv[i])
      }
    ) |>
      unlist()
  }

  pb$terminate()

  tsv
}

filter_vep <- function(v, maf_vep, t) {
  id_include <- switch(v,
    BND = {
      maf_vep |> nothing()
    },
    DEL = {
      maf_vep |> nothing()
    },
    INS = {
      maf_vep |> nothing()
    },
    INV = {
      maf_vep |> nothing()
    },
    DUP = {
      maf_vep |> nothing()
    }
  ) |>
    dplyr::pull(vcf_id) |>
    na.omit() |>
    unique()

  min_qual <- switch(t, cutesv = 0, severus = 0, sniffles = 0, svim = 0, svision = 0, 0)

  id_exclude <- switch(v,
    BND = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomADe_AF) >= 0.01 | as.numeric(gnomADe_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    DEL = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomADe_AF) >= 0.01 | as.numeric(gnomADe_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    INS = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
            as.numeric(gnomADe_AF) >= 0.01 | as.numeric(gnomADe_EAS_AF) >= 0.01 |
            as.numeric(vcf_qual) < min_qual
        )
    },
    INV = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomADe_AF) >= 0.01 | as.numeric(gnomADe_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    DUP = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomADe_AF) >= 0.01 | as.numeric(gnomADe_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    }
  ) |>
    dplyr::pull(vcf_id) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

filter_snpeff <- function(v, tsv_snpeff, t) {
  id_include <- switch(v,
    BND = {
      tsv_snpeff |> nothing()
    },
    DEL = {
      tsv_snpeff |> nothing()
    },
    INS = {
      tsv_snpeff |> nothing()
    },
    INV = {
      tsv_snpeff |> nothing()
    },
    DUP = {
      tsv_snpeff |> nothing()
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  min_qual <- switch(t, cutesv = 0, severus = 0, sniffles = 0, svim = 0, svision = 0, 0)

  id_exclude <- switch(v,
    BND = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    DEL = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    INS = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    INV = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    DUP = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

filter_annotsv <- function(v, tsv_annotsv, terms_relative) {
  id_include <- switch(v,
    BND = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  id_exclude <- switch(v,
    BND = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_gain_AFmax) >= 0.05 & as.numeric(B_loss_AFmax) >= 0.05 & flag_gain & flag_loss
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_loss_AFmax) >= 0.05 & flag_loss) | (!is.na(po_B_loss_allG_coord) & flag_po_loss)
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_ins_AFmax) >= 0.05 & flag_ins) | (!is.na(po_B_gain_allG_coord) & flag_po_gain)
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_inv_AFmax) >= 0.05 & flag_inv
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_gain_AFmax) >= 0.05 & flag_gain) | (!is.na(po_B_gain_allG_coord) & flag_po_gain)
        )
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

merge_svid <- function(list_vep, list_snpeff, list_annotsv, stringent = TRUE) {
  if (stringent) {
    ids_include <- unique(
      c(list_vep$id_include, list_snpeff$id_include, list_annotsv$id_include)
    )
    ids_exclude <- unique(
      c(list_vep$id_exclude, list_snpeff$id_exclude, list_annotsv$id_exclude)
    )
    ids <- setdiff(ids_include, ids_exclude) |>
      union(list_annotsv$id_include) |>
      unique()
  } else {
    ids_vep <- setdiff(list_vep$id_include, list_vep$id_exclude)
    ids_snpeff <- setdiff(list_snpeff$id_include, list_snpeff$id_exclude)
    ids <- setdiff(c(ids_vep, ids_snpeff), list_annotsv$id_exclude) |>
      union(list_annotsv$id_include) |>
      unique()
  }

  ids
}

filter_svid <- function(caller, sv_type, id_tab, annotsv_result, vep_result, snpeff_result, terms_relative, ...) {
  if (nrow(vep_result) == 0) {
    list_vep <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_vep <- filter_vep(v = sv_type, maf_vep = vep_result, t = caller)
  }
  if (nrow(snpeff_result) == 0) {
    list_snpeff <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_snpeff <- filter_snpeff(v = sv_type, tsv_snpeff = snpeff_result, t = caller)
  }
  if (nrow(annotsv_result) == 0) {
    list_annotsv <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_annotsv <- filter_annotsv(v = sv_type, tsv_annotsv = annotsv_result, terms_relative = terms_relative)
  }

  ids <- merge_svid(list_vep = list_vep, list_snpeff = list_snpeff, list_annotsv = list_annotsv, ...)

  ids
}
