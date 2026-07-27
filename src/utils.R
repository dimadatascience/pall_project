
get_sample_id <- function(filepath) {
  tools::file_path_sans_ext(basename(filepath)) %>%
    sub("\\.cnv\\.annotated$", "", .) %>%
    sub("\\.target\\.counts$", "", .)
}

process_cnv_sample <- function(cnv_file) {
  sample_id <- strsplit(basename(cnv_file), "\\.")[[1]][[1]]

  # 1. Leggi CNV
  cnv <- tryCatch({
    read.delim(cnv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    message(sprintf("Errore nella lettura di %s: %s", cnv_file, e$message))
    return(NULL)
  })

  if (is.null(cnv) || !"All.protein.coding.genes" %in% colnames(cnv)) return(NULL)
  cnv <- cnv[cnv$All.protein.coding.genes != "", ]
  cnv$CNV_Status <- ifelse(cnv$Type == "DEL", "Del",
                           ifelse(cnv$Type == "DUP", "Amp", NA))
  cnv <- cnv[!is.na(cnv$CNV_Status), ]
  if (nrow(cnv) == 0) return(NULL)

  cnv_gr <- GRanges(
    seqnames = cnv$Chromosome,
    ranges = IRanges(start = cnv$Start, end = cnv$End)
  )

  # 🔹 File gene-level (solo se ti serve)
  gene_annotation_df <- cnv %>%
    mutate(Sample = sample_id,
           Start_Position = Start,
           End_Position = End) %>%
    separate_rows(All.protein.coding.genes, sep = ",\\s*") %>%
    mutate(Hugo_Symbol = All.protein.coding.genes) %>%
    dplyr::select(Hugo_Symbol, Sample, CNV_Status, Chromosome, Start_Position, End_Position, Classification)

  # Ritorna entrambi
  return(list(
    cnv_list = gene_annotation_df
  ))
}

apply_priority <- function(df, colnames_list, priority, new_colname = "Clinical_Significance_Summary") {

  colnames_vec <- unlist(colnames_list)
  df <- as.data.frame(df)

  # Funzione che prende un vettore di tag e sceglie quello con priorità più alta
  select_highest_priority <- function(tags) {
    tags <- trimws(tags)
    tags <- tags[tags != "" & !is.na(tags)]

    for (p in priority) {
      if (p %in% tags) return(p)
    }
    return(NA_character_)
  }

  df[[new_colname]] <- apply(df[, colnames_vec, drop=FALSE], 1, select_highest_priority)

  return(data.table(df))
}

simplify_clinvar <- function(value, priority, mapping) {
  if (is.na(value) || trimws(value) == "") {
    return("Not_Provided")
  }

  # Split per tutti i separatori
  tags <- unlist(strsplit(value, "[|/;,]"))
  tags <- trimws(tags)

  # Applica il mapping ai tag
  mapped_tags <- sapply(tags, function(x) {
    if (x %in% names(mapping)) {
      return(mapping[[x]])
    } else {
      return(x)
    }
  }, USE.NAMES = FALSE)

  # Trova il primo tag con priorità più alta nel mapping
  for (p in priority) {
    if (p %in% mapped_tags) {
      return(p)
    }
  }

  # Se nessun match nella priorità, restituisce il primo tag mappato
  return(mapped_tags[1])
}

preproc_oncocn <- function(dt,
                           cyto_file,
                           gain = 0.30,
                           loss = -0.30,
                           top_regions = 50) {

  # ---- 1) Preprocessing ----
  dt_test <- dt %>%
    mutate(
      Chromosome = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome)),
      original_call = ifelse(is.na(Type) | Type == "", NA, Type),
      log2R_seg = ifelse(Segment_Mean > 0, log2(Segment_Mean), NA)   # evita -Inf
    )

  # ---- 2) GRanges ----
  cnv_gr <- GRanges(
    seqnames = dt_test$Chromosome,
    ranges = IRanges(start = dt_test$Start, end = dt_test$End),
    sample = dt_test$Sample,
    log2R = dt_test$log2R_seg,
    original_call = dt_test$original_call
  )

  # ---- 3) Cytoband ----
  cyto <- read.delim(cyto_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(cyto) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  cyto_gr <- GRanges(
    seqnames = cyto$chrom,
    ranges = IRanges(start = cyto$chromStart + 1, end = cyto$chromEnd),
    name = cyto$name,
    gieStain = cyto$gieStain
  )

  cyto_gr$region <- paste0(seqnames(cyto_gr), cyto_gr$name)

  # ---- 4) Match seqlevels ----
  seqlevelsStyle(cnv_gr) <- "UCSC"
  common <- intersect(seqlevels(cyto_gr), seqlevels(cnv_gr))
  cyto_gr <- keepSeqlevels(cyto_gr, common, pruning.mode = "coarse")
  cnv_gr  <- keepSeqlevels(cnv_gr, common, pruning.mode = "coarse")

  # ---- 5) Overlap ----
  ov <- findOverlaps(cnv_gr, cyto_gr, ignore.strand = TRUE)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)

  ints <- pintersect(cnv_gr[qh], cyto_gr[sh])

  df_raw <- data.frame(
    sample = cnv_gr$sample[qh],
    region = cyto_gr$region[sh],
    log2R = cnv_gr$log2R[qh],
    original_call = cnv_gr$original_call[qh],
    w = width(ints),
    stringsAsFactors = FALSE
  )

  # ---- 6) Weighted mean per regione ----
  df <- df_raw %>%
    group_by(sample, region) %>%
    summarise(
      log2R = weighted.mean(log2R, w, na.rm = TRUE),
      original_call = original_call[which.max(w)],
      .groups = "drop"
    )

  # ---- 7) Tua classificazione: SOLO DEL/DUP ----
  df <- df %>%
    mutate(
      my_call = case_when(
        log2R <= loss ~ "DEL",
        log2R >= gain ~ "DUP",
        TRUE ~ NA_character_
      )
    )

  # ---- 8) Final call = DRAGEN se disponibile ----
  df <- df %>%
    mutate(
      final_call = ifelse(
        is.na(original_call) | original_call == "",
        my_call,
        original_call
      )
    ) %>%
    filter(!is.na(final_call))

  # ---- 9) Matrice regione × sample ----
  mat <- df %>%
    dplyr::select(sample, region, final_call) %>%
    pivot_wider(names_from = sample, values_from = final_call, values_fill = "") %>%
    group_by(region) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ";")),
              .groups = "drop") %>%
    tibble::column_to_rownames("region") %>%
    as.matrix()

  # ---- 10) Selezione top regioni ----
  region_freq <- apply(mat, 1, function(x) sum(x != ""))
  top_regions <- names(sort(region_freq, decreasing = TRUE))[1:top_regions]
  mat_top <- mat[top_regions, , drop = FALSE]

  # ---- 11) Oncoprint setup ----
  col_map <- c(DUP = "#e41a1c", DEL = "#377eb8")

  alter_fun <- list(
    background = function(x, y, w, h) {
      grid::grid.rect(x, y, w, h,
                      gp = grid::gpar(fill = "#f2f2f2", col = NA))
    },
    DUP = function(x, y, w, h) {
      grid::grid.rect(x, y, w * 0.9, h * 0.9,
                      gp = grid::gpar(fill = col_map["DUP"], col = NA))
    },
    DEL = function(x, y, w, h) {
      grid::grid.rect(x, y, w * 0.9, h * 0.9,
                      gp = grid::gpar(fill = col_map["DEL"], col = NA))
    }
  )

  get_type <- function(x) strsplit(x, ";")[[1]]

  return(list(
    mat_top = mat_top,
    col_map = col_map,
    df = df
  ))
}

save_pdf_plot <- function(plot_expr, pdf_dir) {

  # Recupera opzioni del chunk corrente
  opts <- knitr::opts_current$get()

  # Nome del chunk
  chunk_name <- opts$label
  if (is.null(chunk_name) || chunk_name == "") {
    chunk_name <- "unnamed_chunk"
  }

  # Dimensioni del chunk
  fig_width  <- opts$fig.width  %||% 10
  fig_height <- opts$fig.height %||% 8

  # Assicurati che la cartella esista
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

  # Gestione del nome con indice progressivo
  index <- 1
  pdf_file <- file.path(pdf_dir, paste0(chunk_name, "-", index, ".pdf"))
  while (file.exists(pdf_file)) {
    index <- index + 1
    pdf_file <- file.path(pdf_dir, paste0(chunk_name, "-", index, ".pdf"))
  }

  # Salva il PDF
  pdf(pdf_file, width = fig_width, height = fig_height)
  eval(plot_expr)   # esegue il plot
  dev.off()

  message("PDF salvato in: ", pdf_file)
}

plot_cn_oncoprint <- function(cn_block,
                              title,
                              remove_empty_columns = TRUE) {

  mat_top <- cn_block$mat_top
  col_map <- cn_block$col_map

  alter_fun <- list(
    background = function(x, y, w, h) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(fill = "#f2f2f2", col = NA)
      )
    },
    DUP = function(x, y, w, h) {
      grid::grid.rect(
        x, y, w * 0.9, h * 0.9,
        gp = grid::gpar(fill = col_map["DUP"], col = NA)
      )
    },
    DEL = function(x, y, w, h) {
      grid::grid.rect(
        x, y, w * 0.9, h * 0.9,
        gp = grid::gpar(fill = col_map["DEL"], col = NA)
      )
    }
  )

  missing_samples <- setdiff(all_samples, colnames(mat_top))

  if (length(missing_samples) > 0) {

    mat_add <- matrix(
      "",
      nrow = nrow(mat_top),
      ncol = length(missing_samples),
      dimnames = list(
        rownames(mat_top),
        missing_samples
      )
    )

    mat_top <- cbind(mat_top, mat_add)

  }

  mat_top <- mat_top[, all_samples]

  anno_df <- clinical_sample_info %>%
    filter(Tumor_Sample_Barcode %in% colnames(mat_top)) %>%
    arrange(T1221, Sample_group)

  mat_top <- mat_top[, anno_df$Tumor_Sample_Barcode]

  bottom_anno <- HeatmapAnnotation(
    T1221 = anno_df$T1221,
    Sample_Group = anno_df$Sample_group,
    col = list(
      T1221 = t1221_colors,
      Sample_Group = sample_group_colors
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      T1221 = list(title = "T1221"),
      Sample_Group = list(title = "Sample Group")
    )
  )

  ht <- ComplexHeatmap::oncoPrint(
    mat_top,
    get_type = get_type,
    alter_fun = alter_fun,
    col = col_map,
    top_annotation = HeatmapAnnotation(
      `Alteration burden` = anno_oncoprint_barplot(),
      annotation_name_side = "left"
    ),
    left_annotation = rowAnnotation(
      `Alteration freq` = anno_oncoprint_barplot(),
      annotation_name_rot = 0
    ),
    bottom_annotation = bottom_anno,
    column_order = anno_df$Tumor_Sample_Barcode,
    remove_empty_rows = TRUE,
    remove_empty_columns = remove_empty_columns,
    show_column_names = TRUE,
    column_labels = anno_df$Tumor_Sample_Barcode,
    column_names_gp = grid::gpar(fontsize = 8, rot = 45),
    column_title = title
  )

  draw(ht)

}
