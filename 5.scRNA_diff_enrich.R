CONFIG <- list(
  seurat_rds = "result.seurat.rds",
  comparison_name = "M_3_vs_V_3",
  group_col = "cell",
  group1 = "M_3",
  group2 = "V_3",
  assay = "RNA",
  presto_slot = "data",
  species = "Pig",
  de_software = "presto",
  test_use = "wilcox",
  min_pct = 0.10,
  logfc_threshold = 0.25,
  threshold_type = "p_val_adj",
  threshold_value = 0.05,
  downsample = FALSE,
  max_cells_per_ident = 1000,
  downsample_seed = 123,
  enrich_method = "ORA",
  databases = c("GO", "KEGG", "Reactome", "hallmarker-H"),
  go_ontology = "BP",
  enrich_pvalue_cutoff = 0.05,
  enrich_qvalue_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  outdir = "diff_enrich_out"
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(tolower(x) %in% c("1", "true", "t", "yes", "y"))
  FALSE
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

get_meta_vector <- function(obj, colname) {
  if (!colname %in% colnames(obj[[]])) {
    stop("Group column not found in metadata: ", colname)
  }
  obj[[colname]][, 1, drop = TRUE]
}

get_species_info <- function(species) {
  species <- tolower(species)
  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("org.Hs.eg.db required for human species")
    }
    return(list(
      orgdb = getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
      kegg_code = "hsa",
      reactome_code = "human",
      msigdbr_species = "Homo sapiens"
    ))
  }
  if (species == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("org.Mm.eg.db required for mouse species")
    }
    return(list(
      orgdb = getExportedValue("org.Mm.eg.db", "org.Mm.eg.db"),
      kegg_code = "mmu",
      reactome_code = "mouse",
      msigdbr_species = "Mus musculus"
    ))
  }
  stop("Only human and mouse species are supported")
}

subset_to_groups <- function(obj, group_col, group1, group2) {
  meta_group <- get_meta_vector(obj, group_col)
  keep_cells <- colnames(obj)[meta_group %in% c(group1, group2)]
  if (length(keep_cells) == 0) {
    stop("No cells found for group1/group2 in column ", group_col)
  }
  subset(obj, cells = keep_cells)
}

downsample_groups <- function(obj, group_col, max_cells_per_ident, seed = 123) {
  if (is.null(max_cells_per_ident) || !is.finite(max_cells_per_ident) || max_cells_per_ident <= 0) {
    return(obj)
  }
  set.seed(seed)
  md <- obj[[]]
  keep_cells <- unlist(lapply(split(rownames(md), md[[group_col]]), function(cells) {
    if (length(cells) <= max_cells_per_ident) return(cells)
    sample(cells, size = max_cells_per_ident, replace = FALSE)
  }), use.names = FALSE)
  subset(obj, cells = keep_cells)
}

standardize_de_table <- function(df, software, group1, group2) {
  df$gene <- rownames(df)
  if ("avg_logFC" %in% colnames(df) && !"avg_log2FC" %in% colnames(df)) {
    df$avg_log2FC <- df$avg_logFC
  }
  if (!"p_val" %in% colnames(df)) {
    stop("p_val column not found in DE results")
  }
  if (!"p_val_adj" %in% colnames(df)) {
    if ("padj" %in% colnames(df)) {
      df$p_val_adj <- df$padj
    } else {
      df$p_val_adj <- stats::p.adjust(df$p_val, method = "BH")
    }
  }
  if (!"pct.1" %in% colnames(df) && "pct_in" %in% colnames(df)) {
    df$pct.1 <- df$pct_in
  }
  if (!"pct.2" %in% colnames(df) && "pct_out" %in% colnames(df)) {
    df$pct.2 <- df$pct_out
  }
  if (!"pct.1" %in% colnames(df)) df$pct.1 <- NA_real_
  if (!"pct.2" %in% colnames(df)) df$pct.2 <- NA_real_
  if (!"avg_log2FC" %in% colnames(df)) {
    stop("avg_log2FC/avg_logFC/logFC not found in DE results")
  }
  df$group1 <- group1
  df$group2 <- group2
  df$software <- software
  df
}

run_de_by_presto <- function(obj, config) {
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("presto package required")
  }
  obj <- subset_to_groups(obj, config$group_col, config$group1, config$group2)
  if (as_bool(config$downsample)) {
    obj <- downsample_groups(obj, config$group_col, config$max_cells_per_ident, config$downsample_seed)
  }
  Seurat::DefaultAssay(obj) <- config$assay
  group_vec <- get_meta_vector(obj, config$group_col)
  expr_mat <- Seurat::GetAssayData(obj, assay = config$assay, slot = config$presto_slot)
  res <- presto::wilcoxauc(X = expr_mat, y = group_vec)
  res <- res[res$group == config$group1, , drop = FALSE]
  if ("feature" %in% colnames(res)) {
    rownames(res) <- res$feature
  }
  if ("logFC" %in% colnames(res) && !"avg_log2FC" %in% colnames(res)) {
    res$avg_log2FC <- res$logFC
  }
  standardize_de_table(res, software = "presto", group1 = config$group1, group2 = config$group2)
}

run_de_by_findmarker <- function(obj, config) {
  obj <- subset_to_groups(obj, config$group_col, config$group1, config$group2)
  if (as_bool(config$downsample)) {
    obj <- downsample_groups(obj, config$group_col, config$max_cells_per_ident, config$downsample_seed)
  }
  obj <- Seurat::SetIdent(obj, value = config$group_col)
  res <- Seurat::FindMarkers(
    object = obj,
    ident.1 = config$group1,
    ident.2 = config$group2,
    assay = config$assay,
    test.use = config$test_use,
    min.pct = config$min_pct,
    logfc.threshold = config$logfc_threshold,
    max.cells.per.ident = if (as_bool(config$downsample)) config$max_cells_per_ident else Inf,
    verbose = FALSE
  )
  standardize_de_table(res, software = "findmarker", group1 = config$group1, group2 = config$group2)
}

run_diff_analysis <- function(obj, config) {
  software <- tolower(config$de_software)
  if (software == "presto") {
    res <- run_de_by_presto(obj, config)
    keep <- rep(TRUE, nrow(res))
    keep <- keep & (pmax(res$pct.1, res$pct.2, na.rm = TRUE) >= config$min_pct | is.na(res$pct.1) | is.na(res$pct.2))
    keep <- keep & abs(res$avg_log2FC) >= config$logfc_threshold
    res <- res[keep, , drop = FALSE]
    return(res[order(res$p_val_adj, res$p_val), , drop = FALSE])
  }
  if (software == "findmarker") {
    res <- run_de_by_findmarker(obj, config)
    return(res[order(res$p_val_adj, res$p_val), , drop = FALSE])
  }
  stop("de_software must be presto or findmarker")
}

filter_significant_genes <- function(de_res, config) {
  threshold_col <- config$threshold_type
  if (!threshold_col %in% colnames(de_res)) {
    stop("Threshold column not found in DE results: ", threshold_col)
  }
  keep <- de_res[[threshold_col]] <= config$threshold_value
  keep <- keep & abs(de_res$avg_log2FC) >= config$logfc_threshold
  de_res[keep, , drop = FALSE]
}

convert_symbol_to_entrez <- function(genes, orgdb) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler required for enrichment analysis")
  }
  if (length(genes) == 0) return(data.frame())
  suppressWarnings(clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  ))
}

make_ranked_entrez <- function(de_res, orgdb) {
  conv <- convert_symbol_to_entrez(de_res$gene, orgdb)
  if (nrow(conv) == 0) return(numeric())
  merged <- merge(
    conv,
    de_res[, c("gene", "avg_log2FC")],
    by.x = "SYMBOL",
    by.y = "gene",
    all.x = TRUE,
    all.y = FALSE
  )
  merged <- merged[!is.na(merged$avg_log2FC), , drop = FALSE]
  merged <- merged[order(abs(merged$avg_log2FC), decreasing = TRUE), , drop = FALSE]
  merged <- merged[!duplicated(merged$ENTREZID), , drop = FALSE]
  ranked <- merged$avg_log2FC
  names(ranked) <- merged$ENTREZID
  sort(ranked, decreasing = TRUE)
}

make_hallmark_term2gene <- function(species_info) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr required for hallmark-H")
  }
  hallmark <- msigdbr::msigdbr(species = species_info$msigdbr_species, category = "H")
  unique(hallmark[, c("gs_name", "gene_symbol")])
}

run_ora_single_db <- function(db_name, sig_genes, species_info, config) {
  if (length(sig_genes) == 0) return(NULL)
  if (db_name == "GO") {
    entrez <- convert_symbol_to_entrez(sig_genes, species_info$orgdb)$ENTREZID
    if (length(entrez) == 0) return(NULL)
    return(clusterProfiler::enrichGO(
      gene = unique(entrez),
      OrgDb = species_info$orgdb,
      keyType = "ENTREZID",
      ont = config$go_ontology,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      qvalueCutoff = config$enrich_qvalue_cutoff,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      readable = TRUE
    ))
  }
  if (db_name == "KEGG") {
    entrez <- convert_symbol_to_entrez(sig_genes, species_info$orgdb)$ENTREZID
    if (length(entrez) == 0) return(NULL)
    return(clusterProfiler::enrichKEGG(
      gene = unique(entrez),
      organism = species_info$kegg_code,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      qvalueCutoff = config$enrich_qvalue_cutoff,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size
    ))
  }
  if (db_name == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("ReactomePA required for Reactome enrichment")
    }
    entrez <- convert_symbol_to_entrez(sig_genes, species_info$orgdb)$ENTREZID
    if (length(entrez) == 0) return(NULL)
    return(ReactomePA::enrichPathway(
      gene = unique(entrez),
      organism = species_info$reactome_code,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      qvalueCutoff = config$enrich_qvalue_cutoff,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      readable = TRUE
    ))
  }
  if (db_name == "hallmarker-H") {
    term2gene <- make_hallmark_term2gene(species_info)
    return(clusterProfiler::enricher(
      gene = unique(sig_genes),
      TERM2GENE = term2gene,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      qvalueCutoff = config$enrich_qvalue_cutoff,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size
    ))
  }
  stop("Unsupported database: ", db_name)
}

run_gsea_single_db <- function(db_name, de_res, species_info, config) {
  ranked_entrez <- make_ranked_entrez(de_res, species_info$orgdb)
  if (db_name == "GO") {
    if (length(ranked_entrez) == 0) return(NULL)
    return(clusterProfiler::gseGO(
      geneList = ranked_entrez,
      OrgDb = species_info$orgdb,
      keyType = "ENTREZID",
      ont = config$go_ontology,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      verbose = FALSE
    ))
  }
  if (db_name == "KEGG") {
    if (length(ranked_entrez) == 0) return(NULL)
    return(clusterProfiler::gseKEGG(
      geneList = ranked_entrez,
      organism = species_info$kegg_code,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      verbose = FALSE
    ))
  }
  if (db_name == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("ReactomePA required for Reactome enrichment")
    }
    if (length(ranked_entrez) == 0) return(NULL)
    return(ReactomePA::gsePathway(
      geneList = ranked_entrez,
      organism = species_info$reactome_code,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      verbose = FALSE
    ))
  }
  if (db_name == "hallmarker-H") {
    term2gene <- make_hallmark_term2gene(species_info)
    ranked_symbol <- de_res$avg_log2FC
    names(ranked_symbol) <- de_res$gene
    ranked_symbol <- sort(ranked_symbol[!duplicated(names(ranked_symbol))], decreasing = TRUE)
    return(clusterProfiler::GSEA(
      geneList = ranked_symbol,
      TERM2GENE = term2gene,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      pvalueCutoff = config$enrich_pvalue_cutoff,
      verbose = FALSE
    ))
  }
  stop("Unsupported database: ", db_name)
}

write_table_if_exists <- function(x, file) {
  if (is.null(x)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  utils::write.csv(df, file = file, row.names = FALSE)
}

run_enrichment <- function(de_all, de_sig, config) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler required for enrichment analysis")
  }
  species_info <- get_species_info(config$species)
  dbs <- config$databases
  if (length(dbs) == 0) return(list())
  method <- toupper(config$enrich_method)
  results <- list()
  if (method == "ORA") {
    up_genes <- unique(de_sig$gene[de_sig$avg_log2FC > 0])
    down_genes <- unique(de_sig$gene[de_sig$avg_log2FC < 0])
    for (db in dbs) {
      results[[paste0(db, "_up")]] <- run_ora_single_db(db, up_genes, species_info, config)
      results[[paste0(db, "_down")]] <- run_ora_single_db(db, down_genes, species_info, config)
    }
    return(results)
  }
  if (method == "GSEA") {
    for (db in dbs) {
      results[[db]] <- run_gsea_single_db(db, de_all, species_info, config)
    }
    return(results)
  }
  stop("enrich_method must be ORA or GSEA")
}

save_results <- function(de_all, de_sig, enrich_res, config) {
  outdir <- ensure_dir(file.path(config$outdir, config$comparison_name))
  utils::write.csv(de_all, file = file.path(outdir, paste0(config$comparison_name, ".de_all.csv")), row.names = FALSE)
  utils::write.csv(de_sig, file = file.path(outdir, paste0(config$comparison_name, ".de_sig.csv")), row.names = FALSE)
  if (length(enrich_res) > 0) {
    for (nm in names(enrich_res)) {
      write_table_if_exists(
        enrich_res[[nm]],
        file = file.path(outdir, paste0(config$comparison_name, ".", nm, ".csv"))
      )
    }
  }
  saveRDS(
    list(
      config = config,
      de_all = de_all,
      de_sig = de_sig,
      enrich = enrich_res
    ),
    file = file.path(outdir, paste0(config$comparison_name, ".all_results.rds"))
  )
}

run_diff_enrich_pipeline <- function(config = CONFIG) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required")
  }
  obj <- readRDS(config$seurat_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input file is not a Seurat object: ", config$seurat_rds)
  }
  de_all <- run_diff_analysis(obj, config)
  de_sig <- filter_significant_genes(de_all, config)
  enrich_res <- run_enrichment(de_all, de_sig, config)
  save_results(de_all, de_sig, enrich_res, config)
  invisible(list(de_all = de_all, de_sig = de_sig, enrich = enrich_res))
}

if (sys.nframe() == 0) {
  run_diff_enrich_pipeline(CONFIG)
}
