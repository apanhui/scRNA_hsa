CONFIG <- list(
  inputs = c("M_14_1_SGS2406069","M_14_2_SGS2406069","M_3_1_SGS2406069","M_3_2_SGS2406069","V_14_1_SGS2406069","V_14_2_SGS2406069","V_3_1_SGS2406069","V_3_2_SGS2406069"),
  sample_names = NULL,
  species = "human",
  min_cells_per_gene = 3,
  umi_range = c(3, 30000),
  feature_range = c(300, 30000),
  mito_related_genes = NULL,
  mito_method = "MAD",
  mito_mad_n = 3,
  mito_percent_max = NULL,
  per_sample_mode = "none",
  per_sample_max_cells = NULL,
  integrate_method = "harmony",
  hvg_method = "vst",
  hvg_nfeatures = 2000,
  integrate_pca_dims = 30,
  cluster_resolution = 0.5,
  cluster_pca_dims = 30,
  outdir = "scRNA_out",
  prefix = "result",
  markers_do = FALSE
)

build_params_from_config <- function(config) {
  params <- default_params()
  params$input$species <- config$species
  params$input$sample_names <- config$sample_names
  params$filter$min_cells_per_gene <- config$min_cells_per_gene
  params$filter$umi_range <- config$umi_range
  params$filter$feature_range <- config$feature_range
  params$filter$mito$related_genes <- config$mito_related_genes
  params$filter$mito$method <- config$mito_method
  params$filter$mito$mad_n <- config$mito_mad_n
  params$filter$mito$percent_max <- config$mito_percent_max
  params$filter$per_sample$mode <- config$per_sample_mode
  params$filter$per_sample$max_cells <- config$per_sample_max_cells
  params$integrate$method <- config$integrate_method
  params$integrate$hvg_method <- config$hvg_method
  params$integrate$hvg_nfeatures <- config$hvg_nfeatures
  params$integrate$pca_dims <- config$integrate_pca_dims
  params$cluster$resolution <- config$cluster_resolution
  params$cluster$pca_dims <- config$cluster_pca_dims
  params$output$outdir <- config$outdir
  params$output$prefix <- config$prefix
  params$markers$do <- config$markers_do
  params
}

default_params <- function() {
  list(
    env = list(
      seurat_major = NULL
    ),
    input = list(
      assay = "RNA",
      project = "scRNA",
      species = "human",
      sample_names = NULL
    ),
    filter = list(
      min_cells_per_gene = 3,
      umi_range = c(200, 100000),
      feature_range = c(100, 30000),
      mito = list(
        related_genes = NULL,
        method = "MAD",
        mad_n = 3,
        percent_max = NULL,
        colname = "percent.mito"
      ),
      per_sample = list(
        mode = "none",
        max_cells = NULL,
        seed = 1
      )
    ),
    integrate = list(
      method = "harmony",
      hvg_method = "vst",
      hvg_nfeatures = 2000,
      pca_dims = 30
    ),
    cluster = list(
      resolution = 0.5,
      pca_dims = 30,
      umap = list(min_dist = 0.3)
    ),
    markers = list(
      do = FALSE,
      only_pos = TRUE,
      min_pct = 0.1,
      logfc_threshold = 0.25
    ),
    output = list(
      outdir = "scRNA_out",
      prefix = "result",
      save_rds = TRUE,
      save_metadata_csv = TRUE,
      save_markers_csv = TRUE
    )
  )
}
merge_params <- function(base, override) {
  if (is.null(override)) return(base)
  for (nm in names(override)) {
    if (is.list(base[[nm]]) && is.list(override[[nm]])) {
      base[[nm]] <- merge_params(base[[nm]], override[[nm]])
    } else {
      base[[nm]] <- override[[nm]]
    }
  }
  base
}
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(tolower(x) %in% c("1", "true", "t", "yes", "y"))
  FALSE
}
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}
infer_mito_genes <- function(obj, species = c("human", "mouse"), custom_genes = NULL) {
  if (!is.null(custom_genes) && length(custom_genes) > 0) return(custom_genes)
  species <- match.arg(species)
  feats <- rownames(obj)
  if (is.null(feats)) return(character())
  pat <- if (species == "human") "^MT-" else "^mt-"
  grep(pat, feats, value = TRUE)
}
apply_qc_filters <- function(obj, params) {
  assay <- params$input$assay
  if (!assay %in% names(obj@assays)) {
    assay <- Seurat::DefaultAssay(obj)
  }
  Seurat::DefaultAssay(obj) <- assay
  nCount_col <- paste0("nCount_", assay)
  nFeature_col <- paste0("nFeature_", assay)
  if (!all(c(nCount_col, nFeature_col) %in% colnames(obj[[]]))) {
    stop("QC metrics not found in metadata: ", nCount_col, " / ", nFeature_col)
  }
  umi_rng <- params$filter$umi_range
  feat_rng <- params$filter$feature_range
  if (length(umi_rng) != 2 || length(feat_rng) != 2) stop("umi_range/feature_range must be numeric vectors of length 2")
  md <- obj[[]]
  keep_cells <- rownames(md)[
    md[[nCount_col]] >= umi_rng[[1]] &
      md[[nCount_col]] <= umi_rng[[2]] &
      md[[nFeature_col]] >= feat_rng[[1]] &
      md[[nFeature_col]] <= feat_rng[[2]]
  ]
  obj <- subset(obj, cells = keep_cells)
  mito_col <- params$filter$mito$colname
  mito_genes <- infer_mito_genes(obj, species = params$input$species, custom_genes = params$filter$mito$related_genes)
  if (length(mito_genes) > 0) {
    obj[[mito_col]] <- Seurat::PercentageFeatureSet(obj, features = mito_genes)
    method <- toupper(params$filter$mito$method %||% "MAD")
    if (method == "MAD") {
      v <- obj[[mito_col]][, 1, drop = TRUE]
      thr <- stats::median(v) + as.numeric(params$filter$mito$mad_n) * stats::mad(v)
      keep_cells <- rownames(obj[[]])[obj[[mito_col]][, 1, drop = TRUE] <= thr]
      obj <- subset(obj, cells = keep_cells)
    } else if (method %in% c("THRESHOLD", "CUSTOM")) {
      mx <- params$filter$mito$percent_max
      if (is.null(mx) || !is.finite(mx)) stop("mito.percent_max required when method=threshold/custom")
      keep_cells <- rownames(obj[[]])[obj[[mito_col]][, 1, drop = TRUE] <= mx]
      obj <- subset(obj, cells = keep_cells)
    } else {
      stop("Unknown mito.method: ", params$filter$mito$method)
    }
  }
  obj
}
cap_cells_per_sample <- function(obj, group_col, max_cells, seed = 1) {
  if (is.null(max_cells) || !is.finite(max_cells) || max_cells <= 0) return(obj)
  if (!group_col %in% colnames(obj[[]])) stop("cap_cells_per_sample requires group_col to exist in metadata: ", group_col)
  set.seed(seed)
  md <- obj[[]]
  keep <- unlist(lapply(split(rownames(md), md[[group_col]]), function(cells) {
    if (length(cells) <= max_cells) return(cells)
    sample(cells, size = max_cells, replace = FALSE)
  }), use.names = FALSE)
  subset(obj, cells = keep)
}
normalize_and_hvg <- function(obj, params) {
  Seurat::NormalizeData(obj, verbose = FALSE) |>
    Seurat::FindVariableFeatures(selection.method = params$integrate$hvg_method, nfeatures = params$integrate$hvg_nfeatures, verbose = FALSE)
}
run_cca_or_rpca <- function(objs, params, reduction) {
  objs <- lapply(objs, normalize_and_hvg, params = params)
  features <- Seurat::SelectIntegrationFeatures(object.list = objs, nfeatures = params$integrate$hvg_nfeatures)
  objs <- lapply(objs, function(o) {
    o <- Seurat::ScaleData(o, features = features, verbose = FALSE)
    Seurat::RunPCA(o, features = features, npcs = params$integrate$pca_dims, verbose = FALSE)
  })
  anchors <- Seurat::FindIntegrationAnchors(object.list = objs, anchor.features = features, reduction = reduction, dims = seq_len(params$integrate$pca_dims))
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(params$integrate$pca_dims))
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  integrated <- Seurat::RunPCA(integrated, npcs = params$integrate$pca_dims, verbose = FALSE)
  list(obj = integrated, reduction = "pca")
}
run_harmony <- function(objs, params, group_col) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony package required for method=harmony")
  }
  merged <- Reduce(function(a, b) merge(a, b), objs)
  merged <- normalize_and_hvg(merged, params = params)
  merged <- Seurat::ScaleData(merged, verbose = FALSE)
  merged <- Seurat::RunPCA(merged, npcs = params$integrate$pca_dims, verbose = FALSE)
  merged <- Seurat::RunHarmony(merged, group.by.vars = group_col, reduction = "pca", dims.use = seq_len(params$integrate$pca_dims))
  list(obj = merged, reduction = "harmony")
}
run_merge <- function(objs, params) {
  merged <- Reduce(function(a, b) merge(a, b), objs)
  merged <- normalize_and_hvg(merged, params = params)
  merged <- Seurat::ScaleData(merged, verbose = FALSE)
  merged <- Seurat::RunPCA(merged, npcs = params$integrate$pca_dims, verbose = FALSE)
  list(obj = merged, reduction = "pca")
}
cluster_and_umap <- function(obj, params, reduction_use) {
  dims_use <- seq_len(params$cluster$pca_dims)
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_use, dims = dims_use, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = params$cluster$resolution, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = reduction_use, dims = dims_use, min.dist = params$cluster$umap$min_dist, verbose = FALSE)
  obj
}
write_outputs <- function(obj, params) {
  outdir <- ensure_dir(params$output$outdir)
  prefix <- params$output$prefix
  if (as_bool(params$output$save_rds)) {
    saveRDS(obj, file = file.path(outdir, paste0(prefix, ".seurat.rds")))
  }
  if (as_bool(params$output$save_metadata_csv)) {
    utils::write.csv(obj[[]], file = file.path(outdir, paste0(prefix, ".metadata.csv")), row.names = TRUE)
  }
}
calc_markers <- function(obj, params) {
  if (!as_bool(params$markers$do)) return(NULL)
  Seurat::FindAllMarkers(
    obj,
    only.pos = as_bool(params$markers$only_pos),
    min.pct = params$markers$min_pct,
    logfc.threshold = params$markers$logfc_threshold
  )
}
`%||%` <- function(a, b) if (!is.null(a)) a else b
load_single_input <- function(path, params) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS file does not contain a Seurat object: ", path)
    return(obj)
  }
  if (!dir.exists(path)) stop("Input path does not exist: ", path)
  m <- Seurat::Read10X(path)
  Seurat::CreateSeuratObject(counts = m, min.cells = params$filter$min_cells_per_gene, project = params$input$project, assay = params$input$assay)
}
run_sc_pipeline <- function(inputs, params = default_params()) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package required")
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject package required")
  if (is.null(params$env$seurat_major)) {
    pv <- as.character(utils::packageVersion("Seurat"))
    params$env$seurat_major <- as.integer(strsplit(pv, "\\.")[[1]][[1]])
  }
  if (is.character(inputs) && length(inputs) == 1L && grepl(",", inputs, fixed = TRUE)) {
    inputs <- strsplit(inputs, ",", fixed = TRUE)[[1]]
  }
  if (!is.character(inputs) || length(inputs) < 1) stop("inputs must be a vector of paths")
  objs <- lapply(inputs, load_single_input, params = params)
  sample_names <- params$input$sample_names
  if (is.null(sample_names) || length(sample_names) != length(objs)) {
    sample_names <- paste0("sample", seq_along(objs))
  }
  for (i in seq_along(objs)) {
    objs[[i]]$sample <- sample_names[[i]]
  }
  objs <- lapply(objs, apply_qc_filters, params = params)
  per_sample_mode <- tolower(params$filter$per_sample$mode %||% "none")
  if (per_sample_mode %in% c("max_cells", "cap", "single_sample_max_cells")) {
    objs <- lapply(objs, function(o) cap_cells_per_sample(o, group_col = "sample", max_cells = params$filter$per_sample$max_cells, seed = params$filter$per_sample$seed))
  }
  method <- tolower(params$integrate$method)
  integ <- if (method == "harmony") {
    run_harmony(objs, params = params, group_col = "sample")
  } else if (method == "cca") {
    run_cca_or_rpca(objs, params = params, reduction = "cca")
  } else if (method == "rpca") {
    run_cca_or_rpca(objs, params = params, reduction = "rpca")
  } else if (method == "merge") {
    run_merge(objs, params = params)
  } else {
    stop("Unknown integrate.method: ", params$integrate$method)
  }
  obj <- integ$obj
  obj <- cluster_and_umap(obj, params = params, reduction_use = integ$reduction)
  write_outputs(obj, params = params)
  markers <- calc_markers(obj, params = params)
  if (!is.null(markers) && as_bool(params$output$save_markers_csv)) {
    outdir <- ensure_dir(params$output$outdir)
    utils::write.csv(markers, file = file.path(outdir, paste0(params$output$prefix, ".markers.csv")), row.names = FALSE)
  }
  invisible(obj)
}
if (sys.nframe() == 0) {
  params <- build_params_from_config(CONFIG)
  run_sc_pipeline(inputs = CONFIG$inputs, params = params)
}
