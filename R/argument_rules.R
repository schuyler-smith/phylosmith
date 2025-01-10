#' Function to check inputs to arguments
#'
#' Checks that arguments passed to functions are correct types.
#' @author Schuyler D. Smith
#' @useDynLib phylosmith
#' @usage check_args()
#' @param ... Any argument to a function in the phylosmith package
#' 
check_args <- function(...) {
  args <- list(...)
  if ("phyloseq_obj" %in% names(args)) {
    check_phyloseq_object(args[["phyloseq_obj"]])
  }
  if ("otu_table" %in% names(args)) {
    otu_table <- args[["otu_table"]]
    check_phyloseq(otu_table)
  }
  if ("tax_table" %in% names(args)) {
    tax_table <- args[["tax_table"]]
    check_phyloseq(tax_table)
  }
  if ("sam_data" %in% names(args)) {
    sam_data <- args[["sam_data"]]
    check_phyloseq(sam_data)
  }
  if ("refseq" %in% names(args)) {
    refseq <- args[["refseq"]]
    check_phyloseq(refseq)
  }
  if ("phytree" %in% names(args)) {
    phytree <- args[["phytree"]]
    check_phyloseq(phytree)
  }

  if ("treatment" %in% names(args)) {
    treatment <- args[["treatment"]]
    check_treatment(args[["phyloseq_obj"]], treatment)
  }
  if ("variables" %in% names(args)) {
    variables <- args[["variables"]]
    check_treatment(args[["phyloseq_obj"]], variables)
  }
  if ("replicate_samples" %in% names(args)) {
    replicate_samples <- args[["replicate_samples"]]
    check_treatment(args[["phyloseq_obj"]], replicate_samples)
  }
  if ("subset" %in% names(args)) {
    subset <- args[["subset"]]
    check_subset(args[["phyloseq_obj"]], args[["treatment"]], subset)
  }
  if ("labels" %in% names(args)) {
    labels <- args[["labels"]]
    check_treatment(args[["phyloseq_obj"]], labels)
  }

  if ("classification" %in% names(args)) {
    classification <- args[["classification"]]
    check_classification(args[["phyloseq_obj"]], classification)
  }

  if ("frequency" %in% names(args)) {
    frequency <- args[["frequency"]]
    check_frequency(frequency)
  }
  if ("p" %in% names(args)) {
    p <- args[["p"]]
    check_frequency(p)
  }
  if ("rho" %in% names(args)) {
    rho <- args[["rho"]]
    check_frequency(rho)
  }
  if ("abundance_threshold" %in% names(args)) {
    abundance_threshold <- args[["abundance_threshold"]]
    check_frequency(abundance_threshold)
  }

  if ("by_treatment" %in% names(args)) {
    by_treatment <- args[["by_treatment"]]
    check_boolean(by_treatment)
  }
  if ("hierarchical" %in% names(args)) {
    hierarchical <- args[["hierarchical"]]
    check_boolean(hierarchical)
  }
  if ("use_taxonomic_names" %in% names(args)) {
    use_taxonomic_names <- args[["use_taxonomic_names"]]
    check_boolean(use_taxonomic_names)
  }
  if ("below" %in% names(args)) {
    below <- args[["below"]]
    check_boolean(below)
  }
  if ("drop_samples" %in% names(args)) {
    drop_samples <- args[["drop_samples"]]
    check_boolean(drop_samples)
  }
  if ("na_rm" %in% names(args)) {
    na_rm <- args[["na_rm"]]
    check_boolean(na_rm)
  }
  if ("cluster" %in% names(args)) {
    cluster <- args[["cluster"]]
    check_boolean_or_frequency(cluster)
  }
  if ("relative_abundance" %in% names(args)) {
    relative_abundance <- args[["relative_abundance"]]
    check_boolean(relative_abundance)
  }
  if ("points" %in% names(args)) {
    points <- args[["points"]]
    check_boolean(points)
  }
  if ("verbose" %in% names(args)) {
    verbose <- args[["verbose"]]
    check_boolean(verbose)
  }
  if ("merge" %in% names(args)) {
    merge <- args[["merge"]]
    check_boolean(merge)
  }
  if ("grid" %in% names(args)) {
    grid <- args[["grid"]]
    check_boolean(grid)
  }

  if ("treatment_labels" %in% names(args)) {
    labels <- args[["treatment_labels"]]
    if (!is.null(labels)){check_strings(labels)}
  }
  if ("sample_labels" %in% names(args)) {
    labels <- args[["sample_labels"]]
    if (!is.null(labels)){check_strings(labels)}
  }
  if ("classification_labels" %in% names(args)) {
    labels <- args[["classification_labels"]]
    if (!is.null(labels)){check_strings(labels)}
  }
  if ("nodes_of_interest" %in% names(args)) {
    nodes_of_interest <- args[["nodes_of_interest"]]
    check_strings(nodes_of_interest)
  }
  
  if ("dist_method" %in% names(args)) {
    method <- args[["dist_method"]]
    options <- c("euclidian", "manhattan", "canberra", "bray", "kulczynski",
      "gower", "morisita", "horn", "mountford", "jaccard", "raup", "binomial",
      "chao", "altGower", "cao", "mahalanobis", "clark")
    check_options(method, options)
  }
  if ("corr_method" %in% names(args)) {
    method <- args[["corr_method"]]
    options <- c("pearson", "kendall", "spearman")
    check_options(method, options)
  }
  if ("diversity_index" %in% names(args)) {
    index <- args[["diversity_index"]]
    options <- c("shannon", "simpson", "invsimpson")
    check_options(index, options)
  }
  if ("transformation" %in% names(args)) {
    transformation <- args[["transformation"]]
    options <- c("none", "relative_abundance", "mean", "median", "sd", 
      "log", "log10", "log1p", "log2",
      "asn", "atanh", "boxcox", "exp", "identity", "logit",
      "probability", "probit", "reciprocal", "reverse", "sqrt")
    check_options(transformation, options)
  }

  if ("permutations" %in% names(args)){
    permutations <- args[["permutations"]]
    check_positive_numeric(permutations)
  }  
  if ("buffer" %in% names(args)){
    buffer <- args[["buffer"]]
    check_positive_numeric(buffer)
  }
  if ("x" %in% names(args)){
    x <- args[["x"]]
    check_positive_numeric(x)
  }
  if ("y" %in% names(args)){
    y <- args[["y"]]
    check_positive_numeric(y)
  }
  if ("sig_fig" %in% names(args)){
    sig_fig <- args[["sig_fig"]]
    check_positive_numeric(sig_fig)
  }
  if ("perplexity" %in% names(args)){
    perplexity <- args[["perplexity"]]
    check_positive_numeric(perplexity)
  }

  if ("cores" %in% names(args)) {
    check_cores(args[["cores"]])
  }
  if ("taxa_to_remove" %in% names(args)) {
    taxa_to_remove <- args[["taxa_to_remove"]]
    check_string(taxa_to_remove)
  }
  if ("taxa_to_extract" %in% names(args)) {
    taxa_to_extract <- args[["taxa_to_extract"]]
    check_string(taxa_to_extract)
  }
  if ("nodes_of_interest" %in% names(args)) {
    nodes_of_interest <- args[["nodes_of_interest"]]
    check_string(nodes_of_interest)
  }
  if ("n" %in% names(args)) {
    check_n(args[["n"]])
  }
  if ("co_occurrence_table" %in% names(args)) {
    check_co_occurrence_table(args[["co_occurrence_table"]])
  }

  if ("colors" %in% names(args)) {
    colors <- args[["colors"]]
    check_colors(colors)
  }
  if ("cluster_colors" %in% names(args)) {
    cluster_colors <- args[["cluster_colors"]]
    check_colors(cluster_colors)
  }
  if ("node_colors" %in% names(args)) {
    node_colors <- args[["node_colors"]]
    check_colors(node_colors)
  }
  
}



check_phyloseq_object <- function(phyloseq_obj) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class object", call. = FALSE)
  }
}
check_phyloseq <- function(phyloseq_obj) {
  attr <- deparse(substitute(phyloseq_obj))
  if (is.null(phyloseq::access(phyloseq_obj, attr))) {
    stop("`phyloseq_obj` must contain `", attr, "()` information",
      call. = FALSE)
  }
}

check_treatment <- function(phyloseq_obj, arg) {
  arg_name <- deparse(substitute(arg))
  if (!(is.null(arg)) &
    is.null(phyloseq::access(phyloseq_obj, "sam_data"))) {
  stop("`", arg_name, "` provided but `phyloseq_obj` does not contain 
sample_data()", call. = FALSE)
    }
  columns <- colnames(phyloseq::access(phyloseq_obj, "sam_data"))
  if (any(!(arg %in% columns))) if (any(!(arg %in% 'Sample'))) {
    stop("`", arg_name, "` must be at least one column name from the
phyloseq::sample_data(phyloseq_obj) or specified as 'Sample'", call. = FALSE)
  }
}

check_subset <- function(phyloseq_obj, treatment, subset) {
  check_treatment(phyloseq_obj, treatment)
  if (!is.null(subset)){
    if (is.null(treatment)) {
      stop("`subset` provided but `treatment` is null", call. = FALSE)
    }
    phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
    treatment_name <- paste(treatment, collapse = sep)
    treatment_vals <- data.table::data.table(as(phyloseq::access(phyloseq_obj, "sam_data"), "data.frame"))[[treatment_name]]
    treatment_vals <- unique(treatment_vals)
    if (any(!subset %in% treatment_vals)){
      missing <- paste(subset[!subset %in% treatment_vals], collapse = "', '")
      stop("`subset` must be values contained in `treatment`. c('", 
      paste(missing, collapse = "', '"), "') not found.", call. = FALSE)
    }
  }
}

check_classification <- function(phyloseq_obj, classification) {
  if (!(is.null(classification)) &
    is.null(phyloseq::access(phyloseq_obj, "tax_table"))) {
    stop(
"`classification` provided by `phyloseq_obj` does not contain tax_table()",
      call. = FALSE)
  }
  if (any(!(classification %in% colnames(phyloseq::access(phyloseq_obj,
    "tax_table"))))) {
    stop(
"`classification` must be at least one column name from the tax_table()",
    call. = FALSE)
  }
}

check_boolean <- function(arg) {
    arg_name <- deparse(substitute(arg))
    if (!(is.logical(arg))) {
      stop("`", arg_name, "` must must be either`TRUE`, or `FALSE`",
      call. = FALSE)
    }
}

check_frequency <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (is.null(arg)){return()}
  if (!(is.numeric(arg)) |
    !(arg >= 0 & arg <= 1)) {
      stop("`", arg_name, "` must be a numeric value between and 1",
      call. = FALSE)
  }
}

check_boolean_or_frequency <- function(arg) {
  if (!is.logical(arg) && (!is.numeric(arg) || arg < 0 || arg > 1)) {
    stop("`cluster` must be a boolean or a numeric value between 0 and 1.")
  }
}

check_string <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (!all(is.character(arg))) {
    stop("`", arg_name, "` must be string values", call. = FALSE)
  }
}

check_options <- function(arg, options) {
  arg_name <- deparse(substitute(arg))
  tryCatch(
    match.arg(arg, options),
    error = function(e) {
      stop("Invalid argument '", arg, "' given to '", arg_name, "'. 
  Valid arguments are: '", paste(options, collapse = "', '"))
    }
  )
}

check_positive_numeric <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (!(is.numeric(arg)) | !(arg >= 0)) {
    stop("`", arg_name, "` must be a numeric value greater than 0.", 
      call. = FALSE)
  }
}
check_cores <- function(cores) {
  if (!(is.numeric(cores)) |
      !(cores >= 0 & cores <= (parallel::detectCores() - 1))) {
    stop(
"`cores` must be a numeric value between 0 and ", (parallel::detectCores() - 1),
"\n(upper limit is set by the cores available on machine)",
      call. = FALSE
    )
  }
}
check_permuted_rhos <- function(permuted_rhos) {
  if (!(is.data.frame(permuted_rhos))) {
    stop("`permuted_rhos` must be at data.frame object", call. = FALSE)
  }
}
check_n <- function(n) {
  if (!(is.numeric(n)) & n != "all") {
    stop(
"`n` must be either 'all' or a numeric value less than the number of treatments
being compared", call. = FALSE )
  }
}

check_co_occurrence_table <- function(co_occurrence_table) {
  if (!(is.null(co_occurrence_table)) & !(is.data.frame(co_occurrence_table))) {
    stop("`co_occurrence_table` must be at data.frame object", call. = FALSE)
  }
}

check_strings <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (!is.character(arg) || any(!is.character(arg))) {
    stop("`", arg_name, "` must be a string or vector of strings.", 
      call. = FALSE)
  }
}

check_colors <- function(arg) {
  arg_name <- deparse(substitute(arg))
  allowed_colors <- c(names(RColorBrewer::brewer.pal.info), colors())
  if (is.character(arg)) {
    if (!(arg %in% allowed_colors) && arg != "default") {
      stop("`", arg_name, "` must be a palette name from `RColorBrewer`, a 
vector of `colors` that R accepts, or 'default'.", call. = FALSE)
    }
  } else if (!is.vector(arg) || any(!arg %in% allowed_colors)) {
    stop("`", arg_name, "` must be a palette name from `RColorBrewer`, a vector
of `colors` that R accepts, or 'default'.", call. = FALSE)
  }
}

# check_colors <- function(colors) {
#   if (any(colors == 'default')) {
#   }
#   if (any(colors %in% 'rev') | any(colors %in% 'reverse')) {
#   }
#   if (any(!(colors %in% colors()))) {
#     if (any(colors %in% rownames(RColorBrewer::brewer.pal.info))) {
#     } 
#   }
# }