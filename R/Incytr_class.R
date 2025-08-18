#' The Incytr Class
#'
#' The Incytr object is created based on the following inputs: single-cell transcriptomic data matrix (normalized), meta data with cell cluster and condition labels.
#' When inputting an data matrix, genes should be in rows and cells in columns. rownames and colnames should be included.
#'
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom data.table data.table
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))

#' The key slots used in the Incytr object are described below.
#'
#' @slot data.raw raw count data matrix
#' @slot data normalized data matrix for CellChat analysis (Genes should be in rows and cells in columns)
#' @slot expr.bygroup a list of data frames storing the average gene expression level in each cluster in each condition
#' @slot pathways data frame storing the inferred L-R-EM-T pathways
#' @slot pathways_5steps data frame storing the inferred 5-step pathways for future exploration
#' @slot meta data frame storing the information associated with each cell
#' @slot sender a char for the cluster name used as the sender cell group
#' @slot receiver a char for the cluster name used as the receiver cell group
#' @slot idents factor of the cell cluster labels, with cell barcodes
#' @slot conditions a char for the condition labels (no more than two) to use in the analysis
#' @slot SigProb data frame storing the predicted signaling probability of each L-R-EM-T pathway in each condition
#' @slot p_value data frame storing the p-value of each L-R-EM-T pathway in each condition
#' @slot sc_FC list of data frames storing the gene log2FoldChange values calculated using DESeq2
#' @slot pr_FC data frame storing the fold change values of each cluster between conditions, calculated via the proteomics data
#' @slot ps_FC data frame storing the fold change values of each cluster between conditions, calculated via the ps data
#' @slot py_FC data frame storing the fold change values of each cluster between conditions, calculated via the py data
#' @slot Evaluation data frame storing the evaluation scores of the L-R-EM-T pathways
#' @slot kl data frame of the kinase library
#' @slot kl.pathways data frame storing all the inferred signaling-involved kinases (SiKs) of each pathways
#' @slot EI list storing the gene excluiveness index (EI) in each condition
#' @slot kl.explore data frame storing the inferred singaling-related kinases (SrKs) and their EI
#' @slot options List of miscellaneous data, such as parameters used throughout analysis
#'
#' @exportClass Incytr
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
object <- methods::setClass("Incytr",
                            slots = c(data.raw = 'AnyMatrix',
                                      data = 'AnyMatrix',
                                      expr.bygroup = 'list',
                                      pathways = "data.frame",
                                      pathways_5steps = "data.frame",
                                      meta = "data.frame",
                                      sender = 'character',
                                      receiver = 'character',
                                      idents = "AnyFactor",
                                      conditions = 'character',
                                      SigProb = "data.frame",
                                      p_value = "data.frame",
                                      sc_FC = 'list',
                                      pr_FC = "data.frame",
                                      ps_FC = "data.frame",
                                      py_FC = "data.frame",
                                      Evaluation = "data.frame",
                                      kl = "data.frame",
                                      kl.pathways = "data.frame",
                                      EI = 'list',
                                      kl.explore = "data.frame",
                                      options = "list")
)


#' Create a new Incytr object from a data matrix
#'
#' @param object a normalized (NOT count) data matrix (genes by cells)
#' @param meta a data frame (rows are cells with rownames) consisting of the cell cluster and condition labels
#' @param sender a cluster name used as the sender cell group
#' @param receiver a cluster name used as the receiver cell group
#' @param group.by a char name of the variable in meta data, defining cell groups.
#' If input is a data matrix and group.by is NULL, the input `meta` should contain a column named 'labels',
#' If input is a Seurat or SingleCellExperiment object, USER must provide `group.by` to define the cell groups. e.g, group.by = "ident" for Seurat object
#' @param conditions condition labels (no more than two) to use in the analysis.
#' If the input meta does not have the 'condition' variable, it will be created (as a factor) with one condition label 'condition1',
#' Then if the input conditions is NULL, the factor levels of the variable 'condition' will be used.
#' @param assay Assay to use when the input is a Seurat object
#' @param do.sparse whether use sparse format
#'
#' @return an Incytr object
#' @export
#' @importFrom methods as new
#' @importFrom Matrix t
create_Incytr <- function(object,
                          meta = NULL,
                          sender = NULL,
                          receiver = NULL,
                          group.by = NULL,
                          conditions = NULL,
                          assay = NULL,
                          do.sparse = T){

  # data matrix as input
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    message("Create an Incytr object from a data matrix")
    data <- object
    if (is.null(group.by)) {
      group.by <- "labels"
    }
  }

  if ( is.null(sender) | is.null(receiver) ){
    stop("Please specify the 'sender' and 'receiver'!")
  }else if ( !(sender %in% levels(factor(meta[, group.by]))) | !(receiver %in% levels(factor(meta[, group.by]))) ){
    stop("The inputs 'sender' and 'receiver' must be the clusters!")
  }

  if (!is.null(meta)) {
    meta <- as.data.frame(x = meta)
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame")
    }
    if (!identical(rownames(meta), colnames(data))) {
      cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
      warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }

  if( !"condition" %in% colnames(meta) ){
    meta$condition = "condition1"
  }

  meta$condition <- factor(meta$condition)
  if (is.null(conditions)){
    conditions <- as.character(levels(meta$condition))
  }

  idents = factor(meta[, group.by])
  names(idents) <- colnames(data)

  object <- methods::new(Class = "Incytr",
                         data = data,
                         meta = meta,
                         sender = sender,
                         receiver = receiver,
                         idents = idents,
                         conditions = conditions)
  object@options$mode <- "single"
  return(object)

}








