#' Infer the highly expressed genes of a select cluster
#'
#' @param object Incytr object
#' @param meta data frame storing the information associated with each cell
#' @param group.by variable name used to label the cells
#' @param group.select cluster name for inferring highly expressed genes
#' @param geneuse a list of genes used for the inference, use NULL if all genes are using. The default value is NULL.
#' @param cutoff_percentile the persentile used to identifiy the highly expressed genes. The default value is cutoff_percentile = 0.5, which means if a gene's average expression level is higher than the 50th persentile then it is identified as a highly expressed gene.
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#'
#' @return a data frame
#' @export
Find_highexp_gene <- function(object,
                              meta,
                              group.by,
                              group.select,
                              geneuse = NULL,
                              cutoff_percentile = 0.5,
                              mean_method = NULL){

  # data matrix as input
  data <- object
  if (is.null(group.by)) {
    group.by <- "labels"
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

  ##############################################################################

  # find the cell barcodes
  meta = meta[meta[, group.by] == group.select, ]
  celluse = rownames(meta)

  data <- data[ , celluse]

  # find all the genes that will be used
  if( is.null(geneuse) ){
    geneuse <- rownames(data)
  }else{
    geneuse <- intersect(rownames(data), geneuse)
  }

  if( length(celluse)==0 | length(geneuse)==0 ){
    stop("No cell/gene is used!")
  }

  M1 = as.matrix(data[geneuse, celluse])
  if(length(geneuse)==1){
    df1 <- as.data.frame(M1)
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M1))
  }
  df1$group_label <- group.select

  q = c(.25, .5, .75)

  # for condition 1
  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
  }

  output = as.data.frame(t(yy[, -1]))
  colnames(output) = q1$group_label
  gene_symbol = rownames(output)
  output <- cbind(gene_symbol, output)

  if(is.null(cutoff_percentile)){
    return(output)
  }else{

    object.2 = object[object>0]
    cutoff_exp = quantile(object.2, cutoff_percentile)
    output = output[ output[ , 2] >= cutoff_exp, ]

    if(nrow(output)==0){
      stop("No gene found using current filter!")
    }
  }

  return(output)
}



#' Calculate the average expression level of each cluster in two conditions from the single-cell data
#'
#' @param object Incytr object
#' @param meta data frame storing the information associated with each cell
#' @param group.by variable name used to label the cells
#' @param conditions names of the conditions, or levels of the condition factor
#' @param omics a char to note the data type, only used in the variable names. Must be NULL or one of following: pr, ps, py.
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#'
#' @return a data frame
#' @export
dataprepare_Expr_bygroup <- function(object,
                                     meta = NULL,
                                     group.by = NULL,
                                     conditions = NULL,
                                     omics = NULL,
                                     mean_method = NULL){

  # data matrix as input
  data <- object
  if (is.null(group.by)) {
    group.by <- "labels"
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

  meta$condition <- factor(meta$condition)
  if (is.null(conditions)){
    conditions <- as.character(levels(meta$condition))
  }

  idents = factor(meta[, group.by])
  names(idents) <- colnames(data)

  # gather cell barcodes for each condition
  condition.label <- meta$condition
  names(condition.label) <- colnames(data)

  # select data of the condition 1
  barcodes.con1 <- names(condition.label)[condition.label == conditions[1]]

  # select data of the condition 2
  barcodes.con2 <- names(condition.label)[condition.label == conditions[2]]

  ##############################################################################

  # find the cell barcodes of condition1 and condition2
  celluse <- list(condition1 = barcodes.con1, condition2 = barcodes.con2)

  # find all the genes that will be used
  geneuse <- rownames(data)

  # expression data of condition 1 (row: cell; column: gene; the last column: cell group label)
  M1 = as.matrix(data[ geneuse, celluse$condition1])
  if(length(geneuse)==1){
    df1 <- as.data.frame(M1)
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M1))
  }
  df1$group_label <- fct_drop(idents[names(idents) %in% celluse$condition1])
  df1$group_label <- as.character(df1$group_label)

  # expression data of condition 2 (row: cell; column: gene; the last column: cell group label)
  M2 = as.matrix(data[ geneuse, celluse$condition2])
  if(length(geneuse)==1){
    df2 <- as.data.frame(M2)
  }else if(length(geneuse)>1){
    df2 <- as.data.frame(t(M2))
  }
  df2$group_label <- fct_drop(idents[names(idents) %in% celluse$condition2])
  df2$group_label <- as.character(df2$group_label)

  q = c(.25, .5, .75)

  # for condition 1
  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
  }

  output_con1 = as.data.frame(t(yy[, -1]))
  colnames(output_con1) = q1$group_label
  gene_symbol = rownames(output_con1)
  output_con1 <- cbind(gene_symbol, output_con1)

  # for condition 2
  if(is.null(mean_method)){

    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df2)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df2)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
  } else if(mean_method=="mean"){
    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df2)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
  }

  output_con2 = as.data.frame(t(yy[, -1]))
  colnames(output_con2) = q1$group_label
  gene_symbol = rownames(output_con2)
  output_con2 <- cbind(gene_symbol, output_con2)

  if( !is.null(omics) & length(intersect(omics, c("pr", "ps", "py")))==0 ){
    warning("The input 'omics' must be NULL or one of following: pr, ps, py. Otherwise it will be set to NULL.")
    omics = NULL
  }

  if(!is.null(omics)){
    colnames(output_con1)[2:ncol(output_con1)] = paste0(colnames(output_con1)[2:ncol(output_con1)], "_", omics)
    colnames(output_con2)[2:ncol(output_con2)] = paste0(colnames(output_con2)[2:ncol(output_con2)], "_", omics)
  }

  output = list(condition1 = output_con1, condition2 = output_con2)
  names(output) = conditions

  return(output)
}



#' Gather cell barcodes for each condition from the meta data
#'
#' @param object Incytr object
#'
#' @return a list
#' @export
barcodes_bycondition <- function(object){

  condition.label <- object@meta$condition
  names(condition.label) <- colnames(object@data)

  barcodes <- list()

  # select data of the condition 1
  barcodes.con1 <- names(condition.label)[condition.label == object@conditions[1]]
  barcodes[[1]] = barcodes.con1

  if( length(object@conditions) == 2 ){
    # select data of the condition 2
    barcodes.con2 <- names(condition.label)[condition.label == object@conditions[2]]
    barcodes[[2]] = barcodes.con2
  }else{
    barcodes[[2]] = character()
  }

  names(barcodes) = c("condition1", "condition2")

  return(barcodes)
}


#' Update the Incytr object when the inferred L-T paths have been changed/updated
#'
#' @param object Incytr object
#' @param Path list of the new path names
#'
#' @return an Incytr object
#' @export
object_update <- function(object, Path){

  Path.orignal = object@pathways$Path

  if( nrow(object@pathways)>0  ){
    object@pathways = object@pathways[object@pathways$Path %in% Path, ]
  }

  if( nrow(object@SigProb)>0  ){
    object@SigProb = object@SigProb[object@SigProb$Path %in% Path, ]
  }

  if( nrow(object@p_value)>0  ){
    object@p_value = object@p_value[object@p_value$Path %in% Path, ]
  }

  if( nrow(object@pr_FC)>0  ){
    df = object@pr_FC
    df$Path = Path.orignal
    df = df[df$Path %in% Path, 1:(ncol(df)-1)]
    object@pr_FC = df
  }

  if( nrow(object@ps_FC)>0  ){
    df = object@ps_FC
    df$Path = Path.orignal
    df = df[df$Path %in% Path, 1:(ncol(df)-1)]
    object@ps_FC = df
  }

  if( nrow(object@py_FC)>0  ){
    df = object@py_FC
    df$Path = Path.orignal
    df = df[df$Path %in% Path, 1:(ncol(df)-1)]
    object@py_FC = df
  }

  if( nrow(object@Evaluation)>0  ){
    object@Evaluation = object@Evaluation[object@Evaluation$Path %in% Path, ]
  }

  if( nrow(object@kl.pathways)>0  ){
    object@kl.pathways = object@kl.pathways[object@kl.pathways$Path %in% Path, ]

    if(length(object@EI)==2){
      k1 = as.character(object@kl.pathways$SiK_R_of_EM[!is.na(object@kl.pathways$SiK_R_of_EM)] )
      k2 = as.character(object@kl.pathways$SiK_R_of_T[!is.na(object@kl.pathways$SiK_R_of_T)] )
      k3 = as.character(object@kl.pathways$SiK_EM_of_T[!is.na(object@kl.pathways$SiK_EM_of_T)] )
      k4 = as.character(object@kl.pathways$SiK_EM_of_R[!is.na(object@kl.pathways$SiK_EM_of_R)] )
      k5 = as.character(object@kl.pathways$SiK_T_of_R[!is.na(object@kl.pathways$SiK_T_of_R)] )
      k6 = as.character(object@kl.pathways$SiK_T_of_EM[!is.na(object@kl.pathways$SiK_T_of_EM)] )
      k = unique(c(k1, k2, k3, k4, k5, k6))
      object@EI[[1]] = object@EI[[1]][object@EI[[1]]$Gene %in% k, ]
      object@EI[[2]] = object@EI[[2]][object@EI[[2]]$Gene %in% k, ]
    }
  }

  return(object)

}



#' Inferring signaling pathways based on the user-select genes & database
#'
#' @param object Incytr object
#' @param DB database used for the inference
#' @param gene.use_Sender a list of gene used as the sender genes (ligands). Used all the genes when it is NULL.
#' @param gene.use_Receiver a list of gene used as the receiver genes (receptors, EMs, and targets). Used all the genes when it is NULL.
#' @param ligand a list of gene used to filter the ligand
#' @param receptor a list of gene used to filter the receptor
#' @param em a list of gene used to filter the em
#' @param target a list of gene used to filter the target
#'
#' @return an Incytr object
#' @export
pathway_inference <- function(object,
                              DB,
                              gene.use_Sender = NULL,
                              gene.use_Receiver = NULL,
                              ligand = NULL,
                              receptor = NULL,
                              em = NULL,
                              target = NULL){

  # check if the gene.use is null
  if( is.null(gene.use_Sender) | is.null(gene.use_Receiver)){
    message("The input 'gene.use_Sender' or 'gene.use_Receiver' is null, set it to all the genes in the data.")
    gene_use = rownames(object@data)
  }
  # merge all the gene inputs
  gene_use = unique(c(gene.use_Sender, gene.use_Receiver))

  # filter the database
  Role <- c("Ligand", "Receptor", "EM", "Target")
  for (i in 1:3) {
    DB[[i]] = DB[[i]][ (DB[[i]]$from %in% gene_use)&(DB[[i]]$to %in% gene_use), 1:2 ]
    colnames(DB[[i]]) <- c(Role[i], Role[i+1])
  }

  SigPath <- inner_join(DB[[1]], DB[[2]], by = join_by(Receptor), multiple = "all")
  SigPath <- inner_join(SigPath, DB[[3]], by = join_by(EM), multiple = "all")

  # filter the pathways based on the Sender genes & Receiver genes
  SigPath <- SigPath[ (SigPath$Ligand %in% gene.use_Sender) &
                        (SigPath$Receptor %in% gene.use_Receiver) &
                        (SigPath$EM %in% gene.use_Receiver) &
                        (SigPath$Target %in% gene.use_Receiver), ]

  if(!is.null(ligand)){
    SigPath <- SigPath[ SigPath$Ligand %in% ligand, ]
  }

  if(!is.null(receptor)){
    SigPath <- SigPath[ SigPath$Receptor %in% receptor, ]
  }

  if(!is.null(em)){
    SigPath <- SigPath[ SigPath$EM %in% em, ]
  }

  if(!is.null(target)){
    SigPath <- SigPath[ SigPath$Target %in% target, ]
  }

  if(nrow(SigPath)==0){
    stop("No signal pathway found!")
  }

  # remove pathways with duplicated genes
  boolean_dup <- apply(SigPath, 1, duplicated)
  index = apply(boolean_dup, 2, sum)
  SigPath$index = index
  SigPath = SigPath[SigPath$index == 0, 1:4]

  SigPath$Path <- paste(SigPath$Ligand, SigPath$Receptor, SigPath$EM, SigPath$Target, sep = "*")
  SigPath <- SigPath[!duplicated(SigPath$Path), ]

  object@pathways <- SigPath
  return(object)

}


#' Calculate the average expression level in each cell group for both conditions
#'
#' @param object Incytr object
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#'
#' @return a data frame
#' @export
Expr_bygroup <- function(object,
                         mean_method = NULL){

  # find the cell barcodes of condition1 and condition2
  celluse <- barcodes_bycondition(object)

  # find all the genes that will be used
  geneuse <- unique(c(object@pathways$Ligand,
                      object@pathways$Receptor,
                      object@pathways$EM,
                      object@pathways$Target))

  q = c(.25, .5, .75)

  # expression data of condition 1 (row: cell; column: gene; the last column: cell group label)
  M1 = as.matrix(object@data[ geneuse, celluse$condition1])
  if(length(geneuse)==1){
    df1 <- as.data.frame(M1)
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M1))
  }
  df1$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition1])
  df1$group_label <- as.character(df1$group_label)

  # for condition 1
  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    object@expr.bygroup[[1]] = as.data.frame(t(yy))
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    object@expr.bygroup[[1]] = as.data.frame(t(yy[, -1]))
  }

  # object@expr.bygroup[[1]] = as.data.frame(t(yy[, -1]))
  colnames(object@expr.bygroup[[1]]) <- q1$group_label
  object@expr.bygroup[[1]]$Gene = rownames(object@expr.bygroup[[1]])

  # expression data of condition 2 (row: cell; column: gene; the last column: cell group label)
  if( length(celluse$condition2)>0 ){
    M2 = as.matrix(object@data[ geneuse, celluse$condition2])
    if(length(geneuse)==1){
      df2 <- as.data.frame(M2)
    }else if(length(geneuse)>1){
      df2 <- as.data.frame(t(M2))
    }
    df2$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition2])
    df2$group_label <- as.character(df2$group_label)

    # for condition 2
    if(is.null(mean_method)){

      q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
      q2 <- setDT(df2)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
      q3 <- setDT(df2)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

      yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
      # rownames(yy) <- q1$group_label
      object@expr.bygroup[[2]] = as.data.frame(t(yy))
    } else if(mean_method=="mean"){
      q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
      yy <- setDT(df2)[, lapply(.SD, mean), keyby = group_label]
      # rownames(yy) <- q1$group_label
      object@expr.bygroup[[2]] = as.data.frame(t(yy[, -1]))
    }

    # object@expr.bygroup[[2]] = as.data.frame(t(yy[, -1]))
    colnames(object@expr.bygroup[[2]]) <- q1$group_label
    object@expr.bygroup[[2]]$Gene = rownames(object@expr.bygroup[[2]])
  }else{
    object@expr.bygroup[[2]] = data.frame()
  }

  names(object@expr.bygroup) = c("condition1", "condition2")

  return(object)
}


#' Define hill function
#'
#' @param x variable
#' @param K,N parameters
#'
#' @return a value
#' @export
hill <- function(x,K,N){
  return(x^N/(x^N+K^N))
}


#' Calculate the log2FoldChang & adjlog2FC for the data
#'
#' @param df a data frame with the following variables: "gene_symbol", "condition1", "condition2"
#' @param correction the value added to all the data to avoid division over zero, the default value is 0.0001
#' @param q the percentile used to determine "low expression level" and is used when calculating the adjlog2FC, the default value is 0.75
#'
#' @return a data frame
#' @export
Cal_foldchange <- function(df, correction = 0.0001, q = NULL){

  no.Zero <- sum(df[, 2:3]==0)

  if(is.null(correction)){
    correction = 0.0001
  }

  if(isTRUE(no.Zero!=0) & correction==0 ){
    stop("There is '0' in the data, please consider using 'correction' option.")
  }else if(isTRUE(no.Zero!=0)){
    # print("The input data has been corrected.")
    df[, 2:3] = df[, 2:3] + correction
  }

  # calculate the log2 fold change
  df$log2FC = log2(df$condition1/df$condition2)

  # calculate the adjusted log2FC
  if(is.null(q)){
    q = 0.75
  }

  # transform df into a vector
  df1 = unlist(df[, 2:3])
  # df1 = as.numeric(df1[df1 != "NA"])
  th = quantile(df1, q, na.rm=TRUE)

  N = 2
  Vmax = apply(df[, 2:3], 1, max)
  a = data.frame(V1 = 2*(Vmax^N)/(Vmax^N + th^N), V2 = 1)
  adj = apply(a, 1, min)

  df$aFC = df$log2FC*adj

  return(df)

}


#' Calculate the signaling probability and its fold change for each L-T pathway based on the gene expression data
#'
#' @param object Incytr object
#' @param K parameter for the half occupation level in the hill function, the default value is 0.5
#' @param N Hill coefficient, the defaulty value is 2
#' @param cutoff_SigProb cutoff value for selecting the pathways that pass this value in at least one condition. In default, it is NULL.
#' @param correction the value added to all the data to avoid division over zero, the default value is 0.0001
#' @param q the percentile used to determine "low expression level" and is used when calculating the adjlog2FC, the default value is 0.75
#'
#' @return an Incytr object
#' @export
Cal_SigProb <- function(object, K = 0.5, N = 2, cutoff_SigProb = NULL,
                        correction = 0.0001, q = NULL){

  sender = object@sender
  receiver = object@receiver

  # expression data for the L in the sender group (condition 1)
  scdata0 <- object@expr.bygroup$condition1[ , c(sender,"Gene")]
  scdata.l <- data.frame(Ligand = scdata0$Gene, Ligand.value = scdata0[ , 1])

  ##############################################################################
  # expression data for the R, EM, T in the receiver group (condition 1)
  scdata <- object@expr.bygroup$condition1[ , c(receiver,"Gene")]

  scdata.r <- data.frame(Receptor = scdata$Gene, Receptor.value = scdata[ , 1])
  scdata.em <- data.frame(EM = scdata$Gene, EM.value = scdata[ , 1])
  scdata.t <- data.frame(Target = scdata$Gene, Target.value = scdata[ , 1])

  # filter the pathways
  net.filtered <- object@pathways
  net.filtered$Sender.group <- sender
  net.filtered$Receiver.group <- receiver
  net.filtered <- net.filtered[ net.filtered$Ligand %in% scdata.l$Ligand &
                                  net.filtered$Receptor %in% scdata.r$Receptor &
                                  net.filtered$EM %in% scdata.em$EM &
                                  net.filtered$Target %in% scdata.t$Target , ]

  # combine net.filtered and the data
  net.filtered <- left_join(net.filtered, scdata.l, by = join_by(Ligand))
  net.filtered <- left_join(net.filtered, scdata.r, by = join_by(Receptor))
  net.filtered <- left_join(net.filtered, scdata.em, by = join_by(EM))
  net.filtered <- left_join(net.filtered, scdata.t, by = join_by(Target))

  df1 <- net.filtered

  # calculate the signaling probability
  layer1 = df1$Ligand.value*df1$Receptor.value
  layer2 = df1$Receptor.value*df1$EM.value
  layer3 = df1$EM.value*df1$Target.value
  P1 <- hill(layer1,K,N)*hill(layer2,K,N)*hill(layer3,K,N)
  df1$SigProb <- P1

  ##############################################################################
  # expression data for the L in the sender group (condition 2)
  if( nrow(object@expr.bygroup$condition2)>0 ){
    scdata0 <- object@expr.bygroup$condition2[ , c(sender,"Gene")]
    scdata.l <- data.frame(Ligand = scdata0$Gene, Ligand.value = scdata0[ , 1])

    # expression data for the R, EM, T in the receiver group (condition 2)
    scdata <- object@expr.bygroup$condition2[ , c(receiver,"Gene")]

    scdata.r <- data.frame(Receptor = scdata$Gene, Receptor.value = scdata[ , 1])
    scdata.em <- data.frame(EM = scdata$Gene, EM.value = scdata[ , 1])
    scdata.t <- data.frame(Target = scdata$Gene, Target.value = scdata[ , 1])

    # filter the pathways
    net.filtered <- object@pathways
    net.filtered$Sender.group <- sender
    net.filtered$Receiver.group <- receiver
    net.filtered <- net.filtered[ net.filtered$Ligand %in% scdata.l$Ligand &
                                    net.filtered$Receptor %in% scdata.r$Receptor &
                                    net.filtered$EM %in% scdata.em$EM &
                                    net.filtered$Target %in% scdata.t$Target , ]

    # combine net.filtered and the data
    net.filtered <- left_join(net.filtered, scdata.l, by = join_by(Ligand))
    net.filtered <- left_join(net.filtered, scdata.r, by = join_by(Receptor))
    net.filtered <- left_join(net.filtered, scdata.em, by = join_by(EM))
    net.filtered <- left_join(net.filtered, scdata.t, by = join_by(Target))

    df2 <- net.filtered

    # calculate the signaling probability
    layer1 = df2$Ligand.value*df2$Receptor.value
    layer2 = df2$Receptor.value*df2$EM.value
    layer3 = df2$EM.value*df2$Target.value
    P2 <- hill(layer1,K,N)*hill(layer2,K,N)*hill(layer3,K,N)
    df2$SigProb <- P2
  }else{
    P2 = rep(NA, length(df1$Path))
  }

  ##############################################################################
  object@SigProb= data.frame(Path = df1$Path,
                             SigProb_condition1 = P1,
                             SigProb_condition2 = P2)

  names(object@SigProb) = c("Path",
                            paste0("SigProb_", object@conditions[1]),
                            paste0("SigProb_", object@conditions[2]))
  # save the SigProb info to a new data frame
  if( nrow(object@expr.bygroup$condition2)>0 ){
    # calculate the log2FC & aFC
    df = data.frame(gene_symbol = object@SigProb$Path,
                    condition1 = object@SigProb[ , 2],
                    condition2 = object@SigProb[ , 3])

    df = Cal_foldchange(df, correction = correction, q = q)

    object@SigProb$log2FC = df$log2FC
    object@SigProb$aFC = df$aFC

    # cutoff_SigProb value: select the paths that have high SigProb in AT LEAST one condition
    if( !is.null(cutoff_SigProb) ){
      df.update = object@SigProb[ (object@SigProb[ , 2] > cutoff_SigProb) |
                                    (object@SigProb[ , 3] > cutoff_SigProb),  ]
      if(nrow(df.update)>0){
        Path.update = df.update$Path
        object = object_update(object = object, Path = Path.update)
      }else{
        stop("No signaling pathway inferred using current filter, the filter has been ignored.")
      }
    }
  }else{
    # cutoff_SigProb value: select the paths that have high SigProb in AT LEAST one condition
    if( !is.null(cutoff_SigProb) ){
      df.update = object@SigProb[ (object@SigProb[ , 2] > cutoff_SigProb),  ]
      if(nrow(df.update)>0){
        Path.update = df.update$Path
        object = object_update(object = object, Path = Path.update)
      }else{
        stop("No signaling pathway inferred using current filter, the filter has been ignored.")
      }
    }
  }

  return(object)

}


#' Calculate the fold change (log2FoldChange) for the genes in each L-T pathway based on the scRNA-seq data using DESeq2
#'
#' @param object Incytr object
#' @param count.matrix the count matrix used to calculate the fold change
#' @param pseudocount add 1 pseudocount to avoid division over zero issue. Default is TRUE.
#' @param fitType an option from the DESeq2, in default, fitType = c("parametric", "local", "mean", "glmGamPoi")
#' @param sfType an option from the DESeq2, in default, sfType = c("ratio", "poscounts", "iterate")
#'
#' @return an Incytr object
#' @export
Cal_scFC <- function(object, count.matrix = NULL,
                     pseudocount = TRUE,
                     fitType = c("parametric", "local", "mean", "glmGamPoi"),
                     sfType = c("ratio", "poscounts", "iterate")){

  if(is.null(count.matrix)){
    stop("The input 'count.matrix' cannot be null.")
  }else if(nrow(count.matrix)==0 | ncol(count.matrix)==0){
    stop("The input 'count.matrix' cannot be empty.")
  }else{
    object@data.raw = count.matrix
  }

  # library(DESeq2)
  # gene lists
  gene.sender = unique(object@pathways$Ligand)
  gene.receiver = unique(c(object@pathways$Receptor, object@pathways$EM, object@pathways$Target))
  # cell barcodes
  cell.sender = names(object@idents[object@idents == object@sender])
  cell.receiver = names(object@idents[object@idents == object@receiver])
  # count matrix for sender and receiver
  Data.sender = object@data.raw[rownames(object@data.raw) %in% gene.sender,
                                colnames(object@data.raw) %in% cell.sender]
  Data.receiver = object@data.raw[rownames(object@data.raw) %in% gene.receiver,
                                  colnames(object@data.raw) %in% cell.receiver]
  # add a pseudo-count of 1 to all counts
  if(isTRUE(pseudocount)){
    Data.sender = Data.sender + 1
    Data.receiver = Data.receiver + 1
  }
  # condition factors for sender and receiver cells
  meta = object@meta
  rownames(meta) = colnames(object@data)
  condition.sender = meta[cell.sender, ]$condition
  condition.receiver = meta[cell.receiver, ]$condition
  # set up colData
  colData.sender <- data.frame(row.names = colnames(Data.sender), condition = condition.sender)
  colData.receiver <- data.frame(row.names = colnames(Data.receiver), condition = condition.receiver)
  # DESeq analysis
  dds.sender <- DESeqDataSetFromMatrix(Data.sender, colData.sender, design = ~condition)

  try.sender <- try( dds.sender <- DESeq(dds.sender, fitType = fitType, sfType = sfType) )
  if( "try-error" %in% class(try.sender) ){
    dds.sender <- estimateSizeFactors(dds.sender)
    dds.sender <- estimateDispersionsGeneEst(dds.sender)
    dispersions(dds.sender) <- mcols(dds.sender)$dispGeneEst
    dds.sender <- nbinomWaldTest(dds.sender)
  }

  dds.receiver <- DESeqDataSetFromMatrix(Data.receiver, colData.receiver, design = ~condition)

  try.receiver <- try( dds.receiver <- DESeq(dds.receiver, fitType = fitType, sfType = sfType) )
  if( "try-error" %in% class(try.receiver) ){
    dds.receiver <- estimateSizeFactors(dds.receiver)
    dds.receiver <- estimateDispersionsGeneEst(dds.receiver)
    dispersions(dds.receiver) <- mcols(dds.receiver)$dispGeneEst
    dds.receiver <- nbinomWaldTest(dds.receiver)
  }
  # export the results
  res.sender = results(dds.sender, contrast = c("condition", object@conditions))
  res.receiver = results(dds.receiver, contrast = c("condition", object@conditions))

  object@sc_FC[[1]] = res.sender
  object@sc_FC[[2]] = res.receiver

  names(object@sc_FC) = c(paste0("DESeq2_", object@sender), paste0("DESeq2_", object@receiver))

  return(object)
}


#' Integrate the proteomics and phosphorylation data to the Incytr object, and calculate the log2FC, adjlog2FC
#'
#' @param object Incytr object
#' @param pr.data_condition1 data frame of the proteomics data in condition1, must have the variable "gene_symbol".
#' @param pr.data_condition2 data frame of the proteomics data in condition2, must have the variable "gene_symbol"
#' @param pr.correction the correction used when calculating the fold change to avoid division over zero issue, the default value is NULL
#' @param pr.q the percentile used to determine "low expression level" and is used when calculating the adjlog2FC, the default value is 0.75
#' @param ps.data_condition1 similar to the above
#' @param ps.data_condition2 similar to the above
#' @param ps.correction similar to the above
#' @param ps.q similar to the above
#' @param py.data_condition1 similar to the above
#' @param py.data_condition2 similar to the above
#' @param py.correction similar to the above
#' @param py.q similar to the above
#'
#' @return an Incytr object
#' @export
Integr_multiomics <- function(object,
                              pr.data_condition1 = NULL,
                              pr.data_condition2 = NULL,
                              pr.correction = NULL,
                              pr.q = NULL,
                              ps.data_condition1 = NULL,
                              ps.data_condition2 = NULL,
                              ps.correction = NULL,
                              ps.q = NULL,
                              py.data_condition1 = NULL,
                              py.data_condition2 = NULL,
                              py.correction = NULL,
                              py.q = NULL){

  # option 1: pr data
  if(is.null(pr.data_condition1) | is.null(pr.data_condition2)){
    message("Proteomics data not found for at least one condition, this part has been skipped.")
  }else{
    # pr data from two conditions for sender & receiver
    pr.data_s1 = pr.data_condition1[, c("gene_symbol", paste0(object@sender,"_pr"))]
    pr.data_s2 = pr.data_condition2[, c("gene_symbol", paste0(object@sender,"_pr"))]
    pr.data_r1 = pr.data_condition1[, c("gene_symbol", paste0(object@receiver,"_pr"))]
    pr.data_r2 = pr.data_condition2[, c("gene_symbol", paste0(object@receiver,"_pr"))]
    # merge data of two conditions
    pr.data_s = inner_join(pr.data_s1, pr.data_s2, by = c("gene_symbol"))
    pr.data_r = inner_join(pr.data_r1, pr.data_r2, by = c("gene_symbol"))

    # normaliation between assays via limma
    xs <- as.matrix(pr.data_s[ , 2:3])
    ys <- normalizeBetweenArrays(xs)
    pr.data_s[ , 2:3] = ys
    colnames(pr.data_s) = c("gene_symbol", "condition1", "condition2")

    xr <- as.matrix(pr.data_r[ , 2:3])
    yr <- normalizeBetweenArrays(xr)
    pr.data_r[ , 2:3] = yr
    colnames(pr.data_r) = c("gene_symbol", "condition1", "condition2")

    pr.data_s = Cal_foldchange(pr.data_s, correction = pr.correction, q = pr.q)[ , c("gene_symbol", "log2FC", "aFC")]
    pr.data_r = Cal_foldchange(pr.data_r, correction = pr.correction, q = pr.q)[ , c("gene_symbol", "log2FC", "aFC")]

    df <- object@pathways

    df <- left_join(df, pr.data_s, by = c("Ligand" = "gene_symbol") )
    df <- left_join(df, pr.data_r, by = c("Receptor" = "gene_symbol") )
    df <- left_join(df, pr.data_r, by = c("EM" = "gene_symbol") )
    df <- left_join(df, pr.data_r, by = c("Target" = "gene_symbol") )

    df <- df[ , 6:13]
    colnames(df) <- c("Ligand_pr_log2FC", "Ligand_pr_aFC",
                      "Receptor_pr_log2FC", "Receptor_pr_aFC",
                      "EM_pr_log2FC", "EM_pr_aFC",
                      "Target_pr_log2FC", "Target_pr_aFC")

    object@pr_FC = df
  }

  # option 2: ps data
  if(is.null(ps.data_condition1) | is.null(ps.data_condition2)){
    message("pS fold change data not found for at least one condition, this part has been skipped.")
  }else{
    # ps fold change data
    ps.data_s1 = ps.data_condition1[, c("gene_symbol", paste0(object@sender,"_ps"))]
    ps.data_s2 = ps.data_condition2[, c("gene_symbol", paste0(object@sender,"_ps"))]
    ps.data_r1 = ps.data_condition1[, c("gene_symbol", paste0(object@receiver,"_ps"))]
    ps.data_r2 = ps.data_condition2[, c("gene_symbol", paste0(object@receiver,"_ps"))]
    # merge data of two conditions
    ps.data_s = inner_join(ps.data_s1, ps.data_s2, by = c("gene_symbol"))
    ps.data_r = inner_join(ps.data_r1, ps.data_r2, by = c("gene_symbol"))

    # normaliation between assays via limma
    xs <- as.matrix(ps.data_s[ , 2:3])
    ys <- normalizeBetweenArrays(xs)
    ps.data_s[ , 2:3] = ys
    colnames(ps.data_s) = c("gene_symbol", "condition1", "condition2")

    xr <- as.matrix(ps.data_r[ , 2:3])
    yr <- normalizeBetweenArrays(xr)
    ps.data_r[ , 2:3] = yr
    colnames(ps.data_r) = c("gene_symbol", "condition1", "condition2")

    ps.data_s = Cal_foldchange(ps.data_s, correction = ps.correction, q = ps.q)[ , c("gene_symbol", "log2FC", "aFC")]
    ps.data_r = Cal_foldchange(ps.data_r, correction = ps.correction, q = ps.q)[ , c("gene_symbol", "log2FC", "aFC")]

    df <- object@pathways

    df <- left_join(df, ps.data_s, by = c("Ligand" = "gene_symbol") )
    df <- left_join(df, ps.data_r, by = c("Receptor" = "gene_symbol") )
    df <- left_join(df, ps.data_r, by = c("EM" = "gene_symbol") )
    df <- left_join(df, ps.data_r, by = c("Target" = "gene_symbol") )

    df <- df[ , 6:13]
    colnames(df) <- c("Ligand_ps_log2FC", "Ligand_ps_aFC",
                      "Receptor_ps_log2FC", "Receptor_ps_aFC",
                      "EM_ps_log2FC", "EM_ps_aFC",
                      "Target_ps_log2FC", "Target_ps_aFC")

    object@ps_FC = df
  }

  # option 3: py data
  if(is.null(py.data_condition1) | is.null(py.data_condition2)){
    message("pY fold change data not found for at least one condition, this part has been skipped.")
  }else{
    # py fold change data
    py.data_s1 = py.data_condition1[, c("gene_symbol", paste0(object@sender,"_py"))]
    py.data_s2 = py.data_condition2[, c("gene_symbol", paste0(object@sender,"_py"))]
    py.data_r1 = py.data_condition1[, c("gene_symbol", paste0(object@receiver,"_py"))]
    py.data_r2 = py.data_condition2[, c("gene_symbol", paste0(object@receiver,"_py"))]
    # merge data of two conditions
    py.data_s = inner_join(py.data_s1, py.data_s2, by = c("gene_symbol"))
    py.data_r = inner_join(py.data_r1, py.data_r2, by = c("gene_symbol"))

    # normaliation between assays via limma
    xs <- as.matrix(py.data_s[ , 2:3])
    ys <- normalizeBetweenArrays(xs)
    py.data_s[ , 2:3] = ys
    colnames(py.data_s) = c("gene_symbol", "condition1", "condition2")

    xr <- as.matrix(py.data_r[ , 2:3])
    yr <- normalizeBetweenArrays(xr)
    py.data_r[ , 2:3] = yr
    colnames(py.data_r) = c("gene_symbol", "condition1", "condition2")

    py.data_s = Cal_foldchange(py.data_s, correction = py.correction, q = py.q)[ , c("gene_symbol", "log2FC", "aFC")]
    py.data_r = Cal_foldchange(py.data_r, correction = py.correction, q = py.q)[ , c("gene_symbol", "log2FC", "aFC")]

    df <- object@pathways

    df <- left_join(df, py.data_s, by = c("Ligand" = "gene_symbol") )
    df <- left_join(df, py.data_r, by = c("Receptor" = "gene_symbol") )
    df <- left_join(df, py.data_r, by = c("EM" = "gene_symbol") )
    df <- left_join(df, py.data_r, by = c("Target" = "gene_symbol") )

    df <- df[ , 6:13]
    colnames(df) <- c("Ligand_py_log2FC", "Ligand_py_aFC",
                      "Receptor_py_log2FC", "Receptor_py_aFC",
                      "EM_py_log2FC", "EM_py_aFC",
                      "Target_py_log2FC", "Target_py_aFC")

    object@py_FC = df
  }

  return(object)

}


#' Calculating the signaling probability each L-T pathway using the scRNA-seq data during the permutation test
#'
#' @param object Incytr object
#' @param K parameter for the half occupation level in the hill function, the default value is 0.5 (used when predicting the signaling probability)
#' @param N Hill coefficient, the defaulty value is 2 (used when predicting the signaling probability)
#'
#' @return a Incytr object
#' @export
Cal_SigProb_ptest <- function(object, K = 0.5, N = 2){

  sender = object@sender
  receiver = object@receiver

  # expression data for the L in the sender group (condition 1)
  scdata0 <- object@expr.bygroup$condition1[ , c(sender,"Gene")]
  scdata.l <- data.frame(Ligand = scdata0$Gene, Ligand.value = scdata0[ , 1])

  # expression data for the R, EM, T in the receiver group (condition 1)
  scdata <- object@expr.bygroup$condition1[ , c(receiver,"Gene")]

  scdata.r <- data.frame(Receptor = scdata$Gene, Receptor.value = scdata[ , 1])
  scdata.em <- data.frame(EM = scdata$Gene, EM.value = scdata[ , 1])
  scdata.t <- data.frame(Target = scdata$Gene, Target.value = scdata[ , 1])

  # filter the pathways
  net.filtered <- object@pathways
  net.filtered$Sender.group <- sender
  net.filtered$Receiver.group <- receiver
  net.filtered <- net.filtered[ net.filtered$Ligand %in% scdata.l$Ligand &
                                  net.filtered$Receptor %in% scdata.r$Receptor &
                                  net.filtered$EM %in% scdata.em$EM &
                                  net.filtered$Target %in% scdata.t$Target , ]

  # combine net.filtered and the data
  net.filtered <- left_join(net.filtered, scdata.l, by = join_by(Ligand))
  net.filtered <- left_join(net.filtered, scdata.r, by = join_by(Receptor))
  net.filtered <- left_join(net.filtered, scdata.em, by = join_by(EM))
  net.filtered <- left_join(net.filtered, scdata.t, by = join_by(Target))

  df1 <- net.filtered
  # calculate the communication probability
  layer1 = df1$Ligand.value*df1$Receptor.value
  layer2 = df1$Receptor.value*df1$EM.value
  layer3 = df1$EM.value*df1$Target.value
  P1 <- hill(layer1,K,N)*hill(layer2,K,N)*hill(layer3,K,N)
  df1$SigProb <- P1

  ##############################################################################
  # expression data for the L in the sender group (condition 2)
  if( nrow(object@expr.bygroup$condition2)>0 ){
    scdata0 <- object@expr.bygroup$condition2[ , c(sender,"Gene")]
    scdata.l <- data.frame(Ligand = scdata0$Gene, Ligand.value = scdata0[ , 1])

    # expression data for the R, EM, T in the receiver group (condition 2)
    scdata <- object@expr.bygroup$condition2[ , c(receiver,"Gene")]

    scdata.r <- data.frame(Receptor = scdata$Gene, Receptor.value = scdata[ , 1])
    scdata.em <- data.frame(EM = scdata$Gene, EM.value = scdata[ , 1])
    scdata.t <- data.frame(Target = scdata$Gene, Target.value = scdata[ , 1])

    # filter the pathways
    net.filtered <- object@pathways
    net.filtered$Sender.group <- sender
    net.filtered$Receiver.group <- receiver
    net.filtered <- net.filtered[ net.filtered$Ligand %in% scdata.l$Ligand &
                                    net.filtered$Receptor %in% scdata.r$Receptor &
                                    net.filtered$EM %in% scdata.em$EM &
                                    net.filtered$Target %in% scdata.t$Target , ]

    # combine net.filtered and the data
    net.filtered <- left_join(net.filtered, scdata.l, by = join_by(Ligand))
    net.filtered <- left_join(net.filtered, scdata.r, by = join_by(Receptor))
    net.filtered <- left_join(net.filtered, scdata.em, by = join_by(EM))
    net.filtered <- left_join(net.filtered, scdata.t, by = join_by(Target))

    df2 <- net.filtered
    # calculate the communication probability
    layer1 = df2$Ligand.value*df2$Receptor.value
    layer2 = df2$Receptor.value*df2$EM.value
    layer3 = df2$EM.value*df2$Target.value
    P2 <- hill(layer1,K,N)*hill(layer2,K,N)*hill(layer3,K,N)
    df2$SigProb <- P2
  }else{
    P2 = rep(NA, length(df1$Path))
  }

  ##############################################################################

  # save the SigWeight info to a new data frame
  object@SigProb= data.frame(Path = df1$Path,
                             SigProb_condition1 = P1,
                             SigProb_condition2 = P2)
  names(object@SigProb) = c("Path",
                            paste0("SigProb_", object@conditions[1]),
                            paste0("SigProb_", object@conditions[2]))

  return(object)

}


#' Calculate the p_value based on the scRNA-seq data for pathways in each condition using permutation test
#'an Incytr object
#' @param object Incytr object
#' @param K parameter for the half occupation level in the hill function, the default value is 0.5 (used when predicting the signaling probability)
#' @param N Hill coefficient, the defaulty value is 2 (used when predicting the signaling probability)
#' @param nboot number of permutation test performed, the default value is 10 for a faster running time. However, the reasonable value should be at least 100.
#' @param seed.use the seed use in the random selection
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#' @param FDR_method the method used in the FDR correction (in the function p.adjust), the default selection is "BH"
#' @param cutoff_p_value the cutoff value used to select the pathway that has a lower p-value in at least one condition
#'
#' @return an Incytr object
#' @export
Permutation_test <- function(object,
                             K = 0.5,
                             N = 2,
                             nboot = 10,
                             seed.use = 1L,
                             mean_method = NULL,
                             FDR_method = "BH",
                             cutoff_p_value = NULL){

  nC = ncol(object@data)
  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))

  ptest.SigProb = list()

  for (j in 1:length(object@conditions)) {
    ptest.SigProb[[j]] = matrix(0, nrow(object@SigProb), nboot)
  }

  ptest.object = object

  for (i in 1:nboot) {
    ptest.object@idents = object@idents[permutation[, i]]
    names(ptest.object@idents) = names(object@idents)

    ptest.object <- Expr_bygroup(ptest.object,
                                 mean_method = mean_method)

    # Predict signaling prob. (from transcriptomics data)
    ptest.object <- Cal_SigProb_ptest(ptest.object, K = K, N = N)

    if(length(object@conditions) == 2){
      ptest.SigProb[[1]][ , i] = ptest.object@SigProb[ , (1+1)]
      ptest.SigProb[[2]][ , i] = ptest.object@SigProb[ , (2+1)]
    }else{
      ptest.SigProb[[1]][ , i] = ptest.object@SigProb[ , (1+1)]
    }

  }

  ##############################################################################

  p.adj = list()
  for (j in 1:length(object@conditions)) {
    m <- t(sapply(object@SigProb[ , (j+1)], rep, nboot)) - ptest.SigProb[[j]]
    p = format(round(rowSums(m <= 0)/nboot, 4), nsmall = 4)
    # FDR correction
    p.adj[[j]] = p.adjust(p, method = FDR_method)
  }

  if(length(object@conditions) == 2){
    object@p_value = data.frame(Path = object@pathways$Path, condition1 = p.adj[[1]], condition2 = p.adj[[2]])
  }else{
    object@p_value = data.frame(Path = object@pathways$Path, condition1 = p.adj[[1]], condition2 = NA)
  }

  names(object@p_value) = c("Path",
                            paste0("p_value_", object@conditions[1]),
                            paste0("p_value_", object@conditions[2]))

  # filter the pathways
  if(!is.null(cutoff_p_value) ){
    if(length(object@conditions) == 2){
      df.update = object@p_value[ object@p_value[ , paste0("p_value_", object@conditions[1])] <= cutoff_p_value |
                                    object@p_value[ , paste0("p_value_", object@conditions[2])] <= cutoff_p_value , ]
    }else{
      df.update = object@p_value[ object@p_value[ , paste0("p_value_", object@conditions[1])] <= cutoff_p_value, ]
    }

    if(nrow(df.update)>0){
      Path.update = df.update$Path
      object = object_update(object = object, Path = Path.update)
    }else{
      stop("No signaling pathway inferred using current filter, the filter has been ignored.")
    }
  }

  return(object)

}


#' Infer the 5-step pathways based on the L-T pathway table for fututr exploration
#'
#' @param object Incytr object
#'
#' @return an Incytr object
#' @export
Infer_5step_pathways <- function(object){
  # get the pathway table (4 steps)
  df = object@pathways
  # find the genes which appear as both EM and Target
  multirole <- intersect(df$EM, df$Target)
  # create a list to store the results
  fivesteps.list <- list()

  for (i in 1:length(multirole)) {
    # the select gene is the EM
    df1 <- df[df$EM == multirole[i], c("Ligand", "Receptor", "EM", "Target")]
    colnames(df1) = c("Ligand", "Receptor", "EM_2", "Target")
    # the select gene is the Target
    df2 <- df[df$Target == multirole[i], c("Ligand", "Receptor", "EM", "Target")]
    colnames(df2) = c("Ligand", "Receptor", "EM_1", "EM_2")

    df3 <- inner_join(df1, df2, join_by(Ligand == Ligand, Receptor == Receptor, EM_2 == EM_2), multiple = "all")
    df3 <- df3[!duplicated(df3), ]

    fivesteps.list[[i]] <- df3[ , c("Ligand", "Receptor", "EM_1", "EM_2", "Target")]
  }

  fivesteps.table <- rbindlist(fivesteps.list, use.names = TRUE)

  fivesteps.table = fivesteps.table[!duplicated(fivesteps.table), ]

  object@pathways_5steps = fivesteps.table

  return(object)
}


#' Infer genes with high proteomics fold change
#'
#' @param pr.data_condition1 data frame for proteomics data in condition 1, must have variable "gene_symbol"
#' @param pr.data_condition2 data frame for proteomics data in condition 2, must have variable "gene_symbol"
#' @param cell_group levels for the cell cluster labels
#' @param style select the fold change formula between "log2FC" and "aFC", the default value is "aFC"
#' @param cutoff keep the gene with fold change value passing the cutoff, in default cutoff = 1
#' @param pr.correction the correction used when calculating the fold change to avoid division over zero issue, the default value is NULL
#' @param pr.q the percentile used to determine "low expression level" and is used when calculating the adjlog2FC, the default value is 0.75
#'
#' @return a data frame
#' @export
proteomics_gene <- function(pr.data_condition1,
                            pr.data_condition2,
                            cell_group,
                            style = "aFC",
                            cutoff = 1,
                            pr.correction = 0.00001,
                            pr.q = NULL){

  if(is.null(style) | !style %in% c("log2FC", "aFC")){
    style = "aFC"
  }

  pr.gene = list()

  for (i in 1:length(cell_group)) {
    pr.data_1 = pr.data_condition1[, c("gene_symbol", paste0(cell_group[i], "_pr"))]
    pr.data_2 = pr.data_condition2[, c("gene_symbol", paste0(cell_group[i],"_pr"))]
    # merge data of two conditions
    pr.data = inner_join(pr.data_1, pr.data_2, by = c("gene_symbol"))
    # normalization between assays via limma
    xs <- as.matrix(pr.data[ , 2:3])
    ys <- normalizeBetweenArrays(xs)
    pr.data[ , 2:3] = ys
    colnames(pr.data) = c("gene_symbol", "condition1", "condition2")
    # calculate the fold change value
    pr.data = Cal_foldchange(pr.data, correction = pr.correction, q = pr.q)[ , c("gene_symbol", "log2FC", "aFC")]
    pr.data$cluster = cell_group[i]

    pr.gene[[i]] = pr.data
  }

  pr.gene_output = rbindlist(pr.gene)

  if(style == "aFC"){
    pr.gene_output = pr.gene_output[abs(pr.gene_output$aFC) >= cutoff, ]
  }else{
    pr.gene_output = pr.gene_output[abs(pr.gene_output$log2FC) >= cutoff, ]
  }

  return(pr.gene_output)

}
