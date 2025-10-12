#' Calculate the exclusiveness index (EI) for the genes between clusters based on the scRNA-seq data
#'
#' @param df data frame of the expression level in clusters
#' @param cell_group the group/cluster to calculate the EI
#' @param fold_threshold if the difference between the selected value and the second highest value passes the threshold, then we think it is "highly exclusive", EI = 1. The default setting is fold_threshold = 10.
#'
#' @return a data frame
#' @export
Cal_EI <- function(df,
                   cell_group,
                   fold_threshold = 10){
  # check the input 'cell_group'
  if( length(setdiff(cell_group, colnames(df)))>0 ){
    message("Some cell types are not included in the data.")
    cell_group = intersect(cell_group, colnames(df))
  }
  # check the input 'fold_threshold'
  if( is.null(fold_threshold) ){
    fold_threshold = 10
  }else if( fold_threshold<1 ){
    message("The input 'fold_threshold' must be no less than 1, has been set to 1.")
    fold_threshold = 1
  }

  exclu = df[, cell_group] # use this format to store the EI
  for (i in 1:length(cell_group)) {
    # the highest one for each row
    max.row <- apply(df[, cell_group],1, max)
    # the lowest one for each row
    min.row <- apply(df[, cell_group],1, min)
    # the second highest one for each row
    second <- apply(df[, cell_group],1, sort)[length(cell_group)-1, ]

    select.col <- df[, cell_group[i]]
    # the distance to the lowest value: 1 -- the selected one is the highest, 0 -- is the lowest
    distance1 = select.col - min.row
    distance.max = max.row - min.row
    distance1[distance1 == 0] = 0.00001
    distance.max[distance.max == 0] = 0.00001
    porp.1 <- distance1/distance.max
    # the difference between the select value and the second highest value
    select.col[select.col == 0] = 0.00001
    second[second == 0] = 0.00001
    porp.2 <- select.col/second

    # the index for the selected column
    porp <- vector("numeric", length = length(porp.1))

    for (j in 1:length(porp.1)) {
      if(porp.1[j]==1 & porp.2[j]>fold_threshold & max.row[j] != min.row[j]){
        porp[j] = 1
      }else if(max.row[j] == min.row[j]){
        porp[j] = 0
      }else{
        porp[j] = porp.1[j]*0.99
      }
    }
    # save the EI for the selected column
    exclu[ , cell_group[i]] <- porp
  }


  return(exclu)
}



#' Integrate the kinase library, infer the SiKs and calculate the EI
#'
#' @param object Incytr object
#' @param kldata data frame of the kinase library with variables 'gene', 'site_pos', and 'motif.geneName'
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#' @param cell_group the group/cluster to calculate the EI
#' @param fold_threshold if the difference between the selected value and the second highest value passes the threshold, then we think it is "highly exclusive", EI = 1. The default setting is fold_threshold = 10.
#'
#' @return an Incytr object
#' @export
Integr_kinasedata <- function(object,
                              kldata,
                              mean_method = NULL,
                              cell_group,
                              fold_threshold = 10){

  if(is.null(kldata)){
    stop("The input kinase data cannot be empty.")
  }else if( !is.data.frame(kldata) ){
    stop("The input 'kldata' must be a data frame with variables 'gene', 'site_pos', and 'motif.geneName'.")
  }else if( !all(c('gene', 'site_pos', 'motif.geneName') %in% colnames(kldata)) ){
    stop("The input 'kldata' must be a data frame with variables 'gene', 'site_pos', and 'motif.geneName'.")
  }

  # genes in the pathways
  gene.use = unique(c(object@pathways$Ligand,
                      object@pathways$Receptor,
                      object@pathways$EM,
                      object@pathways$Target))
  ##############################################################################
  # Part 1: integrate the kinase data
  # kinase data with pathway genes: gene -- substrate; motif.gene -- kinase
  kldata = kldata[ kldata$gene %in% gene.use & kldata$motif.geneName %in% rownames(object@data),
                   c('gene', 'site_pos', 'motif.geneName') ]

  object@kl = kldata
  ##############################################################################
  # Part 2: get motif genes for each pathway in each case
  # case1: R is a SiK of EM
  # case2: R is a SiK of T
  # case3: EM is a SiK of T
  # case4: EM is a SiK of R
  # case5: T is a SiK of R
  # case6: T is a SiK of EM
  ncase = 6
  pathdf = object@pathways
  kl.pathways = data.frame(Path = pathdf$Path,
                           case1 = 0,
                           case2 = 0,
                           case3 = 0,
                           case4 = 0,
                           case5 = 0,
                           case6 = 0) # used to store the results, initially they are all zeros
  # specify the kinases and substrates for all 6 cases
  case_key = data.frame(Kinase = c("Receptor", "Receptor", "EM", "EM", "Target", "Target"),
                        Substrate = c("EM", "Target", "Target", "Receptor", "Receptor", "EM"))
  k.gene = c() # used to store all kinase genes
  for (i in 1:ncase) {
    # grab the corresponding col
    x = pathdf[ , c(case_key$Kinase[i], case_key$Substrate[i], "Path") ]
    colnames(x)[1:2] = c('motif.geneName', 'gene')
    y = kldata[ , c('gene', 'motif.geneName') ]
    # find the matched paths
    z = inner_join(x, y,
                   by = c('motif.geneName' = 'motif.geneName', 'gene' = 'gene'),
                   multiple = "all")
    Path.selected = unique(z$Path)
    # record the results in the data frame
    kl.pathways[kl.pathways$Path %in% Path.selected, i+1] = 1
    # add the inferred kinase (motif) genes
    k.gene = c( unique(object@pathways[ object@pathways$Path %in% Path.selected, case_key$Kinase[i] ]), k.gene )
  }

  if(length(k.gene)==0){
    message("No satisfied SiK found!")

    object@kl.pathways = data.frame(Path = object@pathways$Path,
                                    SiK_R_of_EM = as.character(NA),
                                    SiK_R_of_T = as.character(NA),
                                    SiK_EM_of_T = as.character(NA),
                                    SiK_EM_of_R = as.character(NA),
                                    SiK_T_of_R = as.character(NA),
                                    SiK_T_of_EM = as.character(NA))

    SiK_namelist = colnames(object@kl.pathways)[2:(ncase+1)]
    EI.df = data.frame(matrix(ncol = 2*ncase, nrow = length(object@pathways$Path)))
    # name the columns
    for (j in 1:ncase) {
      colnames(EI.df)[2*j-1] = paste0(SiK_namelist[j], "_EI_", object@conditions[1])
      EI.df[ , (2*j-1)] = as.numeric(NA)
      colnames(EI.df)[2*j] = paste0(SiK_namelist[j], "_EI_", object@conditions[2])
      EI.df[ , (2*j)] = as.numeric(NA)
    }
    # merge
    object@kl.pathways = cbind(object@kl.pathways, EI.df)

    object@kl.pathways$s4_1 = 0
    object@kl.pathways$s4_2 = 0
    colnames(object@kl.pathways)[(ncol(object@kl.pathways)-1):ncol(object@kl.pathways)] = c(paste0("SiK_score_", object@conditions[1]), paste0("SiK_score_", object@conditions[2]))
    return(object)
  }

  ##############################################################################
  # Part 3: calculate the motif gene expression in the selected cell_group in two conditions
  # find the cell barcodes of organ1 and organ2
  celluse <- barcodes_bycondition(object)

  # find all the kinase genes that will be used
  geneuse <- unique(k.gene)

  # expression data of condition 1 (row: cell; column: gene; the last column: cell group label)
  M1 = as.matrix(object@data[ geneuse, celluse$condition1])
  if(length(geneuse)==1){
    df1 <- as.data.frame(M1)
    colnames(df1) = geneuse
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M1))
  }
  df1$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition1])
  df1$group_label <- as.character(df1$group_label)

  # expression data of condition 2 (row: cell; column: gene; the last column: cell group label)
  M2 = as.matrix(object@data[ geneuse, celluse$condition2])
  if(length(geneuse)==1){
    df2 <- as.data.frame(M2)
    colnames(df2) = geneuse
  }else if(length(geneuse)>1){
    df2 <- as.data.frame(t(M2))
  }
  df2$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition2])
  df2$group_label <- as.character(df2$group_label)

  q = c(.25, .5, .75)

  # for condition 1
  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_1 = as.data.frame(t(yy)) # row: cells, column: groups
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_1 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  }

  # motifexpr.bygroup_1 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  colnames(motifexpr.bygroup_1) <- q1$group_label

  # for condition 2
  if(is.null(mean_method)){

    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df2)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df2)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_2 = as.data.frame(t(yy)) # row: cells, column: groups
  } else if(mean_method=="mean"){
    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df2)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_2 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  }

  # motifexpr.bygroup_2 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  colnames(motifexpr.bygroup_2) <- q1$group_label

  ##############################################################################
  # Part 4: calculate the motif gene EI in the selected cell_group in two conditions
  df_1 = Cal_EI(df = motifexpr.bygroup_1,
                cell_group = cell_group,
                fold_threshold = fold_threshold)

  df_2 = Cal_EI(df = motifexpr.bygroup_2,
                cell_group = cell_group,
                fold_threshold = fold_threshold)

  df_1$Gene = rownames(df_1)
  df_2$Gene = rownames(df_2)
  object@EI[[1]] = df_1
  object@EI[[2]] = df_2
  names(object@EI) = c("condition1", "condition2")

  ##############################################################################
  # Part 5: get the EI for each kinase genes (which are stored in Part 2)
  kl.info = kl.pathways
  # create two data frames to store the results of two conditions
  kl.pathways_1 = kl.pathways
  kl.pathways_2 = kl.pathways
  rownames(kl.pathways_1) = kl.pathways$Path
  rownames(kl.pathways_2) = kl.pathways$Path
  for (i in 1:ncase) {
    # find the pathways have value 1
    Path.selected = kl.pathways$Path[kl.pathways[, i+1]==1]
    # select the corresponding kinase (motif) genes
    Gene_Path = object@pathways[ object@pathways$Path %in% Path.selected, c(case_key$Kinase[i], "Path")]
    colnames(Gene_Path) = c("Kinase gene", "Path")
    # grab the kinase gene name
    kl.info[ kl.info[, i+1]==1 , i+1] = object@pathways[ kl.info[, i+1]==1, case_key$Kinase[i]]
    kl.info[ kl.info[, i+1]==0 , i+1] = NA
    # condition 1
    df_EI_1 = inner_join(df_1[ , c(object@receiver, "Gene")], Gene_Path,
                         by = c('Gene' = 'Kinase gene'),
                         multiple = "all")
    kl.pathways_1[ df_EI_1$Path , i+1] = df_EI_1[ , object@receiver]
    # condition 2
    df_EI_2 = inner_join(df_2[ , c(object@receiver, "Gene")], Gene_Path,
                         by = c('Gene' = 'Kinase gene'),
                         multiple = "all")
    kl.pathways_2[ df_EI_2$Path , i+1] = df_EI_2[ , object@receiver]
  }
  # calculate the s4 for each condition
  s4_1 = apply(kl.pathways_1[, 2:(ncase+1)], 1, sum)/ncase
  s4_2 = apply(kl.pathways_2[, 2:(ncase+1)], 1, sum)/ncase

  object@kl.pathways = data.frame(Path = names(s4_1),
                                  SiK_R_of_EM = kl.info$case1,
                                  SiK_R_of_T = kl.info$case2,
                                  SiK_EM_of_T = kl.info$case3,
                                  SiK_EM_of_R = kl.info$case4,
                                  SiK_T_of_R = kl.info$case5,
                                  SiK_T_of_EM = kl.info$case6)

  ##############################################################################
  # # Part 6: grab the kinase gene' name's EI
  kl.EI = kl.info
  case.list = c("case1", "case2", "case3", "case4", "case5", "case6")
  SiK_namelist = colnames(object@kl.pathways)[2:(ncase+1)]

  for (j in 1:ncase) {
    if(!all(is.na(kl.EI[ , case.list[j]]))){
      colnames(kl.EI)[j+1] = 'Gene'
      kl.EI = left_join(kl.EI, df_1[ , c("Gene", object@receiver)], by = c('Gene' = 'Gene'))
      kl.EI = left_join(kl.EI, df_2[ , c("Gene", object@receiver)], by = c('Gene' = 'Gene'))
      colnames(kl.EI)[j+1] = case.list[j]
    }else{
      kl.EI$new_1 = NA
      kl.EI$new_2 = NA
    }

    colnames(kl.EI)[ncol(kl.EI)-1] = paste0(SiK_namelist[j], "_EI_", object@conditions[1])
    colnames(kl.EI)[ncol(kl.EI)] = paste0(SiK_namelist[j], "_EI_", object@conditions[2])
  }

  kl.EI <- kl.EI[ , !colnames(kl.EI) %in% case.list]
  object@kl.pathways = left_join(object@kl.pathways, kl.EI, by = join_by(Path))

  # merge the s4 scores
  object@kl.pathways$s4_1 = s4_1
  object@kl.pathways$s4_2 = s4_2
  colnames(object@kl.pathways)[(ncol(object@kl.pathways)-1):ncol(object@kl.pathways)] = c(paste0("SiK_score_", object@conditions[1]), paste0("SiK_score_", object@conditions[2]))

  return(object)
}



#' Identify the SrKs and SiKs, recorded in a data frame in Incytr
#'
#' @param object Incytr object
#' @param kldata data frame of the kinase library with variables 'gene', 'site_pos', and 'motif.geneName'
#' @param mean_method the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".
#' @param cell_group the group/cluster to calculate the EI
#' @param fold_threshold if the difference between the selected value and the second highest value passes the threshold, then we think it is "highly exclusive", EI = 1. The default setting is fold_threshold = 10.
#' @param exp_cutoff keep the SrKs whose expression level is higher than the cutoff. In default, exp_cutoff = 0
#'
#' @return an Incytr object
#' @export
Kinase_exploration <- function(object,
                               kldata,
                               mean_method = NULL,
                               cell_group,
                               fold_threshold = 10,
                               exp_cutoff = 0){

  if(is.null(kldata)){
    stop("The input kinase data cannot be empty.")
  }else if( !is.data.frame(kldata) ){
    stop("The input 'kldata' must be a data frame with variables 'gene', 'site_pos', and 'motif.geneName'.")
  }else if( !all(c('gene', 'site_pos', 'motif.geneName') %in% colnames(kldata)) ){
    stop("The input 'kldata' must be a data frame with variables 'gene', 'site_pos', and 'motif.geneName'.")
  }

  # genes in the pathways (R, EM. and T)
  gene.use = unique(c(object@pathways$Receptor,
                      object@pathways$EM,
                      object@pathways$Target))

  gene.measured = rownames(object@data)
  ##############################################################################
  # Part 1: get kinase genes (can be any measured genes), whose substrate is in gene.use (the R, EM. and T in the pathways)
  # "gene" -- substrate ; "motif" -- kinase
  y = kldata[ (kldata$gene %in% gene.use) & (kldata$motif.geneName %in% gene.measured) , ]

  ##############################################################################
  # Part 2: calculate the kinase gene expression in the selected cell_group in two conditions
  celluse <- barcodes_bycondition(object)

  # find all the kinase genes that will be used
  geneuse <- unique(y$motif.geneName)

  # expression data of condition 1 (row: cell; column: gene; the last column: cell group label)
  M1 = as.matrix(object@data[ geneuse, celluse$condition1])
  if(length(geneuse)==1){
    df1 <- as.data.frame(M1)
    colnames(df1) = geneuse
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M1))
  }
  df1$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition1])
  df1$group_label <- as.character(df1$group_label)

  # expression data of condition 2 (row: cell; column: gene; the last column: cell group label)
  M2 = as.matrix(object@data[ geneuse, celluse$condition2])
  if(length(geneuse)==1){
    df2 <- as.data.frame(M2)
    colnames(df2) = geneuse
  }else if(length(geneuse)>1){
    df2 <- as.data.frame(t(M2))
  }
  df2$group_label <- fct_drop(object@idents[names(object@idents) %in% celluse$condition2])
  df2$group_label <- as.character(df2$group_label)

  q = c(.25, .5, .75)

  # for condition 1
  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_1 = as.data.frame(t(yy)) # row: cells, column: groups
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_1 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  }

  # motifexpr.bygroup_1 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  colnames(motifexpr.bygroup_1) <- q1$group_label

  # for condition 2
  if(is.null(mean_method)){

    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df2)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df2)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_2 = as.data.frame(t(yy)) # row: cells, column: groups
  } else if(mean_method=="mean"){
    q1 <- setDT(df2)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df2)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    motifexpr.bygroup_2 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  }

  # motifexpr.bygroup_2 = as.data.frame(t(yy[, -1])) # row: cells, column: groups
  colnames(motifexpr.bygroup_2) <- q1$group_label

  # remove genes that: do not satisfy exp_cutoff in both conditions
  exp.check_1 = motifexpr.bygroup_1[ , object@receiver]
  names(exp.check_1) = rownames(motifexpr.bygroup_1)
  exp.check_2 = motifexpr.bygroup_2[ , object@receiver]
  names(exp.check_2) = rownames(motifexpr.bygroup_2)

  gene.filtered = unique(c(names(exp.check_1[exp.check_1 > exp_cutoff]) , names(exp.check_2[exp.check_2 > exp_cutoff])))

  if(length(gene.filtered)==0){
    return(object)
  }

  motifexpr.bygroup_1 = motifexpr.bygroup_1[ gene.filtered , ]
  motifexpr.bygroup_2 = motifexpr.bygroup_2[ gene.filtered , ]

  # update y
  y = y[ y$motif.geneName %in% gene.filtered, ]

  ##############################################################################
  # Part 3: calculate the motif gene EI in the selected cell_group in two conditions
  df_1 = Cal_EI(df = motifexpr.bygroup_1,
                cell_group = cell_group,
                fold_threshold = fold_threshold)

  df_2 = Cal_EI(df = motifexpr.bygroup_2,
                cell_group = cell_group,
                fold_threshold = fold_threshold)

  df_1$Gene = rownames(df_1)
  df_2$Gene = rownames(df_2)

  ##############################################################################
  # Part 4: get the EI for each kinase (motif) genes
  df_1 = df_1[ , c(object@receiver, "Gene")]
  df_2 = df_2[ , c(object@receiver, "Gene")]

  y = left_join(y, df_1, by = c('motif.geneName' = 'Gene'), multiple = "all")
  colnames(y)[ncol(y)] = paste0(colnames(y)[ncol(y)], "_", object@conditions[1])

  y = left_join(y, df_2, by = c('motif.geneName' = 'Gene'), multiple = "all")
  colnames(y)[ncol(y)] = paste0(colnames(y)[ncol(y)], "_", object@conditions[2])

  y$Receiver.group = object@receiver
  y$Kinase.type = "SrK"

  SiK.list = y[y$motif.geneName %in% gene.use, ]$motif.geneName
  if(length(SiK.list)>0){
    y[y$motif.geneName %in% SiK.list, ]$Kinase.type = "SiK"
  }

  # add the results
  object@kl.explore = y

  return(object)
}
