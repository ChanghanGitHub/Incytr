#' Define Logistic function
#'
#' @param x input variable
#' @param k parameter controls the overall slope
#'
#' @return a value
#' @export
logi <- function(x,k){
  return( 2/(1+exp(-k*x))-1 )
}



#' Calculate the TPDS, PPDS, PhPDS, and multimodel score for the L-T pathways
#'
#' @param object Incytr object
#' @param score.weight weight for proteomics, ps, and pr data, the default values are all 0.5.
#' @param k_logi parameter for the Logistic funtion, the defaulte value is 2
#' @param style the fold change value used in the evaluation, should be either "log2FC" or "aFC", the default value is "aFC"
#' @param abs.value the values can be "None", "All", "Ligand", "Receptor", "EM", "Target". The default is "None" means the fold change values used in the evaluation can be positive or negative. If it is "All" then the absolute value of the fold changes will be used, ignoring the direction of the change.
#' @param cutoff_TPDS selecting the L-T pathways whose TPDS is higher than the cutoff
#'
#' @return an Incytr object
#' @export
Pathway_evaluation <- function(object,
                               score.weight = NULL,
                               k_logi = NULL,
                               style = NULL,
                               abs.value = "None",
                               cutoff_TPDS = NULL){
  # check the input 'score.weight'
  if(length(score.weight)!=3){
    score.weight = c(0.5, 0.5, 0.5)
  }

  # parameter for logistic function
  if(is.null(k_logi)){
    k_logi = 2
  }

  # check the input 'style', should be either "log2FC" or "aFC"
  if(is.null(style)){
    style = "aFC"
  }else if( !(style %in% c("log2FC", "aFC")) ){
    style = "aFC"
  }

  # check the input 'abs.value'
  if(is.null(abs.value)){ abs.value = "None" }
  if( abs.value == "None"){
    case = c(0, 0, 0, 0)
  }else if( abs.value == "All" ){
    case = c(1, 1, 1, 1)
  }else if( !all(is.element(abs.value, c("Ligand", "Receptor", "EM", "Target"))) ){
    stop("Should specify the component(s) that need to use the absolute value!")
  }else if("Ligand" %in% abs.value){
    case[1] = 1
  }else if("Receptor" %in% abs.value){
    case[2] = 1
  }else if("EM" %in% abs.value){
    case[3] = 1
  }else if("Target" %in% abs.value){
    case[4] = 1
  }

  # part 1 -- scdata (base, s0)
  if(nrow(object@SigProb)==0){
    stop("No SigProb data available, should perform the analysis first.")
  }else{
    # calculate the score
    df = object@SigProb[, style]
    df[is.na(df)] = 0 # 'NA'
    s0 = logi(df, k = k_logi)
  }

  # part 2 -- pr data (add-on, s1)
  if(nrow(object@pr_FC)==0){
    print("No protein fold change data available, skip this part.")
    s1 = 0
  }else{
    # get the data
    df = object@pr_FC[, c(paste0("Ligand_pr_", style),
                          paste0("Receptor_pr_", style),
                          paste0("EM_pr_", style),
                          paste0("Target_pr_", style))]
    df[is.na(df)] = 0 # 'NA'
    # calculate the score
    s1 = 0
    for (i in 1:4) {
      if(case[i]==0){
        s1 = s1 + logi(df[, i], k = k_logi)/4
      }else{
        s1 = s1 + abs(logi(df[, i], k = k_logi))/4
      }
    }
  }

  # part 3 -- ps data  (add-on, s2)
  if(nrow(object@ps_FC)==0){
    print("No pS fold change datan available, skip this part.")
    s2 = 0
  }else{
    # get the data
    df = object@ps_FC[, c(paste0("Ligand_ps_", style),
                          paste0("Receptor_ps_", style),
                          paste0("EM_ps_", style),
                          paste0("Target_ps_", style))]
    df[is.na(df)] = 0 # 'NA'
    # calculate the score
    s2 = 0
    for (i in 1:4) {
      if(case[i]==0){
        s2 = s2 + logi(df[, i], k = k_logi)/4
      }else{
        s2 = s2 + abs(logi(df[, i], k = k_logi))/4
      }
    }
  }

  # part 4 -- py data  (add-on, s3)
  if(nrow(object@py_FC)==0){
    print("No pY fold change datan available, skip this part.")
    s3 = 0
  }else{
    # get the data
    df = object@py_FC[, c(paste0("Ligand_py_", style),
                          paste0("Receptor_py_", style),
                          paste0("EM_py_", style),
                          paste0("Target_py_", style))]
    df[is.na(df)] = 0 # 'NA'
    # calculate the score
    s3 = 0
    for (i in 1:4) {
      if(case[i]==0){
        s3 = s3 + logi(df[, i], k = k_logi)/4
      }else{
        s3 = s3 + abs(logi(df[, i], k = k_logi))/4
      }
    }
  }

  # calculate the final score
  s4 = s0 + score.weight[1]*s1 + score.weight[2]*s2 + score.weight[3]*s3
  object@Evaluation = data.frame(Path = object@pathways$Path,
                                 TPDS = s0,
                                 PPDS = s1,
                                 PhPDS_ps = s2,
                                 PhPDS_py = s3,
                                 multimodel_score = s4)

  # filter the pathways
  if(!is.null(cutoff_TPDS) ){
    df.update = object@Evaluation[abs(object@Evaluation$TPDS) >= cutoff_TPDS,  ]

    if(nrow(df.update)>0){
      Path.update = df.update$Path
      object = object_update(object = object, Path = Path.update)
    }else{
      stop("No signaling pathway inferred using current filter, the filter has been ignored.")
    }
  }

  return(object)
}



#' Recording how many components (in one pathway) is up/down-regulated
#'
#' @param df data frame of the L-T pathway information, usually is selected from the Incytr object
#' @param cutoff the value defining how much change is up/down-regulated, the default value is 0
#'
#' @return a data frame
#' @export
Generate_indicator <- function(df, cutoff = NULL){
  # get the maximum of the df
  df1 = unlist(df)
  df1 = as.numeric(df1[df1 != "NA"])
  m = max(na.omit(df1))

  # check cutoff input
  if(is.null(cutoff)){ cutoff = 0 }

  # set 'NA' value
  M = max(m+1, cutoff+1)
  df[is.na(df)] = M
  # number of up-regulated components
  up = rowSums(df>cutoff & df<M)
  # number of down-regulated components
  down = rowSums(df<cutoff)

  output =  data.frame(up = up, down = down)

  return(output)
}



#' Merge the tables from the pairwise analysis (i.e. multiple rounds of analysis)
#'
#' @param results.list a list of tables for multiple rounds of analysis
#' @param cutoff_PDS selecting pathways whose PDS is higher than the cutoff value, the default value is NULL
#' @param group same pathway may appear in different Sender-Receiver communications, grouping the pathways by "Sender.group", "Receiver.group", or "Both". The default value is "Both".
#' @param select.top select the top pathways based on the PDS value, in default select.top = NULL
#'
#' @return a data frame
#' @export
Merge_results <- function(results.list, cutoff_PDS = NULL, group = "Both", select.top = NULL){

  if(is.null(group)){
    group = "Both"
  }else if( !group %in% c("Sender.group", "Receiver.group", "Both")){
    stop("The input 'group' must be 'Sender.group', 'Receiver.group', or 'Both'.")
  }

  if(length(results.list)==1){
    print("Only one result table found.")
    return(results.list)
  }else if(length(results.list)==0){
    stop("The input 'results.list' cannot be empty.")
  }

  # check the colnames
  var.names = colnames(results.list[[1]])
  if(!all(c("PDS", "Sender.group", "Receiver.group") %in% var.names)){
    stop("The tables must have columns 'PDS', 'Sender.group', and 'Receiver.group'.")
  }

  for (i in 2:length(results.list)) {
    if( !identical(var.names, colnames(results.list[[i]])) ){
      stop("The colnames of all tables must be the same.")
    }
  }

  # merge all tables
  df <- results.list[[1]]
  for (i in 2:length(results.list)) {
    df <- rbind(df, results.list[[i]])
  }

  # apply cutoff
  if(!is.null(cutoff_PDS)){
    df <- df[abs(df$PDS) >= cutoff_PDS, ]
  }

  # apply select.top & group.by
  if(group == "Both"){
    df %>%
      group_by(ID_2) -> df
  }else if(group == "Sender.group"){
    df %>%
      group_by(Sender.group) -> df
  }else{
    df %>%
      group_by(Receiver.group) -> df
  }

  if(!is.null(select.top)){
    df %>%
      arrange(desc(df$PDS)) %>%
      slice_head(n = select.top) -> df.top

    df %>%
      arrange(desc(df$PDS)) %>%
      slice_tail(n = select.top) -> df.bottom

    df <- rbind(df.top, df.bottom)
    df <- df[!duplicated(df), ]
  }

  df %>%
    ungroup() -> df

  return(df)
}



#' Export the results as a table from the Incytr object
#'
#' @param object Incytr object
#' @param indicator determine if we need the indicators to show the number of changed components in the pathways. In default, indicator = FALSE
#'
#' @return a data frame
#' @export
Export_results <- function(object,
                           indicator = FALSE){

  # obtain the pathway table
  df <- object@pathways
  df$Sender.group = object@sender
  df$Receiver.group = object@receiver

  if( length(object@sc_FC)==2 ){
    # export the log2FoldChange of all genes in the sender group between two conditions (Ligand)
    df_sender <- data.frame(Gene = rownames(object@sc_FC[[1]]), log2FoldChange = object@sc_FC[[1]]$log2FoldChange)
    # export the log2FoldChange of all genes in the receiver group between two conditions (R, EM, T)
    df_receiver <- data.frame(Gene = rownames(object@sc_FC[[2]]), log2FoldChange = object@sc_FC[[2]]$log2FoldChange)

    df <- left_join(df, df_sender, by = c("Ligand" = "Gene"))
    df <- left_join(df, df_receiver, by = c("Receptor" = "Gene"))
    df <- left_join(df, df_receiver, by = c("EM" = "Gene"))
    df <- left_join(df, df_receiver, by = c("Target" = "Gene"))
    colnames(df)[ (ncol(df)-3):ncol(df) ] = c("Ligand_sclog2FC", "Receptor_sclog2FC", "EM_sclog2FC", "Target_sclog2FC")
  }


  # check each part and merge the results into the table 'df'
  if(nrow(object@SigProb)>0){
    if(length(object@conditions)==2){
      df <- left_join(df, object@SigProb, join_by(Path))
    }else{
      df <- left_join(df, object@SigProb[ , 1:2], join_by(Path))
    }
  }

  # add the p-value columns and apply the cutoff values
  if(nrow(object@p_value)>0){
    if(length(object@conditions)==2){
      df <- left_join(df, object@p_value, join_by(Path))
    }else{
      df <- left_join(df, object@p_value[ , 1:2], join_by(Path))
    }
  }else{
    message("No p-value information found.")
  }

  if(length(object@conditions)==1){
    if(nrow(df)==0){
      stop("No pathway inferred under current filters.")
    }

    df$ID_1 <- paste0(df$Path, "_", df$Sender.group, "_", df$Receiver.group)
    df$ID_2 <- paste0(df$Sender.group, "_", df$Receiver.group)

    return(df)
  }

  ##############################################################################

  if(nrow(object@pr_FC)>0){
    df <- cbind(df, object@pr_FC)
  }

  if(nrow(object@ps_FC)>0){
    df <- cbind(df, object@ps_FC)
  }

  if(nrow(object@py_FC)>0){
    df <- cbind(df, object@py_FC)
  }

  # check if we need to add the indicators
  if(is.null(indicator)){
    indicator = FALSE
  }
  # add the indicators: scFC_up/down, prFC_up/down, psFC_up/down, pyFC_up/down
  if(isTRUE(indicator)){
    # sc part
    if( length(object@sc_FC)==2 ){
      ind_sc = Generate_indicator(df[, c("Ligand_sclog2FC", "Receptor_sclog2FC", "EM_sclog2FC", "Target_sclog2FC")], cutoff = 0)
      colnames(ind_sc) = c("sc_up", "sc_down")
      df <- cbind(df, ind_sc)
    }
    # pr part
    if(nrow(object@pr_FC)>0){
      ind_pr = Generate_indicator(object@pr_FC[, c(1,3,5,7)], cutoff = 0)
      colnames(ind_pr) = c("pr_up", "pr_down")
      df <- cbind(df, ind_pr)
    }
    # ps part
    if(nrow(object@ps_FC)>0){
      ind_ps = Generate_indicator(object@ps_FC[, c(1,3,5,7)], cutoff = 0)
      colnames(ind_ps) = c("ps_up", "ps_down")
      df <- cbind(df, ind_ps)
    }
    # py part
    if(nrow(object@py_FC)>0){
      ind_py = Generate_indicator(object@py_FC[, c(1,3,5,7)], cutoff = 0)
      colnames(ind_py) = c("py_up", "py_down")
      df <- cbind(df, ind_py)
    }
  }

  # add kinase analysis results
  if(nrow(object@kl.pathways)>0){
    df = inner_join(df, object@kl.pathways, by = join_by(Path))
  }else{
    print("No kinase analysis result found.")
  }

  # add the PDS data frame ("multimodel_score", "PDS") and calculate the difference between two scores
  df <- cbind(df, object@Evaluation[ , 2:ncol(object@Evaluation)])

  if(nrow(df)==0){
    stop("No pathway inferred under current filters.")
  }

  df$ID_1 <- paste0(df$Path, "_", df$Sender.group, "_", df$Receiver.group)
  df$ID_2 <- paste0(df$Sender.group, "_", df$Receiver.group)

  return(df)
}


#' Based on the multi-model score, add the KPDS (if applicable) and calculate the PDS (the final score).
#'
#' @param object Incytr object
#' @param KPDS.weight weight for the KPDS, the default value is 0.5.
#' @param cutoff_PDS selecting pathways whose PDS is higher than the cutoff value, the default value is NULL
#'
#' @return an Incytr object
#' @export
Cal_PDS <- function(object,
                    KPDS.weight = NULL,
                    cutoff_PDS = NULL){

  # check the input 'KPDS.weight'
  if(is.null(KPDS.weight)){
    KPDS.weight = 0.5
  }

  if( length(object@EI)==0 ){
    object@Evaluation$PDS = object@Evaluation$multimodel_score
  }else{
    # Combine the Kinase score with previous calculated "final score"
    df = data.frame(score03 = object@Evaluation[ , "multimodel_score"],
                    s4_1 = object@kl.pathways[, paste0("SiK_score_", object@conditions[1])],
                    s4_2 = object@kl.pathways[, paste0("SiK_score_", object@conditions[2])])
    df$score04 = df$score03

    df$score04[df$score03>0] = df$score04[df$score03>0] + KPDS.weight*df$s4_1[df$score03>0]
    df$score04[df$score03==0] = df$score04[df$score03==0] + KPDS.weight*(df$s4_1[df$score03==0]-df$s4_2[df$score03==0])
    df$score04[df$score03<0] = df$score04[df$score03<0] - KPDS.weight*df$s4_2[df$score03<0]

    object@Evaluation$PDS = df$score04
  }

  # filter the pathways
  if(!is.null(cutoff_PDS) ){
    df.update = object@Evaluation[ abs(object@Evaluation$PDS) >= cutoff_PDS, ]

    if(nrow(df.update)>0){
      Path.update = df.update$Path
      object = object_update(object = object, Path = Path.update)
    }else{
      stop("No signaling pathway inferred using current filter, the filter has been ignored.")
    }
  }

  return(object)

}






