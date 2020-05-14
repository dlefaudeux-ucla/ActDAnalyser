#' Get fit result for a specific regression
#'
#' Get fit result for a specific regression including replicates information
#'
#' @param counts list of numeric vectors. List of log2 normalised counts for a
#' specific experiment with each element correpsonding to one replicate.
#' @param time list of numeric vectors. List of corresponding time after
#' ActinomycinD treatment for a specific experiment, with each element correpsonding to one
#' replicate.
#' @param indices list of integers vectors. List of indices of useable counts,
#' with each element correpsonding to one replicate.
#' @param start integers vector. Index of the starting point for the regression,
#' with each element correpsonding to one replicate.
#' @param end integers vector. Index of the ending point for the regression,
#' with each element correpsonding to one replicate.
#' @param skip integers vector. Index of skipped point for the regression if any,
#' with each element correpsonding to one replicate.
#' @param rep_order integers vector. Replicate number giving the order of the
#' previous parameters
#' @param conf_levels numeric vector. Value of the confindence intervals wanted.
#' Values must be between 0 and 1. If it is a vector, it will calculate the confidence
#' intervals for all the values given.
#' @return A vector containing fit results with the following information:
#' \tabular{ll}{
#'   ReplicatesNumber \tab replicates number from \code{metadata$ExperimentNumber}
#'   used, if multiple replicates were used they are separated by a comma \cr
#'   Start \tab indice of the starting point for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   End \tab indice of the point that was skipped for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   Intercept \tab intercept result for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   Slope \tab slope result of the regression \cr
#'   Adj Rsquare \tab adjusted R square of the regression \cr
#'   DF \tab degrees of freedoms for the regression \cr
#'   CI_X_LB \tab lower bound of the X\% confidence interval of the slope \cr
#'   CI_X_UB \tab upper bound of the X\% confidence interval of the slope \cr
#' }
#' @keywords internal
#'
regression_function  <- function(counts, time, indices, start, end, skip, rep_order, conf_levels){
  nrep <- length(counts)
  read_used <- as.vector(unlist(sapply(1:nrep, function(r){counts[[r]][indices[[r]][ indices[[r]] >= start[r] & indices[[r]] <= end[r] & !(indices[[r]] %in% skip[r])]]})))
  time_used <- as.vector(unlist(sapply(1:nrep, function(r){time[[r]][indices[[r]][ indices[[r]] >= start[r] & indices[[r]] <= end[r] & !(indices[[r]] %in% skip[r])]]})))

  if(nrep >1){
    batch <- as.vector(unlist(sapply(1:nrep,function(r){rep(paste('r', r, sep=''), length(indices[[r]][ indices[[r]] >= start[r] & indices[[r]] <= end[r] & !(indices[[r]] %in% skip[r])]))})))
    fit_all <- stats::lm(read_used ~ batch + time_used)
  }else{
    fit_all <- stats::lm(read_used ~ time_used)
  }

  confint <- lapply(conf_levels, function(cl){
    stats::confint(fit_all, level = cl)
    }
  )

  regression_summary <- summary(fit_all)
  coefs <- stats::coefficients(fit_all)

  res <- c(paste(rep_order, collapse = ','), paste(start,collapse=','), paste(end,collapse=','), paste(skip,collapse=','),
           paste(coefs[names(coefs) != "time_used"] + c(0,rep(coefs[names(coefs) != "time_used"][1], length(coefs[names(coefs) != "time_used"]) - 1 )), collapse=','), # if mulptiple the first is the true intercept but the others are the difference in intercepts compared to the first intercept
           coefs["time_used"],
           regression_summary$adj.r.squared,
           fit_all$df.residual,
           unlist(lapply(confint, function(l){l["time_used",]})))
  return(res)
}

#' Get fit result for a specific gene and experiment
#'
#' Get fit result for a specific gene including replicates information.
#' All possible fit will be tried and the one with the best adjusted R^2
#' will be selected.
#'
#' @param counts numeric vector. log2 normalised counts for a specific experiment.
#' @param time numeric vector. Corresponding time after ActinomycinD treatment
#' for a specific experiment.
#' @param replicates integer vector. Corresponding replicate number as given by
#' \code{metadata$ExperimentReplicate}.
#' @param qc character vector. Correspondinc QC status as given by
#' \code{metadata$QC}.
#' @param threshold numeric. log2 threshold for which a replicate with all counts
#' below it will not be considered.
#' @param conf_levels numeric vector. Value of the confindence intervals wanted.
#' Values must be between 0 and 1. If it is a vector, it will calculate the confidence
#' intervals for all the values given.
#' @inherit regression_function return
#' @inherit inferDecay details
#' @keywords internal
#'
decay_with_replicates_1exp <- function(counts, time, replicates = rep(1,length(counts)), qc, threshold, conf_levels) {
  # replace counts of qc failed by NA
  counts[qc == 'QC failed' | is.infinite(counts)] <- NA

  # check that we can use replicates, i.e. not all data of a single replicate is below the threshold.
  # If we can remove the replicate and treat a no replicates.
  exp_to_consider <- sapply(unique(replicates),function(r){any(counts[replicates == r] > threshold, na.rm=T)})
  enough_points <- sapply(unique(replicates), function(r){sum(!is.na(counts[replicates == r])) >=3})
  if(any(exp_to_consider & enough_points)){
    rep_tmp <- replicates[replicates %in% unique(replicates)[exp_to_consider & enough_points]]
    counts <- lapply(unique(rep_tmp),function(r){counts[replicates == r]})
    time <- lapply(unique(rep_tmp),function(r){time[replicates == r]})
  }else{
    res <- rep(NA, 14)
    return(res)
  }
  replicates <- rep_tmp

  # get number of usable replicate
  nrep <- length(unique(replicates))
  rep_order <- unique(replicates)

  # reoder by time to be sure
  t_order <- lapply(time, order)
  counts <- lapply(1:nrep, function(r){counts[[r]][t_order[[r]]]})
  time <- lapply(1:nrep, function(r){time[[r]][t_order[[r]]]})

  # get indices of usable points
  indices <- lapply(counts, function(x){which(!is.na(x))})

  # get indices of portential start point
  potential_start <- lapply(1:nrep, function(r){indices[[r]][time[[r]][indices[[r]]] <= 60 & (1:length(indices[[r]]) <= length(indices[[r]]) - 2) ]}) # check time and at least 3 points for regression

  # get max slope/h between max (from potential start point) and any other points (with at least 1 intermediate point) to allow for some tolerance to pick up the start point
  index_max <- lapply(1:nrep, function(r){potential_start[[r]][which.max(counts[[r]][potential_start[[r]]])]})
  time_max <- lapply(1:nrep, function(r){time[[r]][index_max[[r]]]})
  max_slope <- lapply(1:nrep, function(r){abs(max((counts[[r]][index_max[[r]]] - counts[[r]][indices[[r]][indices[[r]] > index_max[[r]]+1]])/(time[[r]][indices[[r]][indices[[r]] > index_max[[r]]+1]] - time_max[[r]])*60, na.rm=T))})

  # calculate threhold to pick potential start point
  threshold_start <- lapply(1:nrep, function(r){counts[[r]][index_max[[r]]]-.25 * max_slope[[r]]})

  # pick potential start point based ont the threshold
  potential_start2 <- lapply(1:nrep, function(r){potential_start[[r]][counts[[r]][potential_start[[r]]] >= threshold_start[[r]]]})

  # list all regression to try
  # start/end/skip matrix with row for each replicate, col for the different regression to try
  start <- end <- skip <- c()

  for (r in 1:nrep){
    tmp_skip <- tmp_end <- tmp_start <- c()
    for (s in potential_start2[[r]]){
      potential_end <- indices[[r]][indices[[r]] > indices[[r]][which(indices[[r]] == s) + 1]] # need at least 3 points
      nb_intermediary_removable_pts <- sapply(potential_end, function(e){ l <- length(indices[[r]][indices[[r]] > s & indices[[r]] < e ]); replace(l, l==1,0)})

      tmp_skip <- c(tmp_skip,unlist(sapply(1:length(potential_end), function(i){if(nb_intermediary_removable_pts[i] > 0){c(NA,indices[[r]][indices[[r]] > s & indices[[r]] < potential_end[i]])}else{NA}})))
      tmp_end <- c(tmp_end,unlist(sapply(1:length(potential_end), function(i){rep(potential_end[i], 1+nb_intermediary_removable_pts[i])})))
      tmp_start <- c(tmp_start,rep(s, length(potential_end) + sum(nb_intermediary_removable_pts)))
    }

    if(is.null(skip)){
      skip <- rbind(skip, tmp_skip)
      end <- rbind(end,tmp_end)
      start <- rbind(start, tmp_start)
    }else{
      skip <- rbind(matrix(rep(skip,length(tmp_skip)),nrow=nrow(skip)), rep(tmp_skip, each = ncol(skip)))
      end <- rbind(matrix(rep(end, length(tmp_end)),nrow=nrow(end)), rep(tmp_end, each = ncol(end)))
      start <- rbind(matrix(rep(start,length(tmp_start)),nrow=nrow(start)), rep(tmp_start, each = ncol(start)))
    }
  }


  all_regression <- as.data.frame(t(sapply(1:ncol(start), function(i){regression_function(counts = counts, time = time, indices = indices, start = start[,i], end = end[,i], skip = skip[,i], rep_order = rep_order, conf_levels = conf_levels)})),stringsAsFactors = F)
  colnames(all_regression) <- c("ReplicateNumbers", "Start", "End", "Skip", "Intercept", "Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))
  all_regression[,c("Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))] <- apply(all_regression[,c("Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))],2,as.numeric)

  if( any(all_regression$Slope <= 0)){
    best <- as.vector(unlist(all_regression[all_regression$Slope <= 0,][which.max(all_regression[all_regression$Slope <= 0 ,"Adj Rsquare"]),]))
  }else{
    best <- rep(NA, 14)
  }
  return(best)
}


#' Convert slope result to HL results
#'
#' Convert slope results obtained with \code{\link{inferDecay}} to half-lives
#' results
#'
#' @param result named list of data.frame. Each element of the list correspond
#' to the slope result of one experiment. This is obtained from \code{\link{inferDecay}}.
#' @return A list containing half-life results for all genes and all experiments.
#' Each element corresponds to the result of a single experiment with the following
#' information:
#' \tabular{ll}{
#'   HL \tab estimated half-life in minute \cr
#'   CI_X_LB \tab lower bound of the X\% confidence interval of the half-life \cr
#'   CI_X_UB \tab upper bound of the X\% confidence interval of the half-life \cr
#'   Adj Rsquare \tab adjusted R square of the regression \cr
#' }
#' @details
#' This function will try all the possible linear regression starting at any
#' point within the first hour of ActinomycinD treatment and for which the
#' counts are within the max counts within the first hour and the max within
#' the first hour -1/4 of the max decay rate per hour. Then it will try all
#' possible linear regressions allowing one intermediary point to be excluded
#' and with at least 3 remaining points. It will select the fit with  the best
#' adjusted R^2 and a negative slope as result.
#' @export
#'
convertSlopeToHL <- function(result){
  res_conv <- list()
  for (r in 1:length(result)){
    res_conv[[r]] <- result[[r]][,c(6,9:14,7)]
    res_conv[[r]][,1:7] <- -1/(result[[r]][,c(6,9:14)])
    res_conv[[r]][,1:7][res_conv[[r]][,1:7] < 0] <- Inf
    colnames(res_conv[[r]])[1] <- "HL"
  }
  names(res_conv) <- names(result)
  return(res_conv)
}

#' Infer slope from normalised data
#'
#' Infer slope and its confidence interval from the ActinomycinD normalised data.
#' Depending on the number of genes and experiment it can take some time.
#' Most use to plot the result of single genes.
#'
#' @param norm_counts numeric matrix or data.frame. Normalised data
#' as given by the second element of the result of \code{\link{SpikeInsNormalisation}}.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param threshold numeric. If an experiment as all normalised counts below
#' this threshold no fit are going to be perfomed.
#' @param conf_levels numeric vector. Value of the confindence intervals wanted.
#' Values must be between 0 and 1. If it is a vector, it will calculate the confidence
#' intervals for all the values given.
#' @return A list which elements corresponds to the result of one experiment. \cr
#' Each element is a data.frame with the following columns:
#' \tabular{ll}{
#'   ReplicatesNumber \tab replicates number from \code{metadata$ExperimentNumber}
#'   used, if multiple replicates were used they are separated by a comma \cr
#'   Start \tab indice of the starting point for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   End \tab indice of the point that was skipped for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   Intercept \tab intercept result for the regression , if multiple
#'   replicates were used they are separated by a comma \cr
#'   Slope \tab slope result of the regression \cr
#'   Adj Rsquare \tab adjusted R^2 of the regression \cr
#'   DF \tab degrees of freedoms for the regression \cr
#'   CI_X_LB \tab lower bound of the X\% confidence interval of the slope \cr
#'   CI_X_UB \tab upper bound of the X\% confidence interval of the slope \cr
#' }
#' @details
#' This function will try all the possible linear regression starting at any
#' point within the first hour of ActinomycinD treatment and for which the
#' counts are within the max counts within the first hour and the max within
#' the first hour -1/4 of the max decay rate per hour. Then it will try all
#' possible linear regressions allowing one intermediary point to be excluded
#' and with at least 3 remaining points. It will select the fit with the best
#' adjusted R^2 as result.
#' @examples
#' attach(MEF_dataset)
#' decay <- inferDecay(norm_counts[1:100,], metadata)
#' @export
# wrapper to infer decay for all experiments
inferDecay <- function(norm_counts, metadata, threshold = 32, conf_levels = c(0.95,0.75,0.5)){
  if( threshold <= 0 ){stop("threshold parameter must be > 0")}
  res <- list()
  exp_order <- unique(metadata$ExperimentNumber)
  exp_name <- sapply(exp_order,function(r){metadata$ExperimentName[metadata$ExperimentNumber == r][1]})
  nexp <- length(exp_order)
  for (r in 1:nexp){
    # print(paste('Exp',r))
    reads_exp <- norm_counts[, colnames(norm_counts) %in% metadata$SampleName[metadata$ExperimentNumber == exp_order[r]]]
    time_exp <- metadata$ActDTime[match(colnames(reads_exp), metadata$SampleName)]
    replicates <- metadata$ExperimentReplicate[match(colnames(reads_exp), metadata$SampleName)]
    qc <- metadata$QC[match(colnames(reads_exp), metadata$SampleName)]
    # res[[r]] <- t(sapply(1:nrow(reads_exp),function(i){if(i %% 100 == 0){print(i)}; decay_with_replicates_1exp(reads = log2(reads_exp[i,]) , time = time_exp , replicates = replicates, qc = qc, threshold = log2(threshold))}))
    res[[r]] <- as.data.frame(t(sapply(1:nrow(reads_exp),function(i){decay_with_replicates_1exp(counts = log2(reads_exp[i,]) , time = time_exp , replicates = replicates, qc = qc, threshold = log2(threshold), conf_levels = conf_levels)})))
    colnames(res[[r]]) <- c("ReplicateNumbers", "Start", "End", "Skip", "Intercept", "Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))
    res[[r]][,c("Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))] <- apply(res[[r]][,c("Slope", "Adj Rsquare", "DF", paste('CI', rep(conf_levels*100,each=2), rep(c('LB','UB'), length(conf_levels)), sep='_'))],2,as.numeric)
    row.names(res[[r]]) <- row.names(reads_exp)
  }
  names(res) <- exp_name
  return(res)
}
