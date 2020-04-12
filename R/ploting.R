#' Evaluate ggplot layer added based on control list
#' @param p ggplot2 object. Current plot to add a layer to.
#' @param layer character. gpplot2 function to use.
#' @param args named list. Arguments names and their associated value of the
#' \code{layer} function.
#' @return A gggplot2 object.
#' @keywords internal
# eval layer based on layer name and arguments
evalLayer <- function(p, layer, args){
  p <- eval(parse(text = paste( 'p + ggplot2::', layer, '(', paste(names(args), ' = args[["',names(args), '"]]', sep='',collapse=', '),')', sep='')))
  return(p)
}

#' Plot boxplots of spike-ins counts before and after normalisation
#'
#' This function plots spikeins counts before and after normalisation.
#' It assumes that before and after colnames are found in the column
#' SampleName of the metadata
#'
#' @param before numeric matrix or data.frame. Spike-ins counts matrix before
#' normalisation.
#' @param after numeric matrix or data.frame. Spike-ins counts matrix after
#' normalisation.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param control list. Further control ggplot2 output. This must be a named
#' list with the names corresponding to ggplot2 function and the associated
#' value is also a named list of the ggplot2 function argument name and its
#' associated value.
#' @examples
#' attach(MEF_dataset)
#' plotSpikeIns(before = spikeins, after = norm_spikeins, metadata = metadata)
#'
#' # Add some control over the ggplot2 output
#' plotSpikeIns(before = spikeins, after = norm_spikeins, metadata = metadata,
#'              control = list(labs=list(title= "Spike-ins counts before vs after normalisation")))
#' @export
plotSpikeIns <- function(before, after, metadata, control = NULL){
  dat_g <- rbind(data.frame(Counts = unlist(before),
                            Sample = rep(colnames(before), each=nrow(before)),
                            Type='Before'),
                 data.frame(Counts = as.vector(unlist(after)),
                            Sample = rep(colnames(after), each=nrow(after)),
                            Type='After'))
  dat_g <- cbind(dat_g, metadata[match(dat_g$Sample, metadata$SampleName),])

  # order experiment by ActD time
  dat_g <- dat_g[order(dat_g$ActDTime),]
  dat_g$Sample <- factor(as.character(dat_g$Sample), levels = unique(as.character(dat_g$Sample)))
  dat_g$ActDTime <- factor(as.character(dat_g$ActDTime), levels = as.character(sort(unique(dat_g$ActDTime))))

  # create plot
  p <- ggplot2::ggplot(dat_g, ggplot2::aes(x=Sample, y=Counts, fill = ExperimentName))
  p <- p + ggplot2::theme_bw()

  if('geom_boxplot' %in% names(control)){
    p <- evalLayer(p, 'geom_boxplot', args)
  }else{
    p <- p + ggplot2::geom_boxplot()
  }

  # transparent black to create a gradient
  p <- p + ggplot2::geom_boxplot(mapping = ggplot2::aes(alpha=ActDTime), fill = 'black', color=NA)
  suppressWarnings(p <- p + ggplot2::scale_alpha_discrete(range=c(0.6,0.05)))

  # adjust plotting based on control argument
  if('facet_grid' %in% names(control)){
    args <- control$facet_grid
    if( !('rows' %in% names(args))){
      args <- c(args, list(rows= '~ ExperimentName*ExperimentReplicate*Type'))
    }
    if( !('scales' %in% names(args))){
      args <- c(args, list(scales = 'free_x'))
    }
    if( !('drop' %in% names(args))){
      args <- c(args, list(drop= T))
    }
    p <- evalLayer(p, 'facet_grid', args)
  }else{
    p <- p + ggplot2::facet_grid( rows = . ~ ExperimentName*ExperimentReplicate*Type, scales = 'free_x', drop=T)
  }

  if('scale_y_continuous' %in% names(control)){
    args <- control$scale_y_continuous
    if(!('trans' %in% names(args))){
      args <- c(args, list(trans = 'log10'))
    }
    p <- evalLayer(p, 'scale_y_continuous', args)
  }else{
    p <- p + ggplot2::scale_y_continuous(trans='log10')
  }

  if('theme' %in% names(control)){
    args <- control$theme
    if(! ('axis.text.x' %in% names(args))){
      args <- c(args, list(axis.text.x = ggplot2::element_text(angle=90)))
    }
    p <- evalLayer(p, 'theme', args)
  }else{
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))
  }

  if(!is.null(control)){
    for(i in 1:length(control)){
      if(!(names(control)[i] %in% c("scale_y_continuous", "theme", "geom_boxplot", 'facet_grid'))){
        p <- evalLayer(p, names(control)[i], control[[i]])
      }
    }
  }

  print(p)
}


#' Plot bargraph of library size before and after normalisation
#'
#' This function plots observed library size before normalisation and effective
#' normalised library size.
#' It assumes that before and after colnames are found in the column
#' SampleName of the metadata
#'
#' @param before numeric matrix or data.frame. Gene counts matrix before
#' normalisation.
#' @param after numeric matrix or data.frame. Gene counts matrix after
#' normalisation.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param control list. Further control ggplot2 output. This must be a named
#' list with the names corresponding to ggplot2 function and the associated
#' value is also a named list of the ggplot2 function argument name and its
#' associated value.
#' @examples
#' attach(MEF_dataset)
#' plotLibrarySize(before = counts, after = norm_counts, metadata = metadata)
#'
#' #' # Add some control over the ggplot2 output
#' plotLibrarySize(before = counts, after = norm_counts, metadata = metadata,
#'                 control = list(labs=list(y="Effective Library size")))
#' @export
plotLibrarySize <- function(before, after, metadata, control = NULL){
  dat_g <- rbind(data.frame(Library = colSums(before),
                            Sample = colnames(before),
                            Type='Before'),
                 data.frame(Library = colSums(after),
                            Sample = colnames(after),
                            Type='After'))
  dat_g <- cbind(dat_g, metadata[match(dat_g$Sample, metadata$SampleName),])

  # order experiment by ActD time
  dat_g <- dat_g[order(dat_g$ActDTime),]
  dat_g$Sample <- factor(as.character(dat_g$Sample), levels = unique(as.character(dat_g$Sample)))
  dat_g$ActDTime <- factor(as.character(dat_g$ActDTime), levels = unique(as.character(sort(dat_g$ActDTime))))

  # creat plot
  p <- ggplot2::ggplot(dat_g, ggplot2::aes(x=Sample, y= Library, fill= ExperimentName)) + ggplot2::theme_bw()

  if('geom_bar' %in% names(control)){
    args <- c(stat = 'identity', control$geom_bar)
    p <- evalLayer(p, 'geom_bar', args)
  }else{
    p <- p + ggplot2::geom_bar(stat = 'identity')
  }
 p <- p + ggplot2::geom_bar(mapping = ggplot2::aes(alpha=ActDTime), stat = 'identity', fill = 'black', color=NA)
 suppressWarnings( p <- p  + ggplot2::scale_alpha_discrete(range=c(0.6,0.05)))

  if('facet_grid' %in% names(control)){
    args <- control$facet_grid
    if( !('rows' %in% names(args))){
      args <- c(args, list(rows= '~ ExperimentName*ExperimentReplicate*Type'))
    }
    if( !('scales' %in% names(args))){
      args <- c(args, list(scales = 'free_x'))
    }
    if( !('drop' %in% names(args))){
      args <- c(args, list(drop= T))
    }
    p <- evalLayer(p, 'facet_grid', args)
  }else{
    p <- p + ggplot2::facet_grid( rows = . ~ ExperimentName*ExperimentReplicate*Type, scales = 'free_x', drop=T)
  }

  if('theme' %in% names(control)){
    args <- control$theme
    if(! ('axis.text.x' %in% names(args))){
      args <- c(args, list(axis.text.x = ggplot2::element_text(angle=90)))
    }
    p <- evalLayer(p, 'theme', args)
  }else{
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))
  }

  if(!is.null(control)){
    for(i in 1:length(control)){
      if(!(names(control)[i] %in% c("scale_y_continuous", "theme", "geom_boxplot", 'facet_grid'))){
        p <- evalLayer(p, names(control)[i], control[[i]])
      }
    }
  }

  print(p)
}


#' Plot PCA gene counts before and after normalisation
#'
#' This function plots PCA of gene counts before and after normalisation.
#' It assumes that before and after colnames are found in the column
#' SampleName of the metadata
#'
#' @param before numeric matrix or data.frame. Ggene counts matrix before
#' normalisation.
#' @param after numeric matrix or data.frame. Gene counts matrix after
#' normalisation.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param dim1 integer. Principal component to plot in x axis e.g. 1 will plot
#' the first principal component as x axis.
#' @param dim2 integer. Principal component to plot in y axis e.g. 2 will plot
#' the second principal component as y axis.
#' @param control list. Further control ggplot2 output. This must be a named
#' list with the names corresponding to ggplot2 function and the associated
#' value is also a named list of the ggplot2 function argument name and its
#' associated value.
#' @examples
#' attach(MEF_dataset)
#' plotPCA(before = counts, after = norm_counts, metadata = metadata)
#' @export
plotPCA <- function(before, after, metadata, dim1 = 1, dim2 = 2, control = NULL){
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package \"gridExtra\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package \"ggpubr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  pca_before <- stats::prcomp(t(before), center = T, scale. = T)
  pca_after <- stats::prcomp(t(after), center = T, scale. = T)

  dat_ggplot <- rbind(data.frame(pca_before$x, Type='Before', metadata[match(row.names(pca_before$x),metadata$SampleName),]),
                      data.frame(pca_after$x, Type='After', metadata[match(row.names(pca_after$x),metadata$SampleName),]))

  # calculate % of variance explained
  var_exp_before <- pca_before$sdev^2 / sum(pca_before$sdev^2) * 100
  var_exp_after <- pca_after$sdev^2 / sum(pca_after$sdev^2) * 100

  # create plot for before and after
  for(type in c('Before', 'After')){
    p <- ggplot2::ggplot(dat_ggplot[dat_ggplot$Type == type,],
                         ggplot2::aes(x = get(paste('PC',dim1,sep='')),
                             y = get(paste('PC', dim2,sep='')),
                             col = ExperimentName,
                             size = ActDTime,
                             shape = QC)
                        )
    p <- p + ggplot2::theme_bw()
    var <- switch(type, Before = var_exp_before, After = var_exp_after)
    p <- p + ggplot2::labs(x = paste('PC', dim1, ' (',round(var[dim1],2),'%)', sep = ''),
                           y = paste('PC', dim2, ' (',round(var[dim2],2),'%)', sep = ''))

    if('geom_point' %in% names(control)){
      args <- control$geom_point
      if( !('alpha' %in% names(args))){
        args <- c(args, list(alpha= '0.5'))
      }
      p <- evalLayer(p, 'geom_point', args)
    }else{
      p <- p + ggplot2::geom_point(alpha=0.5)
    }

    if('scale_shape_manual' %in% names(control)){
      args <- control$scale_shape_manual
      if(!('values' %in% names(args))){
        args <- c(args, list(values = c('Included'=16,'QC failed'=4)))
      }
      p <- evalLayer(p, 'scale_shape_manual', args)
    }else{
      p <- p + ggplot2::scale_shape_manual(values = c('Included'=16,'QC failed'=4))
    }

    if('facet_grid' %in% names(control)){
      args <- control$facet_grid
      if( !('cols' %in% names(args))){
        args <- c(args, list(cols = ggplot2::vars(Type)))
      }
      p <- evalLayer(p, 'facet_grid', args)
    }else{
      p <- p + ggplot2::facet_grid(cols=ggplot2::vars(Type))
    }

    if('theme' %in% names(control)){
      args <- control$theme
      if(! ('aspect.ratio' %in% names(args))){
        args <- c(args, list(aspect.ratio = 1))
      }
      p <- evalLayer(p, 'theme', args)
    }else{
      p <- p + ggplot2::theme(aspect.ratio = 1)
    }

    if(!is.null(control)){
      for(i in 1:length(control)){
        if(!(names(control)[i] %in% c("scale_y_continuous", "theme", "geom_boxplot", 'facet_grid'))){
          p <- evalLayer(p, names(control)[i], control[[i]])
        }
      }
    }
    assign(x = paste('p_', type, sep=''), value = p)
  }

  # extract legend
  leg <- ggpubr::get_legend(p_Before)

  p_Before <- p_Before + ggplot2::theme(legend.position = "none")
  p_After <- p_After + ggplot2::theme(legend.position = "none")

  p <- gridExtra::arrangeGrob(grobs = list(p_Before, p_After, leg), nrow=1, widths=c(4,4,2))
  print(ggpubr::as_ggplot(p))
}

#' Plot distribution of half-life across experiments
#'
#' This function plots the overall distribution of half-lives infered.
#'
#' @param HL list. half-life list derived from \code{\link{convertSlopeToHL}}
#' function.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param control list. Further control ggplot2 output. This must be a named
#' list with the names corresponding to ggplot2 function and the associated
#' value is also a named list of the ggplot2 function argument name and its
#' associated value.
#' @examples
#' attach(MEF_dataset)
#' plotDistribution(HL = halflife, metadata = metadata)
#' @export
# plot distribution
plotDistribution <- function(HL, metadata, control = NULL){
  dat_ggplot <- data.frame()
  for(r in 1: length(HL)){
    dat_ggplot <- rbind(dat_ggplot, data.frame(HL = HL[[r]]$HL, metadata[which(metadata$ExperimentName == names(HL)[r])[1],]))
  }

  p <- ggplot2::ggplot(dat_ggplot, ggplot2::aes(x=HL, col=ExperimentName))
  p <- p + ggplot2::theme_bw()

  if('geom_density' %in% names(control)){
    p <- evalLayer(p, 'geom_density', args)
  }else{
    p <- p + ggplot2::geom_density()
  }

  if('scale_x_continuous' %in% names(control)){
    args <- control$scale_x_continuous
    if(!('trans' %in% names(args))){
      args <- c(args, list(trans="log10"))
    }
    p <- evalLayer(p, 'scale_x_continuous', args)
  }else{
    p <- p + ggplot2::scale_x_continuous(trans='log10')
  }

  if(!is.null(control)){
    for(i in 1:length(control)){
      if(!(names(control)[i] %in% c("scale_x_continuous", "geom_density"))){
        p <- evalLayer(p, names(control)[i], control[[i]])
      }
    }
  }

  print(p)
}

#' Plot regression results for specific genes
#'
#' This function plots the fit for specified genes and experiments.
#'
#' @param gene character. Gene id or name. Must be found in \code{gene_infos}
#' column specified by \code{gene_type} argument.
#' @param norm_counts numeric matrix or data.frame. Gene counts matrix after
#' normalisation.
#' @param metadata data.frame. Metadata of the data. Must have specific column
#' names see \code{\link{checkMetadata}} function.
#' @param decay_result list. Decay result as derived from
#' \code{\link{inferDecay}} function.
#' @param gene_infos data.frame. Contain gene information.
#' \code{gene_infos$Geneid} column must correspond to the row names of
#' \code{norm_counts}.
#' @param gene_type character. Column name of \code{gene_infos} in which
#' \code{gene} should be found.
#' @param showRegion logical. Whether to plot or not the confidence interval of
#' the regressions.
#' @param level numeric. Confidence level to show, if \code{showRegion} is
#' \code{TRUE}
#' @param threshold numeric. If all the data from an experiment is below the
#' threshold not fit are going to be shown. Should match the \code{threshold}
#' used in \code{\link{inferDecay}} function.
#' @param experiments character vector. Experiments to plot. Must match names
#' found in the ExperimentName column of \code{metadata}. By default plot all
#' the available experiments.
#' @param control list. Further control ggplot2 output. This must be a named
#' list with the names corresponding to ggplot2 function and the associated
#' value is also a named list of the ggplot2 function argument name and its
#' associated value.
#' @examples
#' attach(MEF_dataset)
#' plotGeneFit(gene = "ENSMUSG00000069020.7", norm_counts = norm_counts,
#'             metadata = metadata, decay_result = decay,
#'             gene_infos = genes_infos)
#' @export
# plot specific gene and experiments
plotGeneFit <- function(gene, norm_counts, metadata, decay_result, gene_infos,
                        gene_type='Geneid', showRegion=T, level=0.95, threshold = 32,
                        experiments = NULL, control=NULL){

  threshold <- log2(threshold)

  # get gene_id and index
  if(gene_type == "Geneid"){
    gene_id <- gene
  }else{
    if(!(gene_type %in% colnames(gene_infos))){
      stop('gene_type must correspond to a column in gene_infos')
    }
    if(!(gene %in% gene_infos[,gene_type])){
      stop(paste0(gene, ' was not found in \"gene_infos$',gene_type,'\"'))
    }
    gene_id <- as.character(gene_infos$Geneid[match(gene,gene_infos[,gene_type])])
  }
  if(!(gene_id %in% row.names(norm_counts))){
    stop(paste0(gene, ' was not found in the data'))
  }
  idx <- which(row.names(norm_counts) == gene_id)

  # create data frame for plot
  dat_ggplot <- cbind(data.frame(log2_norm_counts = log2(norm_counts[idx,]),
                                 Shape = "Excluded", stringsAsFactors = F),
                      metadata[match(colnames(norm_counts), metadata$SampleName),]
                      )
  # calculate RealTime for time series experiments
  dat_ggplot$Time <- dat_ggplot$ActDTime + dat_ggplot$StimulusTime

  # adapt shape based on included, excluded (from regression) or qc failed (from metadata)
  for(exp in unique(dat_ggplot$ExperimentNumber)){
    # redo as in decay_with_replicates_1exp to get indices of included points
    samples <- dat_ggplot$SampleName[dat_ggplot$ExperimentNumber == exp]
    times <- dat_ggplot$ActDTime[dat_ggplot$ExperimentNumber == exp]
    reads <- dat_ggplot$log2_norm_counts[dat_ggplot$ExperimentNumber == exp]
    qc <- dat_ggplot$QC[dat_ggplot$ExperimentNumber == exp]
    replicates <- dat_ggplot$ExperimentReplicate[dat_ggplot$ExperimentNumber == exp]

    reads[qc == 'QC failed' | is.infinite(reads)] <- NA

    exp_to_consider <- sapply(unique(replicates),function(r){any(reads[replicates == r] > threshold, na.rm=T)})
    enough_points <- sapply(unique(replicates), function(r){sum(!is.na(reads[replicates == r])) >=3})
    if(any(exp_to_consider & enough_points)){
      rep_tmp <- replicates[replicates %in% unique(replicates)[exp_to_consider & enough_points]]
      reads <- lapply(unique(rep_tmp),function(r){reads[replicates == r]})
      times <- lapply(unique(rep_tmp),function(r){times[replicates == r]})
      samples <- lapply(unique(rep_tmp),function(r){samples[replicates == r]})
    }else{
      next()
    }
    replicates <- rep_tmp

    # get number of usable replicate
    nrep <- length(unique(replicates))
    rep_order <- unique(replicates)

    # reoder by time to be sure
    t_order <- lapply(times, order)
    reads <- lapply(1:nrep, function(r){reads[[r]][t_order[[r]]]})
    times <- lapply(1:nrep, function(r){times[[r]][t_order[[r]]]})
    samples <- lapply(1:nrep, function(r){samples[[r]][t_order[[r]]]})

    # get indices of usable points
    indices <- lapply(reads, function(x){which(!is.na(x))})

    res <- decay_result[[exp]][row.names(decay_result[[exp]]) == gene_id,]

    cat(paste0('Results for Experiment: ', unique(dat_ggplot$ExperimentName[dat_ggplot$ExperimentNumber == exp]), '\n'))
    print(unlist(convertSlopeToHL(list(decay_result[[exp]][row.names(decay_result[[exp]]) == gene_id,]))))
    cat('\n')

    start <- as.numeric(strsplit(as.character(res$Start), ",")[[1]])
    end <- as.numeric(strsplit(as.character(res$End), ",")[[1]])
    skip <- suppressWarnings(as.numeric(strsplit(as.character(res$Skip), ",")[[1]]))
    samples_included <- unlist(lapply(1:nrep, function(r){samples[[r]][indices[[r]][indices[[r]] >= start[r] & indices[[r]] <= end[r] & !(indices[[r]] %in% skip[r])]]}))

    # change shape for included samples
    dat_ggplot$Shape[dat_ggplot$ExperimentNumber == exp & dat_ggplot$SampleName %in% samples_included ] <-'Included'
  }

  dat_ggplot$Shape[dat_ggplot$QC =='QC failed'] <- "QC failed"

  if(!is.null(experiments)){
    dat_ggplot <- dat_ggplot[dat_ggplot$ExperimentName %in% experiments,]
  }

  suppressWarnings( p <- ggplot2::ggplot(dat_ggplot,ggplot2::aes(x = Time,y=log2_norm_counts,
                                                        colour = ExperimentName,
                                                        group = interaction(ExperimentName, ExperimentReplicate))))
  p <- p + ggplot2::geom_point(ggplot2::aes(shape = Shape))
  p <- p + ggplot2::theme_bw()

  # Recreate fit
  if(length(unique(dat_ggplot$ExperimentReplicate[dat_ggplot$Shape == 'Included'])) > 1){
    fit <- stats::lm(formula = log2_norm_counts ~ ExperimentName:as.character(ExperimentReplicate) + Time * ExperimentName, data = dat_ggplot[dat_ggplot$Shape == 'Included',])
  }else{
    fit <- stats::lm(formula = log2_norm_counts ~ ExperimentName + Time * ExperimentName, data = dat_ggplot[dat_ggplot$Shape == 'Included',])
  }
  fit.df <- ggplot2::fortify(fit)
  fit.df <- cbind(fit.df, dat_ggplot[match(row.names(fit.df),dat_ggplot$SampleName), !(colnames(dat_ggplot) %in% c("log2_norm_counts", "Time", "ExperimentName"))])

  # plot fit
  p <- p + ggplot2::geom_line(data = fit.df, ggplot2::aes(y = .fitted))
  p <- p + ggplot2::coord_cartesian(ylim = ggplot2::layer_scales(p)$y$range$range)

  # calculate error need to do each fit separately
  if(showRegion){
    datnew <- data.frame()
    for(exp in unique(dat_ggplot$ExperimentName)){
      dat_exp <- dat_fit <- data.frame()
      for (r in unique(dat_ggplot$ExperimentReplicate[dat_ggplot$ExperimentName == exp])){
        dat <- dat_ggplot[dat_ggplot$ExperimentName == exp & dat_ggplot$ExperimentReplicate == r & dat_ggplot$Shape == "Included",]
        if(nrow(dat) > 0){
          dat_fit <- suppressWarnings(rbind(dat_fit,dat))
          dat_exp <- suppressWarnings(rbind(dat_exp, cbind(data.frame(Time = seq(min(dat$Time), max(dat$Time), by=1)), dat[1, !(colnames(dat) %in% c("Time", "ActDTime","log2_norm_counts", "SampleName"))])))
        }
      }
      if(nrow(dat_fit) > 0){
        if(length(unique(dat_fit$ExperimentReplicate[dat_fit$Shape == 'Included']))>1){
          fit <- stats::lm(formula = log2_norm_counts ~ as.character(ExperimentReplicate) + Time , data = dat_fit[dat_fit$Shape == 'Included',])
        }else{
          fit <- stats::lm(formula = log2_norm_counts ~ Time , data = dat_fit[dat_fit$Shape == 'Included',])
        }
        conf_interval <- stats::predict(fit, newdata=dat_exp, interval="confidence", level=level)
        datnew <- try(rbind(datnew, cbind(dat_exp, conf_interval)),silent = T)
      }
    }
    if(nrow(datnew) > 0){
      p <- p + ggplot2::geom_ribbon(data = datnew, mapping = ggplot2::aes(x = Time, ymin = lwr, ymax = upr, y = fit, fill = ExperimentName), alpha = 0.1, colour = NA)
    }
  }

  if('labs' %in% names(control)){
    args <- control$labs
    if(!('y' %in% names(args))){
      args <- c(args, list(y='log2(normalized counts)'))
    }
    if(!('x' %in% names(args))){
      args <- c(args, list(x='Time (min)'))
    }
    if(!('title' %in% names(args))){
      args <- c(args, list(title = gene))
    }
    p <- evalLayer(p, 'labs', args)
  }else{
    p <- p + ggplot2::labs(y='log2(normalized counts)', x='Time (min)', title = gene)
  }

  if('scale_color_discrete' %in% names(control)){
    args <- control$scale_color_discrete
    if(!('drop' %in% names(args))){
      args <- c(args, list(drop = F))
    }
    p <- evalLayer(p, 'scale_color_discrete', args)
  }else{
    p <- p + ggplot2::scale_color_discrete(drop=F)
  }

  if('scale_fill_discrete' %in% names(control)){
    args <- control$scale_fill_discrete
    if(!('drop' %in% names(args))){
      args <- c(args, list(drop = F))
    }
    p <- evalLayer(p, 'scale_fill_discrete', args)
  }else{
    p <- p + ggplot2::scale_fill_discrete(drop=F)
  }

  if('scale_shape_manual' %in% names(control)){
    args <- control$scale_shape_manual
    if(!('values' %in% names(args))){
      args <- c(args, list(values = c('Excluded'= 1,'Included'=16,'QC failed'=4)))
    }
    p <- evalLayer(p, 'scale_shape_manual', args)
  }else{
    p <- p + ggplot2::scale_shape_manual(values = c('Excluded'= 1,'Included'=16,'QC failed'=4))
  }

  if(!is.null(control)){
    for(i in 1:length(control)){
      if(!(names(control[i]) %in% c("labs", "scale_shape_manual"))){
        p <- evalLayer(p, names(control)[i], control[[i]])
      }
    }
  }

  suppressWarnings(print(p))
}
