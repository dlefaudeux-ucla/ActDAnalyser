#' Check potential error in the metadata file
#'
#' Check metadata file for obvious discrepancy such SampleName column not
#' matching column names of the counts, not having required columns ...
#'
#' @param metadata data.frame. Metadata data frame to check for consistency.
#' The following columns are expected:
#' \tabular{ll}{
#'        Orgamism \tab organism used for the experiment \cr
#'        Genotype: \tab genotype of the organism used for the experiment \cr
#'        CellType: \tab cell type used for the experiment \cr
#'        Conditionning: \tab conditioning used for the experiment \cr
#'        ExperimentNumber: \tab number of a single ActinomycinD timecourse \cr
#'        ExperimentName: \tab name given of a single ActinomycinD timecourse \cr
#'        ExperimentReplicate: \tab replicate number for the experiment \cr
#'        ActDTime: \tab time after ActinomycinD treatment in min \cr
#'        Stimulus: \tab stimulus used for the experiment \cr
#'        StimulusTime: \tab time between addition of the stimulus and ActinomyinD in min\cr
#'        SampleName: \tab name of the sample. Must correspond to colnames of
#'        full_counts \cr
#'        QC: \tab by default "Included", if for some reason the sample needs to be
#'        removed then can be changed to "QC failed" \cr
#'        GEOid: \tab if data are deposited in GEO, then GEOid of the sample \cr
#' }
#' @param counts data.frame. The counts data  related to the \code{medata}.
#' @param type character. Either 'all' or 'required'. Against which columns to
#' check consistency, all of them or only the required ones.
#' @examples
#' attach(MEF_dataset)
#' checkMetadata(metadata, full_counts, type='all')
#' detach(MEF_dataset)
#' @export
#'
checkMetadata <- function(metadata, counts, type = 'required'){
  # Check metadata has required columns
  if(!(type %in% c('req', 'required', 'all'))){
    stop("type must be 'required' or 'all'")
  }
  req_col <- c('ExperimentNumber', 'ExperimentName', 'ExperimentReplicate','ActDTime','SampleName','QC')
  all_col <- c("Organism", "CellType", "Conditioning", "ExperimentNumber", "ExperimentName", "ExperimentReplicate", "ActDTime", "Stimulus", "StimulusTime", "SampleName", "QC", "GEOid")
  if(type %in% c('req', 'required')){
    if(!all((req_col %in% colnames(metadata)))){
      stop(paste("Metadata must have at least the following columns (case sensitive):\n  ", paste0(req_col,collapse = ' '), "\n\n  Preferably it should have all the following columns: \n  ", paste0(all_col, collapse = ' '), '\n\n  The following columns were not found: \n  ', paste0(req_col[!(req_col %in% colnames(metadata))],collapse=' '),sep= ''))
    }
  }else if( type == 'all'){
    if(!all((all_col %in% colnames(metadata)))){
      stop(paste("Metadata must have all the following columns (case sensitive):\n  ", paste0(all_col,collapse = ' '), '\n\n  The following columns were not found: \n  ', paste0(all_col[!(all_col %in% colnames(metadata))],collapse=' '),sep= ''))
    }
  }

  # Check that ExperimentNames match ExperimentNumber in a 1 to 1 manner
  if(any(sapply(unique(metadata$ExperimentNumber),function(r){length(unique(metadata$ExperimentName[metadata$ExperimentNumber == r]))!=1}))){
    stop('ExperimentNumber and ExperimentName are not a 1 to 1 mapping')
  }

  # Check columns names and sample names match
  if(!all(colnames(counts) %in% metadata$SampleName)){
    stop(paste("Not all columns names of your data are found in the metadata table.\n\n  The following columns are not found in the metadata:\n  ", paste0(colnames(counts)[!(colnames(counts) %in% metadata$SampleNames)], collapse=' ')))
  }

  # check that ExperimentNumber, ExperimentReplicate, and ActDTime are numeric
  if(!is.numeric(metadata$ExperimentNumber)){
    stop('ExperimentNumber column must be numeric.')
  }
  if(!is.numeric(metadata$ExperimentReplicate)){
    stop('ExperimentNumber column must be numeric.')
  }
  if(!is.numeric(metadata$ActDTime)){
    stop('ExperimentNumber column must be numeric.')
  }

  cat('The metadata seems suitable.')
}

#' Filter genes
#'
#' Filter out genes for which no experiment pass the given threshold.
#'
#' @param data data.frame. Raw count data to filter.
#' @param metadata data.frame. Metadata corresponding to the \code{data}.
#' @param threshold integer. Threshold below which gene are going to be removed
#' @return A data.frame which is a subset of \code{data} for which genes that have all
#' counts below the threshold were removed
#' @examples
#' attach(MEF_dataset)
#' filtered <- applyCountThreshold(full_counts, metadata)
#' detach(MEF_dataset)
#' @export
#'
applyCountThreshold <- function(data, metadata, threshold = 32){
  res <- data[,colnames(data) %in% metadata$SampleName]
  keep <- apply(res[,(metadata$QC[match(colnames(res),metadata$SampleName)] != 'QC failed')], 1, function(x){any(x > threshold)})
  res <- res[keep,]
  return(res)
}

#' Separate spike-ins
#'
#' Separate spike-ins from gene data
#'
#' @param data data.frame. filtered count data as the result of
#' \code{\link{applyCountThreshold}}.
#' @param spikein character. Pattern matching in the spike-ins names.
#' @return A list of 2 data.frames. The first on is the a subset of \code{data}
#' with only the spike-ins extracted, the second one corresponds to the genes counts
#' without the spike-ins
#' @examples
#' attach(MEF_dataset)
#' tmp <- extractSpikeIn(filtered, spikein = "ERCC")
#' spikeins <- tmp[[1]]
#' counts <- tmp[[2]]
#' detach(MEF_dataset)
#' @export
#'
extractSpikeIn <- function(data, spikein = "ERCC"){
  if(!is.character(spikein)){
    stop("Spike-ins names need to be character.")
  }
  lines <- unlist(sapply(spikein, function(s){grep(s,row.names(data))}))
  if (length(lines) == 0 ){
    stop("No spike-in found.\nAre you sure you entered the correct name for the spike-in?\nAre you sure the gene id / symbols are in the rownames of the data?")
  }
  res <- list(SpikeIn = data[grep(spikein, row.names(data)),], Genes = data[-grep(spikein, row.names(data)),])
  return(res)
}

#' Normalisation based on spike-ins
#'
#' Normalised counts and spike-ins based on spike-ins counts
#'
#' @param spikein_counts data.frame. Extracted spike-ins count data as the first
#' element of the result from \code{\link{extractSpikeIn}}.
#' @param genes_counts data.frame. Extracted gene count data as the second
#' element of the result from \code{\link{extractSpikeIn}}.
#' @param gene_infos data.frame. Information on the genes. It must include at
#' least a Geneid column which match the row names of the \code{gene_counts}.
#' @param metadata data.frame. Metadata corresponding to the \code{gene_counts}.
#' @param method character. Method to get normalisation factors based on the
#' spike-ins. Must be one method recognized by \code{\link[edgeR]{calcNormFactors}}.
#' By default use "RLE" method.
#' @param p numeric. Percentile (between 0 and 1) of the counts that is aligned
#' when method="upperquartile".
#' @return A list of 2 data.frames. The first on is the normalised version of
#' \code{spikein_counts}, the second one corresponds to the normalised version of
#' \code{genes_counts}.
#' @examples
#' attach(MEF_dataset)
#' tmp <- SpikeInsNormalisation(spikeins, counts, gene_infos, metadata)
#' norm_spikeins <- tmp[[1]]
#' norm_counts <- tmp[[2]]
#' detach(MEF_dataset)
#' @export
#'
SpikeInsNormalisation <- function(spikein_counts, genes_counts, gene_infos, metadata, method="RLE", p = 0.75){
  y_spikeins <- edgeR::DGEList(spikein_counts, genes = gene_infos[match(row.names(spikein_counts), gene_infos$Geneid),], group = metadata$ExperimentNumber[match(colnames(spikein_counts),metadata$SampleName)])
  y_data <- edgeR::DGEList(genes_counts, genes = gene_infos[match(row.names(genes_counts), gene_infos$Geneid),], group = metadata$ExperimentNumber[match(colnames(genes_counts),metadata$SampleName)])

  spikeins_norm_factor <- edgeR::calcNormFactors(y_spikeins, method = "RLE", p = p )

  spikeins_norm_counts <- t(t(spikeins_norm_factor$counts) / (spikeins_norm_factor$samples$lib.size * spikeins_norm_factor$samples$norm.factors / stats::median(spikeins_norm_factor$samples$lib.size * spikeins_norm_factor$samples$norm.factors)))
  y_data_norm_counts <- t(t(y_data$counts) / (spikeins_norm_factor$samples$lib.size * spikeins_norm_factor$samples$norm.factors / stats::median(spikeins_norm_factor$samples$lib.size * spikeins_norm_factor$samples$norm.factors)))

  res <- list(SpikeIn_norm = spikeins_norm_counts, Genes_norm = y_data_norm_counts)
  return(res)
}
