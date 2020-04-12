#' MEF dataset
#'
#' All object created while processing the MEF dataset
#' @format A named list with 10 objects:
#' \describe{
#'    \item{\bold{full_counts} (data.frame)}{ Original count matrix. Columns
#'    represents samples, rows genes. }
#'    \item{\bold{gene_infos} (data.frame)}{ Information on the genes. Must contain a
#'    Geneid column that correspond to row names of full_counts.}
#'    \item{\bold{metadata} (data.frame)}{ A data frame with 49 rows and 13 variables:
#'      \tabular{ll}{
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
#'     }
#'   }
#'   \item{\bold{filtered} (data.frame)}{ Filtered count matrix. Result of
#'   \code{\link{applyCountThreshold}} }
#'   \item{\bold{spikeins} (data.frame)}{ Subset of \bold{filtered} matrix
#'   corresponding to spike-ins data. First element of the result of
#'   \code{\link{extractSpikeIn}} }
#'   \item{\bold{counts} (data.frame)}{ Subset of \bold{filtered} matrix
#'   where spike-ins data was removed. Second element of the result of
#'   \code{\link{extractSpikeIn}}}
#'   \item{\bold{norm_spikeins} (data.frame)}{Normalised matrix
#'   of spike-ins data. First element of the result of
#'   \code{\link{SpikeInsNormalisation}} }
#'   \item{\bold{norm_counts} (data.frame)}{Normalised matrix of the counts.
#'   Second element of the result of \code{\link{SpikeInsNormalisation}} }
#'   \item{\bold{decay} (list)}{Named list of the fit results corresponding
#'   to each experiment. Result from \code{\link{inferDecay}} with the
#'   following columns for each experiment:
#'   \tabular{ll}{
#'       ReplicatesNumber \tab replicates number from \code{metadata$ExperimentNumber}
#'       used, if multiple replicates were used they are separated by a comma \cr
#'       Start \tab indice of the starting point for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       End \tab indice of the point that was skipped for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       Intercept \tab intercept result for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       Slope \tab slope result of the regression \cr
#'       Adj Rsquare \tab adjusted R^2 of the regression \cr
#'       DF \tab degrees of freedoms for the regression \cr
#'       CI_X_LB \tab lower bound of the X\% confidence interval of the slope \cr
#'       CI_X_UB \tab upper bound of the X\% confidence interval of the slope \cr
#'     }
#'   }
#'   \item{\bold{halflife} (list)}{Named list of the fit results corresponding
#'   to each experiment transformed to half-life. Result from \code{\link{convertSlopeToHL}}}
#' }
#'
"MEF_dataset"

#' BMDM_dataset
#'
#' All object created while processing the BMDM dataset
#' @format A named list with 10 objects:
#' \describe{
#'    \item{\bold{full_counts} (data.frame)}{ Original count matrix. Columns
#'    represents samples, rows genes. }
#'    \item{\bold{gene_infos} (data.frame)}{ Information on the genes. Must contain a
#'    Geneid column that correspond to row names of full_counts.}
#'    \item{\bold{metadata} (data.frame)}{ A data frame with 49 rows and 13 variables:
#'      \tabular{ll}{
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
#'     }
#'   }
#'   \item{\bold{filtered} (data.frame)}{ Filtered count matrix. Result of
#'   \code{\link{applyCountThreshold}} }
#'   \item{\bold{spikeins} (data.frame)}{ Subset of \bold{filtered} matrix
#'   corresponding to spike-ins data. First element of the result of
#'   \code{\link{extractSpikeIn}} }
#'   \item{\bold{counts} (data.frame)}{ Subset of \bold{filtered} matrix
#'   where spike-ins data was removed. Second element of the result of
#'   \code{\link{extractSpikeIn}}}
#'   \item{\bold{norm_spikeins} (data.frame)}{Normalised matrix
#'   of spike-ins data. First element of the result of
#'   \code{\link{SpikeInsNormalisation}} }
#'   \item{\bold{norm_counts} (data.frame)}{Normalised matrix of the counts.
#'   Second element of the result of \code{\link{SpikeInsNormalisation}} }
#'   \item{\bold{decay} (list)}{Named list of the fit results corresponding
#'   to each experiment. Result from \code{\link{inferDecay}} with the
#'   following columns for each experiment:
#'   \tabular{ll}{
#'       ReplicatesNumber \tab replicates number from \code{metadata$ExperimentNumber}
#'       used, if multiple replicates were used they are separated by a comma \cr
#'       Start \tab indice of the starting point for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       End \tab indice of the point that was skipped for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       Intercept \tab intercept result for the regression , if multiple
#'       replicates were used they are separated by a comma \cr
#'       Slope \tab slope result of the regression \cr
#'       Adj Rsquare \tab adjusted R^2 of the regression \cr
#'       DF \tab degrees of freedoms for the regression \cr
#'       CI_X_LB \tab lower bound of the X\% confidence interval of the slope \cr
#'       CI_X_UB \tab upper bound of the X\% confidence interval of the slope \cr
#'     }
#'   }
#'   \item{\bold{halflife} (list)}{Named list of the fit results corresponding
#'   to each experiment transformed to half-life. Result from \code{\link{convertSlopeToHL}}}
#' }
#'
"BMDM_dataset"


