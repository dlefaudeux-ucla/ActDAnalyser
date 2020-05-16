library(usethis)
library(biomaRt)

## Read data
MEF_b1 <- read.table(file = 'data-raw/MEF/MEF_ActD_b1_counts.txt', header = T)
MEF_b2 <- read.table(file = 'data-raw/MEF/MEF_ActD_b2_counts.txt', header = T)

# Extract gene_infos and check that they are the same
gene_infos <- MEF_b1[,1:6]
all(gene_infos == MEF_b2[,1:6])

# Add gene symbol in gene_infos
mm9 <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host = 'http://jul2015.archive.ensembl.org')
mm9_gene_table <- getBM(c('ensembl_gene_id', 'external_gene_name'), mart = mm9)

gene_infos$Symbol <-  mm9_gene_table$external_gene_name[match(sapply(gene_infos$Geneid, function(x){substr(x,1,18)}),mm9_gene_table$ensembl_gene_id)]

# Remove gene info from data keep geneid as row.names
MEF_b1 <- MEF_b1[,-c(1:6)]
MEF_b2 <- MEF_b2[,-c(1:6)]
row.names(MEF_b1) <- gene_infos$Geneid
row.names(MEF_b2) <- gene_infos$Geneid

# Add batch to column names
colnames(MEF_b1) <- gsub('.bam', '', x = colnames(MEF_b1), fixed=T)
colnames(MEF_b1) <- paste(colnames(MEF_b1), ".b1", sep="")
colnames(MEF_b2) <- gsub('.bam', '', x = colnames(MEF_b2), fixed=T)
colnames(MEF_b2) <- paste(colnames(MEF_b2), ".b2", sep="")

# Merge b1 and b2
full_counts <- cbind(MEF_b1, MEF_b2)
full_counts <- full_counts[,-c(26,45,46)]

# Create Metadata as user would for App
metadata <- data.frame(
  Organism = 'Mouse',
  Genotype = 'Ifnar-/-',
  CellType = 'MEF',
  Conditioning = 'Naive',
  ExperimentNumber = c(rep(c(1,2,3,4,5),5),rep(c(1,6,3,7), each=6)),
  ExperimentName = c(rep(c('Basal', 'LPS_30min', 'LPS_3h', 'TNF_30min', 'TNF_3h'),5), rep(c('Basal', 'LPS_1h', 'LPS_3h', 'TNF_1h'), each=6)),
  ExperimentReplicate = c(rep(1,5*5),rep(c(2,1,2,1),each=6)),
  ActDTime = c(rep(c(0.5,0,1,3,6),each=5),rep(c(0,0.5,1,2,4,6,0.5,0,1,2,4,6), 2))*60,
  Stimulus = c(rep(c('None', 'LPS', 'LPS', 'TNF', 'TNF'),5), rep(c('None','LPS','LPS','TNF'), each=6)),
  StimulusTime = c(rep(c(0,0.5,3,0.5,3),5),rep(c(0,1,3,1), each=6)) * 60,
  SampleName= colnames(full_counts),
  QC = as.character(rep('Included',49)),
  GEOid = "",
  stringsAsFactors = F)

metadata$ExperimentName <- factor(metadata$ExperimentName, levels=c('Basal', 'LPS_30min', 'LPS_1h','LPS_3h', 'TNF_30min', 'TNF_1h', 'TNF_3h'))
metadata$Stimulus <- factor(metadata$Stimulus, levels=c('None', 'LPS', 'TNF'))

# Load data
checkMetadata(metadata, full_counts, type = 'all')

# remove genes with all reads below threshold by each experiment
filtered <- applyCountThreshold(data = full_counts, metadata = metadata, threshold=32)

# extract Spike-ins
tmp <- extractSpikeIn(data = filtered, spikein = "ERCC")
spikeins <- tmp$SpikeIn
counts <- tmp$Genes
rm(tmp)

# Normalize gene data based on spike-ins
tmp <- SpikeInsNormalisation(spikeins, counts, gene_infos, metadata)
norm_spikeins <- tmp[[1]]
norm_counts <-  tmp[[2]]
rm(tmp)

# plot spike-ins before/after normalisation
plotSpikeIns(before = spikeins, after = norm_spikeins, metadata = metadata, control = list(facet_grid=list(rows = formula('~ Type * ExperimentName * ExperimentReplicate'))))

# plot library size before/after normalisation
plotLibrarySize(before = counts, after = norm_counts, metadata = metadata, control = list(facet_grid=list(rows = formula('~ Type * ExperimentName * ExperimentReplicate'))))

# plot pca of counts before/after normalisation
library(gridExtra)
library(ggpubr)
plotPCA(before = counts, after = norm_counts, metadata = metadata)

# infer decay (including rep informations)
decay <- inferDecay(norm_counts = norm_counts, metadata = metadata, threshold = 32)

halflife <- convertSlopeToHL(decay)

# plot distributions
plotDistribution(HL=halflife, metadata = metadata)

# plot fit for specific gene and experiment(s)
plotGeneFit(gene = "ENSMUSG00000069020.7", norm_counts = norm_counts, metadata = metadata,
            decay_result = decay, gene_infos = gene_infos, gene_type='Geneid', showRegion=T, level=0.95, control=list(facet_grid=list(rows = formula('~Stimulus')), theme = list(aspect.ratio=1)))

plotGeneFit(gene = "Il16", norm_counts = norm_counts, metadata = metadata,
            decay_result = decay, gene_infos = gene_infos, gene_type='Symbol', showRegion=T, level=0.95, control=list(facet_grid=list(rows = formula('~Stimulus')), theme = list(aspect.ratio=1)))

MEF_dataset <- list(full_counts = full_counts, gene_infos = gene_infos,
                    metadata = metadata, filtered = filtered,
                    spikeins = spikeins, counts = counts,
                    norm_spikeins = norm_spikeins, norm_counts = norm_counts,
                    decay = decay, halflife = halflife)

use_data(MEF_dataset, overwrite = T)
