library(usethis)
library(biomaRt)

## Read data
BMDM_1 <- read.table(file = 'data-raw/BMDM/YA009L2_counts.txt', header = T)
BMDM_2 <- read.table(file = 'data-raw/BMDM/YA009L3_counts.txt', header = T)

# Extract gene_infos and check that they are the same
gene_infos <- BMDM_1[,1:6]
all(gene_infos == BMDM_2[,1:6])

# Add gene symbol in gene_infos
mm9 <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host = 'http://jul2015.archive.ensembl.org')
mm9_gene_table <- getBM(c('ensembl_gene_id', 'external_gene_name'), mart = mm9)

sum(sapply(gene_infos$Geneid, function(x){substr(x,1,18)}) %in% mm9_gene_table$ensembl_gene_id) / nrow(gene_infos)

gene_infos$Symbol <- mm9_gene_table$external_gene_name[match(sapply(gene_infos$Geneid, function(x){substr(x,1,18)}),mm9_gene_table$ensembl_gene_id)]

# Remove gene info from data keep geneid as row.names
BMDM_1 <- BMDM_1[,-c(1:6)]
BMDM_2 <- BMDM_2[,-c(1:6)]
row.names(BMDM_1) <- gene_infos$Geneid
row.names(BMDM_2) <- gene_infos$Geneid

# Add batch to column names
colnames(BMDM_1) <- gsub('.bam', '', x = colnames(BMDM_1), fixed=T)
colnames(BMDM_1) <- gsub('X06_filtering.b_deRibo.', '', x = colnames(BMDM_1), fixed=T)

colnames(BMDM_2) <- gsub('.bam', '', x = colnames(BMDM_2), fixed=T)
colnames(BMDM_2) <- gsub('X06_filtering.b_deRibo.', '', x = colnames(BMDM_2), fixed=T)

# Merge b1 and b2
full_counts <- cbind(BMDM_1, BMDM_2)

# Create Metadata as user would for App
metadata <- data.frame(
  Organism = 'Mouse',
  Genotype = 'WT',
  CellType = 'BMDM',
  Conditioning = c(rep('Naive',21),rep('LPA', 21)),
  ExperimentNumber = rep(1:6,each=7),
  ExperimentName = rep(c('Naive_Basal', 'Naive_LPA_1h', 'Naive_LPA_3h', 'LPA_Basal', 'LPA_LPA_1h', 'LPA_LPA_3h'),each=7),
  ExperimentReplicate = rep(1,42),
  ActDTime = rep(c(0,30,60,90,120,240,360),6),
  Stimulus = rep(c('None', 'LPA', 'LPA', 'None', 'LPA', 'LPA'),each=7),
  StimulusTime = rep(c(0,60,180,0,60,180),each=7),
  SampleName= colnames(full_counts),
  QC = as.character(rep('Included',42)),
  GEOid = "",
  stringsAsFactors = F)
metadata$QC[metadata$SampleName %in% c("ActD_Naive_L60A120_b1", "ActD_Naive_L180A0_b1", "ActD_Naive_L180A30_b1", "ActD_LPA_L0A360_b1", "ActD_LPA_L60A30_b1" )] <- "QC failed"

# order important columns
metadata$ExperimentName <- factor(metadata$ExperimentName, levels=paste(rep(c("Naive",'LPA'), each=3),rep(c('Basal', 'LPA_1h','LPA_3h'),2), sep='_'))
metadata$Stimulus <- factor(metadata$Stimulus, levels=c('None', 'LPA'))

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
plotDistribution(HL=halflife)

# plot fit for specific gene and experiment(s)
plotGeneFit(gene = "ENSMUSG00000069020.7", norm_counts = norm_counts, metadata = metadata,
            decay_result = decay, gene_infos = gene_infos, gene_type='Geneid', showRegion=T, level=0.95, control=list(facet_grid=list(rows = formula('~Stimulus')), theme = list(aspect.ratio=1)))

plotGeneFit(gene = "Il16", norm_counts = norm_counts, metadata = metadata,
            decay_result = decay, gene_infos = gene_infos, gene_type='Symbol', showRegion=T, level=0.95, control=list(facet_grid=list(rows = formula('~Stimulus')), theme = list(aspect.ratio=1)))

BMDM_dataset <- list(full_counts = full_counts, gene_infos = gene_infos,
                     metadata = metadata, filtered = filtered,
                     spikeins = spikeins, counts = counts,
                     norm_spikeins = norm_spikeins, norm_counts = norm_counts,
                     decay = decay, halflife = halflife)

usethis::use_data(BMDM_dataset,overwrite=T)
