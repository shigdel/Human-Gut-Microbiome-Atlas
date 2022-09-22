setwd("C:/Users/rsh001/OneDrive - University of Bergen/COST/git_crc/Human-Gut-Microbiome-Atlas/today/Human-Gut-Microbiome-Atlas/data")

"data/HGA"


library(tidyr)
library(mia)
library(miaViz)
library(scater)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(tidyverse)
library(scater)
library(ape)

#Read the data 

samples_df <- read.csv(file ="HGMA.web.metadata.csv", 
                       header = TRUE)
ls (samples_df)
str (samples_df)

samples_df$type <- as.factor(samples_df$type)
samples_df$subtype <- as.factor(samples_df$subtype)
samples_df$Gender <- as.factor(samples_df$Gender)
samples_df$Geography <- as.factor(samples_df$Geography)
samples_df$Sequencer <- as.factor(samples_df$Sequencer)
samples_df$enteroType <- as.factor(samples_df$enteroType)

otu_mat <- read.csv(file ="HGMA.web.MSP.abundance.matrix.csv", 
                    header = TRUE)

tax_mat <- read.table(file='IGC2.1989MSPs.taxo.tsv',sep = '\t', header = TRUE)
ls (tax_mat)

tax_mat <- tax_mat %>% 
  rename(Class = class, Family= family, Genus = genus, Order =order, Phylum =phylum, Species = species)
tree <- read.tree("IGC2.1990MSPs.nwk",text = NULL, tree.names = NULL, skip = 0,comment.char = "", keep.multi = FALSE)


otu_mat <- otu_mat %>% # coordinates
  tibble::column_to_rownames("X")
otu_mat <- as.matrix(otu_mat)
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("msp_name")
samples_df <- samples_df %>%
  tibble::column_to_rownames("sample.ID")

tax_mat <- tax_mat[,c("Phylum", "Class", "Order", "Family", "Genus","Species" )]
se <- TreeSummarizedExperiment(assays = list(counts = otu_mat[,rownames(samples_df)]),
                               colData = samples_df,
                               rowData = tax_mat[rownames(otu_mat),])

# how many taxa and samples the data contains

dim(se)


#TreeSummarizedExperiment object which includes also a rowTree slot

tse <- as(se, "TreeSummarizedExperiment")
rowTree(tse) <- tree

tse


common.nodes <- intersect(rownames(tse), rowTree(tse)$tip.label)
tse <- TreeSummarizedExperiment::subsetByLeaf(x = tse, rowLeaf = common.nodes, updateTree = TRUE)

saveRDS(tse)



colData(tse)

# inspect possible values in metadata variables 
unique(tse$Geography)
unique(tse$type)
unique(tse$Gender)


# inspect possible values for Phylum
unique(rowData(tse)$Phylum)

# show recurrence for each value
rowData(tse)$Phylum %>% table()


#Abundance 
# Add relative abundances
tse <- transformSamples(tse, method = "relabundance")


#Alpha diversity 

# Every index is calculated by default if we don't specify indices.
indices <- c("shannon", "faith")
names <- c("Shannon_diversity", "Faith_diversity")
tse <- estimateDiversity(tse, index = indices, name = names)
# Shows the calculated indices
knitr::kable(head(colData(tse)[names])) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")

# ggplot needs data.frame format as input.
# Here, colData is DataFrame, therefore it needs to be converted to data.frame
shannon_hist <- ggplot(as.data.frame(colData(tse)), 
                       aes(x = Shannon_diversity)) + 
  geom_histogram(bins = 20, fill = "gray", color = "black") +
  labs(x = "Shannon index", y = "Sample frequency")

shannon_hist


cross_plot <- ggplot2::ggplot(as.data.frame(colData(tse)), 
                              aes(x = Shannon_diversity, y = Faith_diversity)) + 
  geom_point() + # Adds points
  geom_smooth(method=lm) + # Adds regression line
  labs(x = "Shannon index", y = "Faith diversity") 

cross_plot

#Visualization
# Creates Shannon boxplot 
shannon_box <- ggplot(as.data.frame(colData(tse)),
                      aes(x = type, 
                          y = Shannon_diversity,
                          fill = cohort)) + 
  geom_boxplot() +
  theme(title = element_text(size = 12)) # makes titles smaller

# Creates Faith boxplot 
faith_box <- ggplot(as.data.frame(colData(tse)), aes(x = type, 
                                                     y = Faith_diversity, 
                                                     fill = cohort)) + 
  geom_boxplot() +
  theme(title = element_text(size = 12)) # makes titles smaller

# Puts them into same picture
gridExtra::grid.arrange(shannon_box, faith_box, nrow = 2)



#Community diversity
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "shannon", 
                              name = "shannon")
head(colData(tse)$shannon)


tse <- mia::estimateFaith(tse,
                          assay_name = "counts")
head(colData(tse)$faith)



plots <- lapply(c("shannon", "faith"),
                plotColData,
                object = tse, colour_by = "Geography")
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")



#Beta diversity 

# Does clr transformation. Pseudocount is added, because data contains zeros. 
tse <- transformCounts(tse, method = "clr", pseudocount = 1)

# Gets clr table
clr_assay <- assays(tse)$clr

# Transposes it to get taxa to columns
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")

# Does principal coordinate analysis
euclidean_pcoa <- ecodist::pco(euclidean_dist)

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])

# Creates the plot
euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "Euclidean PCoA with CLR transformation") +
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_plot



# Adds the variable we later use for coloring to the data frame
euclidean_patient_status_pcoa_df <- cbind(euclidean_pcoa_df,
                                          patient_status = colData(tse)$type)

# Creates a plot
euclidean_patient_status_plot <- ggplot(data = euclidean_patient_status_pcoa_df, 
                                        aes(x=pcoa1, y=pcoa2,
                                            color = patient_status)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "PCoA with Aitchison distances") +
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_patient_status_plot



