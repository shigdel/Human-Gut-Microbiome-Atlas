

"data/HGA"

library("ggplot2"); packageVersion("ggplot2")
library(tidyr)
library(mia)
library(miaViz)
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
colData(tse)
sort(table(colData(tse)$Age))
colData(tse)$age <- colData(tse)$Age
