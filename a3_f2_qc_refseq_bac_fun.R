###############################################################
# Cleaning metatranscriptomic data                            #
# Data: Hiseq - RefSeq Bacteria                               # 
# Mona Parizadeh - 2020-2021                                  #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0’

# Import data #### 
setwd("~/Documents/article3/metatranscriptomics/")
comm = readRDS("aca_rna_refseq_bac_fun.rds")
#comm = as.data.frame(sapply(comm, as.numeric))
comm[1:5,1:5]
#class(comm$gene1)
meta = import_qiime_sample_data("mapping.csv")
meta$year = as.factor(meta$year)
meta$month = as.factor(meta$month)
dim(meta)
fun = read.delim("function_names.tsv")
rownames(fun) = fun$id
dim(fun)
#make phyloseq object
ps = phyloseq(sample_data(meta), otu_table(comm, taxa_are_rows = TRUE))
ps
saveRDS(ps, "ps_ref_bac.rds")

#Explore data ####
summary(taxa_sums(ps))
summary(sample_sums(ps))

#Check for outliers ####
#% NMDS####
#relative abundance
ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
#ordinate
nmds = ordinate(ps, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.ra, nmds, color = "neonic", shape = "year") + 
  theme_bw() + geom_point(size=4) + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3) + 
  geom_point(size = 1) + scale_shape_manual(values = c(19, 1))
#% alpha diversity ####
shn = estimate_richness(ps, split=TRUE, measures="Shannon") 
plot_richness(ps, "neonic","month", measures = "Shannon") +
  geom_text(aes(label = sampleid), size = 3)

# 1.Filter functions with less than 10 reads (for low number of reads) #### 
summary(taxa_sums(ps))
ps.rd10 = prune_taxa(taxa_sums(ps) > 10, ps)
ps.rd10 = prune_samples(sample_sums(ps.rd10)>0,ps.rd10)
ps.rd10
# # No sample filtering (at least 1000 reads) ####
# summary(sample_sums(ps))
# ps.rd10M = prune_samples(sample_sums(ps)>=2000000, ps)
# ps.rd10M = prune_taxa(taxa_sums(ps.rd10M)>0, ps.rd10M)
# ps.rd10M

#Rarefaction ####
#rarefaction curve
source("~/Documents/article3/my_functions.R")
rare_curve = calculate_rarefaction_curves(ps,c('Observed', 'Shannon','Simpson'), 
                                          c(2000000,5000000,10000000,16000000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
# Plot
ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,16000000,200000))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  #labs(title = "Rarefaction curve") +
  labs(x = "\nNumber of sequences", y = "Observed number of functions\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) 
#cutoff
summary(sample_sums(ps))
ps.rare = rarefy_even_depth(ps, sample.size = 2280000, rngseed = 506, trimOTUs = TRUE, replace = TRUE)
ps.rare
rarecurve(otu_table(ps.rare), step=10000,label=FALSE, col = "darkred",
          ylab = "Functions",
          main = "Rarefy to 2,280,000 reads per sample")
saveRDS(ps.rare, "ps_ref_bac_rare.rds")

comm.rare = otu_table(ps.rare)
dim(comm.rare) 
dim(comm)
dim(fun)
fun1 = fun[rownames(comm.rare),]
dim(fun1)
write.table(fun1,"~/Documents/article3/metatranscriptomics/function_names_rare.tsv",sep = "\t", quote = FALSE)

save.image("~/Documents/article3/metatranscriptomics/a3_f2_qc_refseq_bac_fun.Rdata")


