###############################################################
# Cleaning metatranscriptomic data                            #
# Data: Hiseq - seed sybsystem                                # 
# Mona Parizadeh - 2020-2021                                  #
###############################################################


# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0’

# Import data #### 
setwd("~/Documents/article3/metatranscriptomics_dbCor/")
comm = readRDS("aca_rna_subsystem.rds")
comm1 = comm[1:32,]; dim(comm1)
comm1 = as.data.frame(sapply(comm1, as.numeric))
rownames(comm1) = rownames(comm[1:32,])
comm1[1:5,1:5]
#class(comm$gene1)
meta = import_qiime_sample_data("mapping.csv")
meta$year = as.factor(meta$year)
meta$month = as.factor(meta$month)
dim(meta)
hier = read.delim("gene_names.tsv",sep="\t", header=TRUE)
rownames(hier) = hier$gene
dim(hier)

ps = phyloseq(sample_data(meta), otu_table(comm1, taxa_are_rows = FALSE))
ps

#Clean database (level 4) ####  
#for the 100 top functions, verify if there are misannotations, duplicate names w/ a minor difference
#Order based on function abundance
#relative abundance
ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
#melt
ps.tab = psmelt(ps.ra); dim(ps.tab)
#merge
ps.mrg = merge(ps.tab, hier, by.x = "OTU", by.y = "gene", all = TRUE)
head(ps.mrg); dim(ps.mrg)
ps.mrg.ra.ord = ps.mrg %>%
  group_by(level4) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(dplyr::desc(Abundance))  %>% 
  mutate(rel_abund = Abundance/sum(Abundance)*100)

#Top 100 ####
ps.mrg.ra.ord[1:100,]
sum(ps.mrg.ra.ord[1:100,3]) #sum(ps.mrg.ra.ord[1:100,]$rel_abund)
top = ps.mrg.ra.ord[1:100,1]

#Correction ####
#make new table with ids, gene names and counts
comm_new = t(otu_table(ps)); dim(comm_new)
comm_new[1:5,1:5]

View(top[order(top),])

#1.DNA-directed RNA polymerase beta... ####
#find in hier
hier$level4[grep("DNA-directed RNA polymerase beta",hier$level4)] #don't include "h" to avoid upper/lower case issue
hier[which(hier$level4 == "DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)"),]$gene #"gene1390"
hier[which(hier$level4 == "DNA-directed RNA polymerase beta' subunit (EC 2.7.7.6)"),]$gene #"gene1391"
#sum and remove from comm
comm_new[rownames(comm_new) == "gene1390", ]
comm_new[rownames(comm_new) == "gene1391", ]
comm_new1 = comm_new
#remove duplicate row
comm_new1[rownames(comm_new1) == "gene1390", ] <- comm_new1[rownames(comm_new1) == "gene1390", ] + comm_new1[rownames(comm_new1) == "gene1391", ]
comm_new1["gene1390",]
#remove duplicate row
comm_new1 = comm_new1[rownames(comm_new1) != "gene1391", ]; dim(comm_new1)
#remove from hier
hier1 = hier
hier1 = hier1[rownames(comm_new1),];dim(hier1)
ps2 = phyloseq(sample_data(meta), otu_table(comm_new1, taxa_are_rows = TRUE))
ps2
saveRDS(ps2, "ps_sys_dbCor.rds") #for deseq2
write.table(hier1,"~/Documents/article3/metatranscriptomics_dbCor/gene_names_dbCor.tsv",sep = "\t", quote = FALSE)

#Hypothetical proteins removal ####
#verify at all levels
hp.l1 = c(hier1$level1[grep("hypothetical protein",hier1$level1)], hier1$level1[grep("Hypothetical protein",hier1$level1)])
c(hier1$level2[grep("hypothetical protein",hier1$level2)], hier1$level2[grep("Hypothetical protein",hier1$level2)])
c(hier1$level3[grep("hypothetical protein",hier1$level3)], hier1$level3[grep("Hypothetical protein",hier1$level3)])
hp.l4 = c(hier1$level4[grep("hypothetical protein",hier1$level4)], hier1$level4[grep("Hypothetical protein",hier1$level4)])
hier.noHP = hier1[!hier1$level1 %in% hp.l1,];dim(hier.noHP)
hier.noHP = hier.noHP[!hier.noHP$level4 %in% hp.l4,];dim(hier.noHP)
#write.table(hier.noHP,"~/Documents/article3/metatranscriptomics_new/gene_names_dbCor_noHP.tsv",sep = "\t", quote = FALSE)
#accordance
comm.noHP = comm_new1[hier.noHP$gene,]; dim(comm.noHP)
100-(sum(apply(comm.noHP,2,sum))/sum(apply(comm_new1,2,sum))*100) #0.18
100-(ntaxa(ps.noHP)/ntaxa(ps2)*100)
#make phyloseq object
ps.noHP = phyloseq(sample_data(meta), otu_table(comm.noHP, taxa_are_rows = FALSE))
ps.noHP
100-(sum(taxa_sums(ps.noHP))/sum(taxa_sums(ps2))*100)
#saveRDS(ps.noHP, "ps_sys_dbCor_noHP.rds")

#Explore data ####
summary(taxa_sums(ps.noHP))
summary(sample_sums(ps.noHP))

#Check for outliers ####
#% NMDS####
#relative abundance
ps.noHP.ra = transform_sample_counts(ps.noHP, function(otu) otu/sum(otu)) 
#ordinate
nmds = ordinate(ps.noHP, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.noHP.ra, nmds, color = "neonic", shape = "year") + 
  theme_bw() + geom_point(size=4) + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3) + 
  geom_point(size = 1) + scale_shape_manual(values = c(19, 1))
#% alpha diversity ####
shn = estimate_richness(ps.noHP, split=TRUE, measures="Shannon") 
plot_richness(ps.noHP, "neonic","month", measures = "Shannon") +
  geom_text(aes(label = sampleid), size = 3)

# # No gene filtering (for low number of reads) #### 
# summary(taxa_sums(ps.noHP))
# ps.rd10 = prune_taxa(taxa_sums(ps.noHP) > 10, ps.noHP)
# ps.rd10 = prune_samples(sample_sums(ps.rd10)>0,ps.rd10) 
# ps.rd10
# # No sample filtering (at least 1000000 reads) ####
# summary(sample_sums(ps.noHP))
# ps.rd10M = prune_samples(sample_sums(ps.noHP)>=1000000, ps.noHP)
# ps.rd10M = prune_taxa(taxa_sums(ps.rd10M)>0, ps.rd10M)
# ps.rd10M

#Rarefaction ####
#rarefaction curve
source("~/Documents/article3/my_functions.R")
rare_curve = calculate_rarefaction_curves(ps.noHP,c('Observed', 'Shannon','Simpson'), 
                                          c(10000,50000,100000,500000,1000000,2000000,5000000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
#Fig.S1A####
# Plot
rare = ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,5000000,200000))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  #labs(title = "Rarefaction curve") +
  labs(x = "\nNumber of sequences", y = "Observed number of expressed genes\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 20, face = "bold"))
rare
library(cowplot); packageVersion("cowplot") #"‘1.1.1’"
save_plot("Documents/article3/metatranscriptomics_dbCor/graphs/ParizadehM_Fig.S1A_subsys_rare.pdf",rare, ncol = 2, nrow = 2)

#cutoff
summary(sample_sums(ps.noHP))
ps.rare = rarefy_even_depth(ps.noHP, sample.size = 1430000, rngseed = 306, trimOTUs = TRUE, replace = TRUE)
ps.rare
rarecurve(otu_table(ps.rare), step=10000,label=FALSE, col = "darkred",
          ylab = "Functions",
          main = "Rarefy to 1,430,000 reads per sample")
saveRDS(ps.rare, "ps_sys_dbCor_noHP_rare.rds")
100-(sum(taxa_sums(ps.rare))/sum(taxa_sums(ps.noHP))*100)
(dim(otu_table(ps.rare))[1]/dim(otu_table(ps.noHP))[1])*100
ntaxa(ps.rare)/ntaxa(ps.noHP)*100 
comm.rare = otu_table(ps.rare)
dim(comm.rare) #no genes have been removed, no need to update hier table
comm.noHP[1:5,1:5]
comm.rare[1:5,1:5]
#verify
identical(hier.noHP$gene,rownames(comm.rare))
hier.rare = hier.noHP[rownames(comm.rare),];dim(hier.rare)
write.table(hier.rare,"~/Documents/article3/metatranscriptomics_dbCor/gene_names_dbCor_noHP_rare.tsv",sep = "\t", quote = FALSE)

summary(taxa_sums(ps.rare))
gene.rich = estimate_richness(ps.rare, measures = "Observed") #gene per sample (richness)
summary(gene.rich)
sd(gene.rich$Observed, na.rm=TRUE) /  
  sqrt(length(gene.rich$Observed[!is.na(gene.rich$Observed)])) #SE 

save.image("~/Documents/article3/metatranscriptomics_dbCor/a3_s2_qc_subsys.Rdata")



