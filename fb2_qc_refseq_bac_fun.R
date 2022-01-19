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
setwd("~/Documents/article3/metatranscriptomics_dbCor/")
comm = readRDS("aca_rna_refseq_bac_fun.rds")
#comm = as.data.frame(sapply(comm, as.numeric))
comm[1:5,1:5]
#class(comm$gene1)
meta = import_qiime_sample_data("mapping.csv")
meta$year = as.factor(meta$year)
meta$month = as.factor(meta$month)
dim(meta)
fun = read.delim("function_names_bac.tsv")
rownames(fun) = fun$id
dim(fun)
#make phyloseq object
ps = phyloseq(sample_data(meta), otu_table(comm, taxa_are_rows = TRUE))
ps

#Clean database ####  
#for the 100 top functions, verify if there are misannotations, duplicate names w/ a minor difference
#Order based on function abundance
#relative abundance
ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
#melt
ps.tab = psmelt(ps.ra); dim(ps.tab)
#merge
ps.mrg = merge(ps.tab, fun, by.x = "OTU", by.y = "id", all = TRUE)
head(ps.mrg); dim(ps.mrg)
ps.mrg.ra.ord = ps.mrg %>%
  group_by(function.) %>%
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
comm_new = otu_table(ps); dim(comm_new)
comm_new[1:5,1:5]

View(top[order(top),])

#1.DNA-directed RNA polymerase subunit beta... ####
#find in the top functions
top$function.[grep("DNA-directed RNA polymerase subunit beta",top$function.)]
#find in fun
fun[which(fun$function. == "DNA-directed RNA polymerase subunit beta"),] #fun13677
comm_new["fun13677",];sum(comm_new["fun13677",])
fun[which(fun$function. == "DNA-directed RNA polymerase subunit beta'"),] #fun13680
comm_new["fun13680",];sum(comm_new["fun13680",])
comm_new1 = comm_new
comm_new1[rownames(comm_new1) == "fun13677", ] <- comm_new1[rownames(comm_new1) == "fun13677", ] + 
  comm_new1[rownames(comm_new1) == "fun13680", ]
comm_new1["fun13677",]
#remove duplicate row
comm_new1 = comm_new1[rownames(comm_new1) != "fun13680", ]; dim(comm_new1)
#remove from fun
fun1 = fun
fun1 = fun1[rownames(comm_new1),];dim(fun1)
ps2 = phyloseq(sample_data(meta), otu_table(comm_new1, taxa_are_rows = TRUE))
ps2

#Explore data ####
summary(taxa_sums(ps2))
summary(sample_sums(ps2))

#Check for outliers ####
#% NMDS####
#relative abundance
ps2.ra = transform_sample_counts(ps2, function(otu) otu/sum(otu)) 
#ordinate
nmds = ordinate(ps2, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps2.ra, nmds, color = "neonic", shape = "year") + 
  theme_bw() + geom_point(size=4) + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3) + 
  geom_point(size = 1) + scale_shape_manual(values = c(19, 1))
#% alpha diversity ####
shn = estimate_richness(ps2, split=TRUE, measures="Shannon") 
plot_richness(ps2, "neonic","month", measures = "Shannon") +
  geom_text(aes(label = sampleid), size = 3)

# 1.Filter functions with less than 5 reads (for low number of reads) #### 
summary(taxa_sums(ps2))
ps.rdrare = prune_taxa(taxa_sums(ps2) > 5, ps2)
ps.rdrare = prune_samples(sample_sums(ps.rdrare)>0,ps.rdrare)
ps.rdrare
100-ntaxa(ps.rdrare)/ntaxa(ps2)*100 #37.49
100-(sum(taxa_sums(ps.rdrare))/sum(taxa_sums(ps2))*100) #0.03
# # No sample filtering (at least 1000 reads) ####
summary(sample_sums(ps.rdrare))
# ps.rdrareM = prune_samples(sample_sums(ps.rdrare)>=1000, ps.rdrare)
# ps.rdrareM = prune_taxa(taxa_sums(ps.rdrareM)>0, ps.rdrareM)
# ps.rdrareM
saveRDS(ps.rdrare, "ps_ref_bac_dbCor_clean.rds") #for deseq2
comm.rdrare = otu_table(ps.rdrare)
dim(comm.rdrare);dim(comm);dim(fun)
#accordance
fun.rdrare = fun[rownames(comm.rdrare),];dim(fun.rdrare)
write.table(fun.rdrare,"~/Documents/article3/metatranscriptomics_dbCor/function_names_bac_dbCor_clean.tsv",sep = "\t", quote = FALSE)

#Hypothetical proteins removal ####
hp = c(fun.rdrare$function.[grep("hypothetical protein",fun.rdrare$function.)], fun.rdrare$function.[grep("Hypothetical protein",fun.rdrare$function.)])
length(hp)
#(length(hp)/dim(comm)[1])*100
fun.noHP = fun.rdrare[!fun.rdrare$function. %in% hp,];dim(fun.noHP)
#write.table(fun.noHP,"~/Documents/article3/metatranscriptomics_dbCor/function_names_bac_dbCor_clean_noHP.tsv",sep = "\t", quote = FALSE)
#accordance
comm.noHP = comm.rdrare[which(rownames(comm.rdrare) %in% fun.noHP$id),]; dim(comm.noHP)
#100-(sum(apply(comm.noHP,2,sum))/sum(apply(comm.rdrare,2,sum))*100) #19.20%
#make phyloseq object
ps.noHP = phyloseq(sample_data(meta), otu_table(comm.noHP, taxa_are_rows = FALSE))
ps.noHP
#saveRDS(ps.noHP, "ps_sys_bac_noHP.rds")
100-(sum(taxa_sums(ps.noHP))/sum(taxa_sums(ps.rdrare))*100) #19.20
100-(ntaxa(ps.noHP)/ntaxa(ps.rdrare)*100) #0.1

#Rarefaction ####
#rarefaction curve
source("~/Documents/article3/my_functions.R")
summary(sample_sums(ps.noHP))
rare_curve = calculate_rarefaction_curves(ps.noHP,c('Observed', 'Shannon','Simpson'), 
                                          c(50000,100000,200000,500000,
                                            1000000,2000000,5000000,10000000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
#Fig.S2B####
# Plot
rare = ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,10000000,500000), 
                     labels = function(x) format(x, scientific = FALSE)) + #no exponents
  geom_line() + theme_bw() +theme(legend.position = "none") + 
  #labs(title = "Rarefaction curve") +
  labs(x = "\nNumber of sequences", y = "Observed number of expressed genes\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 20, face = "bold"))
rare
library(cowplot); packageVersion("cowplot") #"‘1.1.1’"
save_plot("Documents/article3/metatranscriptomics_dbCor/graphs/ParizadehM_Fig.S1B_bac_fun_rare.pdf",rare, ncol = 2, nrow = 2)

#cutoff
summary(sample_sums(ps.noHP))
ps.rare = rarefy_even_depth(ps.noHP, sample.size = 1800000, rngseed = 506, trimOTUs = TRUE, replace = TRUE)
ps.rare
rarecurve(otu_table(ps.rare), step=10000,label=FALSE, col = "darkred",
          ylab = "Functions",
          main = "Rarefy to 1,800,000 reads per sample")
saveRDS(ps.rare, "ps_ref_bac_dbCor_clean_noHP_rare.rds")
100-(sum(taxa_sums(ps.rare))/sum(taxa_sums(ps.noHP))*100)
(dim(otu_table(ps.rare))[1]/dim(otu_table(ps.noHP))[1])*100
ntaxa(ps.rare)/ntaxa(ps.noHP)*100 
comm.rare = otu_table(ps.rare)
dim(comm.rare);dim(comm);dim(fun)
fun.rare = fun.noHP[rownames(comm.rare),];dim(fun.rare)
write.table(fun.rare,"~/Documents/article3/metatranscriptomics_dbCor/function_names_bac_dbCor_clean_noHP_rare.tsv",sep = "\t", quote = FALSE)

gene.rich = estimate_richness(ps.rare, measures = "Observed") #gene per sample (richness)
summary(gene.rich)
sd(gene.rich$Observed, na.rm=TRUE) /  
  sqrt(length(gene.rich$Observed[!is.na(gene.rich$Observed)])) #SE 

save.image("~/Documents/article3/metatranscriptomics_dbCor/a3_fb2_qc_refseq_bac_fun.Rdata")
