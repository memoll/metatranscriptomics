###############################################################
# Script to import seed sybsystem annotation data             #
# Data: Hiseq - metranscriptomics                             # 
# Mona Parizadeh - 2020-2021                                  #
###############################################################

# import packages ####
library(scales); packageVersion("scales") #‘1.1.1’
library(reshape2); packageVersion("reshape2") #‘1.4.4’
library(vegan); packageVersion("vegan") #‘2.5.7’
library("ggplot2"); packageVersion("ggplot2") #‘3.3.3’
library("data.table"); packageVersion("data.table") #‘1.13.6’
library(dplyr)

# Import data ####
path = "~/Documents/subsys"
list.files(path)
# Get file names ####
files = list.files(path, pattern = "*.reduced", full.names = T, recursive = FALSE)
files
# loading the files 
y <- 0
for (x in files) {
  y <- y + 1
  if (y == 1) {
    data_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(data_table) = c("DELETE", x, "Level4", "Level3", "Level2", "Level1")
    data_table <- data_table[,-1]
    rownames(data_table) <- data_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level2", "Level1")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    data_table <- merge(temp_table, data_table, by = "Level4")  
  }
}

#data_table <- data_table[,-ncol(data_table)]
#View(head(data_table))
colnames(data_table)
head(data_table); dim(data_table)
#table[is.na(data_table)] = 0
rownames(data_table) = data_table$Level4
table_trimmed = data_table[,-1] #remove the extra column contaning rownames (level4)
View(head(table_trimmed))
colnames(table_trimmed)
# Clean sample names
file_names = ""
for (name in colnames(table_trimmed)) {
  file_names = c(file_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
file_names = file_names[-1] #remove ""
file_names_trimmed = ""
for (name in file_names) {
  file_names_trimmed = c(file_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[2])}
file_names_trimmed = file_names_trimmed[-1] #remove ""
file_names_trimmed
dim(table_trimmed);length(file_names_trimmed)
colnames(table_trimmed) = c(file_names_trimmed[1:32],"Level3","Level2","Level1") #replace file names and NAs (levels)

#find if there's any protein with no name at level 4
which(rownames(table_trimmed) == "") 
(sum(table_trimmed[which(rownames(table_trimmed) == ""),][1:32])/sum(table_trimmed[1:32]))*100 #0.19%
#replace it with NA
table_trimmed.l4na = table_trimmed 
rownames(table_trimmed.l4na)[1] = "NA"

#find the protein with no hierarchy
which(rownames(table_trimmed.l4na) == "NO HIERARCHY") #2990
(sum(table_trimmed.l4na[which(rownames(table_trimmed.l4na) == "NO HIERARCHY"),][1:32])/sum(table_trimmed.l4na[1:32]))*100 #9.57%
#keep them

#find if there's any protein with no name at level 1
which(table_trimmed$Level1 == "")
(sum(table_trimmed[which(table_trimmed$Level1 == ""),][1:32])/sum(table_trimmed[1:32]))*100 #9.67%
#replace with NA (it contains the NO HIERARCHY level 4, too)
table_trimmed.l1na = table_trimmed.l4na
table_trimmed.l1na$Level1[which(table_trimmed.l1na$Level1 == "")] = "NA"
which(table_trimmed.l1na$Level1 == "NA") #verify
head(table_trimmed.l1na); dim(table_trimmed.l1na)

#find if there's any protein with no name at level 2
length(which(table_trimmed$Level2 == "")) #a lot!
#find if there's any protein with no name at level 3
which(table_trimmed$Level3 == "") 
#replace with NA (it is the same group of NO HIERARCHY level 4
table_trimmed.l3na = table_trimmed.l1na 
table_trimmed.l3na$Level3[which(table_trimmed.l3na$Level3 == "")] = "NA"

# Export gene names ####
gene_names = rownames(table_trimmed.l3na)
gene_lables = paste0("gene", seq(rownames(table_trimmed.l3na)))
gene_table = as.data.frame(cbind(gene_lables,table_trimmed.l3na$Level1,table_trimmed.l3na$Level2,table_trimmed.l3na$Level3,gene_names))
colnames(gene_table) = c("gene","level1","level2","level3","level4")
head(gene_table);dim(gene_table)
write.table(gene_table,"~/Documents/article3/metatranscriptomics_dbCor/gene_names.tsv",sep = "\t", quote = FALSE)

# Flip table ####
flipped_table = data.frame(t(table_trimmed.l3na), check.names = FALSE) #checknames: checks if the row/colnames are syntactically valid variable names.
identical(colnames(flipped_table),gene_names)
View(head(flipped_table))
View(tail(flipped_table))
# Rename genes
colnames(flipped_table) = gene_lables
flipped_table[1:5,1:5]
rownames(flipped_table); dim(flipped_table)
write.table(flipped_table,"~/Documents/article3/metatranscriptomics_dbCor/aca_rna_subsystem.tsv",sep = "\t", quote = FALSE)
saveRDS(flipped_table, "~/Documents/article3/metatranscriptomics_dbCor/aca_rna_subsystem.rds")

save.image("~/Documents/article3/metatranscriptomics_dbCor/a3_s1_import_subsys.Rdata")
