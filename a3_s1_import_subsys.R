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
table_trimmed = data_table[,-1] #remove the extra column contaning rownames
View(head(table_trimmed))
colnames(table_trimmed)
# Clean sample names
file_names = ""
for (name in colnames(table_trimmed)) {
  file_names = c(file_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
file_names = file_names[-1]
file_names_trimmed = ""
for (name in file_names) {
  file_names_trimmed = c(file_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[2])}
file_names_trimmed = file_names_trimmed[-1]
file_names_trimmed
dim(table_trimmed);length(file_names_trimmed)
colnames(table_trimmed) = c(file_names_trimmed[1:32],"Level3","Level2","Level1") #replace file names and NAs (levels)
which(colnames(table_trimmed) == "") #find if there's any protein with no name
table_trimmed = table_trimmed[-1,] #remove the protein with no hierchy
table_trimmed = table_trimmed[-which(table_trimmed$Level1 == ""),]
head(table_trimmed); dim(table_trimmed)

# Export gene names ####
gene_names = rownames(table_trimmed)
gene_lables = paste0("gene", seq(rownames(table_trimmed)))
gene_table = cbind(gene_lables,table_trimmed$Level1,table_trimmed$Level2,table_trimmed$Level3,gene_names)
colnames(gene_table) = c("gene","level1","level2","level3","level4")
head(gene_table);dim(gene_table)
write.table(gene_table,"~/Documents/article3/metatranscriptomics/gene_names.tsv",sep = "\t", quote = FALSE)

# Flip table ####
flipped_table = data.frame(t(table_trimmed), check.names = FALSE) #checknames: checks if the row/colnames are syntactically valid variable names.
identical(colnames(flipped_table),gene_names)
View(head(flipped_table))
View(tail(flipped_table))
# Rename genes
colnames(flipped_table) = gene_lables
flipped_table[1:5,1:5]
rownames(flipped_table); dim(flipped_table)
write.table(flipped_table,"~/Documents/article3/metatranscriptomics/aca_rna_subsystem.tsv",sep = "\t", quote = FALSE)
saveRDS(flipped_table, "~/Documents/article3/metatranscriptomics/aca_rna_subsystem.rds")

save.image("~/Documents/article3/metatranscriptomics/a3_s1_import_subsys.Rdata")
