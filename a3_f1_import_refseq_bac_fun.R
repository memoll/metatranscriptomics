###############################################################
# Script to import RefSeq Bacteria annotation data            #
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
path = "~/Documents/refseq_fun"
list.files(path)
# Get file names ####
files = list.files(path, pattern = "*.tsv", full.names = T, recursive = FALSE)
files
# loading the files 
y <- 0
for (x in files) {
  y <- y + 1
  if (y == 1) {
    table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(table) = c("DELETE", x, "V3")
    table <- table[,c(3,2)]      }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    table <- merge(table, temp_table, by = "V3", all = T)  }
}
table[is.na(table)] <- 0
rownames(table) = table$V3
View(head(table))
table_trimmed <- data.frame(table[,-1]) #remove V3
# Clean sample names
file_names = ""
for (name in colnames(table_trimmed)) {
  file_names = c(file_names, unlist(strsplit(name, split='_', fixed=TRUE))[3])}
file_names = file_names[-1]
file_names_trimmed = ""
for (name in file_names) {
  file_names_trimmed = c(file_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[2])}
file_names_trimmed = file_names_trimmed[-1]
colnames(table_trimmed) = file_names_trimmed
dim(table_trimmed);length(file_names_trimmed)
table_trimmed = table_trimmed[-1,] #remove unassigned functions
head(table_trimmed); dim(table_trimmed)

# Export function names ####
fun_names = rownames(table_trimmed)
fun_lables = paste0("fun", seq(rownames(table_trimmed)))
fun_table = cbind(fun_lables,fun_names)
colnames(fun_table) = c("id","function")
head(fun_table);dim(fun_table)
write.table(fun_table,"~/Documents/article3/metatranscriptomics/function_names.tsv",sep = "\t", quote = FALSE)
saveRDS(fun_table, "~/Documents/article3/metatranscriptomics/aca_rna_refseq_bac_fun.rds")

save.image("~/Documents/article3/metatranscriptomics/a3_f1_import_refseq_bac_fun.Rdata")


