###############################################################
# Script to import seed sybsystem annotation data
# Data: Hiseq - metranscriptomics
# Mona Parizadeh - 2020-2021                             
###############################################################

# import packages
library(scales); packageVersion("scales") #‘1.1.1’
library(reshape2); packageVersion("reshape2") #‘1.4.4’
library(vegan); packageVersion("vegan") #‘2.5.7’

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="./",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="Subsys_pie_chart.pdf", 
              help="output image name; [default= %default]", metavar="character"),
  make_option(c("-L", "--level"), type="integer", default=1,
              help="level of Subsystems hierarchy for pie chart; [default=%default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print("USAGE: $ Subsystems_pie_charts.R -I working_directory/ -O save.filename -L level (1,2,3,4)")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

cat ("Saving results as ", opt$out, "\n")
save_filename <- opt$out

cat ("Creating pie chart for hierarchy level ", opt$level, "\n")

# import other necessary packages
suppressPackageStartupMessages({
  library("ggplot2")
  library("data.table")
})

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
head(table_trimmed); dim(table_trimmed)

# Flip table ####
flipped_table = data.frame(t(table_trimmed), check.names = FALSE) #checknames: checks if the row/colnames are syntactically valid variable names.
View(head(flipped_table))
View(tail(flipped_table))
rownames(flipped_table); dim(flipped_table)
write.table(flipped_table,"~/Documents/article3/metatranscriptomics/aca_rna_subsystem.tsv",sep = "\t", quote = FALSE)
saveRDS(flipped_table, "~/Documents/article3/metatranscriptomics/aca_rna_subsystem.rds")

# Export gene names ####
gene_names = colnames(flipped_table)
gene_lables = paste0(">gene", seq(colnames(flipped_table)))
gene_table = cbind(gene_lables,gene_names)
head(gene_table);dim(gene_table)
write.table(gene_table,"~/Documents/article3/metatranscriptomics/gene_names.tsv",sep = "\t", quote = FALSE)

identical(gene_names,colnames(flipped_table))
colnames(flipped_table) = paste0(">gene", seq(colnames(flipped_table))) #replace sequence w/ ASV
head(flipped_table)[1:5,1:5];dim(flipped_table)

save.image("~/Documents/article3/metatranscriptomics/a3_s1_import_subsys.Rdata")
