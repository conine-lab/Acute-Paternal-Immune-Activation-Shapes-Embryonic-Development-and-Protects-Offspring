#Code used to process data for
#"Acute Paternal Immune Activation Shapes Embryonic Development and Protects Offspring from Viral Infection"
#load libraries
library(DESeq2)
library(ggplot2)
library(data.table)

#For data processing in DESeq, users should input a csv with the expected counts data
#for each comparison, i.e. control vs PIA in the caput epididymis, there should be its own CSV
#user is also required to have a metadata table containing the sample name in one column and condition in another,  other information to be considered can also be added

#user is asked to provide inputs
#read in expected counts cvs
expected_count_path <- file.choose()
#read in the cvs
expected_counts <- read.csv(expected_count_path)
#read in metadata file
metadata_path <- file.choose()
metadata <- read.csv(metadata_path)
genetotranscript <- expected_counts[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(expected_counts) <- expected_counts[,2]
#now remove column 1 (genes) from data frame
expected_counts <- expected_counts[,-1]
#now remove the additional transcfript row
expected_counts <- expected_counts[,-1]
#convert all entries to integer
expected_counts[] <- sapply(expected_counts, as.integer)
#filter data 
#User needs to specify how many are in each treatment condition, following the metadata sheet
#put in same order as conditions in metadata
#ask for first condition
user_input_str <- readline(prompt="Enter an number of samples in condition 1: ")
# Convert the character input to an integer
cond1 <- as.integer(user_input_str)
#Condition 2
user_input_str <- readline(prompt="Enter an number of samples in condition 2: ")
# Convert the character input to an integer
cond2 <- as.integer(user_input_str)
#use a cut off of 5 tpm in at least half the samples in each group
cond1cut <- cond1/2
cond2cut <- cond2/2
#filter out genes without at least 5 tpm in half the samples in each group 
expected_counts <- expected_counts[which(rowSums(expected_counts[,1:cond1]>=5)>= cond1cut | rowSums(expected_counts[,cond1+1:cond2]>=5)>= cond2cut),]

#create DESeq object
DESeq_dds <- DESeqDataSetFromMatrix(countData = expected_counts, colData = metadata, design = ~condition)
#run DESeq2, 
DESeq_dds <- DESeq(DESeq_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#ask user for the name of the conditions used in metadata
condition1 <- readline(prompt="Enter the treatment conditon from the metadata table: ")
condition2 <- readline(prompt="Enter the control conditon from the metadata table: ")
#do this for each tissue
DESeq_results <- results(DESeq_dds, contrast = c("condition", condition1, condition2))
# Ask the user for the output directory path
out_dir <- readline(prompt = "enter output directory: ")
# Use the input path to set the working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  print("Directory does not exist. Please check the path.")
}
#print out DESeq reults
write.csv(DESeq_results, "DESeqResults.csv", row.names = TRUE)

#reconnect the transcript numbers with the gene name
DESeq_resDF <- as.data.frame(DESeq_results)
DESeq_resDF$geneID <- rownames(DESeq_resDF)
DESeq_resDF <- merge(DESeq_resDF, genetotranscript, by.x = "geneID", by.y = "transcript")
setcolorder(DESeq_resDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
write.csv(DESeq_resDF, "DESeqRes_withGeneName.csv",row.names=TRUE)
