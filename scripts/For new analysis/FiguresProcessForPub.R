#Code used to generate figures for
#"Acute Paternal Immune Activation Shapes Embryonic Development and Protects Offspring from Viral Infection"
#load libraries
library(svglite)
library(ggplot2)
library(data.table)

#code used to make graphs for supplemental figure 2
#this code is used for each graph that is made
#required user input is a sorted the output files with gene names for sig up and down by padj value (smallest to largest)
#the csvs used for the paper analysis were generated using GProfiler
#ask user to choose the csv to read in
UpGoRes_path <- file.choose()
UpGoRes <- read.csv(UpGoRes_path)
#now order the table by the Pvalu and make it so most sig is first (ascending order)
UpGoRes <- UpGoRes[order(UpGoRes$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
UpGoResTop <- head(UpGoRes, n = 10)
UpGoResTop <- UpGoResTop[order(UpGoResTop$negative_log10_of_adjusted_p_value)]
#make a horizontal bar graph with gradient fill, the plum color is for up regulated genes
plot <- ggplot(data = UpGoResTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "thistle1", high = "plum4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,30)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
#ask user for output directory
# Ask the user for the output directory path
out_dir <- readline(prompt = "enter output directory: ")
# Use the input path to set the working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  print("Directory does not exist. Please check the path.")
}
svglite("UpGPro.svg", width=6, height=3.5)
plot 
dev.off()
#for the most down regulated pathways
DnGoRes_path <- file.choose()
DnGoRes <- read.csv(DnGoRes_path)
#now order the table by the Pvalu and make it so most sig is first (ascending order)
DnGoRes <- DnGoRes[order(DnGoResf$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
DnGoResTop <- head(DnGoRes, n = 10)
DnGoResTop <- DnGoResTop[order(DnGoResTop$negative_log10_of_adjusted_p_value)]
#the blue color is for down regulated
plot<- ggplot(data = DnGoResTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "lightsteelblue1", high = "darkslategray4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,12)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
out_dir <- readline(prompt = "enter output directory: ")
# Use the input path to set the working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  print("Directory does not exist. Please check the path.")
}
svglite("DnGPro.svg", width=6, height=3.5)
plot 
dev.off()

#code used to make figures for IPA analysis - Figure 3
library(ggplot2)
#make new graphs from the IPA output
#take in a table of the IPA output for canonical pathways
#want to graph the -log(pvalue) and ZScore
#color bars by zScore and bar length by p-value
#only graph things that have a Z score of absval 1 or greater and greater than 0 pvalue
#read in the file, user file should have input from IPA analysis with pathway names, p-value and zscore
#ask the users for the file
IPA_path <- file.choose()
IPA <- read.csv(IPA_path, header = TRUE)
#subset on things that have a pvalue that isn't 0
IPAFilt <- subset(IPA, X.log.p.value. > 0 & abs(zScore) > .99)
#order the table by pvalue, so highest is first
IPAFilt <- IPAFilt[order(IPAFilt$X.log.p.value., decreasing = TRUE),]
#make the plot
plot <- (ggplot(data = IPAFilt, aes(x = reorder(Ingenuity.Canonical.Pathways, + X.log.p.value.), y = X.log.p.value.))+ 
           geom_col(aes(fill = zScore), width = .8)+
           scale_fill_gradient2(low = "aquamarine4", high = "violetred4")+
           coord_flip()+
           labs(x="Pathway", y = "-log(pvalue)", fill = "zScore")+
           ylim(0,5)+
           theme(aspect.ratio = 1.5/2)+
           theme(axis.title=element_text(size =10)))
#ask user to define the output directory
out_dir <- readline(prompt = "enter output directory: ")
# Use the input path to set the working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  print("Directory does not exist. Please check the path.")
}
library(svglite)
svglite("IPA.svg", width=7, height=5)
plot 
dev.off()


#code to make heatmaps seen in figure 2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
#required input is a gene expression tpm file for your samples of interest
#ask user for path to input file
tpm_file_path <- file.choose()
tpm_file <- read.csv(tpm_file_path, header = TRUE, sep ="\t")
# user will need to select genes of interest to make their list
cat("Enter genes you wish to used in the heatmap (one per line). Press Enter on an empty line when finished:\n")
word_vector_scan <- scan(what = character(), quiet = TRUE)

# Print the resulting vector (optional)
print(word_vector_scan)
#below is the gene list from the paper
#ShortIntersectList <- c("Cxcl10","Oasl1","Ccl5","Ccl2","Rsad2","Mx1","Ccl4","Lag3",
#                      "Ccl11","Oas3","Irf7","Gbp5","Ifit2","Ccl7","Gbp6",
#                     "Isg15","Oas2","Zbp1","Oas1g","Ifit3",
#                      "Ifit1","Gbp3","Cfb","Mnda","Bst2")
#make Rownames The gene names (make first column row names)
#had a duplicated gene but removed it assuming it is not a gene we want
tpm_update <- subset(tpm_file, gene!= "1-Mar")
tpm_update <- subset(tpm_update, gene!= "2-Mar")

#making gene names, row names
rownames(tpm_update) <- tpm_update[,1]

tpm_file2 <- tpm_update[,3:30]

mat <- as.matrix(tpm_file2[rownames(tpm_file2) %in% word_vector_scan, ])
#now to make heatmap
#Individual groups that will be included in the heatmap
#broken up by the tissue category and ordered how the user wants them to appear
#this is the order appearing in the figures
#the data in the paper is 4 different tissues with two conditions and organized with controls first
#this list can be changed by the user based on their input
Caput <- c("cap5","cap6","cap7","cap1","cap2","cap3","cap4")
Corpus <- c("cor5","cor6","cor7","cor1","cor2","cor3","cor4")
Cauda <- c("caud5","caud6","caud7","caud1","caud2","caud3","caud4")
Testes <- c("t5","t6","t7","t1","t2","t3","t4")
#vector to subset
#again vector should be changed by user based on input
KeepSamples <- c(Testes,Caput,Corpus,Cauda)
#subset matrix to only have poly samples
mat2 <- mat[, KeepSamples]
#change to follow user inpu
gene_groups <- factor(rep(c("Testes","Caput",  "Corpus", "Cauda"), 
                          times = c(length(Caput), length(Corpus), length(Cauda), length(Testes))),
                      levels = c("Testes","Caput",  "Corpus", "Cauda"))
#add in a vector to make vector for treatment
# Ensure matrix rows are in the same order as identity
mat2 <- mat2[match(word_vector_scan, rownames(mat2)), ]

myCol <- colorRampPalette(c('grey90', 'red2'))(50)
myBreaks <- seq(-10, 20, length.out = 50)
COLORS = list()
#user can change the colors and the names based upon their data
COLORS[["Tissue"]] = c("Testes" = "#489FA7",
                       "Caput" = "#584482",
                       "Corpus" = "#A4B0D1",
                       "Cauda" = "#883268")

#install.packages("colorRamp2")
library(colorRamp2)

#log transform data to make this look better
#make the poly and contol conditions using col and heatmap annotation
#user should edit to reflect the conditions that they want to show
Cols = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = COLORS[["Tissue"]])), Condition = c(rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4)),
                         col = list(Condition = c("Poly" = "darkred", "Cont" = "snow3")))
mat2 <- log2(mat2+1)
myBreaks <- seq(0, 12, length.out = 50)
heatmap <- Heatmap(mat2, name = "Expression in 
Log2(TPM)",  
                   column_split = factor(gene_groups),
                   cluster_columns = FALSE,
                   show_column_dend = FALSE,
                   cluster_column_slices = TRUE,
                   column_title_gp = gpar(fontsize = 10),
                   column_gap = unit(0.4, "mm"),
                   cluster_rows = TRUE,
                   show_row_dend = FALSE,
                   col = colorRamp2(myBreaks, myCol),
                   row_names_gp = gpar(fontsize = 9),
                   column_title_rot = 0,
                   top_annotation = Cols,
                   show_column_names = FALSE,
                   use_raster = TRUE)
#ask user to define the output directory
out_dir <- readline(prompt = "enter output directory: ")
# Use the input path to set the working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  print("Directory does not exist. Please check the path.")
}
library(svglite)
svglite("Heatmap.svg", width=5.5, height=5.5)
heatmap 
dev.off()



