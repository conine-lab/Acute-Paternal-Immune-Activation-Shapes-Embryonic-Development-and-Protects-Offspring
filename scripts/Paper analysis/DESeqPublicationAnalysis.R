#DESeq analysis used in Acute Paternal Immune Activation Shapes... TME et al. 2026 publication
#this code is all hardcoded and reflects the analysis that was used in the paper

#code for figure 2 analysis 
library(DESeq2)
library(ggplot2)
library(data.table)
#read csv file into R [this has raw count data for all replicates and gene name in first column]
#so made a file with the removed embryos and ordered by condition and sex usind the expected counts output, not counts tpm
#actually made another file that is ordered by condition and by sex becaus eI can't figure it out with subsetting
TME62_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/TME62_Filt_order_Gene_expression_expected.csv', header = TRUE, sep =",")
#preserve gene to transcript match
genetotranscript <- TME62_unfiltered[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(TME62_unfiltered) <- TME62_unfiltered[,2]
#now remove column 1 (genes) from data frame
TME62_unfiltered <- TME62_unfiltered[,-1]
#now remove the additional transcfript row
TME62_unfiltered <- TME62_unfiltered[,-1]
#name each column of your data as replicates i.e. con1, con2 etc.
names(TME62_unfiltered) <- c(colnames(TME62_unfiltered))
# I think this step was unneccessary and they wer actually already named
#change all values to integers - do this for each samples
TME62_unfiltered$TME51Cont_4Cell_1<- as.integer(TME62_unfiltered$TME51Cont_4Cell_1)
TME62_unfiltered$TME51Cont_4Cell_2<- as.integer(TME62_unfiltered$TME51Cont_4Cell_2)
TME62_unfiltered$TME51Cont_4Cell_4<- as.integer(TME62_unfiltered$TME51Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_2<- as.integer(TME62_unfiltered$TME52Cont_4Cell_2)
TME62_unfiltered$TME52Cont_4Cell_4<- as.integer(TME62_unfiltered$TME52Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_5<- as.integer(TME62_unfiltered$TME52Cont_4Cell_5)
TME62_unfiltered$TME52Cont_4Cell_6<- as.integer(TME62_unfiltered$TME52Cont_4Cell_6)
TME62_unfiltered$TME52Cont_4Cell_8<- as.integer(TME62_unfiltered$TME52Cont_4Cell_8)
TME62_unfiltered$TME55Cont_4Cell_2<- as.integer(TME62_unfiltered$TME55Cont_4Cell_2)
TME62_unfiltered$TME55Cont_4Cell_3<- as.integer(TME62_unfiltered$TME55Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_1<- as.integer(TME62_unfiltered$TME59Cont_4Cell_1)
TME62_unfiltered$TME59Cont_4Cell_2<- as.integer(TME62_unfiltered$TME59Cont_4Cell_2)
TME62_unfiltered$TME59Cont_4Cell_3<- as.integer(TME62_unfiltered$TME59Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_4<- as.integer(TME62_unfiltered$TME59Cont_4Cell_4)
TME62_unfiltered$TME59Cont_4Cell_5<- as.integer(TME62_unfiltered$TME59Cont_4Cell_5)
TME62_unfiltered$TME59Cont_4Cell_6<- as.integer(TME62_unfiltered$TME59Cont_4Cell_6)
TME62_unfiltered$TME59Cont_4Cell_7<- as.integer(TME62_unfiltered$TME59Cont_4Cell_7)
TME62_unfiltered$TME59Cont_4Cell_8<- as.integer(TME62_unfiltered$TME59Cont_4Cell_8)
TME62_unfiltered$TME51Poly_4Cell_1<- as.integer(TME62_unfiltered$TME51Poly_4Cell_1)
TME62_unfiltered$TME51Poly_4Cell_2<- as.integer(TME62_unfiltered$TME51Poly_4Cell_2)
TME62_unfiltered$TME51Poly_4Cell_3<- as.integer(TME62_unfiltered$TME51Poly_4Cell_3)
TME62_unfiltered$TME51Poly_4Cell_4<- as.integer(TME62_unfiltered$TME51Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_1<- as.integer(TME62_unfiltered$TME52Poly_4Cell_1)
TME62_unfiltered$TME52Poly_4Cell_2<- as.integer(TME62_unfiltered$TME52Poly_4Cell_2)
TME62_unfiltered$TME52Poly_4Cell_3<- as.integer(TME62_unfiltered$TME52Poly_4Cell_3)
TME62_unfiltered$TME52Poly_4Cell_4<- as.integer(TME62_unfiltered$TME52Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_5<- as.integer(TME62_unfiltered$TME52Poly_4Cell_5)
TME62_unfiltered$TME52Poly_4Cell_6<- as.integer(TME62_unfiltered$TME52Poly_4Cell_6)
TME62_unfiltered$TME52Poly_4Cell_7<- as.integer(TME62_unfiltered$TME52Poly_4Cell_7)
TME62_unfiltered$TME52Poly_4Cell_8<- as.integer(TME62_unfiltered$TME52Poly_4Cell_8)
TME62_unfiltered$TME55Poly_4Cell_1<- as.integer(TME62_unfiltered$TME55Poly_4Cell_1)
TME62_unfiltered$TME55Poly_4Cell_2<- as.integer(TME62_unfiltered$TME55Poly_4Cell_2)
TME62_unfiltered$TME55Poly_4Cell_3<- as.integer(TME62_unfiltered$TME55Poly_4Cell_3)
TME62_unfiltered$TME55Poly_4Cell_4<- as.integer(TME62_unfiltered$TME55Poly_4Cell_4)
TME62_unfiltered$TME55Poly_4Cell_5<- as.integer(TME62_unfiltered$TME55Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_1<- as.integer(TME62_unfiltered$TME59Poly_4Cell_1)
TME62_unfiltered$TME59Poly_4Cell_2<- as.integer(TME62_unfiltered$TME59Poly_4Cell_2)
TME62_unfiltered$TME59Poly_4Cell_3<- as.integer(TME62_unfiltered$TME59Poly_4Cell_3)
TME62_unfiltered$TME59Poly_4Cell_5<- as.integer(TME62_unfiltered$TME59Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_6<- as.integer(TME62_unfiltered$TME59Poly_4Cell_6)
TME62_unfiltered$TME59Poly_4Cell_7<- as.integer(TME62_unfiltered$TME59Poly_4Cell_7)
TME62_unfiltered$TME59Poly_4Cell_8<- as.integer(TME62_unfiltered$TME59Poly_4Cell_8)
TME62_unfiltered$TME52Cont_Morula_1<- as.integer(TME62_unfiltered$TME52Cont_Morula_1)
TME62_unfiltered$TME52Cont_Morula_2<- as.integer(TME62_unfiltered$TME52Cont_Morula_2)
TME62_unfiltered$TME52Cont_Morula_3<- as.integer(TME62_unfiltered$TME52Cont_Morula_3)
TME62_unfiltered$TME52Cont_Morula_4<- as.integer(TME62_unfiltered$TME52Cont_Morula_4)
TME62_unfiltered$TME52Cont_Morula_5<- as.integer(TME62_unfiltered$TME52Cont_Morula_5)
TME62_unfiltered$TME52Cont_Morula_6<- as.integer(TME62_unfiltered$TME52Cont_Morula_6)
TME62_unfiltered$TME52Cont_Morula_8<- as.integer(TME62_unfiltered$TME52Cont_Morula_8)
TME62_unfiltered$TME55Cont_Morula_1<- as.integer(TME62_unfiltered$TME55Cont_Morula_1)
TME62_unfiltered$TME55Cont_Morula_2<- as.integer(TME62_unfiltered$TME55Cont_Morula_2)
TME62_unfiltered$TME55Cont_Morula_3<- as.integer(TME62_unfiltered$TME55Cont_Morula_3)
TME62_unfiltered$TME55Cont_Morula_4<- as.integer(TME62_unfiltered$TME55Cont_Morula_4)
TME62_unfiltered$TME55Cont_Morula_5<- as.integer(TME62_unfiltered$TME55Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_1<- as.integer(TME62_unfiltered$TME59Cont_Morula_1)
TME62_unfiltered$TME59Cont_Morula_2<- as.integer(TME62_unfiltered$TME59Cont_Morula_2)
TME62_unfiltered$TME59Cont_Morula_3<- as.integer(TME62_unfiltered$TME59Cont_Morula_3)
TME62_unfiltered$TME59Cont_Morula_4<- as.integer(TME62_unfiltered$TME59Cont_Morula_4)
TME62_unfiltered$TME59Cont_Morula_5<- as.integer(TME62_unfiltered$TME59Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_6<- as.integer(TME62_unfiltered$TME59Cont_Morula_6)
TME62_unfiltered$TME59Cont_Morula_7<- as.integer(TME62_unfiltered$TME59Cont_Morula_7)
TME62_unfiltered$TME59Cont_Morula_8<- as.integer(TME62_unfiltered$TME59Cont_Morula_8)
TME62_unfiltered$TME52Poly_Morula_1<- as.integer(TME62_unfiltered$TME52Poly_Morula_1)
TME62_unfiltered$TME52Poly_Morula_2<- as.integer(TME62_unfiltered$TME52Poly_Morula_2)
TME62_unfiltered$TME52Poly_Morula_3<- as.integer(TME62_unfiltered$TME52Poly_Morula_3)
TME62_unfiltered$TME52Poly_Morula_4<- as.integer(TME62_unfiltered$TME52Poly_Morula_4)
TME62_unfiltered$TME52Poly_Morula_5<- as.integer(TME62_unfiltered$TME52Poly_Morula_5)
TME62_unfiltered$TME52Poly_Morula_6<- as.integer(TME62_unfiltered$TME52Poly_Morula_6)
TME62_unfiltered$TME52Poly_Morula_8<- as.integer(TME62_unfiltered$TME52Poly_Morula_8)
TME62_unfiltered$TME55Poly_Morula_1<- as.integer(TME62_unfiltered$TME55Poly_Morula_1)
TME62_unfiltered$TME55Poly_Morula_2<- as.integer(TME62_unfiltered$TME55Poly_Morula_2)
TME62_unfiltered$TME55Poly_Morula_3<- as.integer(TME62_unfiltered$TME55Poly_Morula_3)
TME62_unfiltered$TME55Poly_Morula_4<- as.integer(TME62_unfiltered$TME55Poly_Morula_4)
TME62_unfiltered$TME55Poly_Morula_5<- as.integer(TME62_unfiltered$TME55Poly_Morula_5)
TME62_unfiltered$TME55Poly_Morula_6<- as.integer(TME62_unfiltered$TME55Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_3<- as.integer(TME62_unfiltered$TME59Poly_Morula_3)
TME62_unfiltered$TME59Poly_Morula_4<- as.integer(TME62_unfiltered$TME59Poly_Morula_4)
TME62_unfiltered$TME59Poly_Morula_5<- as.integer(TME62_unfiltered$TME59Poly_Morula_5)
TME62_unfiltered$TME59Poly_Morula_6<- as.integer(TME62_unfiltered$TME59Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_8<- as.integer(TME62_unfiltered$TME59Poly_Morula_8)

#filter data - here I use a cutoff of 81 across all replicates, this could be higher
TME62_filtered <- TME62_unfiltered
#split things up by dev stage, will process each stage individually since I'm not comparing across

#split things up by dev stage , will process each stage individually since I'm not comparing across
TME62_filt_4Cell <- TME62_filtered[,1:42]
TME62_filt_Morula <- TME62_filtered[,43:80]

#create a frame with condition names, this can be used for each tissue since they are in the same order
condition_4Cell <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly", "cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly")

#create a data frame for use with DESeq2, make this for each tissue type
colData4Cell <- data.frame(row.names =colnames(TME62_filt_4Cell), condition_4Cell)
View(colData4Cell)

colDataMorula <- data.frame(row.names =colnames(TME62_filt_Morula), condition_Morula)
View(colDataMorula)


#run DESeq2, do this for each condition and sex
TME62_4Cell_dds <- DESeqDataSetFromMatrix(countData = TME62_filt_4Cell, colData = colData4Cell, design = ~condition_4Cell)

TME62_Morula_dds <-  DESeqDataSetFromMatrix(countData = TME62_filt_Morula, colData = colDataMorula, design = ~condition_Morula)

#for each condition in each group, we need to figure out the number
#so number of female control 4 cell, female poly 4 cell, female cont morula, etc
colnames(TME62_4Cell_dds)
colnames(TME62_Morula_dds)
#used this to determine the cut offs (number in each group divided by 2)
cell4_min_Poly <- 12
cell4_min_Cont <- 9

morula_min_Cont <- 9
morula_min_Poly <- 10


#now move on to determining which rows we will keep because the counts are great enough
#determining if the at least half of the samples in either condition have counts for the gene, 5 tpm cut off
keep <- rowSums(counts(TME62_4Cell_dds[,c(1:11, 19:25)])>=5)>= cell4_min_Poly | rowSums(counts(TME62_4Cell_dds[,c(12:18,26:42)])>=5)>= cell4_min_Cont
TME62_4Cell_dds <- TME62_4Cell_dds[keep,]
keep <- rowSums(counts(TME62_Morula_dds[,c(1:13,24:30)])>=5)>= morula_min_Cont | rowSums(counts(TME62_Morula_dds[,c(14:23,31:38)])>=5)>= morula_min_Poly
TME62_Morula_dds <- TME62_Morula_dds[keep,]

#now can continue with same code as before
#run DESeq2, do this for each sample set
TME62_4Cell_dds <- DESeq(TME62_4Cell_dds)
TME62_Morula_dds <- DESeq(TME62_Morula_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each tissue
TME62_4Cell_results <- results(TME62_4Cell_dds, contrast = c("condition_4Cell", "poly", "cont"))

TME62_Morula_results <- results(TME62_Morula_dds, contrast = c("condition_Morula", "poly", "cont"))

#run a summary on each set of samples
summary(TME62_4Cell_results)
summary(TME62_Morula_results)

#write res to CSV do for each tissue
write.csv(TME62_4Cell_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_4Cell_res_DESeq2.csv', row.names = TRUE)

write.csv(TME62_Morula_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_Morula_res_DESeq2.csv', row.names = TRUE)

#pull out significant up and down res Changing to use pvalue, 
Cell_4Sig <- subset(TME62_4Cell_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
Cell_4Sigup <- subset (Cell_4Sig,log2FoldChange > 0)                                                          
Cell_4Sigdn <- subset (Cell_4Sig, log2FoldChange < 0)
write.csv(Cell_4Sigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_Sig_res_up.csv', row.names = TRUE)
write.csv(Cell_4Sigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_Sig_res_dn.csv', row.names = TRUE)

Morula_Sig <- subset(TME62_Morula_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
MorulaSigup <- subset (Morula_Sig,log2FoldChange > 0)                                                          
MorulaSigdn <- subset (Morula_Sig, log2FoldChange < 0)
write.csv(MorulaSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_Sig_res_up.csv', row.names = TRUE)
write.csv(MorulaSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_Sig_res_dn.csv', row.names = TRUE)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
Cell_4SigupDF <- as.data.frame(Cell_4Sigup)
Cell_4SigdnDF <- as.data.frame(Cell_4Sigdn)
MorulaSigupDF <- as.data.frame(MorulaSigup)
MorulaSigdnDF <- as.data.frame(MorulaSigdn)

#now make the geneIDs which is a row name a column in the df so we can match with the gene name
Cell_4SigupDF$geneID <- rownames(Cell_4SigupDF)
Cell_4SigdnDF$geneID <- rownames(Cell_4SigdnDF)
MorulaSigupDF$geneID <- rownames(MorulaSigupDF)
MorulaSigdnDF$geneID <- rownames(MorulaSigdnDF)

#trying to figure out how to get the gene name back
Cell_4SigupDF <- merge(Cell_4SigupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigdnDF <- merge(Cell_4SigdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigupDF <- merge(MorulaSigupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigdnDF <- merge(MorulaSigdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")


#try to switch column order so the gene name isn't at the end
setcolorder(Cell_4SigupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

setcolorder(MorulaSigupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(Cell_4SigupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_Sig_res_up_withGene_ByPvalu01.csv', row.names = TRUE)
write.csv(Cell_4SigdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_Sig_res_dn_withGene_ByPvalu01.csv', row.names = TRUE)

write.csv(MorulaSigupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_Sig_res_up_withGene_ByPvalue01.csv', row.names = TRUE)
write.csv(MorulaSigdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_Sig_res_dn_withGene_ByPvalue01.csv', row.names = TRUE)

#analysis for DEGS by embryo sex
library(DESeq2)
library(ggplot2)
library(data.table)
#read csv file into R [this has raw count data for all replicates and gene name in first column]
#so made a file with the removed embryos and ordered by condition and sex usind the expected counts output, not counts tpm
#actually made another file that is ordered by condition and by sex becaus eI can't figure it out with subsetting
TME62_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/TME62_Filt_order_Gene_expression_expected.csv', header = TRUE, sep =",")
#preserve gene to transcript match
genetotranscript <- TME62_unfiltered[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(TME62_unfiltered) <- TME62_unfiltered[,2]
#now remove column 1 (genes) from data frame
TME62_unfiltered <- TME62_unfiltered[,-1]
#now remove the additional transcfript row
TME62_unfiltered <- TME62_unfiltered[,-1]
#name each column of your data as replicates i.e. con1, con2 etc.
names(TME62_unfiltered) <- c(colnames(TME62_unfiltered))
# I think this step was unneccessary and they wer actually already named
#change all values to integers - do this for each samples
TME62_unfiltered$TME51Cont_4Cell_1<- as.integer(TME62_unfiltered$TME51Cont_4Cell_1)
TME62_unfiltered$TME51Cont_4Cell_2<- as.integer(TME62_unfiltered$TME51Cont_4Cell_2)
TME62_unfiltered$TME51Cont_4Cell_4<- as.integer(TME62_unfiltered$TME51Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_2<- as.integer(TME62_unfiltered$TME52Cont_4Cell_2)
TME62_unfiltered$TME52Cont_4Cell_4<- as.integer(TME62_unfiltered$TME52Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_5<- as.integer(TME62_unfiltered$TME52Cont_4Cell_5)
TME62_unfiltered$TME52Cont_4Cell_6<- as.integer(TME62_unfiltered$TME52Cont_4Cell_6)
TME62_unfiltered$TME52Cont_4Cell_8<- as.integer(TME62_unfiltered$TME52Cont_4Cell_8)
TME62_unfiltered$TME55Cont_4Cell_2<- as.integer(TME62_unfiltered$TME55Cont_4Cell_2)
TME62_unfiltered$TME55Cont_4Cell_3<- as.integer(TME62_unfiltered$TME55Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_1<- as.integer(TME62_unfiltered$TME59Cont_4Cell_1)
TME62_unfiltered$TME59Cont_4Cell_2<- as.integer(TME62_unfiltered$TME59Cont_4Cell_2)
TME62_unfiltered$TME59Cont_4Cell_3<- as.integer(TME62_unfiltered$TME59Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_4<- as.integer(TME62_unfiltered$TME59Cont_4Cell_4)
TME62_unfiltered$TME59Cont_4Cell_5<- as.integer(TME62_unfiltered$TME59Cont_4Cell_5)
TME62_unfiltered$TME59Cont_4Cell_6<- as.integer(TME62_unfiltered$TME59Cont_4Cell_6)
TME62_unfiltered$TME59Cont_4Cell_7<- as.integer(TME62_unfiltered$TME59Cont_4Cell_7)
TME62_unfiltered$TME59Cont_4Cell_8<- as.integer(TME62_unfiltered$TME59Cont_4Cell_8)
TME62_unfiltered$TME51Poly_4Cell_1<- as.integer(TME62_unfiltered$TME51Poly_4Cell_1)
TME62_unfiltered$TME51Poly_4Cell_2<- as.integer(TME62_unfiltered$TME51Poly_4Cell_2)
TME62_unfiltered$TME51Poly_4Cell_3<- as.integer(TME62_unfiltered$TME51Poly_4Cell_3)
TME62_unfiltered$TME51Poly_4Cell_4<- as.integer(TME62_unfiltered$TME51Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_1<- as.integer(TME62_unfiltered$TME52Poly_4Cell_1)
TME62_unfiltered$TME52Poly_4Cell_2<- as.integer(TME62_unfiltered$TME52Poly_4Cell_2)
TME62_unfiltered$TME52Poly_4Cell_3<- as.integer(TME62_unfiltered$TME52Poly_4Cell_3)
TME62_unfiltered$TME52Poly_4Cell_4<- as.integer(TME62_unfiltered$TME52Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_5<- as.integer(TME62_unfiltered$TME52Poly_4Cell_5)
TME62_unfiltered$TME52Poly_4Cell_6<- as.integer(TME62_unfiltered$TME52Poly_4Cell_6)
TME62_unfiltered$TME52Poly_4Cell_7<- as.integer(TME62_unfiltered$TME52Poly_4Cell_7)
TME62_unfiltered$TME52Poly_4Cell_8<- as.integer(TME62_unfiltered$TME52Poly_4Cell_8)
TME62_unfiltered$TME55Poly_4Cell_1<- as.integer(TME62_unfiltered$TME55Poly_4Cell_1)
TME62_unfiltered$TME55Poly_4Cell_2<- as.integer(TME62_unfiltered$TME55Poly_4Cell_2)
TME62_unfiltered$TME55Poly_4Cell_3<- as.integer(TME62_unfiltered$TME55Poly_4Cell_3)
TME62_unfiltered$TME55Poly_4Cell_4<- as.integer(TME62_unfiltered$TME55Poly_4Cell_4)
TME62_unfiltered$TME55Poly_4Cell_5<- as.integer(TME62_unfiltered$TME55Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_1<- as.integer(TME62_unfiltered$TME59Poly_4Cell_1)
TME62_unfiltered$TME59Poly_4Cell_2<- as.integer(TME62_unfiltered$TME59Poly_4Cell_2)
TME62_unfiltered$TME59Poly_4Cell_3<- as.integer(TME62_unfiltered$TME59Poly_4Cell_3)
TME62_unfiltered$TME59Poly_4Cell_5<- as.integer(TME62_unfiltered$TME59Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_6<- as.integer(TME62_unfiltered$TME59Poly_4Cell_6)
TME62_unfiltered$TME59Poly_4Cell_7<- as.integer(TME62_unfiltered$TME59Poly_4Cell_7)
TME62_unfiltered$TME59Poly_4Cell_8<- as.integer(TME62_unfiltered$TME59Poly_4Cell_8)
TME62_unfiltered$TME52Cont_Morula_1<- as.integer(TME62_unfiltered$TME52Cont_Morula_1)
TME62_unfiltered$TME52Cont_Morula_2<- as.integer(TME62_unfiltered$TME52Cont_Morula_2)
TME62_unfiltered$TME52Cont_Morula_3<- as.integer(TME62_unfiltered$TME52Cont_Morula_3)
TME62_unfiltered$TME52Cont_Morula_4<- as.integer(TME62_unfiltered$TME52Cont_Morula_4)
TME62_unfiltered$TME52Cont_Morula_5<- as.integer(TME62_unfiltered$TME52Cont_Morula_5)
TME62_unfiltered$TME52Cont_Morula_6<- as.integer(TME62_unfiltered$TME52Cont_Morula_6)
TME62_unfiltered$TME52Cont_Morula_8<- as.integer(TME62_unfiltered$TME52Cont_Morula_8)
TME62_unfiltered$TME55Cont_Morula_1<- as.integer(TME62_unfiltered$TME55Cont_Morula_1)
TME62_unfiltered$TME55Cont_Morula_2<- as.integer(TME62_unfiltered$TME55Cont_Morula_2)
TME62_unfiltered$TME55Cont_Morula_3<- as.integer(TME62_unfiltered$TME55Cont_Morula_3)
TME62_unfiltered$TME55Cont_Morula_4<- as.integer(TME62_unfiltered$TME55Cont_Morula_4)
TME62_unfiltered$TME55Cont_Morula_5<- as.integer(TME62_unfiltered$TME55Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_1<- as.integer(TME62_unfiltered$TME59Cont_Morula_1)
TME62_unfiltered$TME59Cont_Morula_2<- as.integer(TME62_unfiltered$TME59Cont_Morula_2)
TME62_unfiltered$TME59Cont_Morula_3<- as.integer(TME62_unfiltered$TME59Cont_Morula_3)
TME62_unfiltered$TME59Cont_Morula_4<- as.integer(TME62_unfiltered$TME59Cont_Morula_4)
TME62_unfiltered$TME59Cont_Morula_5<- as.integer(TME62_unfiltered$TME59Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_6<- as.integer(TME62_unfiltered$TME59Cont_Morula_6)
TME62_unfiltered$TME59Cont_Morula_7<- as.integer(TME62_unfiltered$TME59Cont_Morula_7)
TME62_unfiltered$TME59Cont_Morula_8<- as.integer(TME62_unfiltered$TME59Cont_Morula_8)
TME62_unfiltered$TME52Poly_Morula_1<- as.integer(TME62_unfiltered$TME52Poly_Morula_1)
TME62_unfiltered$TME52Poly_Morula_2<- as.integer(TME62_unfiltered$TME52Poly_Morula_2)
TME62_unfiltered$TME52Poly_Morula_3<- as.integer(TME62_unfiltered$TME52Poly_Morula_3)
TME62_unfiltered$TME52Poly_Morula_4<- as.integer(TME62_unfiltered$TME52Poly_Morula_4)
TME62_unfiltered$TME52Poly_Morula_5<- as.integer(TME62_unfiltered$TME52Poly_Morula_5)
TME62_unfiltered$TME52Poly_Morula_6<- as.integer(TME62_unfiltered$TME52Poly_Morula_6)
TME62_unfiltered$TME52Poly_Morula_8<- as.integer(TME62_unfiltered$TME52Poly_Morula_8)
TME62_unfiltered$TME55Poly_Morula_1<- as.integer(TME62_unfiltered$TME55Poly_Morula_1)
TME62_unfiltered$TME55Poly_Morula_2<- as.integer(TME62_unfiltered$TME55Poly_Morula_2)
TME62_unfiltered$TME55Poly_Morula_3<- as.integer(TME62_unfiltered$TME55Poly_Morula_3)
TME62_unfiltered$TME55Poly_Morula_4<- as.integer(TME62_unfiltered$TME55Poly_Morula_4)
TME62_unfiltered$TME55Poly_Morula_5<- as.integer(TME62_unfiltered$TME55Poly_Morula_5)
TME62_unfiltered$TME55Poly_Morula_6<- as.integer(TME62_unfiltered$TME55Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_3<- as.integer(TME62_unfiltered$TME59Poly_Morula_3)
TME62_unfiltered$TME59Poly_Morula_4<- as.integer(TME62_unfiltered$TME59Poly_Morula_4)
TME62_unfiltered$TME59Poly_Morula_5<- as.integer(TME62_unfiltered$TME59Poly_Morula_5)
TME62_unfiltered$TME59Poly_Morula_6<- as.integer(TME62_unfiltered$TME59Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_8<- as.integer(TME62_unfiltered$TME59Poly_Morula_8)

#filter data - here I use a cutoff of 81 across all replicates, this could be higher
TME62_filtered <- TME62_unfiltered
#split things up by dev stage, will process each stage individually since I'm not comparing across

#split things up by dev stage , will process each stage individually since I'm not comparing across
TME62_filt_4Cell <- TME62_filtered[,1:42]
TME62_filt_Morula <- TME62_filtered[,43:80]
#based on xist gene and Eif2s3y gene in EmbryoNewOrderFiltered pull out each based on sex
#filter into diff sex
TME62_Female_4Cell <- TME62_filt_4Cell[,1:18]
TME62_Male_4Cell <- TME62_filt_4Cell[,19:42]
TME62_Female_Morula <- TME62_filt_Morula[,1:23]
TME62_Male_Morula <- TME62_filt_Morula[,24:38]

#create a frame with condition names, this can be used for each tissue since they are in the same order
condition_4Cell_Female <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly")
condition_4Cell_Male <- c("cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula_Female <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula_Male <- c("cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly")

#condition_morula <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
#create a data frame for use with DESeq2, make this for each tissue type
colData4CellFemale <- data.frame(row.names =colnames(TME62_Female_4Cell), condition_4Cell_Female)
View(colData4CellFemale)
colData4CellMale <- data.frame(row.names =colnames(TME62_Male_4Cell), condition_4Cell_Male)
View(colData4CellFemale)
#try to add to dataframe and separate by sex, see embryoNewOrder sheet for sex determinataion
#Cell4_sex <- c("F","M","F","M","F","F","F","F","F","F","F","M","M","M","F","M","F","M","M","F","M","M","M","M","M","M","F","M","M","M","F","M","M","F","M","M","F","M","M","F","F","M")
#MoreColData4Cell <- data.frame(row.names = colnames(TME62_filt_4Cell), condition_4Cell, Cell4_sex)
#can use the above dataframe to subset based on sex solumn
#Male_4Cell <- subset(MoreColData4Cell, Cell4_sex == 'M')

colDataMorulaFemale <- data.frame(row.names =colnames(TME62_Female_Morula), condition_Morula_Female)
View(colDataMorulaFemale)
colDataMorulaMale <- data.frame(row.names =colnames(TME62_Male_Morula), condition_Morula_Male)
View(colDataMorulaMale)

#run DESeq2, do this for each condition and sex
TME62_4Cell_Female_dds <- DESeqDataSetFromMatrix(countData = TME62_Female_4Cell, colData = colData4CellFemale, design = ~condition_4Cell_Female)
TME62_4Cell_Male_dds <- DESeqDataSetFromMatrix(countData = TME62_Male_4Cell, colData = colData4CellMale, design = ~condition_4Cell_Male)

TME62_Morula_Female_dds <-  DESeqDataSetFromMatrix(countData = TME62_Female_Morula, colData = colDataMorulaFemale, design = ~condition_Morula_Female)
TME62_Morula_Male_dds <-  DESeqDataSetFromMatrix(countData = TME62_Male_Morula, colData = colDataMorulaMale, design = ~condition_Morula_Male)

#for each condition in each group, we need to figure out the number
#so number of female control 4 cell, female poly 4 cell, female cont morula, etc
colnames(TME62_4Cell_Female_dds)
#used this to determine the cut offs (number in each group divided by 2)
cell4_Female_min_Poly <- 3.5
cell4_Female_min_Cont <- 5.5
colnames(TME62_4Cell_Male_dds)
cell4_Male_min_Cont <- 3.5
cell4_Male_min_Poly <- 8.5
colnames(TME62_Morula_Female_dds)
morula_Female_min_Cont <- 6.5
morula_Female_min_Poly <- 5
colnames(TME62_Morula_Male_dds)
morula_Male_min_Cont <- 3.5
morula_Male_min_Poly <- 4

#now move on to determining which rows we will keep because the counts are great enough
#determining if the at least half of the samples in either condition have counts for the gene, 5 tpm cut off
keep <- rowSums(counts(TME62_4Cell_Female_dds[,12:18])>=5)>= cell4_Female_min_Poly | rowSums(counts(TME62_4Cell_Female_dds[,1:11])>=5)>= cell4_Female_min_Cont
TME62_4Cell_Female_dds <- TME62_4Cell_Female_dds[keep,]
keep <- rowSums(counts(TME62_4Cell_Male_dds[,1:7])>=5)>= cell4_Male_min_Cont | rowSums(counts(TME62_4Cell_Male_dds[,8:24])>=5)>= cell4_Male_min_Poly
TME62_4Cell_Male_dds <- TME62_4Cell_Male_dds[keep,]
keep <- rowSums(counts(TME62_Morula_Female_dds[,1:13])>=5)>= morula_Female_min_Cont | rowSums(counts(TME62_Morula_Female_dds[,14:23])>=5)>= morula_Female_min_Poly
TME62_Morula_Female_dds <- TME62_Morula_Female_dds[keep,]
keep <- rowSums(counts(TME62_Morula_Male_dds[,1:7])>=5)>= morula_Male_min_Cont | rowSums(counts(TME62_Morula_Male_dds[,8:15])>=5)>= morula_Male_min_Poly
TME62_Morula_Male_dds <- TME62_Morula_Male_dds[keep,]

#now can continue with same code as before
#run DESeq2, do this for each sample set
TME62_4Cell_Female_dds <- DESeq(TME62_4Cell_Female_dds)
TME62_4Cell_Male_dds <- DESeq(TME62_4Cell_Male_dds)
TME62_Morula_Female_dds <- DESeq(TME62_Morula_Female_dds)
TME62_Morula_Male_dds <- DESeq(TME62_Morula_Male_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each tissue
TME62_4Cell_Female_results <- results(TME62_4Cell_Female_dds, contrast = c("condition_4Cell_Female", "poly", "cont"))
TME62_4Cell_Male_results <- results(TME62_4Cell_Male_dds, contrast = c("condition_4Cell_Male", "poly", "cont"))

TME62_Morula_Female_results <- results(TME62_Morula_Female_dds, contrast = c("condition_Morula_Female", "poly", "cont"))
TME62_Morula_MAle_results <- results(TME62_Morula_Male_dds, contrast = c("condition_Morula_Male", "poly", "cont"))

#run a summary on each set of samples
summary(TME62_4Cell_Male_results)
summary(TME62_Morula_MAle_results)

#write res to CSV do for each tissue
write.csv(TME62_4Cell_Female_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_4Cell_FemOnly_res_DESeq2.csv', row.names = TRUE)
write.csv(TME62_4Cell_Male_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_4Cell_MaleOnly_res_DESeq2.csv', row.names = TRUE)

write.csv(TME62_Morula_Female_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_Morula_FemOnly_res_DESeq2.csv', row.names = TRUE)
write.csv(TME62_Morula_MAle_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/TME62_Morula_MaleOnly_res_DESeq2.csv', row.names = TRUE)

#pull out significant up and down res Changing to use pvalue, 
Cell_4SigFem <- subset(TME62_4Cell_Female_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
Cell_4SigFemup <- subset (Cell_4SigFem,log2FoldChange > 0)                                                          
Cell_4SigFemdn <- subset (Cell_4SigFem, log2FoldChange < 0)
write.csv(Cell_4SigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_up.csv', row.names = TRUE)
write.csv(Cell_4SigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_dn.csv', row.names = TRUE)

#pull out significant up and down res Changing to use pvalue
Cell_4SigMale <- subset(TME62_4Cell_Male_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
Cell_4SigMaleup <- subset (Cell_4SigMale,log2FoldChange > 0)                                                          
Cell_4SigMaledn <- subset (Cell_4SigMale, log2FoldChange < 0)
write.csv(Cell_4SigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_up.csv', row.names = TRUE)
write.csv(Cell_4SigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_dn.csv', row.names = TRUE)

#pull out significant up and down res, 1.5 FC or greater, changing ot pvalue
Morula_SigMale <- subset(TME62_Morula_MAle_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
MorulaSigMaleup <- subset (Morula_SigMale,log2FoldChange > 0)                                                          
MorulaSigMaledn <- subset (Morula_SigMale, log2FoldChange < 0)
write.csv(MorulaSigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up.csv', row.names = TRUE)
write.csv(MorulaSigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_dn.csv', row.names = TRUE)

Morula_SigFemale <- subset(TME62_Morula_Female_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
MorulaSigFemup <- subset (Morula_SigFemale,log2FoldChange > 0)                                                          
MorulaSigFemdn <- subset (Morula_SigFemale, log2FoldChange < 0)
write.csv(MorulaSigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up.csv', row.names = TRUE)
write.csv(MorulaSigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_dn.csv', row.names = TRUE)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
Cell_4SigFemupDF <- as.data.frame(Cell_4SigFemup)
Cell_4SigMaleupDF <- as.data.frame(Cell_4SigMaleup)
Cell_4SigFemdnDF <- as.data.frame(Cell_4SigFemdn)
Cell_4SigMalednDF <- as.data.frame(Cell_4SigMaledn)
MorulaSigFemupDF <- as.data.frame(MorulaSigFemup)
MorulaSigMaleupDF <- as.data.frame(MorulaSigMaleup)
MorulaSigFemdnDF <- as.data.frame(MorulaSigFemdn)
MorulaSigMalednDF <- as.data.frame(MorulaSigMaledn)

#now make the geneIDs which is a row name a column in the df so we can match with the gene name
Cell_4SigFemupDF$geneID <- rownames(Cell_4SigFemupDF)
Cell_4SigMaleupDF$geneID <- rownames(Cell_4SigMaleupDF)
Cell_4SigFemdnDF$geneID <- rownames(Cell_4SigFemdnDF)
Cell_4SigMalednDF$geneID <- rownames(Cell_4SigMalednDF)
MorulaSigFemupDF$geneID <- rownames(MorulaSigFemupDF)
MorulaSigMaleupDF$geneID <- rownames(MorulaSigMaleupDF)
MorulaSigFemdnDF$geneID <- rownames(MorulaSigFemdnDF)
MorulaSigMalednDF$geneID <- rownames(MorulaSigMalednDF)

#trying to figure out how to get the gene name back
Cell_4SigFemupDF <- merge(Cell_4SigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMaleupDF <- merge(Cell_4SigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigFemdnDF <- merge(Cell_4SigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMalednDF <- merge(Cell_4SigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemupDF <- merge(MorulaSigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemdnDF <- merge(MorulaSigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMaleupDF <- merge(MorulaSigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMalednDF <- merge(MorulaSigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")


#try to switch column order so the gene name isn't at the end
setcolorder(Cell_4SigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

setcolorder(MorulaSigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(Cell_4SigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(Cell_4SigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(Cell_4SigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(Cell_4SigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)

write.csv(MorulaSigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(MorulaSigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__FemOnly_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(MorulaSigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(MorulaSigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__MaleOnly_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)

#make volcano plots
alpha <- 0.05
TME62_4Cell_Female_results$sig <- -log10(TME62_4Cell_Female_results$padj)
cols <- densCols(TME62_4Cell_Female_results$log2FoldChange, TME62_4Cell_Female_results$sig)
plot(TME62_4Cell_Female_results$log2FoldChange, TME62_4Cell_Female_results$sig,
     col= cols, main = "4 Cell Female", xlab = "log2FC", ylab = "-log10(padj)")
TME62_4Cell_Male_results$sig <- -log10(TME62_4Cell_Male_results$padj)
cols <- densCols(TME62_4Cell_Male_results$log2FoldChange, TME62_4Cell_Male_results$sig)
plot(TME62_4Cell_Male_results$log2FoldChange, TME62_4Cell_Male_results$sig,
     col= cols, main = "4 Cell Male", xlab = "log2FC", ylab = "padj")

TME62_Morula_Female_results$sig <- -log10(TME62_Morula_Female_results$padj)
cols <- densCols(TME62_Morula_Female_results$log2FoldChange, TME62_Morula_Female_results$sig)
plot(TME62_Morula_Female_results$log2FoldChange, TME62_Morula_Female_results$sig,
     col= cols, main = "Morula Female", xlab = "log2FC", ylab = "padj") 

TME62_Morula_MAle_results$sig <- -log10(TME62_Morula_MAle_results$padj)
cols <- densCols(TME62_Morula_MAle_results$log2FoldChange, TME62_Morula_MAle_results$sigPvalue)
plot(TME62_Morula_MAle_results$log2FoldChange, TME62_Morula_MAle_results$sig,
     col= cols, main = "Morula Male", xlab = "log2FC", ylab = "padj") 

#pull out significant up and down res Changing to use padj <.05 and 1.5 FC
Cell_4SigFem <- subset(TME62_4Cell_Female_results, padj <0.05 & abs(log2FoldChange) > 0.58)
Cell_4SigFemup <- subset (Cell_4SigFem,log2FoldChange > 0)                                                          
Cell_4SigFemdn <- subset (Cell_4SigFem, log2FoldChange < 0)
write.csv(Cell_4SigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_up_ByPadj.csv', row.names = TRUE)
write.csv(Cell_4SigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_dn_ByPadj.csv', row.names = TRUE)


Cell_4SigMale <- subset(TME62_4Cell_Male_results, padj <0.05 & abs(log2FoldChange) > 0.58)
Cell_4SigMaleup <- subset (Cell_4SigMale,log2FoldChange > 0)                                                          
Cell_4SigMaledn <- subset (Cell_4SigMale, log2FoldChange < 0)
write.csv(Cell_4SigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_up_ByPadj.csv', row.names = TRUE)
write.csv(Cell_4SigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_dn_ByPadj.csv', row.names = TRUE)

#pull out significant up and down res, 1.5 FC or greater
Morula_SigMale <- subset(TME62_Morula_MAle_results, padj <0.05 & abs(log2FoldChange) > 0.58)
MorulaSigMaleup <- subset (Morula_SigMale,log2FoldChange > 0)                                                          
MorulaSigMaledn <- subset (Morula_SigMale, log2FoldChange < 0)
write.csv(MorulaSigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up_ByPadj.csv', row.names = TRUE)
write.csv(MorulaSigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_dn_ByPadj.csv', row.names = TRUE)

Morula_SigFemale <- subset(TME62_Morula_Female_results, padj <0.05 & abs(log2FoldChange) > 0.58)
MorulaSigFemup <- subset (Morula_SigFemale,log2FoldChange > 0)                                                          
MorulaSigFemdn <- subset (Morula_SigFemale, log2FoldChange < 0)
write.csv(MorulaSigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up_ByPadj.csv', row.names = TRUE)
write.csv(MorulaSigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_dn_ByPadj.csv', row.names = TRUE)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
Cell_4SigFemupDF <- as.data.frame(Cell_4SigFemup)
Cell_4SigMaleupDF <- as.data.frame(Cell_4SigMaleup)
Cell_4SigFemdnDF <- as.data.frame(Cell_4SigFemdn)
Cell_4SigMalednDF <- as.data.frame(Cell_4SigMaledn)
MorulaSigFemupDF <- as.data.frame(MorulaSigFemup)
MorulaSigMaleupDF <- as.data.frame(MorulaSigMaleup)
MorulaSigFemdnDF <- as.data.frame(MorulaSigFemdn)
MorulaSigMalednDF <- as.data.frame(MorulaSigMaledn)

#now make the geneIDs which is a row name a column in the df so we can match with the gene name
Cell_4SigFemupDF$geneID <- rownames(Cell_4SigFemupDF)
Cell_4SigMaleupDF$geneID <- rownames(Cell_4SigMaleupDF)
Cell_4SigFemdnDF$geneID <- rownames(Cell_4SigFemdnDF)
Cell_4SigMalednDF$geneID <- rownames(Cell_4SigMalednDF)
MorulaSigFemupDF$geneID <- rownames(MorulaSigFemupDF)
MorulaSigMaleupDF$geneID <- rownames(MorulaSigMaleupDF)
MorulaSigFemdnDF$geneID <- rownames(MorulaSigFemdnDF)
MorulaSigMalednDF$geneID <- rownames(MorulaSigMalednDF)

#trying to figure out how to get the gene name back
Cell_4SigFemupDF <- merge(Cell_4SigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMaleupDF <- merge(Cell_4SigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigFemdnDF <- merge(Cell_4SigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMalednDF <- merge(Cell_4SigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemupDF <- merge(MorulaSigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemdnDF <- merge(MorulaSigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMaleupDF <- merge(MorulaSigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMalednDF <- merge(MorulaSigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")


#try to switch column order so the gene name isn't at the end
setcolorder(Cell_4SigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

setcolorder(MorulaSigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(Cell_4SigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_up_withGene_ByPadj.csv', row.names = TRUE)
write.csv(Cell_4SigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_dn_withGene_ByPadj.csv', row.names = TRUE)
write.csv(Cell_4SigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_up_withGene_ByPadj.csv', row.names = TRUE)
write.csv(Cell_4SigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_dn_withGene_ByPadj.csv', row.names = TRUE)

write.csv(MorulaSigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up_withGene_ByPadj.csv', row.names = TRUE)
write.csv(MorulaSigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__FemOnly_Sig_res_dn_withGene_ByPadj.csv', row.names = TRUE)
write.csv(MorulaSigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up_withGene_ByPadj.csv', row.names = TRUE)
write.csv(MorulaSigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__MaleOnly_Sig_res_dn_withGene_ByPadj.csv', row.names = TRUE)

#find the overlaps between the sexes at each stages with pvalue of 0.01 cutoff
#pull out significant up and down res Changing to use padj <.05 and 1.5 FC
Cell_4SigFem <- subset(TME62_4Cell_Female_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
Cell_4SigFemup <- subset (Cell_4SigFem,log2FoldChange > 0)                                                          
Cell_4SigFemdn <- subset (Cell_4SigFem, log2FoldChange < 0)
write.csv(Cell_4SigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_up_ByLowerPval01.csv', row.names = TRUE)
write.csv(Cell_4SigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_FemOnly_Sig_res_dn_ByLowerPval01.csv', row.names = TRUE)


Cell_4SigMale <- subset(TME62_4Cell_Male_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
Cell_4SigMaleup <- subset (Cell_4SigMale,log2FoldChange > 0)                                                          
Cell_4SigMaledn <- subset (Cell_4SigMale, log2FoldChange < 0)
write.csv(Cell_4SigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_up_ByLowerPval01.csv', row.names = TRUE)
write.csv(Cell_4SigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4Cell_MaleOnly_Sig_res_dn_ByLowerPval01.csv', row.names = TRUE)

#pull out significant up and down res, 1.5 FC or greater
Morula_SigMale <- subset(TME62_Morula_MAle_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
MorulaSigMaleup <- subset (Morula_SigMale,log2FoldChange > 0)                                                          
MorulaSigMaledn <- subset (Morula_SigMale, log2FoldChange < 0)
write.csv(MorulaSigMaleup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up_ByLowerPval01.csv', row.names = TRUE)
write.csv(MorulaSigMaledn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_dn_ByLowerPval01.csv', row.names = TRUE)

Morula_SigFemale <- subset(TME62_Morula_Female_results, pvalue <0.01 & abs(log2FoldChange) > 0.58)
MorulaSigFemup <- subset (Morula_SigFemale,log2FoldChange > 0)                                                          
MorulaSigFemdn <- subset (Morula_SigFemale, log2FoldChange < 0)
write.csv(MorulaSigFemup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up_ByLowerPval01.csv', row.names = TRUE)
write.csv(MorulaSigFemdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_dn_ByLowerPval01.csv', row.names = TRUE)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
Cell_4SigFemupDF <- as.data.frame(Cell_4SigFemup)
Cell_4SigMaleupDF <- as.data.frame(Cell_4SigMaleup)
Cell_4SigFemdnDF <- as.data.frame(Cell_4SigFemdn)
Cell_4SigMalednDF <- as.data.frame(Cell_4SigMaledn)
MorulaSigFemupDF <- as.data.frame(MorulaSigFemup)
MorulaSigMaleupDF <- as.data.frame(MorulaSigMaleup)
MorulaSigFemdnDF <- as.data.frame(MorulaSigFemdn)
MorulaSigMalednDF <- as.data.frame(MorulaSigMaledn)

#now make the geneIDs which is a row name a column in the df so we can match with the gene name
Cell_4SigFemupDF$geneID <- rownames(Cell_4SigFemupDF)
Cell_4SigMaleupDF$geneID <- rownames(Cell_4SigMaleupDF)
Cell_4SigFemdnDF$geneID <- rownames(Cell_4SigFemdnDF)
Cell_4SigMalednDF$geneID <- rownames(Cell_4SigMalednDF)
MorulaSigFemupDF$geneID <- rownames(MorulaSigFemupDF)
MorulaSigMaleupDF$geneID <- rownames(MorulaSigMaleupDF)
MorulaSigFemdnDF$geneID <- rownames(MorulaSigFemdnDF)
MorulaSigMalednDF$geneID <- rownames(MorulaSigMalednDF)

#trying to figure out how to get the gene name back
Cell_4SigFemupDF <- merge(Cell_4SigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMaleupDF <- merge(Cell_4SigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigFemdnDF <- merge(Cell_4SigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell_4SigMalednDF <- merge(Cell_4SigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemupDF <- merge(MorulaSigFemupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigFemdnDF <- merge(MorulaSigFemdnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMaleupDF <- merge(MorulaSigMaleupDF, genetotranscript, by.x = "geneID", by.y = "transcript")
MorulaSigMalednDF <- merge(MorulaSigMalednDF, genetotranscript, by.x = "geneID", by.y = "transcript")


#try to switch column order so the gene name isn't at the end
setcolorder(Cell_4SigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(Cell_4SigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

setcolorder(MorulaSigFemupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigFemdnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMaleupDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(MorulaSigMalednDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(Cell_4SigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_up_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(Cell_4SigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_FemOnly_Sig_res_dn_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(Cell_4SigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_up_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(Cell_4SigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/4_Cell_MaleOnly_Sig_res_dn_withGene_ByLowerPval01.csv', row.names = TRUE)

write.csv(MorulaSigFemupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_FemOnly_Sig_res_up_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(MorulaSigFemdnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__FemOnly_Sig_res_dn_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(MorulaSigMaleupDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula_MaleOnly_Sig_res_up_withGene_ByLowerPval01.csv', row.names = TRUE)
write.csv(MorulaSigMalednDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/Morula__MaleOnly_Sig_res_dn_withGene_ByLowerPval01.csv', row.names = TRUE)



#this code explains how to make the volcano plots
#load DESeq library
library(DESeq2)
#read csv file into R [this has raw count data for all replicates and gene name in first column]
TME33_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/TME33_BulkRNA_genes_expression_expected_count.tsv', header = TRUE, sep ="\t")
#preserve gene to transcript match
genetotranscript <- TME33_unfiltered[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(TME33_unfiltered) <- TME33_unfiltered[,2]
#now remove column 1 (genes) from data frame
#change all values to integers - do this for each samples
TME33_unfiltered <- TME33_unfiltered[,-1]
#now remove the additional transcfript row
TME33_unfiltered <- TME33_unfiltered[,-1]
#name each column of your data as replicates i.e. con1, con2 etc.
names(TME33_unfiltered) <- c("cap1","cap2","cap3","cap4","cap5","cap6","cap7","caud1","caud2","caud3","caud4","caud5","caud6","caud7","cor1","cor2","cor3","cor4","cor5","cor6","cor7","test1","test2","test3","test4","test5","test6","test7")
#change all values to integers - do this for each samples
TME33_unfiltered$cap1 <- as.integer(TME33_unfiltered$cap1)
TME33_unfiltered$cap2 <- as.integer(TME33_unfiltered$cap2)
TME33_unfiltered$cap3 <- as.integer(TME33_unfiltered$cap3)
TME33_unfiltered$cap4 <- as.integer(TME33_unfiltered$cap4)
TME33_unfiltered$cap5 <- as.integer(TME33_unfiltered$cap5)
TME33_unfiltered$cap6 <- as.integer(TME33_unfiltered$cap6)
TME33_unfiltered$cap7 <- as.integer(TME33_unfiltered$cap7)
TME33_unfiltered$caud1 <- as.integer(TME33_unfiltered$caud1)
TME33_unfiltered$caud2 <- as.integer(TME33_unfiltered$caud2)
TME33_unfiltered$caud3 <- as.integer(TME33_unfiltered$caud3)
TME33_unfiltered$caud4 <- as.integer(TME33_unfiltered$caud4)
TME33_unfiltered$caud5 <- as.integer(TME33_unfiltered$caud5)
TME33_unfiltered$caud6 <- as.integer(TME33_unfiltered$caud6)
TME33_unfiltered$caud7 <- as.integer(TME33_unfiltered$caud7)
TME33_unfiltered$cor1 <- as.integer(TME33_unfiltered$cor1)
TME33_unfiltered$cor2 <- as.integer(TME33_unfiltered$cor2)
TME33_unfiltered$cor3 <- as.integer(TME33_unfiltered$cor3)
TME33_unfiltered$cor4 <- as.integer(TME33_unfiltered$cor4)
TME33_unfiltered$cor5 <- as.integer(TME33_unfiltered$cor5)
TME33_unfiltered$cor6 <- as.integer(TME33_unfiltered$cor6)
TME33_unfiltered$cor7 <- as.integer(TME33_unfiltered$cor7)
TME33_unfiltered$test1 <- as.integer(TME33_unfiltered$test1)
TME33_unfiltered$test2 <- as.integer(TME33_unfiltered$test2)
TME33_unfiltered$test3 <- as.integer(TME33_unfiltered$test3)
TME33_unfiltered$test4 <- as.integer(TME33_unfiltered$test4)
TME33_unfiltered$test5 <- as.integer(TME33_unfiltered$test5)
TME33_unfiltered$test6 <- as.integer(TME33_unfiltered$test6)
TME33_unfiltered$test7 <- as.integer(TME33_unfiltered$test7)

#filter data 
#use a cut off of 5 tpm in at least half the samples in each group
poly <- 2
cont <- 1.5
#keep <- rowSums(counts(TME33_unfiltered[,1:4])>=5)>= poly | rowSums(counts(TME33_unfiltered[,5:7])>=5)>= cont | rowSums(counts(TME33_unfiltered[,8:11])>=5)>= poly| rowSums(counts(TME33_unfiltered[,12:14])>=5)>= cont| rowSums(counts(TME33_unfiltered[,15:18])>=5)>= poly | rowSums(counts(TME33_unfiltered[,19:21])>=5)>= cont | rowSums(counts(TME33_unfiltered[,22:25])>=5)>= poly | rowSums(counts(TME33_unfiltered[,26:28])>=5)>= cont
TME33_filtered <- TME33_unfiltered
#split things up by tissue, will process each tissue individually since I'm not comparing across
TME33_filt_cap <- TME33_filtered[,1:7]
TME33_filt_caud <- TME33_filtered[,8:14]
TME33_filt_cor <- TME33_filtered[,15:21]
TME33_filt_test <- TME33_filtered[,22:28]
#create a frame with condition names, this can be used for each tissue since they are in the same order
tissue_condition <- c("poly","poly","poly","poly","cont","cont","cont")
#create a data frame for use with DESeq2, make this for each tissue type
colDataCap <- data.frame(row.names =colnames(TME33_filt_cap), tissue_condition)
View(colDataCap)
colDataCaud <- data.frame(row.names =colnames(TME33_filt_caud), tissue_condition)
colDataCor <- data.frame(row.names =colnames(TME33_filt_cor), tissue_condition)
colDataTest <- data.frame(row.names =colnames(TME33_filt_test), tissue_condition)
#run DESeq2, do this for each tissue type
Cap_dds <- DESeqDataSetFromMatrix(countData = TME33_filt_cap, colData = colDataCap, design = ~tissue_condition)
Caud_dds <- DESeqDataSetFromMatrix(countData = TME33_filt_caud, colData = colDataCaud, design = ~tissue_condition)
Cor_dds <- DESeqDataSetFromMatrix(countData = TME33_filt_cor, colData = colDataCor, design = ~tissue_condition)
Test_dds <- DESeqDataSetFromMatrix(countData = TME33_filt_test, colData = colDataTest, design = ~tissue_condition)

keep <- rowSums(counts(Cap_dds[,1:4])>=5)>= poly | rowSums(counts(Cap_dds[,5:7])>=5)>= cont 
Cap_dds <- Cap_dds[keep,]
keep <- rowSums(counts(Caud_dds[,1:4])>=5)>= poly | rowSums(counts(Caud_dds[,5:7])>=5)>= cont 
Caud_dds <- Caud_dds[keep,]
keep <- rowSums(counts(Cor_dds[,1:4])>=5)>= poly | rowSums(counts(Cor_dds[,5:7])>=5)>= cont 
Cor_dds <- Cor_dds[keep,]
keep <- rowSums(counts(Test_dds[,1:4])>=5)>= poly | rowSums(counts(Test_dds[,5:7])>=5)>= cont 
Test_dds <- Test_dds[keep,]

#run DESeq2, do this for each tissue
Cap_dds <- DESeq(Cap_dds)
Caud_dds <- DESeq(Caud_dds)
Cor_dds <- DESeq(Cor_dds)
Test_dds <- DESeq(Test_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each tissue
Cap_results <- results(Cap_dds, contrast = c("tissue_condition", "poly", "cont"))
Caud_results <- results(Caud_dds, contrast = c("tissue_condition", "poly", "cont"))
Cor_results <- results(Cor_dds, contrast = c("tissue_condition", "poly", "cont"))
Test_results <- results(Test_dds, contrast = c("tissue_condition", "poly", "cont"))
#run a summary on each set of samples
summary(Cap_results)
summary(Caud_results)
summary(Cor_results)
summary(Test_results)
#write res to CSV do for each tissue
write.csv(Cap_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cap_res_DESeq2.csv', row.names = TRUE)
write.csv(Caud_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Caud_res_DESeq2.csv', row.names = TRUE)
write.csv(Cor_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cor_res_DESeq2.csv', row.names = TRUE)
write.csv(Test_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Test_res_DESeq2.csv', row.names = TRUE)

#pull out significant up and down res, 1.5 FC or greater
CaudResSig <- subset(Caud_results, padj <0.05 & abs(log2FoldChange) > 0.58)
CaudSigup <- subset (CaudResSig,log2FoldChange > 0)                                                          
CaudSigdn <- subset (CaudResSig, log2FoldChange < 0)
write.csv(CaudSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Caud_Sig_res_up.csv', row.names = TRUE)
write.csv(CaudSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Caud_Sig_res_dn.csv', row.names = TRUE)

CapResSig <- subset(Cap_results, padj <0.05 & abs(log2FoldChange) > 0.58)
CapSigup <- subset (CapResSig,log2FoldChange > 0)                                                          
CapSigdn <- subset (CapResSig, log2FoldChange < 0)
write.csv(CapSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cap_Sig_res_up.csv', row.names = TRUE)
write.csv(CapSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cap_Sig_res_dn.csv', row.names = TRUE)

CorResSig <- subset(Cor_results, padj <0.05 & abs(log2FoldChange) > 0.58)
CorSigup <- subset (CorResSig,log2FoldChange > 0)                                                          
CorSigdn <- subset (CorResSig, log2FoldChange < 0)
write.csv(CorSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cor_Sig_res_up.csv', row.names = TRUE)
write.csv(CorSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cor_Sig_res_dn.csv', row.names = TRUE)

TestResSig <- subset(Test_results, padj <0.05 & abs(log2FoldChange) > 0.58)
TestSigup <- subset (TestResSig,log2FoldChange > 0)                                                          
TestSigdn <- subset (TestResSig, log2FoldChange < 0)
write.csv(TestSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Test_Sig_res_up.csv', row.names = TRUE)
write.csv(TestSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Test_Sig_res_dn.csv', row.names = TRUE)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
CapSigUpDF <- as.data.frame(CapSigup)
CapSigDnDF <- as.data.frame(CapSigdn)
CorSigUpDF <- as.data.frame(CorSigup)
CorSigDnDF <- as.data.frame(CorSigdn)
CaudSigUpDF <- as.data.frame(CaudSigup)
CaudSigDnDF <- as.data.frame(CaudSigdn)
TestSigUpDF <- as.data.frame(TestSigup)
TestSigDnDF <- as.data.frame(TestSigdn)
#now make the geneIDs which is a row name a column in the df so we can match with the gene name
CapSigUpDF$geneID <- rownames(CapSigUpDF)
CapSigDnDF$geneID <- rownames(CapSigDnDF)
CorSigUpDF$geneID <- rownames(CorSigUpDF)
CorSigDnDF$geneID <- rownames(CorSigDnDF)
CaudSigUpDF$geneID <- rownames(CaudSigUpDF)
CaudSigDnDF$geneID <- rownames(CaudSigDnDF)
TestSigUpDF$geneID <- rownames(TestSigUpDF)
TestSigDnDF$geneID <- rownames(TestSigDnDF)
#trying to figure out how to get the gene name back
CapSigUpDF <- merge(CapSigUpDF, genetotranscript, by.x = "geneID", by.y = "transcript")
CapSigDnDF <- merge(CapSigDnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
CorSigUpDF <- merge(CorSigUpDF, genetotranscript, by.x = "geneID", by.y = "transcript")
CorSigDnDF <- merge(CorSigDnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
CaudSigUpDF <- merge(CaudSigUpDF, genetotranscript, by.x = "geneID", by.y = "transcript")
CaudSigDnDF <- merge(CaudSigDnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
TestSigUpDF <- merge(TestSigUpDF, genetotranscript, by.x = "geneID", by.y = "transcript")
TestSigDnDF <- merge(TestSigDnDF, genetotranscript, by.x = "geneID", by.y = "transcript")
#try to switch column order so the gene name isn't at the end
library(data.table)
setcolorder(CapSigUpDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(CapSigDnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(CorSigUpDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(CorSigDnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(CaudSigUpDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(CaudSigDnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(TestSigUpDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(TestSigDnDF,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(CaudSigUpDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Caud_Sig_res_up_withGene.csv', row.names = TRUE)
write.csv(CaudSigDnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Caud_Sig_res_dn_withGene.csv', row.names = TRUE)
write.csv(CapSigUpDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cap_Sig_res_up_withGene.csv', row.names = TRUE)
write.csv(CapSigDnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cap_Sig_res_dn_withGene.csv', row.names = TRUE)
write.csv(CorSigUpDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cor_Sig_res_up_withGene.csv', row.names = TRUE)
write.csv(CorSigDnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Cor_Sig_res_dn_withGene.csv', row.names = TRUE)
write.csv(TestSigUpDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Test_Sig_res_up_withGene.csv', row.names = TRUE)
write.csv(TestSigDnDF, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/Test_Sig_res_dn_withGene.csv', row.names = TRUE)

#code used for embryo analysis in figure 3 and 4
TME62_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/TME62_Filt_order_Gene_expression_expected.csv', header = TRUE, sep =",")
#preserve gene to transcript match
genetotranscript <- TME62_unfiltered[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(TME62_unfiltered) <- TME62_unfiltered[,2]
#now remove column 1 (genes) from data frame
TME62_unfiltered <- TME62_unfiltered[,-1]
#now remove the additional transcfript row
TME62_unfiltered <- TME62_unfiltered[,-1]
#name each column of your data as replicates i.e. con1, con2 etc.
names(TME62_unfiltered) <- c(colnames(TME62_unfiltered))
# I think this step was unneccessary and they wer actually already named
#change all values to integers - do this for each samples
TME62_unfiltered$TME51Cont_4Cell_1<- as.integer(TME62_unfiltered$TME51Cont_4Cell_1)
TME62_unfiltered$TME51Cont_4Cell_2<- as.integer(TME62_unfiltered$TME51Cont_4Cell_2)
TME62_unfiltered$TME51Cont_4Cell_4<- as.integer(TME62_unfiltered$TME51Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_2<- as.integer(TME62_unfiltered$TME52Cont_4Cell_2)
TME62_unfiltered$TME52Cont_4Cell_4<- as.integer(TME62_unfiltered$TME52Cont_4Cell_4)
TME62_unfiltered$TME52Cont_4Cell_5<- as.integer(TME62_unfiltered$TME52Cont_4Cell_5)
TME62_unfiltered$TME52Cont_4Cell_6<- as.integer(TME62_unfiltered$TME52Cont_4Cell_6)
TME62_unfiltered$TME52Cont_4Cell_8<- as.integer(TME62_unfiltered$TME52Cont_4Cell_8)
TME62_unfiltered$TME55Cont_4Cell_2<- as.integer(TME62_unfiltered$TME55Cont_4Cell_2)
TME62_unfiltered$TME55Cont_4Cell_3<- as.integer(TME62_unfiltered$TME55Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_1<- as.integer(TME62_unfiltered$TME59Cont_4Cell_1)
TME62_unfiltered$TME59Cont_4Cell_2<- as.integer(TME62_unfiltered$TME59Cont_4Cell_2)
TME62_unfiltered$TME59Cont_4Cell_3<- as.integer(TME62_unfiltered$TME59Cont_4Cell_3)
TME62_unfiltered$TME59Cont_4Cell_4<- as.integer(TME62_unfiltered$TME59Cont_4Cell_4)
TME62_unfiltered$TME59Cont_4Cell_5<- as.integer(TME62_unfiltered$TME59Cont_4Cell_5)
TME62_unfiltered$TME59Cont_4Cell_6<- as.integer(TME62_unfiltered$TME59Cont_4Cell_6)
TME62_unfiltered$TME59Cont_4Cell_7<- as.integer(TME62_unfiltered$TME59Cont_4Cell_7)
TME62_unfiltered$TME59Cont_4Cell_8<- as.integer(TME62_unfiltered$TME59Cont_4Cell_8)
TME62_unfiltered$TME51Poly_4Cell_1<- as.integer(TME62_unfiltered$TME51Poly_4Cell_1)
TME62_unfiltered$TME51Poly_4Cell_2<- as.integer(TME62_unfiltered$TME51Poly_4Cell_2)
TME62_unfiltered$TME51Poly_4Cell_3<- as.integer(TME62_unfiltered$TME51Poly_4Cell_3)
TME62_unfiltered$TME51Poly_4Cell_4<- as.integer(TME62_unfiltered$TME51Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_1<- as.integer(TME62_unfiltered$TME52Poly_4Cell_1)
TME62_unfiltered$TME52Poly_4Cell_2<- as.integer(TME62_unfiltered$TME52Poly_4Cell_2)
TME62_unfiltered$TME52Poly_4Cell_3<- as.integer(TME62_unfiltered$TME52Poly_4Cell_3)
TME62_unfiltered$TME52Poly_4Cell_4<- as.integer(TME62_unfiltered$TME52Poly_4Cell_4)
TME62_unfiltered$TME52Poly_4Cell_5<- as.integer(TME62_unfiltered$TME52Poly_4Cell_5)
TME62_unfiltered$TME52Poly_4Cell_6<- as.integer(TME62_unfiltered$TME52Poly_4Cell_6)
TME62_unfiltered$TME52Poly_4Cell_7<- as.integer(TME62_unfiltered$TME52Poly_4Cell_7)
TME62_unfiltered$TME52Poly_4Cell_8<- as.integer(TME62_unfiltered$TME52Poly_4Cell_8)
TME62_unfiltered$TME55Poly_4Cell_1<- as.integer(TME62_unfiltered$TME55Poly_4Cell_1)
TME62_unfiltered$TME55Poly_4Cell_2<- as.integer(TME62_unfiltered$TME55Poly_4Cell_2)
TME62_unfiltered$TME55Poly_4Cell_3<- as.integer(TME62_unfiltered$TME55Poly_4Cell_3)
TME62_unfiltered$TME55Poly_4Cell_4<- as.integer(TME62_unfiltered$TME55Poly_4Cell_4)
TME62_unfiltered$TME55Poly_4Cell_5<- as.integer(TME62_unfiltered$TME55Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_1<- as.integer(TME62_unfiltered$TME59Poly_4Cell_1)
TME62_unfiltered$TME59Poly_4Cell_2<- as.integer(TME62_unfiltered$TME59Poly_4Cell_2)
TME62_unfiltered$TME59Poly_4Cell_3<- as.integer(TME62_unfiltered$TME59Poly_4Cell_3)
TME62_unfiltered$TME59Poly_4Cell_5<- as.integer(TME62_unfiltered$TME59Poly_4Cell_5)
TME62_unfiltered$TME59Poly_4Cell_6<- as.integer(TME62_unfiltered$TME59Poly_4Cell_6)
TME62_unfiltered$TME59Poly_4Cell_7<- as.integer(TME62_unfiltered$TME59Poly_4Cell_7)
TME62_unfiltered$TME59Poly_4Cell_8<- as.integer(TME62_unfiltered$TME59Poly_4Cell_8)
TME62_unfiltered$TME52Cont_Morula_1<- as.integer(TME62_unfiltered$TME52Cont_Morula_1)
TME62_unfiltered$TME52Cont_Morula_2<- as.integer(TME62_unfiltered$TME52Cont_Morula_2)
TME62_unfiltered$TME52Cont_Morula_3<- as.integer(TME62_unfiltered$TME52Cont_Morula_3)
TME62_unfiltered$TME52Cont_Morula_4<- as.integer(TME62_unfiltered$TME52Cont_Morula_4)
TME62_unfiltered$TME52Cont_Morula_5<- as.integer(TME62_unfiltered$TME52Cont_Morula_5)
TME62_unfiltered$TME52Cont_Morula_6<- as.integer(TME62_unfiltered$TME52Cont_Morula_6)
TME62_unfiltered$TME52Cont_Morula_8<- as.integer(TME62_unfiltered$TME52Cont_Morula_8)
TME62_unfiltered$TME55Cont_Morula_1<- as.integer(TME62_unfiltered$TME55Cont_Morula_1)
TME62_unfiltered$TME55Cont_Morula_2<- as.integer(TME62_unfiltered$TME55Cont_Morula_2)
TME62_unfiltered$TME55Cont_Morula_3<- as.integer(TME62_unfiltered$TME55Cont_Morula_3)
TME62_unfiltered$TME55Cont_Morula_4<- as.integer(TME62_unfiltered$TME55Cont_Morula_4)
TME62_unfiltered$TME55Cont_Morula_5<- as.integer(TME62_unfiltered$TME55Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_1<- as.integer(TME62_unfiltered$TME59Cont_Morula_1)
TME62_unfiltered$TME59Cont_Morula_2<- as.integer(TME62_unfiltered$TME59Cont_Morula_2)
TME62_unfiltered$TME59Cont_Morula_3<- as.integer(TME62_unfiltered$TME59Cont_Morula_3)
TME62_unfiltered$TME59Cont_Morula_4<- as.integer(TME62_unfiltered$TME59Cont_Morula_4)
TME62_unfiltered$TME59Cont_Morula_5<- as.integer(TME62_unfiltered$TME59Cont_Morula_5)
TME62_unfiltered$TME59Cont_Morula_6<- as.integer(TME62_unfiltered$TME59Cont_Morula_6)
TME62_unfiltered$TME59Cont_Morula_7<- as.integer(TME62_unfiltered$TME59Cont_Morula_7)
TME62_unfiltered$TME59Cont_Morula_8<- as.integer(TME62_unfiltered$TME59Cont_Morula_8)
TME62_unfiltered$TME52Poly_Morula_1<- as.integer(TME62_unfiltered$TME52Poly_Morula_1)
TME62_unfiltered$TME52Poly_Morula_2<- as.integer(TME62_unfiltered$TME52Poly_Morula_2)
TME62_unfiltered$TME52Poly_Morula_3<- as.integer(TME62_unfiltered$TME52Poly_Morula_3)
TME62_unfiltered$TME52Poly_Morula_4<- as.integer(TME62_unfiltered$TME52Poly_Morula_4)
TME62_unfiltered$TME52Poly_Morula_5<- as.integer(TME62_unfiltered$TME52Poly_Morula_5)
TME62_unfiltered$TME52Poly_Morula_6<- as.integer(TME62_unfiltered$TME52Poly_Morula_6)
TME62_unfiltered$TME52Poly_Morula_8<- as.integer(TME62_unfiltered$TME52Poly_Morula_8)
TME62_unfiltered$TME55Poly_Morula_1<- as.integer(TME62_unfiltered$TME55Poly_Morula_1)
TME62_unfiltered$TME55Poly_Morula_2<- as.integer(TME62_unfiltered$TME55Poly_Morula_2)
TME62_unfiltered$TME55Poly_Morula_3<- as.integer(TME62_unfiltered$TME55Poly_Morula_3)
TME62_unfiltered$TME55Poly_Morula_4<- as.integer(TME62_unfiltered$TME55Poly_Morula_4)
TME62_unfiltered$TME55Poly_Morula_5<- as.integer(TME62_unfiltered$TME55Poly_Morula_5)
TME62_unfiltered$TME55Poly_Morula_6<- as.integer(TME62_unfiltered$TME55Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_3<- as.integer(TME62_unfiltered$TME59Poly_Morula_3)
TME62_unfiltered$TME59Poly_Morula_4<- as.integer(TME62_unfiltered$TME59Poly_Morula_4)
TME62_unfiltered$TME59Poly_Morula_5<- as.integer(TME62_unfiltered$TME59Poly_Morula_5)
TME62_unfiltered$TME59Poly_Morula_6<- as.integer(TME62_unfiltered$TME59Poly_Morula_6)
TME62_unfiltered$TME59Poly_Morula_8<- as.integer(TME62_unfiltered$TME59Poly_Morula_8)

#filter data - here I use a cutoff of 81 across all replicates, this could be higher
TME62_filtered <- TME62_unfiltered
#split things up by dev stage, will process each stage individually since I'm not comparing across

#split things up by dev stage , will process each stage individually since I'm not comparing across
TME62_filt_4Cell <- TME62_filtered[,1:42]
TME62_filt_Morula <- TME62_filtered[,43:80]
#based on xist gene and Eif2s3y gene in EmbryoNewOrderFiltered pull out each based on sex
#filter into diff sex
#create a frame with condition names, this can be used for each tissue since they are in the same order
condition_4Cell <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly", "cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly")

#create a data frame for use with DESeq2, make this for each tissue type
colData4Cell <- data.frame(row.names =colnames(TME62_filt_4Cell), condition_4Cell)
View(colData4Cell)

colDataMorula <- data.frame(row.names =colnames(TME62_filt_Morula), condition_Morula)
View(colDataMorula)



#run DESeq2, do this for each condition and sex
library(DESeq2)
#run DESeq2, do this for each condition and sex
TME62_4Cell_dds <- DESeqDataSetFromMatrix(countData = TME62_filt_4Cell, colData = colData4Cell, design = ~condition_4Cell)

TME62_Morula_dds <-  DESeqDataSetFromMatrix(countData = TME62_filt_Morula, colData = colDataMorula, design = ~condition_Morula)

#for each condition in each group, we need to figure out the number
#so number of female control 4 cell, female poly 4 cell, female cont morula, etc
colnames(TME62_4Cell_dds)
colnames(TME62_Morula_dds)
#used this to determine the cut offs (number in each group divided by 2)
cell4_min_Poly <- 12
cell4_min_Cont <- 9

morula_min_Cont <- 9
morula_min_Poly <- 10


#now move on to determining which rows we will keep because the counts are great enough
#determining if the at least half of the samples in either condition have counts for the gene, 5 tpm cut off
keep <- rowSums(counts(TME62_4Cell_dds[,c(1:11, 19:25)])>=5)>= cell4_min_Poly | rowSums(counts(TME62_4Cell_dds[,c(12:18,26:42)])>=5)>= cell4_min_Cont
TME62_4Cell_dds <- TME62_4Cell_dds[keep,]
keep <- rowSums(counts(TME62_Morula_dds[,c(1:13,24:30)])>=5)>= morula_min_Cont | rowSums(counts(TME62_Morula_dds[,c(14:23,31:38)])>=5)>= morula_min_Poly
TME62_Morula_dds <- TME62_Morula_dds[keep,]

#now can continue with same code as before
#run DESeq2, do this for each sample set
TME62_4Cell_dds <- DESeq(TME62_4Cell_dds)
TME62_Morula_dds <- DESeq(TME62_Morula_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each tissue
TME62_4Cell_results <- results(TME62_4Cell_dds, contrast = c("condition_4Cell", "poly", "cont"))

TME62_Morula_results <- results(TME62_Morula_dds, contrast = c("condition_Morula", "poly", "cont"))

Cell4X <- c("BC061195","Eda","Gla","Gripap1","Rbm10","Siah1b")
Cell4Y <- c(rep("",6))
Cell4X <- as.data.frame(Cell4X)
Cell4Y <- as.data.frame(Cell4Y)
library(ggplot2)
# making new volcano plots that are zoomed in, and color dots
mycolors <- c("NO" = "black","UP" = "firebrick3","DOWN" = "deepskyblue3","X Chromosome"= "darkolivegreen4","Y Chromosome"="green")
#new volcano plot that have titles, same scale, and use pvalue .01 for cut off
Cell4dfTest <- as.data.frame(TME62_4Cell_results)
Cell4dfTest$geneID <- rownames(Cell4dfTest)
Cell4dfTest <- merge(Cell4dfTest, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell4dfTest$Expression <- "NO"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange > 0.6 & Cell4dfTest$pvalue < 0.01] <- "UP"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange < -0.6 & Cell4dfTest$pvalue < 0.01] <- "DOWN"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange < -0.6 & Cell4dfTest$pvalue < 0.01 & Cell4dfTest$gene %in% Cell4X$Cell4X] <- "X Chromosome"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange > 0.6 & Cell4dfTest$pvalue < 0.01 & Cell4dfTest$gene %in% Cell4X$Cell4X] <- "X Chromosome"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange < -0.6 & Cell4dfTest$pvalue < 0.01 & Cell4dfTest$gene %in% Cell4Y$Cell4Y] <- "Y Chromosome"
Cell4dfTest$Expression[Cell4dfTest$log2FoldChange > 0.6 & Cell4dfTest$pvalue < 0.01 & Cell4dfTest$gene %in% Cell4Y$Cell4Y] <- "Y Chromosome"

p <- ggplot(data = Cell4dfTest, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-8,8)+
  ylim(0,7)+
  ggtitle("4 Cell Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("4CellVolcano.svg", width=5, height=5)
p2
dev.off()

MorulaX <- c("A830080D01Rik","Magea2","Rbmx","Tfe3","Tmem29","Wdr44","Xlr4c")
MorulaX <- as.data.frame(MorulaX)
Morula_results <- as.data.frame(TME62_Morula_results)
Morula_results$geneID <- rownames(Morula_results)
Morula_results <- merge(Morula_results, genetotranscript, by.x = "geneID", by.y = "transcript")
Morula_results$Expression <- "NO"
Morula_results$Expression[Morula_results$log2FoldChange > 0.6 & Morula_results$pvalue < 0.01] <- "UP"
Morula_results$Expression[Morula_results$log2FoldChange < -0.6 & Morula_results$pvalue < 0.01] <- "DOWN"
Morula_results$Expression[Morula_results$log2FoldChange < -0.6 & Morula_results$pvalue < 0.01 & Morula_results$gene %in% MorulaX$MorulaX] <- "X Chromosome"
Morula_results$Expression[Morula_results$log2FoldChange > 0.6 & Morula_results$pvalue < 0.01 & Morula_results$gene %in% MorulaX$MorulaX] <- "X Chromosome"

p <- ggplot(data = Morula_results, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-8,8)+
  ylim(0,7)+
  ggtitle("Morula Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("MorulaVolcano.svg", width=5, height=5)
p2
dev.off()

#by embryo sex
TME62_filtered <- TME62_unfiltered
#split things up by dev stage, will process each stage individually since I'm not comparing across

#split things up by dev stage , will process each stage individually since I'm not comparing across
TME62_filt_4Cell <- TME62_filtered[,1:42]
TME62_filt_Morula <- TME62_filtered[,43:80]
#based on xist gene and Eif2s3y gene in EmbryoNewOrderFiltered pull out each based on sex
#filter into diff sex
TME62_Female_4Cell <- TME62_filt_4Cell[,1:18]
TME62_Male_4Cell <- TME62_filt_4Cell[,19:42]
TME62_Female_Morula <- TME62_filt_Morula[,1:23]
TME62_Male_Morula <- TME62_filt_Morula[,24:38]

#create a frame with condition names, this can be used for each tissue since they are in the same order
condition_4Cell_Female <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly")
condition_4Cell_Male <- c("cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula_Female <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
condition_Morula_Male <- c("cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly")

#condition_morula <- c("cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","cont","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly","poly")
#create a data frame for use with DESeq2, make this for each tissue type
colData4CellFemale <- data.frame(row.names =colnames(TME62_Female_4Cell), condition_4Cell_Female)
View(colData4CellFemale)
colData4CellMale <- data.frame(row.names =colnames(TME62_Male_4Cell), condition_4Cell_Male)
View(colData4CellFemale)
#try to add to dataframe and separate by sex, see embryoNewOrder sheet for sex determinataion
#Cell4_sex <- c("F","M","F","M","F","F","F","F","F","F","F","M","M","M","F","M","F","M","M","F","M","M","M","M","M","M","F","M","M","M","F","M","M","F","M","M","F","M","M","F","F","M")
#MoreColData4Cell <- data.frame(row.names = colnames(TME62_filt_4Cell), condition_4Cell, Cell4_sex)
#can use the above dataframe to subset based on sex solumn
#Male_4Cell <- subset(MoreColData4Cell, Cell4_sex == 'M')

colDataMorulaFemale <- data.frame(row.names =colnames(TME62_Female_Morula), condition_Morula_Female)
View(colDataMorulaFemale)
colDataMorulaMale <- data.frame(row.names =colnames(TME62_Male_Morula), condition_Morula_Male)
View(colDataMorulaMale)

#run DESeq2, do this for each condition and sex
TME62_4Cell_Female_dds <- DESeqDataSetFromMatrix(countData = TME62_Female_4Cell, colData = colData4CellFemale, design = ~condition_4Cell_Female)
TME62_4Cell_Male_dds <- DESeqDataSetFromMatrix(countData = TME62_Male_4Cell, colData = colData4CellMale, design = ~condition_4Cell_Male)

TME62_Morula_Female_dds <-  DESeqDataSetFromMatrix(countData = TME62_Female_Morula, colData = colDataMorulaFemale, design = ~condition_Morula_Female)
TME62_Morula_Male_dds <-  DESeqDataSetFromMatrix(countData = TME62_Male_Morula, colData = colDataMorulaMale, design = ~condition_Morula_Male)

#for each condition in each group, we need to figure out the number
#so number of female control 4 cell, female poly 4 cell, female cont morula, etc
colnames(TME62_4Cell_Female_dds)
#used this to determine the cut offs (number in each group divided by 2)
cell4_Female_min_Poly <- 3.5
cell4_Female_min_Cont <- 5.5
colnames(TME62_4Cell_Male_dds)
cell4_Male_min_Cont <- 3.5
cell4_Male_min_Poly <- 8.5
colnames(TME62_Morula_Female_dds)
morula_Female_min_Cont <- 6.5
morula_Female_min_Poly <- 5
colnames(TME62_Morula_Male_dds)
morula_Male_min_Cont <- 3.5
morula_Male_min_Poly <- 4

#now move on to determining which rows we will keep because the counts are great enough
#determining if the at least half of the samples in either condition have counts for the gene, 5 tpm cut off
keep <- rowSums(counts(TME62_4Cell_Female_dds[,12:18])>=5)>= cell4_Female_min_Poly | rowSums(counts(TME62_4Cell_Female_dds[,1:11])>=5)>= cell4_Female_min_Cont
TME62_4Cell_Female_dds <- TME62_4Cell_Female_dds[keep,]
keep <- rowSums(counts(TME62_4Cell_Male_dds[,1:7])>=5)>= cell4_Male_min_Cont | rowSums(counts(TME62_4Cell_Male_dds[,8:24])>=5)>= cell4_Male_min_Poly
TME62_4Cell_Male_dds <- TME62_4Cell_Male_dds[keep,]
keep <- rowSums(counts(TME62_Morula_Female_dds[,1:13])>=5)>= morula_Female_min_Cont | rowSums(counts(TME62_Morula_Female_dds[,14:23])>=5)>= morula_Female_min_Poly
TME62_Morula_Female_dds <- TME62_Morula_Female_dds[keep,]
keep <- rowSums(counts(TME62_Morula_Male_dds[,1:7])>=5)>= morula_Male_min_Cont | rowSums(counts(TME62_Morula_Male_dds[,8:15])>=5)>= morula_Male_min_Poly
TME62_Morula_Male_dds <- TME62_Morula_Male_dds[keep,]

#now can continue with same code as before
#run DESeq2, do this for each sample set
TME62_4Cell_Female_dds <- DESeq(TME62_4Cell_Female_dds)
TME62_4Cell_Male_dds <- DESeq(TME62_4Cell_Male_dds)
TME62_Morula_Female_dds <- DESeq(TME62_Morula_Female_dds)
TME62_Morula_Male_dds <- DESeq(TME62_Morula_Male_dds)
#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each tissue
TME62_4Cell_Female_results <- results(TME62_4Cell_Female_dds, contrast = c("condition_4Cell_Female", "poly", "cont"))
TME62_4Cell_Male_results <- results(TME62_4Cell_Male_dds, contrast = c("condition_4Cell_Male", "poly", "cont"))

TME62_Morula_Female_results <- results(TME62_Morula_Female_dds, contrast = c("condition_Morula_Female", "poly", "cont"))
TME62_Morula_MAle_results <- results(TME62_Morula_Male_dds, contrast = c("condition_Morula_Male", "poly", "cont"))




Cell4FemaleX <- c("Gm10439","Bex3")
Cell4FemaleX <- as.data.frame(Cell4FemaleX)
# making new volcano plots that are zoomed in, and color dots
mycolors <- c("NO" = "black","UP" = "red","DOWN" = "blue","X Chromosome"="mediumorchid1","Y Chromosome"="green")
#new volcano plot that have titles, same scale, and use pvalue .01 for cut off

#new volcano plot that have titles, same scale, and use pvalue .01 for cut off
Cell4dfFemale <- as.data.frame(TME62_4Cell_Female_results)
Cell4dfFemale$geneID <- rownames(Cell4dfFemale)
Cell4dfFemale <- merge(Cell4dfFemale, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell4dfFemale$Expression <- "NO"
Cell4dfFemale$Expression[Cell4dfFemale$log2FoldChange > 0.6 & Cell4dfFemale$pvalue < 0.01] <- "UP"
Cell4dfFemale$Expression[Cell4dfFemale$log2FoldChange < -0.6 & Cell4dfFemale$pvalue < 0.01] <- "DOWN"
Cell4dfFemale$Expression[Cell4dfFemale$log2FoldChange < -0.6 & Cell4dfFemale$pvalue < 0.01 & Cell4dfFemale$gene %in% Cell4FemaleX$Cell4FemaleX] <- "X Chromosome"
Cell4dfFemale$Expression[Cell4dfFemale$log2FoldChange > 0.6 & Cell4dfFemale$pvalue < 0.01 & Cell4dfFemale$gene %in% Cell4FemaleX$Cell4FemaleX] <- "X Chromosome"

p <- ggplot(data = Cell4dfFemale, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-10,10)+
  ylim(0,10)+
  ggtitle("Female 4Cell Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("4cellFemaleVolcano.svg", width=5, height=5)
p2
dev.off()


Cell4MaleX <- c("Tsx")
Cell4MaleX <- as.data.frame(Cell4MaleX)
Cell4dfMale <- as.data.frame(TME62_4Cell_Male_results)
Cell4dfMale$geneID <- rownames(Cell4dfMale)
Cell4dfMale <- merge(Cell4dfMale, genetotranscript, by.x = "geneID", by.y = "transcript")
Cell4dfMale$Expression <- "NO"
Cell4dfMale$Expression[Cell4dfMale$log2FoldChange > 0.6 & Cell4dfMale$pvalue < 0.01] <- "UP"
Cell4dfMale$Expression[Cell4dfMale$log2FoldChange < -0.6 & Cell4dfMale$pvalue < 0.01] <- "DOWN"
Cell4dfMale$Expression[Cell4dfMale$log2FoldChange < -0.6 & Cell4dfMale$pvalue < 0.01 & Cell4dfMale$gene %in% Cell4MaleX$Cell4MaleX] <- "X Chromosome"
Cell4dfMale$Expression[Cell4dfMale$log2FoldChange > 0.6 & Cell4dfMale$pvalue < 0.01 & Cell4dfMale$gene %in% Cell4MaleX$Cell4MaleX] <- "X Chromosome"

p <- ggplot(data = Cell4dfMale, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-10,10)+
  ylim(0,10)+
  ggtitle("Male 4Cell Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("4cellMaleVolcano.svg", width=5, height=5)
p2
dev.off()

MorulaFemaleX <- c("A830080D01Rik","Gpc4","Tsga8","Tspyl2","Wdr44")
MorulaFemaleX <- as.data.frame(MorulaFemaleX)
MoruladfFemale <- as.data.frame(TME62_Morula_Female_results)
MoruladfFemale$geneID <- rownames(MoruladfFemale)
MoruladfFemale <- merge(MoruladfFemale, genetotranscript, by.x = "geneID", by.y = "transcript")
MoruladfFemale$Expression <- "NO"
MoruladfFemale$Expression[MoruladfFemale$log2FoldChange > 0.6 & MoruladfFemale$pvalue < 0.01] <- "UP"
MoruladfFemale$Expression[MoruladfFemale$log2FoldChange < -0.6 & MoruladfFemale$pvalue < 0.01] <- "DOWN"
MoruladfFemale$Expression[MoruladfFemale$log2FoldChange < -0.6 & MoruladfFemale$pvalue < 0.01 & MoruladfFemale$gene %in% MorulaFemaleX$MorulaFemaleX] <- "X Chromosome"
MoruladfFemale$Expression[MoruladfFemale$log2FoldChange > 0.6 & MoruladfFemale$pvalue < 0.01 & MoruladfFemale$gene %in% MorulaFemaleX$MorulaFemaleX] <- "X Chromosome"

p <- ggplot(data = MoruladfFemale, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-8,8)+
  ylim(0,7)+
  ggtitle("Female Morula Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("MorulaFemaleVolcano.svg", width=5, height=5)
p2
dev.off()

MorulaMaleX <- c("Rbmx","Rpl10","Tfe3","A230072C01Rik","Chic1","Magea3","Nudt10","Tmem29","Xlr3b")
MorulaMaleX <- as.data.frame(MorulaMaleX)
MoruladfMale <- as.data.frame(TME62_Morula_MAle_results)
MoruladfMale$geneID <- rownames(MoruladfMale)
MoruladfMale <- merge(MoruladfMale, genetotranscript, by.x = "geneID", by.y = "transcript")
MoruladfMale$Expression <- "NO"
MoruladfMale$Expression [MoruladfMale$log2FoldChange > 0.6 & MoruladfMale$pvalue < 0.01] <- "UP"
MoruladfMale$Expression [MoruladfMale$log2FoldChange < -0.6 & MoruladfMale$pvalue < 0.01] <- "DOWN"
MoruladfMale$Expression [MoruladfMale$log2FoldChange < -0.6 & MoruladfMale$pvalue < 0.01 & MoruladfMale$gene %in% MorulaMaleX$MorulaMaleX] <- "X Chromosome"
MoruladfMale$Expression [MoruladfMale$log2FoldChange > 0.6 & MoruladfMale$pvalue < 0.01 & MoruladfMale$gene %in% MorulaMaleX$MorulaMaleX] <- "X Chromosome"

p <- ggplot(data = MoruladfMale, aes(x=log2FoldChange, y = -log10(pvalue),col=Expression))+geom_point() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  xlim(-8,8)+
  ylim(0,7)+
  ggtitle("Male Morula Stage Embryos")
p2 <- p + scale_colour_manual(values = mycolors)

p2
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("MorulaMaleVolcano.svg", width=5, height=5)
p2
dev.off()

#code used to do analysis for figure 6
library(DESeq2)
library(ggplot2)
library(data.table)
#will primarily follow the work flow used in TME62 on 9/10/24 to analyze this data
#first we will compare overall conditions the compare by sex
#some conditons have as little as 6 embryos in one sex so ideally more need to be added in the future
#read csv file into R [this has raw count data for all replicates and gene name in first column]
TME123_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/TME123_Expected_counts_prefilt.csv', header = TRUE, sep =",")
#preserve gene to transcript match
genetotranscript <- TME123_unfiltered[,1:2]
#label each row by first column, the will be labeled by transcript
rownames(TME123_unfiltered) <- TME123_unfiltered[,2]
#now remove column 1 (genes) from data frame
TME123_unfiltered <- TME123_unfiltered[,-1]
#now remove the additional transcfript row
TME123_unfiltered <- TME123_unfiltered[,-1]
#name each column of your data as replicates i.e. con1, con2 etc.
names(TME123_unfiltered) <- c(colnames(TME123_unfiltered))
#change all values to integers - do this for each samples
TME123_unfiltered$TME101Cont_NFWMorula1<- as.integer(TME123_unfiltered$TME101Cont_NFWMorula1)
TME123_unfiltered$TME102Cont_NFWMorula1<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula1)
TME123_unfiltered$TME102Cont_NFWMorula6<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula6)
TME123_unfiltered$TME118Cont_NFWMorula2<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula2)
TME123_unfiltered$TME118Cont_NFWMorula4<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula4)
TME123_unfiltered$TME118Cont_NFWMorula4<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula4)
TME123_unfiltered$TME118Cont_NFWMorula5<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula5)
TME123_unfiltered$TME101Cont_NFWMorula2<- as.integer(TME123_unfiltered$TME101Cont_NFWMorula2)
TME123_unfiltered$TME102Cont_NFWMorula3<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula3)
TME123_unfiltered$TME102Cont_NFWMorula2<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula2)
TME123_unfiltered$TME102Cont_NFWMorula4<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula4)
TME123_unfiltered$TME102Cont_NFWMorula5<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula5)
TME123_unfiltered$TME102Cont_NFWMorula7<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula7)
TME123_unfiltered$TME102Cont_NFWMorula8<- as.integer(TME123_unfiltered$TME102Cont_NFWMorula8)
TME123_unfiltered$TME104Cont_NFWMorula2<- as.integer(TME123_unfiltered$TME104Cont_NFWMorula2)
TME123_unfiltered$TME104Cont_NFWMorula3<- as.integer(TME123_unfiltered$TME104Cont_NFWMorula3)
TME123_unfiltered$TME118Cont_NFWMorula1<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula1)
TME123_unfiltered$TME118Cont_NFWMorula3<- as.integer(TME123_unfiltered$TME118Cont_NFWMorula3)
TME123_unfiltered$TME102Cont_TotalMorula1<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula1)
TME123_unfiltered$TME102Cont_TotalMorula3<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula3)
TME123_unfiltered$TME104Cont_TotalMorula2<- as.integer(TME123_unfiltered$TME104Cont_TotalMorula2)
TME123_unfiltered$TME112Cont_TotalMorula1<- as.integer(TME123_unfiltered$TME112Cont_TotalMorula1)
TME123_unfiltered$TME112Cont_TotalMorula2<- as.integer(TME123_unfiltered$TME112Cont_TotalMorula2)
TME123_unfiltered$TME112Cont_TotalMorula3<- as.integer(TME123_unfiltered$TME112Cont_TotalMorula3)
TME123_unfiltered$TME102Cont_TotalMorula2<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula2)
TME123_unfiltered$TME102Cont_TotalMorula4<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula4)
TME123_unfiltered$TME102Cont_TotalMorula5<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula5)
TME123_unfiltered$TME102Cont_TotalMorula6<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula6)
TME123_unfiltered$TME102Cont_TotalMorula7<- as.integer(TME123_unfiltered$TME102Cont_TotalMorula7)
TME123_unfiltered$TME104Cont_TotalMorula1<- as.integer(TME123_unfiltered$TME104Cont_TotalMorula1)
TME123_unfiltered$TME104Cont_TotalMorula3<- as.integer(TME123_unfiltered$TME104Cont_TotalMorula3)
TME123_unfiltered$TME118Cont_TotalMorula1<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula1)
TME123_unfiltered$TME118Cont_TotalMorula3<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula3)
TME123_unfiltered$TME118Cont_TotalMorula4<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula4)
TME123_unfiltered$TME118Cont_TotalMorula5<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula5)
TME123_unfiltered$TME118Cont_TotalMorula6<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula6)
TME123_unfiltered$TME118Cont_TotalMorula7<- as.integer(TME123_unfiltered$TME118Cont_TotalMorula7)
TME123_unfiltered$TME101Cont_SmallMorula1<- as.integer(TME123_unfiltered$TME101Cont_SmallMorula1)
TME123_unfiltered$TME102Cont_SmallMorula2<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula2)
TME123_unfiltered$TME102Cont_SmallMorula4<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula4)
TME123_unfiltered$TME102Cont_SmallMorula6<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula6)
TME123_unfiltered$TME102Cont_SmallMorula8<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula8)
TME123_unfiltered$TME104Cont_SmallMorula2<- as.integer(TME123_unfiltered$TME104Cont_SmallMorula2)
TME123_unfiltered$TME112Cont_smallMorula1<- as.integer(TME123_unfiltered$TME112Cont_smallMorula1)
TME123_unfiltered$TME118Cont_SmallMorula1<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula1)
TME123_unfiltered$TME118Cont_SmallMorula2<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula2)
TME123_unfiltered$TME101Cont_SmallMorula2<- as.integer(TME123_unfiltered$TME101Cont_SmallMorula2)
TME123_unfiltered$TME102Cont_SmallMorula3<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula3)
TME123_unfiltered$TME102Cont_SmallMorula5<- as.integer(TME123_unfiltered$TME102Cont_SmallMorula5)
TME123_unfiltered$TME104Cont_SmallMorula1<- as.integer(TME123_unfiltered$TME104Cont_SmallMorula1)
TME123_unfiltered$TME118Cont_SmallMorula3<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula3)
TME123_unfiltered$TME118Cont_SmallMorula4<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula4)
TME123_unfiltered$TME118Cont_SmallMorula5<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula5)
TME123_unfiltered$TME118Cont_SmallMorula6<- as.integer(TME123_unfiltered$TME118Cont_SmallMorula6)
TME123_unfiltered$TME101PolyIVFMorula1<- as.integer(TME123_unfiltered$TME101PolyIVFMorula1)
TME123_unfiltered$TME101PolyIVFMorula3<- as.integer(TME123_unfiltered$TME101PolyIVFMorula3)
TME123_unfiltered$TME101PolyIVFMorula4<- as.integer(TME123_unfiltered$TME101PolyIVFMorula4)
TME123_unfiltered$TME102PolyIVFMorula2<- as.integer(TME123_unfiltered$TME102PolyIVFMorula2)
TME123_unfiltered$TME102PolyIVFMorula3<- as.integer(TME123_unfiltered$TME102PolyIVFMorula3)
TME123_unfiltered$TME102PolyIVFMorula4<- as.integer(TME123_unfiltered$TME102PolyIVFMorula4)
TME123_unfiltered$TME102PolyIVFMorula5<- as.integer(TME123_unfiltered$TME102PolyIVFMorula5)
TME123_unfiltered$TME102PolyIVFMorula7<- as.integer(TME123_unfiltered$TME102PolyIVFMorula7)
TME123_unfiltered$TME102PolyIVFMorula8<- as.integer(TME123_unfiltered$TME102PolyIVFMorula8)
TME123_unfiltered$TME112PolyIVFMorula1<- as.integer(TME123_unfiltered$TME112PolyIVFMorula1)
TME123_unfiltered$TME112PolyIVFMorula2<- as.integer(TME123_unfiltered$TME112PolyIVFMorula2)
TME123_unfiltered$TME112PolyIVFMorula3<- as.integer(TME123_unfiltered$TME112PolyIVFMorula3)
TME123_unfiltered$TME112PolyIVFMorula5<- as.integer(TME123_unfiltered$TME112PolyIVFMorula5)
TME123_unfiltered$TME118PolyIVFMorula1<- as.integer(TME123_unfiltered$TME118PolyIVFMorula1)
TME123_unfiltered$TME118PolyIVFMorula2<- as.integer(TME123_unfiltered$TME118PolyIVFMorula2)
TME123_unfiltered$TME118PolyIVFMorula3<- as.integer(TME123_unfiltered$TME118PolyIVFMorula3)
TME123_unfiltered$TME101PolyIVFMorula2<- as.integer(TME123_unfiltered$TME101PolyIVFMorula2)
TME123_unfiltered$TME101PolyIVFMorula5<- as.integer(TME123_unfiltered$TME101PolyIVFMorula5)
TME123_unfiltered$TME102PolyIVFMorula6<- as.integer(TME123_unfiltered$TME102PolyIVFMorula6)
TME123_unfiltered$TME112PolyIVFMorula4<- as.integer(TME123_unfiltered$TME112PolyIVFMorula4)
TME123_unfiltered$TME112PolyIVFMorula6<- as.integer(TME123_unfiltered$TME112PolyIVFMorula6)

#filter data - will use cutoff to avoid jackpotting
TME123_filtered <- TME123_unfiltered
#from preprocessing, we know that there are 17 vont+NFW, 19 cont+total RNA, 17 cont+small RNA, 21 poly embryos
#want to remove any genes that do not appear in at least half of one of these groups
#create a frame with the condition for each embryo
condition <- c(rep("cont+NFW",17),rep("cont+total",19), rep("cont+small",17), rep("poly",21))
#create a data frame for use with DESeq2, make this for each tissue type
colData123Embryo <- data.frame(row.names =colnames(TME123_filtered), condition)
View(colData123Embryo)

#run DESeq2, make the dds matrix
TME123Embryo_dds <- DESeqDataSetFromMatrix(countData = TME123_filtered, colData = colData123Embryo, design = ~condition)

#we already know the number in each condition
#determine the number in each group divided by 2 to set cut off
contNFWFilt <- 8.5
contTotFilt <- 9.5
contSmallFilt <- 8.5
PolyFilt <- 10.5
#now move on to determining which rows we will keep because the counts are great enough
#determining if the at least half of the samples in either condition have counts for the gene, 5 cut off
keep <- rowSums(counts(TME123Embryo_dds[,1:17])>=5)>= contNFWFilt | rowSums(counts(TME123Embryo_dds[,18:36])>=5)>= contTotFilt | rowSums(counts(TME123Embryo_dds[,37:53])>=5)>= contSmallFilt | rowSums(counts(TME123Embryo_dds[,54:74])>=5)>= PolyFilt
TME123Embryo <- TME123Embryo_dds[keep,]

#continue to run DESeq2
TME123Embryo_dds <- DESeq(TME123Embryo)

#create matrix with results - the order of conditions here is important. In this example it will calculate log2fc heat/con so increase in heat sample would be positive log2fc
#do this for each combination against NFW cont
TME123Poly_results <- results(TME123Embryo_dds, contrast = c("condition", "poly", "cont+NFW"))
TME123Total_results <- results(TME123Embryo_dds, contrast = c("condition", "cont+total", "cont+NFW"))
TME123Small_results <- results(TME123Embryo_dds, contrast = c("condition", "cont+small", "cont+NFW"))

#run a summary on each set of samples
summary(TME123Poly_results)
summary(TME123Total_results)
summary(TME123Small_results)

#write res to CSV do for each tissue
write.csv(TME123Poly_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/TME123PolyvsNFW.csv', row.names = TRUE)
write.csv(TME123Total_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/TME123TotalvsNFW.csv', row.names = TRUE)

write.csv(TME123Small_results, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/TME123SmallvsNFW.csv', row.names = TRUE)

#pull out significant up and down res Changing to use pvalue, 
PolySig <- subset(TME123Poly_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
PolySigup <- subset (PolySig,log2FoldChange > 0)                                                          
PolySigdn <- subset (PolySig, log2FoldChange < 0)

TotalSig <- subset(TME123Total_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
TotalSigup <- subset (TotalSig,log2FoldChange > 0)                                                          
TotalSigdn <- subset (TotalSig, log2FoldChange < 0)

SmallSig <- subset(TME123Small_results, pvalue <0.05 & abs(log2FoldChange) > 0.58)
SmallSigup <- subset (SmallSig,log2FoldChange > 0)                                                          
SmallSigdn <- subset (SmallSig, log2FoldChange < 0)

#getting a df so we can connect the gene name back to the ID, do this for sig up and down each sample
PolySigup <- as.data.frame(PolySigup)
PolySigdn <- as.data.frame(PolySigdn)
TotalSigup <- as.data.frame(TotalSigup)
SmallSigup <- as.data.frame(SmallSigup)
SmallSigdn<- as.data.frame(SmallSigdn)
TotalSigdn <- as.data.frame(TotalSigdn)

#now make the geneIDs which is a row name a column in the df so we can match with the gene name
PolySigup$geneID <- rownames(PolySigup)
PolySigdn$geneID <- rownames(PolySigdn)
TotalSigup$geneID <- rownames(TotalSigup)
TotalSigdn$geneID <- rownames(TotalSigdn)
SmallSigup$geneID <- rownames(SmallSigup)
SmallSigdn$geneID <- rownames(SmallSigdn)


#trying to figure out how to get the gene name back
PolySigup <- merge(PolySigup, genetotranscript, by.x = "geneID", by.y = "transcript")
PolySigdn <- merge(PolySigdn, genetotranscript, by.x = "geneID", by.y = "transcript")
TotalSigup <- merge(TotalSigup, genetotranscript, by.x = "geneID", by.y = "transcript")
TotalSigdn <- merge(TotalSigdn, genetotranscript, by.x = "geneID", by.y = "transcript")
SmallSigup <- merge(SmallSigup, genetotranscript, by.x = "geneID", by.y = "transcript")
SmallSigdn <- merge(SmallSigdn, genetotranscript, by.x = "geneID", by.y = "transcript")

#try to switch column order so the gene name isn't at the end
setcolorder(PolySigup,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(PolySigdn,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(TotalSigup,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(TotalSigdn,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

setcolorder(SmallSigup,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
setcolorder(SmallSigdn,c("gene","geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

#write out the files with the gene name
write.csv(PolySigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/Poly_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(PolySigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/Poly_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(TotalSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/Total_Sig_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(TotalSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/Total_Sig_res_dn_withGene_diffCutOff.csv', row.names = TRUE)

write.csv(SmallSigup, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/Small_res_up_withGene_diffCutOff.csv', row.names = TRUE)
write.csv(SmallSigdn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/TME123Redo8-5-25/Small_res_dn_withGene_diffCutOff.csv', row.names = TRUE)

#find what is share overlapping changes
PolyvsTotUp <- Reduce(intersect,list(PolySigup$gene,TotalSigup$gene))
PolyvsTotDn <- Reduce(intersect,list(PolySigdn$gene,TotalSigdn$gene))

PolyvsSmallUp <- Reduce(intersect,list(PolySigup$gene,SmallSigup$gene))
PolyvsSmallDn <- Reduce(intersect,list(PolySigdn$gene,SmallSigdn$gene))

#all three
AllThreeUp <- Reduce(intersect,list(PolySigup$gene,SmallSigup$gene,TotalSigup$gene))
AllThreeDn<- Reduce(intersect,list(PolySigdn$gene,SmallSigdn$gene,TotalSigdn$gene))

#write out the overlaps
write.csv(AllThreeUp, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/AllThreeSharedUp.csv', row.names = TRUE)
write.csv(AllThreeDn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/AllThreeSharedDn.csv', row.names = TRUE)
write.csv(PolyvsTotUp, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/PolyandTotSharedUp.csv', row.names = TRUE)
write.csv(PolyvsTotDn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/PolyandTotSharedDn.csv', row.names = TRUE)
write.csv(PolyvsSmallUp, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/PolyandSmallSharedUp.csv', row.names = TRUE)
write.csv(PolyvsSmallDn, '/Users/taymil/Documents/TMEJuradoExpsAndData/TME123/PolyandSmallSharedDn.csv', row.names = TRUE)



