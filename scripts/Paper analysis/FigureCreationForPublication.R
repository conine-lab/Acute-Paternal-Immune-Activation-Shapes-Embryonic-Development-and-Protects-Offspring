#this script is for the figures for publication seen in "Acute Paternal Immune Activation Shapes ..." TME et al 2026
#this script is hardcoded and shows the exact code used to generate the figures seen in the publication

#script for figure 2 heatmpap
TME33_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/TME33BulkSeqDolphinFiles/genes_expression_tpm.tsv', header = TRUE, sep ="\t")
#pull out the specific genes we want
#make heatmap with less genes, removing UBD, SLC15a3, H2-QY
ShortIntersectList <- c("Cxcl10","Oasl1","Ccl5","Ccl2","Rsad2","Mx1","Ccl4","Lag3",
                        "Ccl11","Oas3","Irf7","Gbp5","Ifit2","Ccl7","Gbp6",
                        "Isg15","Oas2","Zbp1","Oas1g","Ifit3",
                        "Ifit1","Gbp3","Cfb","Mnda","Bst2")
#make Rownames The gene names (make first column row names)
#had a duplicated gene but removed it assuming it is not a gene we want
TME33_update <- subset(TME33_unfiltered, gene!= "1-Mar")
TME33_update <- subset(TME33_update, gene!= "2-Mar")
#checked to make sure duplicated were removed
TME33_update[duplicated(TME33_update$gene),]
#making gene names, row names
rownames(TME33_update) <- TME33_update[,1]
TME33_update2 <- TME33_update[,3:30]

mat <- as.matrix(TME33_update2[rownames(TME33_update2) %in% ShortIntersectList, ])
#work on making heatmap
Caput <- c("cap1","cap2","cap3","cap4", "cap5","cap6","cap7")
Corpus <- c("cor1","cor2","cor3","cor4","cor5","cor6","cor7")
Cauda <- c("caud1","caud2","caud3","caud4","caud5","caud6","caud7")
Testes <- c("t1","t2","t3","t4","t5","t6","t7")
#Poly <- c("cap1","cap2","cap3","cap4","cor1","cor2","cor3","cor4","caud1","caud2","caud3","caud4","t1","t2","t3","t4")
#Cont <- c("cap5","cap6","cap7","cor5","cor6","cor7","caud5","caud6","caud7","t5","t6","t7")
#vector to subset
KeepSamples <- c(Testes,Caput,Corpus,Cauda)
#subset matrix to only have poly samples
mat2 <- mat[, KeepSamples]
gene_groups <- factor(rep(c("Testes","Caput",  "Corpus", "Cauda"), 
                          times = c(length(Caput), length(Corpus), length(Cauda), length(Testes))),
                      levels = c("Testes","Caput",  "Corpus", "Cauda"))
#add in a vector to make vector for treatment
#Condition_tissue <- factor(c("Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont"),
#levels = c("Poly","Cont"))
# Ensure matrix rows are in the same order as identity
mat2 <- mat2[match(ShortIntersectList, rownames(mat2)), ]

myCol <- colorRampPalette(c('grey90', 'red2'))(50)
myBreaks <- seq(-10, 20, length.out = 50)
COLORS = list()
COLORS[["Tissue"]] = c("Testes" = "#489FA7",
                       "Caput" = "#584482",
                       "Corpus" = "#A4B0D1",
                       "Cauda" = "#883268")

#install.packages("colorRamp2")
library(colorRamp2)

#log transform data to make this look better
#make the poly and contol conditions using col and heatmap annotation
Cols = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = COLORS[["Tissue"]])), Condition = c(rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3)),
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
                   #bottom_annotation = Cols,
                   #bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = COLORS[["Condition"]]))),
                   show_column_names = FALSE,
                   use_raster = TRUE)

setwd("/Users/taymil/Documents")
svglite("TME33HeatmapBetterIntersectLessGenes.svg", width=5.5, height=5.5)
heatmap 
dev.off()
#make Cytokine and PRR for all samples all conditions
TME33_unfiltered <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/TME33BulkSeqDolphinFiles/genes_expression_tpm.tsv', header = TRUE, sep ="\t")
Cytokine <- c("Tnfrsf1a","Ifnar2","Ifngr1","Il6st","Ifnar1","Ifnlr1","Il6ra")
PRR <- c("Ddx58","Tlr3","Ifih1","Tlr7","Tlr8")
CytoAndPRR <- c(Cytokine, PRR)
TME33_update <- subset(TME33_unfiltered, gene!= "1-Mar")
TME33_update <- subset(TME33_update, gene!= "2-Mar")
#checked to make sure duplicated were removed
TME33_update[duplicated(TME33_update$gene),]
#making gene names, row names
rownames(TME33_update) <- TME33_update[,1]
TME33_update2 <- TME33_update[,3:30]
mat <- as.matrix(TME33_update2[rownames(TME33_update2) %in% CytoAndPRR, ])
#first going to do expression in all conditions
#make indidivual vectors for each tissues
Caput <- c("cap5","cap6","cap7","cap1","cap2","cap3","cap4")
Corpus <- c("cor5","cor6","cor7","cor1","cor2","cor3","cor4")
Cauda <- c("caud5","caud6","caud7","caud1","caud2","caud3","caud4")
Testes <- c("t5","t6","t7","t1","t2","t3","t4")
KeepSamples <- c(Testes,Caput,Corpus,Cauda)
#subset matrix to only have poly samples
mat2 <- mat[, KeepSamples]
gene_groups <- factor(rep(c("Testes","Caput",  "Corpus", "Cauda"), 
                          times = c(length(Caput), length(Corpus), length(Cauda), length(Testes))),
                      levels = c("Testes","Caput",  "Corpus", "Cauda"))
#add in a vector to make vector for treatment
#Condition_tissue <- factor(c("Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont","Poly","Poly", "Poly", "Poly" ,"Cont", "Cont","Cont"),
#levels = c("Poly","Cont"))
# Ensure matrix rows are in the same order as identity

mat2 <- mat2[match(CytoAndPRR, rownames(mat2)), ]

myCol <- colorRampPalette(c('grey95', 'black'))(50)
myBreaks <- seq(-10, 20, length.out = 50)
COLORS = list()
COLORS[["Tissue"]] = c("Testes" = "#489FA7",
                       "Caput" = "#584482",
                       "Corpus" = "#A4B0D1",
                       "Cauda" = "#883268")
#try to make something so we can separate into cytokine vs PRR
gene_type <- factor(rep(c("Cytokine Receptor","PRR"), 
                        times = c(length(Cytokine), length(PRR))),
                    levels = c("Cytokine Receptor","PRR"))
#COLORS[["Type"]] = c("Cytokine"= "blue",
#"PRR" = "purple")
#install.packages("colorRamp2")
library(colorRamp2)

#log transform data to make this look better
#make the poly and contol conditions using col and heatmap annotation
Cols = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = COLORS[["Tissue"]])), Condition = c(rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4), rep("Cont", 3),rep("Poly", 4)),
                         col = list(Condition = c("Poly" = "white", "Cont" = "snow3")))
mat2 <- log2(mat2+1)
myBreaks <- seq(0, 7, length.out = 50)
heatmap <- Heatmap(mat2, name = "Expression in 
Log2(TPM)",  
                   column_split = factor(gene_groups),
                   cluster_columns = FALSE,
                   show_column_dend = FALSE,
                   cluster_column_slices = TRUE,
                   column_title_gp = gpar(fontsize = 10),
                   column_gap = unit(0.4, "mm"),
                   row_split = factor(gene_type),
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   col = colorRamp2(myBreaks, myCol),
                   row_names_gp = gpar(fontsize = 9),
                   column_title_rot = 0,
                   top_annotation = Cols,
                   #left_annotation = rowAnnotation(foo=anno_block(gp=gpar(fill=COLORS[["Type"]]))),
                   #bottom_annotation = Cols,
                   #bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = COLORS[["Condition"]]))),
                   show_column_names = FALSE,
                   use_raster = TRUE)
setwd("/Users/taymil/Documents")
#export as svg so it can be edited in illustrator
library(svglite)
svglite("TME33HeatmapCytoPRR09172024AllConditions.svg", width = 5.5, height = 5.5)
heatmap 
dev.off()

#code used for supplemental figure 2
#sorted the output files with gene names for sig up and down by padj value (smallest to largest)
#took top 50 most sig altered up or down genes and put into gprofiler with default mouse settings
#saved CSV
#read in the results from gprofiler
TestUpGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025TestUp.csv')
#now order the table by the Pvalu and make it so most sig is first (ascending order)
TestUpGProf <- TestUpGProf[order(TestUpGProf$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
TestUpGProfTop <- head(TestUpGProf, n = 10)
TestUpGProfTop <- TestUpGProfTop[order(TestUpGProfTop$negative_log10_of_adjusted_p_value)]
#load ggplot
library(ggplot2)
#make a horizontal bar graph with gradient fill, the plum color is for up regulated genes
plot <- ggplot(data = TestUpGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "thistle1", high = "plum4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,30)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("TestUpGPro.svg", width=6, height=3.5)
plot 
dev.off()
#no we will do this for the Top Down regulated genes in testes
TestDnGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025TestDn.csv')
TestDnGProf <- TestDnGProf[order(TestDnGProf$adjusted_p_value),]
TestDnGProfTop <- head(TestDnGProf, n = 10)
TestDnGProfTop <- TestDnGProfTop[order(TestDnGProfTop$negative_log10_of_adjusted_p_value),]
#the blue color is for down regulated
plot<- ggplot(data = TestDnGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "lightsteelblue1", high = "darkslategray4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,12)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("TestDnGPro.svg", width=6, height=3.5)
plot 
dev.off()
#do this for the up go terms in caput
CapUpGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025CapUp.csv')
#now order the table by the Pvalu and make it so most sig is first (ascending order)
CapUpGProf <- CapUpGProf[order(CapUpGProf$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
CapUpGProfTop <- head(CapUpGProf, n = 10)
CapUpGProfTop <- CapUpGProfTop[order(CapUpGProfTop$negative_log10_of_adjusted_p_value),]
#load ggplot
library(ggplot2)
#make a horizontal bar graph with gradient fill, the plum color is for up regulated genes
plot<- ggplot(data = CapUpGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "thistle1", high = "plum4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,30)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CapUpGPro.svg", width=6, height=3.5)
plot 
dev.off()
#no we will do this for the Top Down regulated genes in Cap
CapDnGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025CapDn.csv')
CapDnGProf <- CapDnGProf[order(CapDnGProf$adjusted_p_value),]
CapDnGProfTop <- head(CapDnGProf, n = 10)
CapDnGProfTop <- CapDnGProfTop[order(CapDnGProfTop$negative_log10_of_adjusted_p_value)]
#the blue color is for down regulated
plot<- ggplot(data = CapDnGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "lightsteelblue1", high = "darkslategray4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+ 
  ylim(0,12)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CapDnGPro.svg", width=6, height=3.5)
plot 
dev.off()
#do this for the up go terms in corpus
CorUpGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025CorUp.csv')
#now order the table by the Pvalu and make it so most sig is first (ascending order)
CorUpGProf <- CorUpGProf[order(CorUpGProf$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
CorUpGProfTop <- head(CorUpGProf, n = 10)
CorUpGProfTop <- CorUpGProfTop[order(CorUpGProfTop$negative_log10_of_adjusted_p_value)]
#load ggplot
library(ggplot2)
#make a horizontal bar graph with gradient fill, the plum color is for up regulated genes
plot<- ggplot(data = CorUpGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "thistle1", high = "plum4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+ 
  ylim(0,30)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CorUpGPro.svg", width=6, height=3.5)
plot 
dev.off()

#no we will do this for the Top Down regulated genes in Corpus
CorDnGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025CorDn.csv')
CorDnGProf <- CorDnGProf[order(CorDnGProf$adjusted_p_value),]
CorDnGProfTop <- head(CorDnGProf, n = 10)
CorDnGProfTop <- CorDnGProfTop[order(CorDnGProfTop$negative_log10_of_adjusted_p_value)]
#the blue color is for down regulated
plot<- ggplot(data = CorDnGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "lightsteelblue1", high = "darkslategray4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,12)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CorDnGPro.svg", width=6, height=3.5)
plot 
dev.off()

#do this for the up go terms in cauda
CaudUpGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025CaudaUp.csv')
#now order the table by the Pvalu and make it so most sig is first (ascending order)
CaudUpGProf <- CaudUpGProf[order(CaudUpGProf$adjusted_p_value),]
#get out the top 10 results (don't want to graph everything)
CaudUpGProfTop <- head(CaudUpGProf, n = 10)
CaudUpGProfTop <- CaudUpGProfTop[order(CaudUpGProfTop$negative_log10_of_adjusted_p_value)]
#load ggplot
library(ggplot2)
#make a horizontal bar graph with gradient fill, the plum color is for up regulated genes
plot<- ggplot(data = CaudUpGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_fill_gradient2(low = "thistle1", high = "plum4")+
  coord_flip()+
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,30)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CaudUpGPro.svg", width=6, height=3.5)
plot 
dev.off()
#no we will do this for the Top Down regulated genes in Cap
CaudDnGProf <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25/gProfiler_mmusculus_9-3-2025_5CaudaDn.csv')
CaudDnGProf <- CaudDnGProf[order(CaudDnGProf$adjusted_p_value),]
CaudDnGProfTop <- head(CaudDnGProf, n = 10)
CaudDnGProfTop <- CaudDnGProfTop[order(CaudDnGProfTop$negative_log10_of_adjusted_p_value)]
#the blue color is for down regulated
plot<- ggplot(data = CaudDnGProfTop, aes(x = reorder(term_name, +negative_log10_of_adjusted_p_value), y = negative_log10_of_adjusted_p_value))+ 
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), width = .8)+
  scale_x_discrete(expand = expansion(add = c(.2,.2))) +
  scale_fill_gradient2(low = "lightsteelblue1", high = "darkslategray4")+
  coord_flip()+ 
  labs(x="Go Analysis Term", y = "-log10(padj)", fill = "")+
  ylim(0,12)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData/TME33and34/BulkSeqRedo9-2-25")
svglite("CaudDnGPro.svg", width=6, height=3.5)
plot 
dev.off()

#code for 3 C and E
#make new graphs from the IPA output
#take in a table of the IPA output for canonical pathways
#want to graph the -log(pvalue) and ZScore
#color bars by zScore and bar length by p-value
#only graph things that have a Z score of absval 1 or greater and greater than 0 pvalue
#try first with the 4 cell overall
#read in the file, this was hand made to get out the stuff I need
Cell4IPA <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/IPAAnalysis102725/TMEForIPA/IPAOutputs_TME62/4CellPathonly.csv', header = TRUE)
#subset on things that have a pvalue that isn't 0
Cell4Filt <- subset(Cell4IPA, X.log.p.value. > 0 & abs(zScore) > .99)
#order the table by pvalue, so highest is first
Cell4Filt <- Cell4Filt[order(Cell4Filt$X.log.p.value., decreasing = TRUE),]
library(ggplot2)
#make the plot
plot <- (ggplot(data = Cell4Filt, aes(x = reorder(Ingenuity.Canonical.Pathways, + X.log.p.value.), y = X.log.p.value.))+ 
           geom_col(aes(fill = zScore), width = .8)+
           scale_fill_gradient2(low = "aquamarine4", high = "violetred4")+
           coord_flip()+
           labs(x="Pathway", y = "-log(pvalue)", fill = "zScore")+
           ylim(0,5)+
           theme(aspect.ratio = 1.5/2)+
           theme(axis.title=element_text(size =10)))
#write out to see if this is better for illustrator
library(svglite)
setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("TME624CellIPA.svg", width=7, height=5)
plot 
dev.off()

#same thing for morula overall
MorulaIPA <- read.csv('/Users/taymil/Documents/TMEJuradoExpsAndData/TME62Files/Redo7-25-25/IPAAnalysis102725/TMEForIPA/IPAOutputs_TME62/MorulaPathonly.csv', header = TRUE)
#subset on things that have a pvalue that isn't 0
MorulaFilt <- subset(MorulaIPA, X.log.p.value. > 0 & abs(zScore) > .99)
#order the table by pvalue, so highest is first
MorulaFilt <- MorulaFilt[order(MorulaFilt$X.log.p.value., decreasing = TRUE),]

#make the plot
plot2 <- ggplot(data = MorulaFilt, aes(x = reorder(Ingenuity.Canonical.Pathways, + X.log.p.value.), y = X.log.p.value.))+ 
  geom_col(aes(fill = zScore), width = .8)+
  scale_fill_gradient2(low = "aquamarine4", high = "violetred4")+
  coord_flip()+
  labs(x="Pathway", y = "-log(pvalue)", fill = "zScore")+
  ylim(0,5)+
  theme(aspect.ratio = 1.5/2)+
  theme(axis.title=element_text(size =10))

setwd("/Users/taymil/Documents/TMEJuradoExpsAndData")
svglite("TME62MorulaIPA.svg", width=6, height=3.5)
plot2 
dev.off()
