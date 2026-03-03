# Acute-Paternal-Immune-Activation-Shapes-Embryonic-Development-and-Protects-Offspring
Code used in the analysis for the paper "Acute Paternal Immune Activation Shapes Embryonic Development..." Taylor Miller-Ensminger et al. 2026. Examples of DESeq2 code for sequencing analysis and code used to generate publication figures are provided

# Scripts
## For new analysis
This folder contains scripts that can be downloaded by the user, ask for user input, and will generate analysis outputs/ figures using the same processes used for this manuscript. These are "ready to run" documents and can easily be used right away for user ease. Each R script contains comments explaining the necessary user input and the corresponding output. It is recommended to run all documents within R studio.
DESeqProcessForPub.R:
  * This documents allows the user to input an expected counts file and return DESeq analysis.
  * Users must have and expected counts file save as a csv that contains only the samples for a specific comparison
  * Users must have the samples ordered by treatment group
  * Users must have a metadata table saved as a csv that contains the sample names and treatment conditions. The columns should be named sample and condition
  * Users will recieve files containing the DESeq results with and with out associated gene names

FiguresProcessForPub.R:
  - This Document allows users to recreate Figure 2D-E, 3C and E, and Supllemental Figure 2 with their own data set
  - Within each section is a comment containing the required csv to the user
  - Users will recieve svg files containing the generated figure(s)


## Paper analysis
This folder contains hard coded (defined paths, sample conditions, sample names) used in the analysis of this paper. The code here can be modified by a user for their own use as well as to replicate the results within this paper (will require path changes when done on a different users computer). This code should primarily serve as examples for other users in the multiple ways the above code can be used. This code also provides the full view of the analysis run for this publication.

DESeqPublicationAnalysis.R:
  + This document contains code for the different DESeq analyses used in this publicaton
  + A comment at the beginning of each section clarifies what analysis is being done


FigureCreationForPublicaton.R:
  * This document contains the code used to generate the figures within this publication
