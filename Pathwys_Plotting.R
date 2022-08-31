########## ANNOTATION (CLEMENT_MSC) ##############
###### OMIXBOX ANNOTATION #######
setwd("/Users/cley/Desktop/CID")
annot <- read.delim("annotation_second_run.txt")
anno.annot <- read.delim("blast2go_annotation.annot")
UpCvE <- read.csv("UPCvE_ID.txt")
colnames(UpCvE) <- "trans"
colnames(anno.annot) <- c("trans","GO","Blast_hit")

upCvE_ann <- merge(anno.annot,UpCvE, "trans")
upCvE_ann2 <- upCvE_ann
colnames(upCvE_ann2) <- upCvE_ann2[1,]
class(upCvE_ann2)
write.table(upCvE_ann2, file = "upCvE_ann.annot")


##### PATHWAYS ANALYSIS 
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/blast2go_real/pathways_analysis")
cve <- read.csv("contrl_vs_exp")
names(cve)
cve1 <- subset(cve, select = c(X,log2FoldChange,lfcSE,pvalue,padj))
names(cve1) <- c("Name","logFC","logCPM","PValue","FDR")
write.table(cve1, file = "cve.tsv")


######## PLOTTING TOP 32 PATHWAYS EXPRESSED ##############
library(ggplot2)
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/blast2go_real/pathways_analysis/results/total pathways")
pathways <- read.csv("pathways.txt", sep = "\t")
pathways1 <- subset(pathways, select = c(Pathway,X.Seqs))
pathways2 <- subset(pathways1, pathways1$X.Seqs >= 100)
write.csv(pathways1, file = "total_expressed_pathways")



pathways3 <- pathways2
pathways3$Pathway <- factor(pathways3$Pathway,
                            levels = pathways3$Pathway[order(pathways3$X.Seqs, decreasing = FALSE)])
g <- ggplot(pathways3, aes(Pathway,X.Seqs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line()) + geom_bar(stat = "identity", fill ="steelblue") + coord_flip() + 
  labs(title = "Expressed pathways", y="Sequences", x="Pathways")
g
ggsave("pathways.pdf", g, "pdf")

#library(RColorBrewer)
#head(brewer.pal.info, 100)

######### CONTROL VS EXPLOITED UP REGULATED PATHWAYS #######
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/blast2go_real/pathways_analysis/results/upCvE")
CvE_pathways <- read.csv("up_pathways_table.txt", sep = "\t")
names(CvE_pathways) <- c("Tags","ID","FDR","P.Value","Nr.Test","Nr.Reference","Non.Annot.Test","Non.Annot.Reference")
CvE_withnames <- merge(CvE_pathways,pathways, by = "ID")
CvE_withnames_subset <- subset(CvE_withnames, select = c(Pathway, Nr.Test))

pathwaysCvE <- CvE_withnames_subset
pathwaysCvE$Pathway <- factor(pathwaysCvE$Pathway,
                              levels = pathwaysCvE$Pathway[order(pathwaysCvE$Nr.Test, decreasing = FALSE)])
gCvE <- ggplot(pathwaysCvE, aes(Pathway,Nr.Test)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line()) + geom_bar(stat = "identity", fill ="steelblue") + coord_flip() + 
  labs(title = "Up-regulated pathways in Exploited", y="Sequences", x="Pathways")
gCvE
ggsave("UpCvEpathways.pdf", gCvE, "pdf")

##########  CONTROL VS POLLUTES UP REGULATED PATHWAYS #######
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/blast2go_real/pathways_analysis/results/upCvP")
CvP_pathways <- read.csv("omicsbox_table.txt", sep = "\t")
names(CvP_pathways) <- c("Tags","ID","FDR","P.Value","Nr.Test","Nr.Reference","Non.Annot.Test","Non.Annot.Reference")
CvP_withnames <- merge(CvP_pathways, pathways, by = "ID")
CvP_withnames_subset <- subset(CvP_withnames, select = c(Pathway, Nr.Test))

pathwaysCvP <- CvP_withnames_subset
pathwaysCvP$Pathway <- factor(pathwaysCvP$Pathway,
                              levels = pathwaysCvP$Pathway[order(pathwaysCvP$Nr.Test, decreasing = FALSE)])
gCvP <- ggplot(pathwaysCvP, aes(Pathway,Nr.Test)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line()) + geom_bar(stat = "identity", fill ="steelblue") + coord_flip() + 
  labs(title = "Up-regulated pathways in Poluted", y="Sequences", x="Pathways")
gCvP
ggsave("UpCvPpathways.pdf", gCvP, "pdf")




