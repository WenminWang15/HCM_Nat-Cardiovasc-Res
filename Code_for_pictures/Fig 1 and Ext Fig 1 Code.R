library(ggplot2)
library(xlsx)
library(xlsxjars)
library(ggdendro)
options(stringsAsFactors = F)


####Fig 1b-----------

mut <- read.xlsx("Source Data Fig. 1.xlsx", 
                 sheetIndex = 1, header = T, startRow = 1, as.data.frame = T, check.names = F)

mut.df <- data.frame()
gene <- c("MYH7", "MYBPC3", "TNNT2", "TNNI3", "MYL2", "TPM1", "ACTC1", "MYL3")
for (i in gene) {
  mut.df <- rbind(mut.df, 
                  data.frame("patient" = unique(mut$`Patient`), 
                             "mutation" = rep(i, length(unique(mut$`Patient`))), 
                             "alter" = mut$Alternations[
                               mut$Gene == i][match(unique(mut$`Patient`), 
                                                   mut$`Patient`[mut$Gene== i])]))}

mut.df <- mut.df[mut.df$patient %in% unique(mut.df$patient[!is.na(mut.df$alter)]),]
mut.df$alter[is.na(mut.df$alter)] <- "na"

MYH7 <- mut.df$patient[mut.df$mutation == "MYH7" & mut.df$alter != "na"]
MYH7.na <- mut.df$patient[mut.df$mutation == "MYH7" & mut.df$alter == "na"]
MYBPC3 <- mut.df$patient[mut.df$mutation == "MYBPC3" & mut.df$alter != "na"]
TNNT2 <- mut.df$patient[mut.df$mutation == "TNNT2" & mut.df$alter != "na"]
TNNI3 <- mut.df$patient[mut.df$mutation == "TNNI3" & mut.df$alter != "na"]
MYL2 <- mut.df$patient[mut.df$mutation == "MYL2" & mut.df$alter != "na"]
TPM1 <- mut.df$patient[mut.df$mutation == "TPM1" & mut.df$alter != "na"]
ACTC1 <- mut.df$patient[mut.df$mutation == "ACTC1" & mut.df$alter != "na"]
MYL3 <- mut.df$patient[mut.df$mutation == "MYL3" & mut.df$alter != "na"]

patient1 <- c(intersect(MYH7, MYBPC3), intersect(MYH7, TNNT2), intersect(MYH7, TNNI3), 
             intersect(MYH7, MYL2), intersect(MYH7, TPM1), intersect(MYH7, ACTC1), 
             intersect(MYH7, MYL3)) 
patient2 <- c(intersect(MYH7.na, MYBPC3), intersect(MYH7.na, TNNT2), intersect(MYH7.na, TNNI3), 
             intersect(MYH7.na, MYL2), intersect(MYH7.na, TPM1), intersect(MYH7.na, ACTC1), 
             intersect(MYH7.na, MYL3))
patient <- c(patient1, unique(MYH7[!MYH7 %in% patient1]), patient2)

mut.df$patient <- factor(mut.df$patient, levels = unique(patient))
mut.df$mutation <- factor(mut.df$mutation, levels = rev(gene))
mut.df$alter <- factor(mut.df$alter, levels = c("Missense Mutation", "Splice Site", "Frame Shift InDel", 
                                                "In Frame InDel", "Nonsense Mutation", "na"))

ggplot(mut.df, aes(patient, mutation)) +
  geom_tile(aes(fill = alter), size = 0.3, color = "white") + 
  scale_fill_manual(values = c("Missense Mutation"  = "#2276B2", 
                               "Splice Site" = "#afd03e", 
                               "Frame Shift InDel" = "#5bb757", 
                               "In Frame InDel" = "#FD4B2B", 
                               "Nonsense Mutation" = "#6f3d90", 
                               "na" = "grey95")) + 
  theme_dendro() + 
  theme(plot.margin = margin(0.1, 0.5, 0.1, 0.5, "cm")) +
  coord_fixed(ratio = 10) +
  labs(x = "", y = "", title = "", fill = "Alteration") +
  theme(plot.title = element_text(vjust = 6, hjust = 0.5,  size = 8, face="plain"))  +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size = 8, face="plain"),
        axis.title.y = element_text(vjust = 5, size = 8, face="plain")) +
  theme(axis.text.x = element_text(size = 0, face="plain", colour = "black", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 8, face="plain", colour = "black")) +
  theme(legend.title = element_text(size = 8, face="plain", colour = "black"), 
        legend.text = element_text(size = 8, face="plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), 
        legend.key.width = unit(0.3, "cm"))

mut.df <- mut.df[mut.df$alter!="na",]
mut.df$value <- 1/349
mut.df$mutation <- factor(mut.df$mutation, levels = gene)
mut.df$alter <- factor(mut.df$alter, levels = c("na", "Missense Mutation", "Splice Site", "Frame Shift InDel", 
                                                "In Frame InDel", "Nonsense Mutation"))

ggplot(mut.df) +
  geom_bar(aes(mutation, value, fill = alter, color = alter),
          width = 0.4, stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c("Missense Mutation"  = "#2276B2", 
                               "Splice Site" = "#afd03e", 
                               "Frame Shift InDel" = "#5bb757", 
                               "In Frame InDel" = "#FD4B2B", 
                               "Nonsense Mutation" = "#6f3d90", 
                               "na" = "white")) +
  scale_color_manual(values = c("Missense Mutation"  = "#2276B2", 
                                "Splice Site" = "#afd03e", 
                                "Frame Shift InDel" = "#5bb757", 
                                "In Frame InDel" = "#FD4B2B", 
                                "Nonsense Mutation" = "#6f3d90", 
                                "na" = "white")) +
  labs(x = '', y = 'Alteration frequency', title = '', fill = "Alteration", color = "Alteration") +
  
  scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3)) +
  theme_classic() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.5) +
  # coord_fixed(ratio = 4) +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size = 8, face = "plain"),
        axis.title.y = element_text(vjust = 5, size = 8, face = "plain", colour = "black")) +
  theme(legend.title = element_text(size = 8, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 8, face = "plain", colour = "black", angle = 45, vjust = 1, hjust = 1), 
        legend.text = element_text(size = 8, face = "plain", colour = "black"),
        legend.spacing = unit(0.3, "cm"), 
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"))



####Ext. Fig 1b-----------

df <- read.xlsx("Source Data Extended Data Fig. 1.xlsx", 
                 sheetIndex = 1, header = T, startRow = 1, as.data.frame = T, check.names = F)
rownames(df) <- df$Metabolites
df <- df[,-1]
df_pca <- prcomp( log2(df), scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = df$Batch) 
percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=Species))+ 
  geom_point(size = 0.2, alpha = 1)+  
  xlab(percentage[1]) +
  ylab(percentage[2])+
  stat_ellipse(level = 0.95, show.legend = F)+
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", size= 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line= element_line(colour = "black"))+
  
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.8)+
  
  theme(plot.title = element_text( size = 8, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 8, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 8, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 8, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black"),
        legend.title = element_text(size = 8, face = "plain", colour = "black"),
        legend.text = element_text(size = 8), 
        legend.key.size = unit(0.4, "cm"), 
        legend.spacing = unit(0.4, "cm"))

