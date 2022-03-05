library(xlsx)
library(xlsxjars)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(ggdendro)
library(pheatmap)
library(clusterProfiler)
options(stringsAsFactors = F)

#Input files
metab_cpd <- read.xlsx("Hu lab.metabolites.KEGG.HMBD.xlsx",
                       sheetIndex = 1, startRow = 1, header = T, as.data.frame = T)

path_metab <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 1, startRow = 1, as.data.frame = T, header = F, check.names = T)

path_intro <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 2, startRow = 1, as.data.frame = T, header = F, check.names = T)

hcm.metab.plasma <- read.xlsx("Source Data Extended Data Fig. 4.xlsx", 
                       sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)

hcm.metab <- read.xlsx("Source Data Extended Data Fig. 2.xlsx", 
                       sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)


rownames(hcm.metab.plasma) <- hcm.metab.plasma$Metabolites
hcm.metab.plasma <- hcm.metab.plasma[,-1]
hcm.metab.plasma.t <- t(hcm.metab.plasma)
hcm.metab.plasma.group <- data.frame(hcm.metab.plasma.t,check.names = F)
hcm.metab.plasma.group$group <-c(rep(c("Ctrl","HCM"),c(60,143)))

hcm.metab.plasma.relative <- apply(hcm.metab.plasma.t, 2, function(x) x/mean(x[1:60]))
hcm.metab.plasma.relative.log <- log2(hcm.metab.plasma.relative[,1:142])

test_results.plasma <- data.frame(Metabolites = colnames(hcm.metab.plasma.t)[1:142])
# wilcox
test_results.plasma$wilcox.test <- apply(hcm.metab.plasma.relative.log[,1:142],2,
                                  function(x) unlist(wilcox.test(as.numeric(x) ~ hcm.metab.plasma.group$group, 
                                                                 data = hcm.metab.plasma.relative.log)[3]))
# adjust p value using BH method
test_results.plasma$wilcox.test_BH <- p.adjust(test_results.plasma$wilcox.test, method = "fdr")

#fold-change
test_results.plasma$FC <- apply(hcm.metab.plasma.group[,1:142], 2, 
                         function(x) mean(x[which(hcm.metab.plasma.group$group == 'HCM')])/mean(x[which(hcm.metab.plasma.group$group == 'Ctrl')]))
test_results.plasma$LOG2FC <- log2(test_results.plasma$FC)

test_results.plasma$cpd <- metab_cpd$KEGG[match(test_results.plasma$Metabolite, metab_cpd$Metabolite)]

# differential altered metabolites FC 1.5 
test_results_diff.plasma <- test_results.plasma[test_results.plasma$wilcox.test_BH < 0.05 & 
                                    abs(test_results.plasma$LOG2FC) > 0.584963, ]

test_results_diff.plasma$path <- path_metab$X1[match(test_results_diff.plasma$cpd, path_metab$X2)]
test_results_diff.plasma$pathway <- path_intro$X2[match(test_results_diff.plasma$path, path_intro$X1)]


rownames(hcm.metab.tissue) <- hcm.metab.tissue$Metabolites
hcm.metab.tissue <- hcm.metab.tissue[,-1]
hcm.metab.tissue <-hcm.metab.tissue[,c(1:365)]

hcm.metab.tissue.t <- t(hcm.metab.tissue)
hcm.metab.tissue.group <- data.frame(hcm.metab.tissue.t,check.names = F)
hcm.metab.tissue.group$group <-c(rep(c("Ctrl","HCM"),c(16,349)))
hcm.metab.tissue.relative <- apply(hcm.metab.tissue.t, 2, function(x) x/mean(x[1:16]))
hcm.metab.tissue.relative.log <- log2(hcm.metab.tissue.relative[,1:154])

test_results.tissue <- data.frame(Metabolites = colnames(hcm.metab.tissue.t)[1:154])

# wilcox
test_results.tissue$wilcox.test <- apply(hcm.metab.tissue.relative.log[,1:154],2,
                                         function(x) unlist(wilcox.test(as.numeric(x) ~ hcm.metab.tissue.group$group, 
                                                                        data = hcm.metab.tissue.relative.log)[3]))
# adjust p value using BH method
test_results.tissue$wilcox.test_BH <- p.adjust(test_results.tissue$wilcox.test, method = "fdr")

#fold-change
test_results.tissue$FC <- apply(hcm.metab.tissue.group[,1:154], 2, 
                                function(x) mean(x[which(hcm.metab.tissue.group$group == 'HCM')])/mean(x[which(hcm.metab.tissue.group$group == 'Ctrl')]))
test_results.tissue$LOG2FC <- log2(test_results.tissue$FC)

test_results.tissue$cpd <- metab_cpd$KEGG[match(test_results.tissue$Metabolite, metab_cpd$Metabolite)]
# differential altered metabolites FC 1.5 
test_results_diff.tissue <- test_results.tissue[test_results.tissue$wilcox.test_BH < 0.05 & 
                                                  abs(test_results.tissue$LOG2FC) > 0.584963, ]

test_results_diff.tissue$path <- path_metab$X1[match(test_results_diff.tissue$cpd, path_metab$X2)]
test_results_diff.tissue$pathway <- path_intro$X2[match(test_results_diff.tissue$path, path_intro$X1)]


#####Ext. Fig 4a#####
hcm.hm <- data.frame(t(hcm.metab.plasma.relative.log), check.names = F)
annotation_col <- data.frame(row.names = colnames(hcm.hm))
annotation_col$group <-rep(c("Ctrl","HCM"),c(60,143))
annotation_colors <- list(group = c("HCM"="#F8766D", "Ctrl"="#00BFC4"))
rg <- c(rep('#104E8B', 100), colorRampPalette(c("#104E8B", "white"))(50)[1:47], 
        colorRampPalette(c("white", "#B22222"))(50)[4:50], rep('#B22222', 100))

p <- pheatmap(hcm.hm, legend = T, scale = "row", 
              color = rg, border_color = 'grey90',
              cellwidth = 1.25, cellheight = 3.25, #cex = 2.5,
              annotation_col = rev(annotation_col),
              annotation_colors = annotation_colors,
              fontsize = 6, fontsize_row = 3, fontsize_col = 0.1, treeheight_row = 0, treeheight_col = 15,
              annotation_names_row = T, annotation_names_col = F,
              clustering_distance_rows = "correlation", clustering_method = "complete", 
              cluster_cols = F, cluster_rows = T,
              main = '')


#####Ext. Fig 4b#####
df_pca <- prcomp( hcm.metab.plasma.relative.log, scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = hcm.metab.plasma.group$group) 

percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,fill=Species))+ 
  scale_fill_manual(values = c("#6FE7DD", "#3490DE"))+
  geom_point(size = 1, alpha = 1, shape = 21, color = "black")+  
  xlab(percentage[1]) +
  ylab(percentage[2])+
  stat_ellipse(level = 0.95, show.legend = F, linetype = 2)+
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1)+
  
  theme(panel.border = element_rect(colour = "black", size= 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line= element_line(colour = "black"))+ 
  
  theme(plot.title = element_text( size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 6, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.4, "cm"), 
        legend.spacing = unit(0.4, "cm"))

######## Fig 4a #######
test_results <- read.xlsx("Source Data Fig. 4.xlsx", 
                          sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)

ggplot(test_results, aes(x = test_results$LOG2FC, y = test_results$log10wilcox.test_BH, colour=test_results$change)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_manual(values=c("#F65A00", "grey","#AA32DF"))+
  geom_text_repel(data = test_results, aes(x = test_results$LOG2FC, 
                                           y = test_results$log10wilcox.test_BH, label = label),
                  size = 1.5,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.4, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)+
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +
  labs(x="Log2(Fold change, HCM/Ctrl)", y="-Log10 (FDR)", color = "")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  xlim(c(-5,5))+
  ylim(c(0,27))+
  theme(axis.title.y = element_text(size = 6, face = "plain"),
        axis.title.x = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        axis.text.y = element_text(size = 6, face = "plain",))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right") +
  theme(legend.text = element_text(size = 6), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.3, "cm"))


######## Fig 4b #######
PT_FC <- read.xlsx("Source Data Fig. 4.xlsx", 
                   sheetIndex = 2, startRow = 1, 
                   as.data.frame = T, header = T, check.names = F)

ggplot(PT_FC, aes(x = PT_FC$`T-LOG2FC`,
                  y = PT_FC$`P-LOG2FC`,color=PT_FC$pathway)) +
  theme_classic() +
  geom_point(size = 4, alpha = 0.8) +
  #scale_colour_manual(values=alpha(c("pink3","darkseagreen4",
  #                                   "darkgoldenrod2", "lightsteelblue4"), 0.6)) + 
  labs(x="T-H/C", y="P-H/C", title="") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(plot.title = element_text(vjust = 6, hjust = 0.5,  size=6, face="bold"))  +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size=6, face="bold"),
        axis.title.y=element_text(vjust = 5, size=6, face="bold")) +
  theme(axis.text.y = element_text(size = 6,face="bold", colour = "black"),
        axis.text.x = element_text(size = 6,face="bold", colour = "black")) +
  geom_vline(xintercept=c(0,0),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = c(0,0),lty=4,col="grey",lwd=0.4)+
  scale_x_continuous(limits=c(-6, 6)) + 
  scale_y_continuous(limits=c(-3, 5))+
  geom_text_repel(data = PT_FC, aes(x = PT_FC$`T-LOG2FC`,
                                    y = PT_FC$`P-LOG2FC`,color=PT_FC$pathway, label = label),
                  size = 2,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.4, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)


######## Fig 4c #######

x_enrich <- enricher(test_results_diff.tissue$cpd, 
                     TERM2GENE = path_metab, TERM2NAME = path_intro,
                     minGSSize = 2, pvalueCutoff = 0.05, pAdjustMethod = "fdr") 

kegg_table <- as.data.frame(x_enrich)  
kegg_table <- na.omit(kegg_table)

kegg_table$FoldEnrich <- apply(kegg_table, 1,
                               function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                 as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
kegg_table$Description <- factor(kegg_table$Description, levels = rev(kegg_table$Description))


# calculate the absolute mean log2(FC) of metabolites in each pathway
path_avelog2FC <- function(x_enrich){
  metabs <- unlist(strsplit(x_enrich[8], "[/]"))
  mean(abs(test_results_diff.tissue[match(metabs, test_results_diff.tissue$cpd),5]))
}

kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC)

TP_kegg <- read.xlsx("Source Data Fig. 4.xlsx", 
                       sheetIndex = 3, startRow = 1, as.data.frame = T, header = T, check.names = F)

ggplot(TP_kegg, aes(Group, Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(1, 3), breaks = c(2,4,6)) +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(1, 1.8, 2.5)) +
  theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.4) +
  labs(x = "", y = "", title = "", 
       fill = "abs. Log2 (Fold change)", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
  ) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))

