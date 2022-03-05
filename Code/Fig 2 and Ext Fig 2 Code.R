library(xlsx)
library(xlsxjars)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(ggdendro)
library(pheatmap)
library(ggcorrplot)
library(clusterProfiler)
options(stringsAsFactors = F)

#Input files
metab_cpd <- read.xlsx("Hu lab.metabolites.KEGG.HMBD.xlsx",
                       sheetIndex = 1, startRow = 1, header = T, as.data.frame = T)

path_metab <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 1, startRow = 1, as.data.frame = T, header = F, check.names = T)

path_intro <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 2, startRow = 1, as.data.frame = T, header = F, check.names = T)

hcm.metab <- read.xlsx("Source Data Extended Data Fig. 2.xlsx", 
                       sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)

rownames(hcm.metab) <- hcm.metab$Metabolites
hcm.metab <- hcm.metab[,-1]
hcm.metab.t <- t(hcm.metab)
hcm.metab.group <- data.frame(hcm.metab.t,check.names = F)
hcm.metab.group$group <-c(rep(c("Ctrl","HCM","DCM"),c(16,349,46)))

hcm.metab.relative <- apply(hcm.metab.t, 2, function(x) x/mean(x[1:16]))
hcm.metab.relative.log <- log2(hcm.metab.relative[,1:154])

HCM.metab.group <-hcm.metab.group[c(1:16,17:365),]
HCM.metab.relative.log <-hcm.metab.relative.log[c(1:16,17:365),]

HCM.test_results <- data.frame(Metabolites = colnames(hcm.metab.t)[1:154])

# wilcox
HCM.test_results$wilcox.test <- apply(HCM.metab.relative.log[,1:154],2,
                                  function(x) unlist(wilcox.test(as.numeric(x) ~ HCM.metab.group$group, 
                                                                 data = HCM.metab.relative.log)[3]))
# adjust p value using BH method
HCM.test_results$wilcox.test_BH <- p.adjust(HCM.test_results$wilcox.test, method = "fdr")

#fold-change
HCM.test_results$FC <- apply(HCM.metab.group[,1:154], 2, 
                         function(x) mean(x[which(HCM.metab.group$group == 'HCM')])/mean(x[which(HCM.metab.group$group == 'Ctrl')]))
HCM.test_results$LOG2FC <- log2(HCM.test_results$FC)


HCM.test_results$cpd <- metab_cpd$KEGG[match(HCM.test_results$Metabolite, metab_cpd$Metabolite)]

# differential altered metabolites FC 1.5 
HCM.test_results_diff <- HCM.test_results[HCM.test_results$wilcox.test_BH < 0.05 & 
                                    abs(HCM.test_results$LOG2FC) > 0.584963, ]


HCM.test_results$change <- ifelse(HCM.test_results$wilcox.test_BH < 0.05 & HCM.test_results$FC > 1.5, "Up", 
                              ifelse(HCM.test_results$wilcox.test_BH < 0.05 & HCM.test_results$FC < 0.67, 
                                     "Down","Not changed"))

DCM.metab.group <-hcm.metab.group[c(1:16,366:411),]
DCM.metab.relative.log <-hcm.metab.relative.log[c(1:16,366:411),]

DCM.test_results <- data.frame(Metabolites = colnames(hcm.metab.t)[1:154])

# wilcox
DCM.test_results$wilcox.test <- apply(DCM.metab.relative.log[,1:154],2,
                                      function(x) unlist(wilcox.test(as.numeric(x) ~ DCM.metab.group$group, 
                                                                     data = DCM.metab.relative.log)[3]))
# adjust p value using BH method
DCM.test_results$wilcox.test_BH <- p.adjust(DCM.test_results$wilcox.test, method = "fdr")

#fold-change
DCM.test_results$FC <- apply(DCM.metab.group[,1:154], 2, 
                             function(x) mean(x[which(DCM.metab.group$group == 'DCM')])/mean(x[which(DCM.metab.group$group == 'Ctrl')]))
DCM.test_results$LOG2FC <- log2(DCM.test_results$FC)
DCM.test_results$cpd <- metab_cpd$KEGG[match(DCM.test_results$Metabolite, metab_cpd$Metabolite)]

# differential altered metabolites FC 1.5 
DCM.test_results_diff <- DCM.test_results[DCM.test_results$wilcox.test_BH < 0.05 & 
                                            abs(DCM.test_results$LOG2FC) > 0.584963, ]


DCM.test_results$change <- ifelse(DCM.test_results$wilcox.test_BH < 0.05 & DCM.test_results$FC > 1.5, "Up", 
                                  ifelse(DCM.test_results$wilcox.test_BH < 0.05 & DCM.test_results$FC < 0.67, 
                                         "Down","Not changed"))

write.xlsx(DCM.test_results, "D://desktop/HCM final/HCM ctrl16-final/Metabolomics DCM VS Ctrl-FC.xlsx",
           row.names = T, col.names = T, append = T)


######## Fig 2a #######
test_results <- read.xlsx("Source Data Fig. 2.xlsx", 
                          sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)

test_results$log10wilcox.test_BH <- -log10(test_results$wilcox.test_BH)

ggplot(test_results, aes(x = test_results$LOG2FC, y = test_results$log10wilcox.test_BH, colour=test_results$change)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("#546de5", "grey","#ff4757"))+
  
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +
  labs(x="Log2(Fold change, HCM/Ctrl)", y="-Log10 (FDR)", color = "")+
  geom_text_repel(data = test_results, aes(x = test_results$LOG2FC, 
                                           y = test_results$log10wilcox.test_BH, label = label),
                  size = 2,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.4, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  xlim(c(-6,5))+
  ylim(c(0,10))+
  theme(axis.title.y = element_text(size = 6, face = "plain"),
        axis.title.x = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        axis.text.y = element_text(size = 6, face = "plain",))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right") +
  theme(legend.text = element_text(size = 6), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.3, "cm"))

######## Fig 2b f Ext. Fig 2e######
hcm.metab.t <-data.frame(hcm.metab.t,check.names = F)
#carnitine
carn <- hcm.metab.t[grep("arnitine", colnames(hcm.metab.t))]
aim.metab <- c("GSH","GSSG")
aim.metab.TCA <- c("cis-Aconitic acid","Citrate", "Isocitrate", "alpha-Ketoglutaric acid")
carn <- hcm.metab.t[colnames(hcm.metab.t) %in% aim.metab]
carn$Gr <- carn$GSSG/carn$GSH

carn <- apply(carn, 2, function(x) x/mean(x[1:16]))
carn <- log2(carn)

group = c("Ctrl","HCM", "DCM")
carn.box <- data.frame()
for (i in c(1:ncol(carn))) {
  carn.box <- rbind(carn.box, 
                    data.frame(value = carn[,i], 
                               metab = rep(colnames(carn)[i], 411), 
                               group = rep(group, c(16,349,46))))
}
carn.box$group <- factor(carn.box$group, level = c("Ctrl", "HCM","DCM"))

my_comparisons <- list(group)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
ggplot(carn.box, aes(x = group, y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.3) + 
  geom_boxplot(aes(fill = group), width = 0.7, outlier.shape = NA, size = 0.3, fatten = 1) +
  facet_wrap(~ metab, scales = "free", strip.position = "top", ncol = 4) + 
  scale_fill_manual(values = c("DCM"="#00BFC4","HCM"="#F8766D", "Ctrl"="#3490DE")) + 
  theme_bw() + # ggtheme
  theme(aspect.ratio = 1.2, 
        panel.border = element_rect(colour = "black", size= 0.3),
        panel.grid = element_line(size = 0)) +
  labs(x = "", y = "", color = "") +
  theme(plot.title = element_text(size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        strip.text = element_text(size = 6, face = "plain", hjust = 0.5, vjust = 0.9), 
        strip.background = element_rect(fill = NA, colour = NA))  +
  theme(axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 6, vjust = 2, face ="plain")) +
  theme(axis.text.x = element_text(size = 6, face ="plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face ="plain", colour = "black")) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     symnum.args = symnum.args, hide.ns = F, label = "p.format")

######## Fig 2d e #######
HCM.metab.t <-hcm.metab.t[c(1:16,17:365),]
HCM.metab.t <-data.frame(HCM.metab.t,check.names = F)

#SAM/SAH,VB6
aim.metab <- c("S-Adenosylmethionine", "S-Adenosylhomocysteine", "Pyridoxal 5'-phosphate", "Pyridoxamine")
carn <- HCM.metab.t[colnames(HCM.metab.t) %in% aim.metab]
carn$Sr <- carn$`S-Adenosylmethionine`/carn$`S-Adenosylhomocysteine`

carn <- apply(carn, 2, function(x) x/mean(x[1:16]))
carn <- log2(carn)

carn.group <- data.frame(carn,check.names = F)
carn.group$group <-c(rep(c("Ctrl","HCM"),c(16,349)))

GSH.test_results <- data.frame(Metabolites = colnames(carn)[1:7])

# wilcox
GSH.test_results$wilcox.test <- apply(carn[,1:7],2,
                                      function(x) unlist(wilcox.test(as.numeric(x) ~ carn.group$group, 
                                                                     data = carn)[3]))


group = c("Ctrl","HCM")
carn.box <- data.frame()
for (i in c(1:ncol(carn))) {
  carn.box <- rbind(carn.box, 
                    data.frame(value = carn[,i], 
                               metab = rep(colnames(carn)[i], 365), 
                               group = rep(group, c(16,349))))
}
carn.box$group <- factor(carn.box$group, level = c("Ctrl", "HCM"))

my_comparisons <- list(group)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
ggplot(carn.box, aes(x = group, y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.3) + 
  geom_boxplot(aes(fill = group), width = 0.7, outlier.shape = NA, size = 0.3, fatten = 1) +
  facet_wrap(~ metab, scales = "free", strip.position = "top", ncol = 4) + 
  scale_fill_manual(values = c("DCM"="#00BFC4","HCM"="#F8766D", "Ctrl"="#3490DE")) + 
  theme_bw() + # ggtheme
  theme(aspect.ratio = 1.4, 
        panel.border = element_rect(colour = "black", size= 0.3),
        panel.grid = element_line(size = 0)) +
  labs(x = "", y = "", color = "") +
  theme(plot.title = element_text(size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        strip.text = element_text(size = 6, face = "plain", hjust = 0.5, vjust = 0.9), 
        strip.background = element_rect(fill = NA, colour = NA))  +
  theme(axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 6, vjust = 2, face ="plain")) +
  theme(axis.text.x = element_text(size = 6, face ="plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face ="plain", colour = "black")) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     symnum.args = symnum.args, hide.ns = F, label = "p.format")

######## Fig 2c #######
x_enrich <- enricher(HCM.test_results_diff$cpd, 
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

ggplot(kegg_table, aes(x = -log10(p.adjust), y = Description)) + 
  geom_point(aes(fill = -log10(p.adjust),  size = FoldEnrich), shape = 21, color = "grey40") +
  scale_fill_viridis_c(direction = -1, end = 0.9, option = "C", limit = c(1, 10), breaks = c(3,6,9)) + 
  scale_size(range = c(1, 6), breaks = c(10, 20, 30)) + 
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 2) +
  labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (FDR)") +
  theme(panel.grid.major.x = element_line(size = 0), 
        panel.grid.minor.x = element_line(size = 0)) +
  theme(plot.title = element_text( size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 6, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.4, "cm"), 
        legend.spacing = unit(0.4, "cm"))+
  xlim(c(1,8))

#######Fig 2g###
G2.Heatmap <- read.xlsx("Source Data Fig. 2.xlsx", 
                        sheetIndex = 5, startRow = 1, as.data.frame = T, header = T, check.names = F)


G2.Heatmap <- as.matrix(G2.Heatmap)
rg <- c(rep('#05C0C5', 100), colorRampPalette(c("#05C0C5", "white"))(50)[1:50], 
        colorRampPalette(c("white", "#FF6C6C"))(50)[4:50], rep('#FF6C6C', 100))

pheatmap(G2.Heatmap,cluster_cols=F,cluster_rows = F,clustering_distance_row="euclidean",
         clustering_distance_cols="euclidean",clustering_method="ward.D",
         color= rg,
         margins=c(5,10),fontsize_row=6,cellheight=5, cellwidth=15)



#####Ext. Fig 2a#####
hcm.hm <- data.frame(t(hcm.metab.relative.log), check.names = F)
annotation_col <- data.frame(row.names = colnames(hcm.hm))
annotation_col$group <-rep(c("Ctrl","HCM","DCM"),c(16,349,46))
annotation_colors <- list(group = c("DCM"="#00BFC4","HCM"="#F8766D", "Ctrl"="#3490DE"))
rg <- c(rep('#3d67a3', 100), colorRampPalette(c("#3d67a3", "white"))(50)[1:47], 
        colorRampPalette(c("white", "#ce1020"))(50)[4:50], rep('#ce1020', 100))

p <- pheatmap(hcm.hm, legend = T, scale = "row", 
              color = rg, border_color = 'grey90',
              cellwidth = 0.7, cellheight = 3, #cex = 2.5,
              annotation_col = rev(annotation_col),
              annotation_colors = annotation_colors,
              fontsize = 6, fontsize_row = 3, fontsize_col = 0.1, treeheight_row = 0, treeheight_col = 15,
              annotation_names_row = T, annotation_names_col = F,
              clustering_distance_rows = "correlation", clustering_method = "complete", 
              cluster_cols = F, cluster_rows = T,
              main = '')

#####Ext. Fig 2b#####
df_pca <- prcomp(hcm.metab.relative.log, scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = hcm.metab.group$group) 

percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,fill=Species))+ 
  scale_fill_manual(values = c("#3490DE","#00BFC4", "#F8766D"))+
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


df_pca <- prcomp(DCM.metab.relative.log, scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = DCM.metab.group$group) 

percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,fill=Species))+ 
  scale_fill_manual(values = c("#3490DE","#00BFC4"))+
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

######## Ext. Fig 2d #######
DH_FC <- read.xlsx("Source Data Extended Data Fig. 2.xlsx", 
                   sheetIndex = 2, startRow = 1, as.data.frame = T, header = T, check.names = F)

ggplot(DH_FC, aes(x = DH_FC$`HCM-LOG2FC`,
                  y = DH_FC$`DCM-LOG2FC`,color=DH_FC$pathway)) +
  theme_classic() +
  geom_point(size = 1.5, alpha = 0.8) +
  #scale_colour_manual(values=alpha(c("darkgoldenrod2","pink3","darkseagreen4",
   #                                   "lightsteelblue4"), 0.6)) + 
  labs(x="HCM-H/C", y="DCM-H/C", title="") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(plot.title = element_text(vjust = 6, hjust = 0.5,  size=6, face="bold"))  +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size=6, face="bold"),
        axis.title.y=element_text(vjust = 5, size=6, face="bold")) +
  theme(axis.text.y = element_text(size = 6,face="bold", colour = "black"),
        axis.text.x = element_text(size = 6,face="bold", colour = "black")) +
  geom_vline(xintercept=c(0,0),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = c(0,0),lty=4,col="grey",lwd=0.4)+
  scale_x_continuous(limits=c(-6, 4)) + 
  scale_y_continuous(limits=c(-5, 3))+
  geom_text_repel(data = DH_FC, aes(x = DH_FC$`HCM-LOG2FC`,
                                    y = DH_FC$`DCM-LOG2FC`,color=DH_FC$pathway, label = label),
                  size = 2,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.4, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)

######## Ext. Fig 2f #######
rownames(HCM.test_results_diff) <- HCM.test_results_diff$Metabolites
HCM.test_results_diff <- HCM.test_results_diff[,-1]

HCM.metab.relative.log.t <- data.frame(t(HCM.metab.relative.log), check.names = F)

metab.diff.corr <- HCM.metab.relative.log.t[match(rownames(HCM.test_results_diff), 
                                                  rownames(HCM.metab.relative.log.t)), ]
metab.diff.corr <- metab.diff.corr[,17:365]

corr <-round(cor(t(metab.diff.corr), method = "spearman"), 3)
p.mat <-round(cor_pmat(t(metab.diff.corr), method = "spearman"), 4)

ggcorrplot(corr, method = "square",p.mat = p.mat,insig = "blank",type = "lower",
           hc.order = TRUE, hc.method = "ward.D2",
           tl.cex= 3,
           #outline.color = "white", 
           ggtheme = theme( panel.background = element_rect(fill="white"), 
                            panel.grid.major.x=element_blank(),
                            panel.grid.major.y=element_blank()),
           lab = T, lab_size = 0.5)+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1)

