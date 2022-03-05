library(xlsx)
library(xlsxjars)
library(ggplot2)
library(pheatmap)


#####Fig 7a###########
MSEA <- read.xlsx("Source Data Fig. 7.xlsx",
                   sheetIndex = 1, startRow = 1, header = T, as.data.frame = T,check.names = F)

MSEA$Term<-factor(MSEA$Term,levels = rev(MSEA$Term)) 

ggplot(MSEA,aes(x=Term,y=MSEA$`Fold Enrichment`,fill=MSEA$FDR)) +
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#F4D6C9', high='#A98FC3',breaks = seq(0,25,5)) +  
  scale_x_discrete(expand = c(0,0.8))+ 
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
        legend.spacing = unit(0.4, "cm"))

#####Fig 7b###########
G2.Heatmap <- read.xlsx("Source Data Fig. 7.xlsx",
                    sheetIndex = 2, startRow = 1, header = T, as.data.frame = T,check.names = F)

rownames(G2.Heatmap) <- G2.Heatmap$`Gene name`
G2.Heatmap <- G2.Heatmap[,-1]

rg = colorRampPalette(c("#8989FF", "white", "#FF6C6C"))(100)

p <- pheatmap(G2.Heatmap, legend = T, scale = "none", 
              color = rg, border_color = '#999999',
              cellwidth = 7, cellheight = 7, #cex = 2.5,
              fontsize = 6, fontsize_row = 6, fontsize_col = 0.1, treeheight_row = 0, treeheight_col = 15,
              annotation_names_row = T, annotation_names_col = F,
              clustering_distance_rows = "correlation", clustering_method = "complete", 
              cluster_cols = F, cluster_rows = F,
              main = '')


