library(xlsx)
library(ConsensusClusterPlus)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggdendro)
library(pheatmap)
library(survival)
library(survminer)
library(dplyr)
library(cluster)
library(clusterProfiler)
library(hrbrthemes)
library(viridis)
library(dunn.test)
options(stringsAsFactors = F)

#Input files
metab_cpd <- read.xlsx("Hu lab.metabolites.KEGG.HMBD.xlsx",
                       sheetIndex = 1, startRow = 1, header = T, as.data.frame = T)

path_metab <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 1, startRow = 1, as.data.frame = T, header = F, check.names = T)

path_intro <- read.xlsx("Metabolic pathway.xlsx",
                        sheetIndex = 2, startRow = 1, as.data.frame = T, header = F, check.names = T)
#metabolomics
hcm.metab <- read.xlsx("Source Data Extended Data Fig. 2.xlsx", 
                       sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)
rownames(hcm.metab) <- hcm.metab$Metabolites
hcm.metab <- hcm.metab[,-1]
hcm.metab <- data.frame((hcm.metab[,17:365]), check.names = F)

hcm.metab.t <- t(hcm.metab)
metab_name <- data.frame(ori = colnames(hcm.metab.t), check = make.names(colnames(hcm.metab.t)))
colnames(hcm.metab.t) <- make.names(colnames(hcm.metab.t))
hcm.metab.t <- data.frame(hcm.metab.t, check.names = F)

#lipidomics
hcm.lipo <- read.xlsx("Source Data Fig. 6.xlsx", 
                      sheetIndex = 6, startRow = 1, as.data.frame = T, header = T, check.names = F)

rownames(hcm.lipo) <- hcm.lipo$Metabolites
hcm.lipo <- hcm.lipo[,-1]
hcm.lipo <- data.frame((hcm.lipo[,17:365]), check.names = F)

hcm.lipo.t <- t(hcm.lipo)
lipo_name <- data.frame(ori = colnames(hcm.lipo.t), check = make.names(colnames(hcm.lipo.t)))
colnames(hcm.lipo.t) <- make.names(colnames(hcm.lipo.t))
hcm.lipo.t <- data.frame(hcm.lipo.t, check.names = F)


#combing
hcm.com <- cbind(hcm.metab.t, hcm.lipo.t)

# clinical information
clinical <- read.xlsx("HCM clinical information.xlsx", 
                      sheetIndex = 1, header = T, startRow = 1, 
                      as.data.frame = T, check.names = F,encoding= "UTF-8")
rownames(clinical) <- clinical$Patient
clinical <- clinical[,-1]
cli_name <- data.frame(ori = colnames(clinical), check = make.names(colnames(clinical)))
colnames(clinical) <- make.names(colnames(clinical))

# Fig 6a Ext. Fig 5a-d ------------------------------
top <- nrow(hcm.metab) * 1
mads <- order(apply(hcm.metab, 1, mad), decreasing = T)[1:top]
df.filtered <- as.matrix(hcm.metab[mads, ])
df.filtered.z <- t(scale(t(df.filtered), scale = T, center = T))

results <- ConsensusClusterPlus(df.filtered.z, 
                                maxK = 6, 
                                reps = 1000, pItem = 0.8, pFeature = 0.8, 
                                seed = 999999999, distance = "pearson", clusterAlg = "kmdist", 
                                title = "Fig 6-metab", plot = "pdf") 

clustRes <- data.frame(Cluster = results[[3]][["consensusClass"]])
clustRes$Value <- rep(1,349)
clustRes <- clustRes[order(clustRes$Cluster, decreasing = F),]

# Fig 6b ------------------------------
rownames(clustRes) <- clustRes$Patient
clustRes <- clustRes[,-1]
hcm.heatmap <-hcm.metab
hcm.heatmap <- hcm.heatmap[match(rownames(clustRes), colnames(hcm.heatmap))]

annotation_col <- clustRes[,]
annotation_col$Cluster <- as.factor(annotation_col$Cluster)
annotation_colors <- list(Cluster = c("1" = "#A5CDE2", "2" = "#1E78B3", "3"="#B2DF8A"))
rg = c(rep('dodgerblue4', times = 50),
       colorRampPalette(c("dodgerblue4", "white", "firebrick"))(45),
       rep("firebrick", times = 50))
p <- pheatmap(hcm.heatmap, legend = T, scale = "row", 
              color = rg, border_color = 'grey90',
              cellwidth = 0.4, cellheight = 0.8, #cex = 2.5, 
              annotation_col = rev(annotation_col),
              annotation_colors = annotation_colors,
              fontsize = 7, fontsize_row = 0.5, fontsize_col = 0.1,  
              annotation_names_row = T, treeheight_row = 10,
              annotation_names_col = F,
              clustering_distance_rows = "correlation", clustering_method = "average", 
              cluster_cols = F, cluster_rows = T,
              main = '')



# Ext.Fig 5e f------------
clust <- data.frame(Sample = rownames(clustRes), Cluster = clustRes[,1], Value = clustRes[,2])
clust$Cluster <- as.factor(clust$Cluster)

clust$NYHA2 <- factor(clinical$NYHA..class.I.II..III.IV.[match(clust$Sample, rownames(clinical))],c('1', '2'))
clust$MWT <- clinical$MWT..mm.[match(clust$Sample, rownames(clinical))]

clust$MWT_1 <- ifelse(clust$MWT < 30, 1, 2)
clust$MWT_1 <- factor(clust$MWT_1, levels = c("1", "2"), 
                      labels = c("( < 30)", "(>= 30)"))

if(T){
  NYHA2_color <- c("#8D8DFF", "#FFA0A0")
  MWT_color <- c("#93ADC3", "#F2DCDD")
}

ggplot(clust) +
  
#geom_bar(mapping = aes(Cluster, Value, fill = NYHA2, color = NYHA2),
  #         width = 0.4, stat='identity', position='fill') +
  #  scale_fill_manual(values = NYHA2_color) +
  #  scale_color_manual(values = NYHA2_color) +
  # labs(x = 'Subtype', y = 'Proportion', title = '', fill = "NYHA level", color = "NYHA level")+
  

    geom_bar(mapping = aes(Cluster, Value, fill = MWT_1, color = MWT_1),
                   width = 0.4, stat='identity', position='fill') +
           scale_fill_manual(values = MWT_color) +
     scale_color_manual(values = MWT_color) +
       labs(x = 'Subtype', y = 'Proportion', title = '', fill = "MWT (mm)", color = "MWT (mm)")+

theme_classic() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 3) +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size = 6, face = "plain"),
        axis.title.y=element_text(vjust = 5, size = 6, face = "plain")) +
  theme(legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        legend.spacing = unit(0.3, "cm"), 
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.25),
        axis.line.y=element_line(linetype=1,color="black",size=0.25),
        axis.ticks.x=element_line(color="black",size=0.25,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.25,lineend = 1))

fisher.test(x = clust$Cluster, y = clust$MWT_1)

# Fig 6c ---------------------------

clinical.surv <- clinical[,c(1:3,22:24)]
clinical.surv <- na.omit(clinical.surv)
clinical.surv$cluster <- clustRes$Cluster[match(rownames(clinical.surv), rownames(clustRes))]
clinical.surv$cluster <- factor(clinical.surv$cluster, levels = c(1,2,3))
clinical.surv$OS <- clinical.surv$Survival.time.after.surgery..year.[]
clinical.surv$OS_status <- ifelse(clinical.surv$cardiovascular.death > 0, 1, 0 )

coxph(Surv(OS, OS_status) ~ cluster, data = clinical.surv)
logrank <- survdiff(Surv(OS, OS_status) ~ cluster, data = clinical.surv)
p.value <- 1 - pchisq(logrank$chisq, length(logrank$n) -1)
tcga_sfit <- survfit(Surv(OS, OS_status) ~ cluster, data = clinical.surv)

pairwise_survdiff(Surv(OS, OS_status) ~ cluster, data = clinical.surv)

ggsurvplot(tcga_sfit, 
           ggtheme = theme_bw(), palette = c("#1668B7", "#67B166", "#EC4000"),
           linetype = c(1,1,1), size = c(0.4, 0.4,0.4), censor.shape = 3, censor.size = 1,
           font.title = c(7,"black"), font.x = c(7,"black"), font.y = c(7, "black"), 
           font.tickslab = c(7, "black"), font.legend = c(7,"black"),
           
           title = "", xlab = 'Time (years)', ylab = 'Survival probability', 
           # legend
           legend = "right", legend.title = "", legend.labs = c("S-I", "S-II", "S-III"), 
           
           ylim = c(0.84, 1), # xscale = 'd_m', break.x.by = ((365.25/12)), 
           
           conf.int = F, pval.method = F, pval = T, pval.method.size = 3.6, pval.size = 2.5, 
           pval.method.coord = c(0,0.85), pval.coord = c(0,0.87),
           
           risk.table = F)



# Fig 6d ---------------------------

df.test <- data.frame(t(log2(hcm.metab)), check.names = F)
df.test <- df.test[match(rownames(clustRes), rownames(df.test)), ]
df.test$Cluster <- clustRes$Cluster
df.test$Cluster <- as.factor(df.test$Cluster)

#Kruskal.test with Dunn comparion
KW.test <- data.frame()
Dunn.test <- data.frame()
KW.Dunn <-data.frame()
for (i in c(1:154)) {
  KWRes<-kruskal.test(df.test[[i]]~Cluster,data=df.test)
  KW.test <- rbind(KW.test, unlist(KWRes[3]))
  DunnRes <-dunn.test(df.test[[i]], df.test$Cluster, method="BH", altp = T)
  Dunn.test <- rbind(Dunn.test, unlist(DunnRes[4]))
  KW.Dunn <- cbind(KW.test, Dunn.test)
}
colnames(KW.Dunn) <- c("p",unlist(DunnRes[5]))
rownames(KW.Dunn) <- colnames(df.test)[1:154]  

comparison <- combn(3, 2)
df.test.fc <- data.frame(t(hcm.metab), check.names = F)
df.test.fc <- df.test.fc[match(rownames(clustRes), rownames(df.test.fc)), ]
df.test.fc$Cluster <- clustRes$Cluster
df.test.fc$Cluster <- as.factor(df.test.fc$Cluster)

for (i in c(1:3)) {
  FCname <- paste0("FC", comparison[,i][1], "-", comparison[,i][2])
  KW.Dunn[[FCname]] <- apply(df.test.fc[, 1:154], 2, 
                                 function(x) 
                                   mean(x[which(df.test.fc$Cluster == comparison[,i][1])])/
                                   mean(x[which(df.test.fc$Cluster == comparison[,i][2])]))
  LOG2FCname <- paste0("LOG2FC", comparison[,i][1], "-", comparison[,i][2])
  KW.Dunn[[LOG2FCname]] <- log2(KW.Dunn[[FCname]])
}

dm.df <- data.frame()
for (i in c(2:4)) {
  dm <- KW.Dunn[KW.Dunn[,i] < 0.05 & abs(KW.Dunn[, i+i+2]) > log2(1.25), ]
  dm$group <- colnames(KW.Dunn)[i]
  dm$metab <- rownames(dm)
  dm.df <- rbind(dm.df, dm)
}

# point plot of KEGG enrichment results
C123_kegg<-read.xlsx("Source Data Fig. 6.xlsx", 
                     sheetIndex = 4, startRow = 1, as.data.frame = T, header = T, check.names = F)

ggplot(C123_kegg, aes(Group, Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(1, 3), breaks = c(2,7,12)) +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.5, 0.8, 1.1)) +
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


# Fig 6f Ext. Fig5 i-l---------------------------
top <- nrow(hcm.lipo) * 1
mads <- order(apply(hcm.lipo, 1, mad), decreasing = T)[1:top]
df.filtered <- as.matrix(hcm.lipo[mads, ])
df.filtered.z <- t(scale(t(df.filtered), scale = T, center = T))

results <- ConsensusClusterPlus(df.filtered.z, 
                                maxK = 6, 
                                reps = 1000, pItem = 0.8, pFeature = 0.8, 
                                seed = 999999999, distance = "pearson", clusterAlg = "kmdist", 
                                title = "Fig 6-lipo", plot = "pdf") 

clustRes <- data.frame(Cluster = results[[4]][["consensusClass"]])
clustRes$Value <- rep(1,349)


# Ext. Fig5 m n ---------------------------

clust <- data.frame(Sample = rownames(clustRes), Cluster = clustRes[,1], Value = clustRes[,2])
clust$Cluster <- as.factor(clust$Cluster)

clust$NYHA2 <- factor(clinical$NYHA..class.I.II..III.IV.[match(clust$Sample, rownames(clinical))],c('1', '2'))
clust$MWT <- clinical$MWT..mm.[match(clust$Sample, rownames(clinical))]

clust$MWT_1 <- ifelse(clust$MWT < 30, 1, 2)
clust$MWT_1 <- factor(clust$MWT_1, levels = c("1", "2"), 
                      labels = c("( < 30)", "(>= 30)"))

if(T){
  NYHA2_color <- c("#8D8DFF", "#FFA0A0")
  MWT_color <- c("#93ADC3", "#F2DCDD")
}

ggplot(clust) +
    geom_bar(mapping = aes(Cluster, Value, fill = NYHA2, color = NYHA2),
           width = 0.4, stat='identity', position='fill') +
    scale_fill_manual(values = NYHA2_color) +
   scale_color_manual(values = NYHA2_color) +
    labs(x = 'Subtype', y = 'Proportion', title = '', fill = "NYHA level", color = "NYHA level")+
  
 #   geom_bar(mapping = aes(Cluster, Value, fill = MWT_1, color = MWT_1),
  #           width = 0.4, stat='identity', position='fill') +
  #   scale_fill_manual(values = MWT_color) +
# scale_color_manual(values = MWT_color) +
# labs(x = 'Subtype', y = 'Proportion', title = '', fill = "MWT (mm)", color = "MWT (mm)")+ 
  theme_classic() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 3) +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size = 6, face = "plain"),
        axis.title.y=element_text(vjust = 5, size = 6, face = "plain")) +
  theme(legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        legend.spacing = unit(0.3, "cm"), 
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.25),
        axis.line.y=element_line(linetype=1,color="black",size=0.25),
        axis.ticks.x=element_line(color="black",size=0.25,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.25,lineend = 1))

fisher.test(x = clust$Cluster, y = clust$MWT_1, workspace = 2e7)


# Fig 6g ---------------------------
clinical.surv <- clinical[,c(1:3,22:24)]
clinical.surv <- na.omit(clinical.surv)
clinical.surv$cluster <- clustRes$Cluster[match(rownames(clinical.surv), rownames(clustRes))]
clinical.surv$cluster <- factor(clinical.surv$cluster, levels = c(1,2,3))
clinical.surv$OS <- clinical.surv$Survival.time.after.surgery..year.[]
clinical.surv$OS_status <- ifelse(clinical.surv$cardiovascular.death > 0, 1, 0 )

coxph(Surv(OS, OS_status) ~ cluster, data = clinical.surv)
logrank <- survdiff(Surv(OS, OS_status) ~ cluster, data = clinical.surv)
p.value <- 1 - pchisq(logrank$chisq, length(logrank$n) -1)
tcga_sfit <- survfit(Surv(OS, OS_status) ~ cluster, data = clinical.surv)

pairwise_survdiff(Surv(OS, OS_status) ~ cluster, data = clinical.surv)

ggsurvplot(tcga_sfit, 
           ggtheme = theme_bw(), palette = c("#1668B7", "#67B166", "#EC4000"),
           linetype = c(1,1,1), size = c(0.4, 0.4,0.4), censor.shape = 3, censor.size = 1,
           font.title = c(7,"black"), font.x = c(7,"black"), font.y = c(7, "black"), 
           font.tickslab = c(7, "black"), font.legend = c(7,"black"),
           
           title = "", xlab = 'Time (years)', ylab = 'Survival probability', 
           # legend
           legend = "right", legend.title = "", legend.labs = c("S-I", "S-II", "S-III"), 
           
           ylim = c(0.84, 1), # xscale = 'd_m', break.x.by = ((365.25/12)), 
           
           conf.int = F, pval.method = F, pval = T, pval.method.size = 3.6, pval.size = 2.5, 
           pval.method.coord = c(0,0.85), pval.coord = c(0,0.87),
           
           risk.table = F)


# Fig 6h ---------------------------
df.test <- data.frame(t(log2(hcm.lipo)), check.names = F)
df.test <- df.test[match(rownames(clustRes), rownames(df.test)), ]
df.test$Cluster <- clustRes$Cluster
df.test$Cluster <- as.factor(df.test$Cluster)


#Kruskal.test with Dunn comparion
KW.test <- data.frame()
Dunn.test <- data.frame()
KW.Dunn <-data.frame()
for (i in c(1:768)) {
  KWRes<-kruskal.test(df.test[[i]]~Cluster,data=df.test)
  KW.test <- rbind(KW.test, unlist(KWRes[3]))
  DunnRes <-dunn.test(df.test[[i]], df.test$Cluster, method="BH", altp = T)
  Dunn.test <- rbind(Dunn.test, unlist(DunnRes[4]))
  KW.Dunn <- cbind(KW.test, Dunn.test)
}
colnames(KW.Dunn) <- c("p",unlist(DunnRes[5]))
rownames(KW.Dunn) <- colnames(df.test)[1:768]  


comparison <- combn(3, 2)
df.test.fc <- data.frame(t(hcm.lipo), check.names = F)
df.test.fc <- df.test.fc[match(rownames(clustRes), rownames(df.test.fc)), ]
df.test.fc$Cluster <- clustRes$Cluster
df.test.fc$Cluster <- as.factor(df.test.fc$Cluster)

for (i in c(1:3)) {
  FCname <- paste0("FC", comparison[,i][1], "-", comparison[,i][2])
  KW.Dunn[[FCname]] <- apply(df.test.fc[, 1:768], 2, 
                             function(x) 
                               mean(x[which(df.test.fc$Cluster == comparison[,i][1])])/
                               mean(x[which(df.test.fc$Cluster == comparison[,i][2])]))
  LOG2FCname <- paste0("LOG2FC", comparison[,i][1], "-", comparison[,i][2])
  KW.Dunn[[LOG2FCname]] <- log2(KW.Dunn[[FCname]])
}

dm.df <- data.frame()
for (i in c(2:4)) {
  dm <- KW.Dunn[KW.Dunn[,i] < 0.05 & abs(KW.Dunn[, i+i+2]) > log2(1.25), ]
  dm$group <- colnames(KW.Dunn)[i]
  dm$metab <- rownames(dm)
  dm.df <- rbind(dm.df, dm)
}

####plot####
cluster.lipo.bubble <- read.xlsx("Source Data Fig. 6.xlsx", 
                                 sheetIndex = 8, startRow = 1, as.data.frame = T, header = T, check.names = F)
cluster.lipo.bubble %>%
  arrange(desc(`LOG10 P`)) %>%
  ggplot(aes(x=number, y=`LOG2FC2-1`, size=`LOG10 P`, color=`lipid species`)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(0.3, 3), name="Population (M)" +
               scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A"))+
  theme(plot.title = element_text( size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 6, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.text = element_text(size = 2), 
        legend.key.size = unit(0.1, "cm"), 
        legend.spacing = unit(0.1, "cm"))+
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.45)+
  ylim(-1,3)+
  geom_hline(yintercept=c(-0.32,0.32),lty=4,col="grey",lwd=0.4)




#### Ext. Fig 5g h p ####

clinical.surv <- clinical[,c(1:3,22:24)]
clinical.surv <- na.omit(clinical.surv)
hcm.regres <- cbind(clinical, 
                    hcm.metab.t[match(rownames(clinical), rownames(hcm.metab.t)),],
                    hcm.lipo.t[match(rownames(clinical), rownames(hcm.lipo.t)), ])

hcm.regres <- hcm.regres[rownames(hcm.regres) %in% rownames(clinical.surv),]

if(T){
  for (i in c(25:ncol(hcm.regres))) {
    hcm.regres$level <- ifelse(hcm.regres[,i] > quantile(hcm.regres[,i], 0.5), "high", "low")
    clinical.surv[paste(colnames(hcm.regres)[i])] <- hcm.regres$level[match(rownames(clinical.surv), rownames(hcm.regres))]
  }
  clinical.surv <- na.omit(clinical.surv)
}

options(warn = -1)
p.value.total <- data.frame(colnames(clinical.surv)[7:928])
pl = c()

tcga_sfit <- coxph(Surv(Survival.time.after.surgery..year., cardiovascular.death) 
                   ~ clinical.surv$TAG48.4.C16.0., data = clinical.surv)
summary(tcga_sfit)


for (x in c(7:ncol(clinical.surv))) {
  
  tcga_sfit <- survfit(Surv(Survival.time.after.surgery..year., cardiovascular.death) ~ clinical.surv[,x], data = clinical.surv)
  logrank <- survdiff(Surv(Survival.time.after.surgery..year., cardiovascular.death) ~ clinical.surv[,x], data = clinical.surv)
  p.value <- 1 - pchisq(logrank$chisq, length(logrank$n) -1)
  p.value.total$p <-p.value
  if (p.value < 0.05) {
    
    pdf(paste("Ext. Fig 5",colnames(clinical.surv)[x] , ".pdf", sep = ""), 
        width = 2.82, height = 1.76, onefile = FALSE) 
    plot <- ggsurvplot(tcga_sfit, 
                       ggtheme = theme_bw(), 
                       # palette = c("firebrick", "dodgerblue4"), 
                       linetype = c(1,1), size = c(0.4, 0.4),
                       censor.shape = 3, censor.size = 1,
                       
                       title = comb_name$ori[comb_name$check == colnames(clinical.surv)[x]], xlab = 'Time (years)', ylab = 'Survival probability', 
                       
                       risk.table = F, tables.height = 0.3, fontsize = 2.6,
                       risk.table.col = "strata", risk.table.title = '', risk.table.y.text = T, 
                       tables.theme = clean_theme(), surv.plot.height = 0.5,
                       
                       legend = "right", legend.title = "", legend.labs = c("High", "Low"), 
                       
                       ylim = c(0.83, 1), # xscale = 'd_m', break.x.by = ((365.25/12)), 
                       font.title = c(6,"black"), font.x = c(6,"black"), font.y = c(6, "black"), 
                       font.tickslab = c(6, "black"), font.legend = c(6,"black"),
                       
                       conf.int = F, pval.method = F, pval = T, pval.method.size = 3.6, pval.size = 2.5, 
                       pval.method.coord = c(0,0.85), pval.coord = c(0,0.85))
    
    
    print(plot)
    dev.off()
  }
}

p.value_diff <- p.value.total[p.value.total$pval < 0.05,]


#### Ext. Fig 5o ####
data <-read.xlsx("Source Data Extended Data Fig. 5.xlsx", 
                 sheetIndex = 3, startRow = 1, as.data.frame = T, header = T, check.names = F)
data %>%
  arrange(desc(data$`Log10 (FDR)`)) %>%
  ggplot(aes(x=data$number, y=data$`Log10 (FDR)`, size=data$N, color=data$`lipid species`)) +
  geom_point(alpha=.5) +
  scale_size(range = c(2, 5), name="Population (M)" +
               scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A"))+
  theme(plot.title = element_text( size = 7, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 7, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 7, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 7, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 7, face = "plain", colour = "black"),
        legend.title = element_text(size = 7, face = "plain", colour = "black"),
        legend.text = element_text(size = 4), 
        legend.key.size = unit(0.4, "cm"), 
        legend.spacing = unit(0.4, "cm"))+
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.5)+
  
  geom_hline(yintercept = 0 ,lty=4,col="grey",lwd=0.2) 



