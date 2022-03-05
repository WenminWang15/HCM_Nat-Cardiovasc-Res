library(xlsx)
library(xlsxjars)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(survival)
library(dplyr)
library(survminer)
library(hrbrthemes)
library(viridis)
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
hcm.metab <- hcm.metab[,c(1:365)]
hcm.metab.t <- t(hcm.metab)

hcm.metab.relative <- apply(hcm.metab.t, 2, function(x) x/mean(x[1:16]))
hcm.metab.relative.log <- log2(hcm.metab.relative[,1:154])


metab_name <- data.frame(ori = colnames(hcm.metab.t), check = make.names(colnames(hcm.metab.t)))
colnames(hcm.metab.t) <- make.names(colnames(hcm.metab.t))
hcm.metab.t <- data.frame((hcm.metab.t[17:365,]), check.names = F)
hcm.metab.relative <- data.frame((hcm.metab.relative[17:365,]), check.names = F)
hcm.metab.relative.log <- data.frame((hcm.metab.relative.log[17:365,]), check.names = F)



#lipidomics
hcm.lipo <- read.xlsx("Source Data Extended Data Fig. 3.xlsx", 
                      sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)

rownames(hcm.lipo) <- hcm.lipo$Metabolites
hcm.lipo <- hcm.lipo[,-1]
hcm.lipo.t <- t(hcm.lipo)
hcm.lipo.relative <- apply(hcm.lipo.t, 2, function(x) x/mean(x[1:16]))
hcm.lipo.relative.log <- log2(hcm.lipo.relative[,1:768])

lipo_name <- data.frame(ori = colnames(hcm.lipo.t), check = make.names(colnames(hcm.lipo.t)))
colnames(hcm.lipo.t) <- make.names(colnames(hcm.lipo.t))
hcm.lipo.t <- data.frame((hcm.lipo.t[17:365,]), check.names = F)
hcm.lipo.relative <- data.frame((hcm.lipo.relative[17:365,]), check.names = F)
hcm.lipo.relative.log <- data.frame((hcm.lipo.relative.log[17:365,]), check.names = F)

#combing
hcm.com <- cbind(hcm.metab.relative.log, hcm.lipo.relative.log)

# clinical information
clinical <- read.xlsx("HCM clinical information.xlsx", 
                      sheetIndex = 1, header = T, startRow = 1, 
                      as.data.frame = T, check.names = F,encoding= "UTF-8")
rownames(clinical) <- clinical$Patient
clinical <- clinical[,-1]
cli_name <- data.frame(ori = colnames(clinical), check = make.names(colnames(clinical)))
colnames(clinical) <- make.names(colnames(clinical))




########Fig 5a#####
hcm.com$Mutation <- clinical$Mutation[match(rownames(clinical), rownames(hcm.com))]
hcm.com <- na.omit(hcm.com)
hcm.com$Mutation <- ifelse(hcm.com$Mutation > 0, "Mutation","No mutation")

hcm.comb.del <- hcm.com[,-923]
df_pca <- prcomp(hcm.comb.del, scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = hcm.com$Mutation) 

percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,fill=Species))+ 
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
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

########Fig 5b#####
hcm.com <- cbind(hcm.metab.relative.log, hcm.lipo.relative.log)
cli.muta <- read.xlsx("Source Data Fig. 5.xlsx", 
                      sheetIndex = 1, header = T, startRow = 1, 
                      as.data.frame = T, check.names = F,encoding= "UTF-8")
rownames(cli.muta) <- cli.muta$Patient
hcm.com$muta.type <- cli.muta$Gene[match(rownames(cli.muta), rownames(hcm.com))]
hcm.com.cli <- na.omit(hcm.com)

hcm.comb.del <- hcm.com.cli[,-2]
df_pca <- prcomp(hcm.comb.del, scale.=TRUE)
df_pcs <-data.frame(df_pca$x, Species = hcm.com.cli$muta.type) 

percentage<-round((df_pca$sdev)^2/sum((df_pca$sdev)^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,fill=Species))+ 
  scale_fill_manual(values = c("#00BFC4", "#F8766D","#D7AB33","#33CB85","#D296FF","#96BE33","#FF81D6","#33BAFF"))+
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



########Fig 5c#####
hcm.com <- cbind(hcm.metab.t, hcm.lipo.t)
rownames(cli.muta) <- cli.muta$Patient
hcm.com$muta.type <- cli.muta$Gene[match(rownames(cli.muta), rownames(hcm.com))]
hcm.com.cli <- na.omit(hcm.com)

hcm.com.MYH7 <- hcm.com.cli[which(hcm.com.cli$muta.type == 'MYH7'),]
hcm.com.MYBPC3 <- hcm.com.cli[which(hcm.com.cli$muta.type == 'MYBPC3'),]
hcm.com.doub <- rbind(hcm.com.MYH7,hcm.com.MYBPC3)

test_results <- data.frame(Metabolites = colnames(hcm.com)[2:923])

# wilcox
test_results$wilcox.test <- apply(hcm.com.doub[,2:923],2,
                                  function(x) unlist(wilcox.test(as.numeric(x) ~ hcm.com.doub$muta.type, 
                                                                 data = hcm.com.doub)[3]))
# adjust p value using BH method
test_results$wilcox.test_BH <- p.adjust(test_results$wilcox.test, method = "fdr")

#fold-change
test_results$FC <- apply(hcm.com[,2:923], 2, 
                         function(x) mean(x[which(hcm.com.doub$muta.type == 'MYBPC3')])
                         /mean(x[which(hcm.com.doub$muta.type == 'MYH7')]))
test_results$LOG2FC <- log2(test_results$FC)

test_results_diff <- test_results[test_results$wilcox.test_BH < 0.05 & 
                                    abs(test_results$LOG2FC) > 0.321928, ]

test_results$log10wilcox.test <- -log10(test_results$wilcox.test)

test_results$change <- ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC > 1.25, "Up", 
                              ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC < 0.8, 
                                     "Down","Not changed"))


test_results <- read.xlsx("Source Data Fig. 5.xlsx", 
                       sheetIndex = 2, startRow = 1, as.data.frame = T, header = T, check.names = F)

ggplot(test_results, aes(x = test_results$LOG2FC, y = test_results$log10wilcox.test, colour=test_results$change)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("grey","#ff4757"))+
 
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +
  

  labs(x="Log2(Fold change, MYBPC3/MYH7)", y="-Log10 (FDR)", color = "")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ylim(0,4)+
  theme(axis.title.y = element_text(size = 6, face = "plain"),
        axis.title.x = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        axis.text.y = element_text(size = 6, face = "plain",))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right") +
  theme(legend.text = element_text(size = 6), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.3, "cm"))

############Fig 5d##########
hcm.regres <- cbind(clinical, 
                    hcm.metab.t[match(rownames(clinical), rownames(hcm.metab.t)),])

hcm.or <- hcm.regres
m = "NYHA..class.I.II..III.IV."
hcm.or[[m]] <- ifelse(hcm.or[[m]] == 1, 0, 1)
subgroup <- levels(factor(hcm.or[[m]]))
subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])

hcm.or$Gender.0M1F. <- as.factor(hcm.or$Gender.0M1F.)
hcm.or[,25:178] <- log10(hcm.or[,25:178])
hcm.or[,25:178] <- apply(hcm.or[,25:178], 2, scale)

ORRes.df <- data.frame()
for (x in colnames(hcm.or)[25:178]) {
  fmla <- as.formula(paste(m," ~ ", paste(c("Gender.0M1F.", "Age..year.", "BMI..m2.kg.", x), collapse= " + ")))
  ORRes <- glm(fmla, data = hcm.or, family = binomial())
  # unisum <- summary(ORRes)
  OR <- exp(coef(ORRes))[5]
  OR <- round(OR, 2)
  ci <- exp(confint(ORRes))[5,]
  ci <- round(ci, 2)
  cito <- paste0(ci[1], " - ", ci[2])
  
  ORRes.s <- data.frame(coef(summary(ORRes)))
  ORRes.s$BH <- p.adjust(ORRes.s$Pr...z.., method = "BH")
  ORRes.df <- rbind(ORRes.df, 
                    data.frame("Metab" = x, "NYHA" = subgroup1, "OR" = OR, 
                               "95%CI" = cito, "p.value" = ORRes.s$Pr...z..[5], "BH" = ORRes.s$BH[5]))
}

ORRes.df <- ORRes.df[ORRes.df$BH < 0.05, ]
ORRes.df$name <- lipo_name$ori[match(ORRes.df$Metab, lipo_name$check)]

ORRes.df$cpd <- metab_cpd$KEGG[match(ORRes.df$name, metab_cpd$Metabolite)]
ORRes.df$path <- path_metab$X1[match(ORRes.df$cpd, path_metab$X2)]
ORRes.df$pathway <- path_intro$X2[match(ORRes.df$path, path_intro$X1)]

############Fig 5e##########

#lipidomics
hcm.lipo <- read.xlsx("Source Data Fig. 6.xlsx", 
                      sheetIndex = 6, startRow = 1, as.data.frame = T, header = T, check.names = F)

rownames(hcm.lipo) <- hcm.lipo$Metabolites
hcm.lipo <- hcm.lipo[,-1]
hcm.lipo.t <- t(hcm.lipo)
hcm.lipo.relative <- apply(hcm.lipo.t, 2, function(x) x/mean(x[1:16]))
hcm.lipo.relative.log <- log2(hcm.lipo.relative[,1:768])

lipo_name <- data.frame(ori = colnames(hcm.lipo.t), check = make.names(colnames(hcm.lipo.t)))
colnames(hcm.lipo.t) <- make.names(colnames(hcm.lipo.t))
hcm.lipo.t <- data.frame((hcm.lipo.t[17:365,]), check.names = F)
hcm.lipo.relative <- data.frame((hcm.lipo.relative[17:365,]), check.names = F)
hcm.lipo.relative.log <- data.frame((hcm.lipo.relative.log[17:365,]), check.names = F)


hcm.regres <- cbind(clinical, 
                    hcm.lipo.t[match(rownames(clinical), rownames(hcm.lipo.t)),])
                     
comb_name <- rbind (metab_name,lipo_name)
hcm.or <- hcm.regres
m = "NYHA..class.I.II..III.IV."
hcm.or[[m]] <- ifelse(hcm.or[[m]] == 1, 0, 1)
subgroup <- levels(factor(hcm.or[[m]]))
subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])

hcm.or$Gender.0M1F. <- as.factor(hcm.or$Gender.0M1F.)
hcm.or[,25:792] <- log10(hcm.or[,25:792])
hcm.or[,25:792] <- apply(hcm.or[,25:792], 2, scale)

ORRes.df <- data.frame()
for (x in colnames(hcm.or)[25:792]) {
  fmla <- as.formula(paste(m," ~ ", paste(c("Gender.0M1F.", "Age..year.", "BMI..m2.kg.", x), collapse= " + ")))
  ORRes <- glm(fmla, data = hcm.or, family = binomial())
  # unisum <- summary(ORRes)
  OR <- exp(coef(ORRes))[5]
  OR <- round(OR, 2)
  ci <- exp(confint(ORRes))[5,]
  ci <- round(ci, 2)
  cito <- paste0(ci[1], " - ", ci[2])
  
  ORRes.s <- data.frame(coef(summary(ORRes)))
  ORRes.s$BH <- p.adjust(ORRes.s$Pr...z.., method = "BH")
  ORRes.df <- rbind(ORRes.df, 
                    data.frame("Metab" = x, "NYHA" = subgroup1, "OR" = OR, 
                               "95%CI" = cito, "p.value" = ORRes.s$Pr...z..[5], "BH" = ORRes.s$BH[5]))
}

ORRes.df <- ORRes.df[ORRes.df$BH < 0.05, ]
ORRes.df$name <- lipo_name$ori[match(ORRes.df$Metab, lipo_name$check)]

########Fig 5f#####
hcm.regres <- cbind(clinical, 
                    hcm.metab.t[match(rownames(clinical), rownames(hcm.metab.t)),],
                    hcm.lipo.t[match(rownames(clinical), rownames(hcm.lipo.t)), ])

hcm.or <- hcm.regres
m = "NYHA..class.I.II..III.IV."
# y values must be 0 <= y <= 1
hcm.or[[m]] <- ifelse(hcm.or[[m]] == 1, 0, 1)
subgroup <- levels(factor(hcm.or[[m]]))
subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])

hcm.or$Gender.0M1F. <- as.factor(hcm.or$Gender.0M1F.)
hcm.or[,25:946] <- log10(hcm.or[,25:946])
hcm.or[,25:946] <- apply(hcm.or[,25:946], 2, scale)

fmla <- as.formula(paste(m," ~ ", paste(c(
                                            "Galactose.1.phosphate", "X6.Phosphogluconic.acid",
                                            "Ribose.5.phosphate", "Lactic.acid", 
                                            "G6P.F6P","Fructose.1.6.bisphosphate",
                                            "Inosine","Hypoxanthine","Uridine",
                                            "AMP","dCMP","GMP","Guanosine",
                                            "dTMP", "Xanthosine","Uracil",
                                            "Deoxyinosine","IMP","Xanthine",
                                            "X3.0.Carnitine","gamma.Glu.Cys",
                                            "GSSG","ADP.ribose","X16.0.Carnitine",
                                            "X18.0.Carnitine","X12.0.Carnitine",
                                            "X14.0.Carnitine","UDP.galactose",
                                            "TAG40.0.C16.0.","TAG40.0.C8.0.","TAG40.1.C8.0.",
                                            "TAG40.2.C14.0.","TAG40.2.C18.2.","TAG42.0.C18.0.",
                                            "TAG42.0.C8.0.","TAG42.1.C16.0.","TAG42.1.C8.0.",
                                            "TAG42.2.C16.0.","TAG42.2.C18.0.","TAG42.2.C18.2.",
                                            "TAG42.2.C8.0.","TAG44.0.C18.0.","TAG44.1.C18.0.",
                                            "TAG44.1.C8.0.","TAG44.2.C18.2.","TAG44.2.C8.0.",
                                            "TAG46.0.C18.0.","TAG46.2.C18.0.","TAG46.2.C18.2.",
                                            "TAG46.3.C18.0.","TAG46.4.C20.4.","TAG48.4.C18.0.",
                                            "TAG48.4.C20.4."), collapse= " + ")))

  ORRes <- glm(fmla, data = hcm.or, family = binomial())
  summary(ORRes)
  OR <- exp(coef(ORRes))
  OR <- round(OR, 2)
  ci <- exp(confint(ORRes))
  ci <- round(ci, 2)
  cito <- paste0(ci[1], " - ", ci[2])
  ORRes.s <- data.frame(coef(summary(ORRes)))
  ORRes.df <- data.frame()
  ORRes.df <- rbind(ORRes.df, 
                    data.frame( "NYHA" = subgroup1, "OR" = OR, 
                               "p.value" = ORRes.s$Pr...z.. ))
  

