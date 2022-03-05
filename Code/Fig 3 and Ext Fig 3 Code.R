library(xlsx)
library(xlsxjars)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(ggdendro)
library(mixOmics)
library(viridis)
library(hrbrthemes)
library(dplyr)
options(stringsAsFactors = F)


hcm.lipo <- read.xlsx("Source Data Extended Data Fig. 3.xlsx", 
                      sheetIndex = 1, startRow = 1, as.data.frame = T, header = T, check.names = F)
rownames(hcm.lipo) <- hcm.lipo$Metabolites
hcm.lipo <- hcm.lipo[,-1]
hcm.lipo.t <- t(hcm.lipo)
hcm.lipo.group <- data.frame(hcm.lipo.t,check.names = F)
hcm.lipo.group$group <-c(rep(c("Ctrl","HCM"),c(16,349)))


hcm.lipo.relative <- apply(hcm.lipo.t, 2, function(x) x/mean(x[1:16]))
hcm.lipo.relative.log <- log2(hcm.lipo.relative[,1:768])

test_results <- data.frame(Metabolites = colnames(hcm.lipo.t)[1:768])

# wilcox.test
test_results$wilcox.test <- apply(hcm.lipo.relative.log[,1:768],2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ hcm.lipo.group$group, 
                             data = hcm.lipo.relative.log)[3]))
# adjust p value using BH method
test_results$wilcox.test_BH <- p.adjust(test_results$wilcox.test, method = "fdr")

#fold-change
test_results$FC <- apply(hcm.lipo.group[,1:768], 2, 
                  function(x) mean(x[which(hcm.lipo.group$group == 'HCM')])/mean(x[which(hcm.lipo.group$group == 'Ctrl')]))
test_results$LOG2FC <- log2(test_results$FC)

# differential altered lipoolites FC 1.5 
test_results_diff <- test_results[test_results$wilcox.test_BH < 0.05 & 
                                    abs(test_results$LOG2FC) > 0.584963, ]

test_results$change <- ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC > 1.5, "Up", 
                                  ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC < 0.67, 
                                         "Down","Not changed"))

#####Fig 3a#####
lipo.bubble <- read.xlsx("Source Data Fig. 3.xlsx",
                         sheetIndex = 1, startRow = 1, header = T, as.data.frame = T,check.names = F)

lipo.bubble %>%
  arrange(desc(`FDR`)) %>%
  ggplot(aes(x=number, y=`Log2 (FC) relative to Ctrl`, size=`FDR`, color=`lipid species`)) +
  geom_point(alpha=0.4) +
  scale_size(range = c(0.1, 3.5), name="Population (M)" +
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
  geom_hline(yintercept = c(-0.58,0.58),lty=4,col="black",lwd=0.4)

#####Fig 3b Ext. Fig 3e#####
lipo.total <-data.frame(hcm.lipo.relative, check.names = F)
lipo.total$TAG_total <- apply(lipo.total[,grep("TAG", colnames(lipo.total))], 1, mean)
lipo.total$LPC_total <- apply(lipo.total[,grep("LPC", colnames(lipo.total))], 1, mean)
lipo.total$LPE_total <- apply(lipo.total[,grep("LPE", colnames(lipo.total))], 1, mean)
lipo.total$PC_total <- apply(lipo.total[,grep("^PC", colnames(lipo.total))], 1, mean)
lipo.total$PE_total <- apply(lipo.total[,grep("^PE", colnames(lipo.total))], 1, mean)
lipo.total$PG_total <- apply(lipo.total[,grep("^PG", colnames(lipo.total))], 1, mean)
lipo.total$PI_total <- apply(lipo.total[,grep("^PI", colnames(lipo.total))], 1, mean)
lipo.total$DAG_total <- apply(lipo.total[,grep("^DAG", colnames(lipo.total))], 1, mean)
lipo.total$PS_total <- apply(lipo.total[,grep("^PS", colnames(lipo.total))], 1, mean)
lipo.total$FFA_total <- apply(lipo.total[,grep("^FFA", colnames(lipo.total))], 1, mean)
lipo.total$SM_total <- apply(lipo.total[,grep("^SM", colnames(lipo.total))], 1, mean)
lipo.total$Cer_total <- apply(lipo.total[,grep("^Cer", colnames(lipo.total))], 1, mean)
lipo.total$HexCer_total <- apply(lipo.total[,grep("^HexCer", colnames(lipo.total))], 1, mean)

#Total
carn <- lipo.total[grep("total", colnames(lipo.total))]
carn <- apply(carn, 2, function(x) x/mean(x[1:16]))
carn <-log2(carn)
colnames(carn) <- gsub("_total","",colnames(carn))

group = c("Ctrl", "HCM")
carn.box <- data.frame()
for (i in c(1:ncol(carn))) {
  carn.box <- rbind(carn.box, 
                    data.frame(value = carn[,i], 
                               lipo = rep(colnames(carn)[i], 365), 
                               group = rep(group, c(16,349))))
}
carn.box$group <- factor(carn.box$group, level = c("Ctrl", "HCM"))

my_comparisons <- list(group)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

ggplot(carn.box, aes(x = group, y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, size = 0.3) + 
  geom_boxplot(aes(fill = group), width = 0.7, outlier.shape = NA, size = 0, fatten = 1,outlier.colour = NA) +
  facet_wrap(~ lipo, scales = "free", strip.position = "top", ncol = 4) + 
  geom_point(position="jitter",col=1,pch=16,cex=0.02)+
  # geom_point(aes(color = group), size = 1.6, shape = 1) +
  scale_fill_manual(values = c("HCM" = "#f8766d", "Ctrl" = "#00bfc4")) + 
  theme_bw() + # ggtheme
  theme(aspect.ratio = 1.5, 
        panel.border = element_rect(colour = "black", size= 0.5),
        panel.grid = element_line(size = 0.5)) +
  labs(x = "", y = "", color = "") +
  theme(plot.title = element_text(size = 8, vjust = 6, hjust = 0.5, face = "plain"), 
        strip.text = element_text(size = 8, face = "plain", hjust = 0.5, vjust = 0.9), 
        strip.background = element_rect(fill = NA, colour = NA))  +
  theme(axis.title.x = element_text(size = 8, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 8, vjust = 2, face ="plain")) +
  theme(axis.text.x = element_text(size = 8, face ="plain", colour = "black"), 
        axis.text.y = element_text(size = 8, face ="plain", colour = "black")) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm")) +
  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     symnum.args = symnum.args, hide.ns = F, label = "p.format")


#####Fig 3d Ext. Fig 3f#####
hcm.lipo.t <-data.frame(hcm.lipo.t, check.names = F)
hcm.TAG.bond <- hcm.lipo.t[grep("TAG", colnames(hcm.lipo.t))]
colnames(hcm.TAG.bond) <- substr(colnames(hcm.TAG.bond),7,7)

hcm.TAG.group <- data.frame(hcm.TAG.bond,check.names = F)
hcm.TAG.group$group <-c(rep(c("Ctrl","HCM"),c(16,349)))


hcm.TAG.relative <- apply(hcm.TAG.bond, 2, function(x) x/mean(x[1:16]))
hcm.TAG.relative.log <- log2(hcm.TAG.relative[,1:453])

test_results <- data.frame(Metabolites = colnames(hcm.TAG.bond)[1:453])

# wilcox.test
test_results$wilcox.test <- apply(hcm.TAG.relative.log[,1:453],2,
                                  function(x) unlist(wilcox.test(as.numeric(x) ~ hcm.TAG.group$group, 
                                                                 data = hcm.TAG.relative.log)[3]))
# adjust p value using BH method
test_results$wilcox.test_BH <- p.adjust(test_results$wilcox.test, method = "fdr")

#fold-change
test_results$FC <- apply(hcm.TAG.group[,1:453], 2, 
                         function(x) mean(x[which(hcm.TAG.group$group == 'HCM')])/mean(x[which(hcm.TAG.group$group == 'Ctrl')]))
test_results$LOG2FC <- log2(test_results$FC)

test_results$change <- ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC > 1.5, "Up", 
                              ifelse(test_results$wilcox.test_BH < 0.05 & test_results$FC < 0.67, 
                                     "Down","Not changed"))


TAG.bond.mean <- data.frame(matrix(ncol=9, nrow = 365))
rownames(TAG.bond.mean) <- rownames(hcm.TAG.bond)[1:365]
colnames(TAG.bond.mean) <-c(0:8)
TAG.bond.mean[,1] <- apply(hcm.TAG.bond[,grep("0", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,2] <- apply(hcm.TAG.bond[,grep("1", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,3] <- apply(hcm.TAG.bond[,grep("2", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,4] <- apply(hcm.TAG.bond[,grep("3", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,5] <- apply(hcm.TAG.bond[,grep("4", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,6] <- apply(hcm.TAG.bond[,grep("5", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,7] <- apply(hcm.TAG.bond[,grep("6", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,8] <- apply(hcm.TAG.bond[,grep("7", colnames(hcm.TAG.bond))], 1, mean)
TAG.bond.mean[,9] <- apply(hcm.TAG.bond[,grep("8", colnames(hcm.TAG.bond))], 1, mean)

TAG.bond.mean$group <-c(rep(c("Ctrl","HCM"),c(16,349)))

TAG.bond.FC <- data.frame(bondnumber = colnames(TAG.bond.mean)[1:9])
TAG.bond.FC$FC <- apply(TAG.bond.mean[,1:9], 2, 
                         function(x) mean(x[which(TAG.bond.mean$group == 'HCM')])
                                    /mean(x[which(TAG.bond.mean$group == 'Ctrl')]))
TAG.bond.FC$wilcox.test <- apply(TAG.bond.mean[,1:9],2,
                                  function(x) unlist(wilcox.test(as.numeric(x) ~ TAG.bond.mean$group, 
                                                                 data = TAG.bond.mean)[3]))


hcm.TAG.carbon  <- hcm.lipo.t[grep("TAG", colnames(hcm.lipo.t))]
colnames(hcm.TAG.carbon) <- substr(colnames(hcm.TAG.carbon),4,5)
TAG.carbon.mean <- data.frame(matrix(ncol=11, nrow = 365))
rownames(TAG.carbon.mean) <- rownames(hcm.TAG.carbon)[1:365]
colnames(TAG.carbon.mean) <-c("36","38","40","42","44","46","48","50","52","54","56")
TAG.carbon.mean[,1] <- apply(hcm.TAG.carbon[,grep("36", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,2] <- apply(hcm.TAG.carbon[,grep("38", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,3] <- apply(hcm.TAG.carbon[,grep("40", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,4] <- apply(hcm.TAG.carbon[,grep("42", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,5] <- apply(hcm.TAG.carbon[,grep("44", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,6] <- apply(hcm.TAG.carbon[,grep("46", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,7] <- apply(hcm.TAG.carbon[,grep("48", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,8] <- apply(hcm.TAG.carbon[,grep("50", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,9] <- apply(hcm.TAG.carbon[,grep("52", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,10] <- apply(hcm.TAG.carbon[,grep("54", colnames(hcm.TAG.carbon))], 1, mean)
TAG.carbon.mean[,11] <- apply(hcm.TAG.carbon[,grep("56", colnames(hcm.TAG.carbon))], 1, mean)

TAG.carbon.mean$group <-c(rep(c("Ctrl","HCM"),c(16,349)))
TAG.carbon.FC <- data.frame(carbonnumber = colnames(TAG.carbon.mean)[1:11])
TAG.carbon.FC$FC <- apply(TAG.carbon.mean[,1:11], 2, 
                        function(x) mean(x[which(TAG.carbon.mean$group == 'HCM')])
                        /mean(x[which(TAG.carbon.mean$group == 'Ctrl')]))
TAG.carbon.FC$log2FC <-log2(TAG.carbon.FC$FC)


#####Ext. Fig 3a#####

Y <-c(rep(c("Ctrl","HCM"),c(16,349)))
df_plsda <- plsda(hcm.lipo.relative.log,Y, ncomp = 4)
plsda_result_eig <- {df_plsda$explained_variance$X}[1:2]
sample_site <- data.frame(df_plsda$variates)[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('plsda1', 'plsda2')
sample_site$group <- c(rep(c("Ctrl","HCM"),c(16,349)))

ggplot(sample_site,aes(x=plsda1,y=plsda2,fill=group))+ 
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_point(size = 1, alpha = 1, shape = 21, color = "black")+  
  stat_ellipse(level = 0.95, show.legend = F, linetype = 2)+
  labs(x = paste('PLS-DA axis1 ( ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), 
       y = paste('PLS-DA axis2 ( ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))+
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

set.seed(365) 
perf.plsda.srbct <- perf(df_plsda, validation = "Mfold", folds = 5, 
                         progressBar = FALSE, auc = TRUE, nrepeat = 10)
plot(perf.plsda.srbct, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(df_plsda, roc.comp = 3)

######Ext. Fig 3d-e#####
test_results <- read.xlsx("Source Data Extended Data Fig. 3.xlsx",
                          sheetIndex = 2, startRow = 1, header = T, as.data.frame = T,check.names = F)

ggplot(test_results, aes(x = test_results$LOG2FC, y = test_results$log10wilcox.test_BH, colour=test_results$change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("#546de5", "grey","#ff4757"))+
  geom_text_repel(data = test_results, aes(x = test_results$LOG2FC, y = test_results$log10wilcox.test_BH, label = label),
                  size = 1.5,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.4, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)+
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +

  labs(x="Log2(Fold change, HCM/Ctrl)", y="-Log10 (FDR)", color = "")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  xlim(c(-4,5))+
  ylim(c(0,8))+
  theme(axis.title.y = element_text(size = 6, face = "plain"),
        axis.title.x = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        axis.text.y = element_text(size = 6, face = "plain",))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right") +
  theme(legend.text = element_text(size = 6), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.3, "cm"))


