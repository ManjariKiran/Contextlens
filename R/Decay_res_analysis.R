library(ggplot2)
library(introdataviz)
library(forcats)
x <- read.table("/Users/mkiran/Work/DP/deseq_analysis1/try/file", sep="\t",header = T)
x1 <- data.frame(x)

a <- x1[which(x1$padj1 > 0.1),]
a$sig <- "unique_4-6"
b <- x1[which(x1$padj2 > 0.1),]
b$sig <- "unique_no_norm"
c <- x1[which((x1$padj1 < 0.1) & (x1$padj2 < 0.1)),]
c$sig <- "common"
data <- rbind(a,b,c)
head(data)
str(data)
data$sig

cbPalette <- c("#c2a5cf","#0571b0","#ca0020")
p <- ggplot(data, aes(x = lfc1, y = lfc2)) + geom_point(aes(colour = factor(sig)))+ 
    xlab("LFC (No normalization)") + ylab("LFC (4-6 radial distance normalization)") + 
    scale_color_manual(values=cbPalette) +
      theme_classic()

#       scale_y_continuous(limits = c(-1, 5)) + scale_x_continuous(limits = c(-1,5))
p
p + theme(
  axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),
  axis.title.x = element_text(color="darkgrey", size=16, face="bold"),
  axis.title.y = element_text(color="darkgrey", size=16, face="bold"),
  legend.position="right",legend.text = element_text(size=16),legend.title = element_text(size=16)
) + geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") + 
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey")


a <- x1[which((x1$padj1 > 0.1) & (x1$lfc2 > 0)),]
dim(a)
a$sig <- "unique_4-6"
b <- x1[which((x1$padj2 > 0.1) & (x1$lfc1 > 0)),]
dim(b)
b$sig <- "unique_no_norm"
c <- x1[which((x1$padj1 < 0.1) & (x1$padj2 < 0.1)),]
d <- c[which((c$lfc1 > 0) & (c$lfc2 > 0)),]
dim(d)
d$sig <- "common"
data1 <- rbind(b,d,a)
head(data1)
cbPalette <- c("#c2a5cf","#0571b0","#ca0020")
p <- ggplot(data1, aes(x = fct_inorder(sig), y= basemean2)) + 
  geom_boxplot(aes(colour = factor(sig))) + 
  xlab("Down") + ylab("count") + 
  scale_color_manual(values=cbPalette) +
  theme_classic()
p + theme(
  axis.text.y=element_text(size=16),axis.text.x=element_text(size=12),
  axis.title.x = element_text(color="darkgrey", size=16, face="bold"),
  axis.title.y = element_text(color="darkgrey", size=16, face="bold"),
  legend.position="right",legend.text = element_text(size=16),legend.title = element_text(size=16)
) 
