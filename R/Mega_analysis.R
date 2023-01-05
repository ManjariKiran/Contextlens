library(ggplot2)
library(introdataviz)
library(forcats)
x <- read.table("MEGA/MEGA_file", sep="\t",header = T)
x1 <- data.frame(x)
x2 <- x1[complete.cases(x1), ]
dim(x1)
x2
a <- x2[which((x2$padj1 < 0.05) & (abs(x2$log2FoldChange1) > 0.585)),]
b <- x2[which((x2$padj < 0.05) & (abs(x2$log2FoldChange) > 0.585)),]
a1 <- which((x2$padj1 < 0.05) & (abs(x2$log2FoldChange1) > 0.585))
b1 <- which((x2$padj < 0.05) & (abs(x2$log2FoldChange) > 0.585))
common <- intersect(a1,b1)
c <- x2[(common),]
a <- a[!(row.names(a) %in% common),]
b <- b[!(row.names(b) %in% common),]
dim(a)
dim(b)
dim(c)
a$sig <- "unique_8-10"
b$sig <- "unique_no_norm"
c$sig <- "common"
data <- rbind(a,b,c)
d <- as.numeric(rownames(data))
d <- x2[d,]
d$sig <- "base"
data <- rbind(a,b,c,d)
head(data)
str(data)
data$sig
c
cbPalette <- c("lightgrey","#c2a5cf","#0571b0","#ca0020")
p <- ggplot(data, aes(x = log2FoldChange, y = log2FoldChange1)) + geom_point(aes(colour = factor(sig)))+
  xlab("LFC (No normalization)") + ylab("LFC (8-10 radial distance normalization)") +
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


a <- x2[which((x2$padj1 < 0.05) & ((x2$log2FoldChange1) < -0.585)),]
b <- x2[which((x2$padj < 0.05) & ((x2$log2FoldChange) < -0.585)),]
a1 <- which((x2$padj1 < 0.05) & ((x2$log2FoldChange1) < -0.585))
b1 <- which((x2$padj < 0.05) & ((x2$log2FoldChange) < -0.585))
common <- intersect(a1,b1)
c <- x2[(common),]
a <- a[!(row.names(a) %in% common),]
b <- b[!(row.names(b) %in% common),]
dim(a)
dim(b)
dim(c)
a$sig <- "unique_8-10"
b$sig <- "unique_no_norm"
c$sig <- "common"
data1 <- rbind(b,c,a)
head(data1)
cbPalette <- c("#c2a5cf","#0571b0","#ca0020")
p <- ggplot(data1, aes(x = fct_inorder(sig), y= baseMean1)) +
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

a[which(a[order(a$padj1),]$padj > 0.1),]
