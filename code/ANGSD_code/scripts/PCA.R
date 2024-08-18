# Load metadata
library(tidyverse)
library(stringr)

rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
name <- args[2]
meta.path <- args[3]
ids.path <- args[4]
setwd(out.dir)

# Load metadata
meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample"
ids <- read.table(ids.path)[1]
names(ids) <- "BioSample"
PCA <- read.table("PCA.cov")
e<-eigen(PCA )

eigen.values <- data.frame("values" = c(e$values), "Principal_Component" = c(1:length(e$values)))
print(eigen.values)

eigen.values.plot <- ggplot(eigen.values,aes(x=as.factor(Principal_Component),y=values)) +
    geom_bar(stat="identity") +
    xlab("Principal Component") + ylab("Eigenvalue")+
    geom_text(aes(label=format(round(values,2), nsmall = 2)), position=position_dodge(width=0.9), vjust=-0.25,size=3)+
  theme_bw() 

eigen.values_plot.name <- paste0(out.dir,"/",name,"_eigen_value.png")

ggsave(eigen.values_plot.name,eigen.values.plot)


eigen.vectors <- as.data.frame(e$vectors)
names(eigen.vectors) <- paste("PC",1:ncol(eigen.vectors),sep="")
pca.ids <- cbind(ids,eigen.vectors)
head(pca.ids )
all_data.merged <- merge(pca.ids,meta_data,by="BioSample")
plots.bar <- all_data.merged  %>%
  ggplot(.,aes(x=PC1,PC2,color=Location,shape=Sex)) + 
  geom_point() +
  xlab("PCA1") + ylab("PCA2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
PCA1_2_plot.name <- paste0(out.dir,"/",name,"PC1_2.png")

ggsave(PCA1_2_plot.name ,plots.bar)


plots.bar <- all_data.merged  %>%
  ggplot(.,aes(x=PC1,PC2,color=Depth.mean,shape=Sex)) + 
  geom_point() +
  xlab("PCA1") + ylab("PCA2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
PCA1_2_plot.name <- paste0(out.dir,"/",name,"PC1_2_depth.png")

ggsave(PCA1_2_plot.name ,plots.bar)
