# Install and load necessary packages if not already installed
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(vcfR)
library(adegenet)
library(parallel)
library(ggplot2)
library(ape)
rm(list=ls())
# Set project directory and paths
proj.dir <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/variants/", sep="")

# Load metadata
meta.path <- paste(data.dir, "/Samples.Metadata.csv", sep="")
meta_data <- read.csv(meta.path)

# Load genlight object
rdata.geno.file <- paste(data.dir, "/geno.Rdata", sep="")
load(rdata.geno.file)


# Perform PCA
pca_obj <- glPca(genlight_obj, parallel = TRUE, n.cores = 4,,nf = 17,alleleAsUnit =TRUE)

snpgl.fclust <- find.clusters(genlight_obj, n.pca = 18, stat = "BIC", method = "ward",
                              max.n.clust = 8, 
                              criterion = "smoothNgoesup", glPca = pca_obj , n.iter = 10000) #w




scatter(pca_obj, posi="bottomright")

# Prepare data for plotting
pca_scores <- as.data.frame(pca_obj$scores)
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
pca_scores$Sample <- rownames(pca_scores)
pca_scores <- merge(pca_scores, meta_data, by.x = "Sample", by.y = "BioSample")

# Create PCA plot with ggplot2
pc12_location <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Location , shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "PCA of Genlight Object", x = "Principal Component 1", y = "Principal Component 2") +
  scale_shape_manual(values = c(17, 16), name = "Sex") + 
# Add shape legend for ploidy
  theme_minimal() + theme(legend.position = "top")

print(pc12_location )


pc13_location <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = Location ,  shape = as.factor(Sex))) +
  geom_point(size = 3) + 
  labs(title = "PCA of Genlight Object", x = "Principal Component 1", y = "Principal Component 3")+
  scale_shape_manual(values = c(17, 16), name = "Sex") + # Add shape legend for ploidy
  theme_minimal()   + theme(legend.position = "top")

print(pc13_location)


pc23_location <- ggplot(pca_scores, aes(x = PC2, y = PC3, color =  Location, shape = as.factor(Sex))) +
  geom_point(size = 3) + 
  labs(title = "PCA of Genlight Object", x = "Principal Component 2", y = "Principal Component 3")+
  scale_shape_manual(values = c(17, 16), name = "Sex") + # Add shape legend for ploidy
  theme_minimal()  + theme(legend.position = "top")

print(pc23_location)

pca_scores$grps <- snpgl.fclust$grp
pc12_cluster <- ggplot(pca_scores, aes(x = PC1, y = PC2, color =  grps, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "PCA of Genlight Object", x = "Principal Component 1", y = "Principal Component 2") +
  scale_shape_manual(values = c(17, 16), name = "Sex") + # Add shape legend for ploidy
  theme_minimal()  + theme(legend.position = "top")

print(pc12_cluster)
library(ggpubr)
ggarrange(pc12_location,pc12_cluster)

# Run the DAPC using disease status to group samples
location.dapc <- dapc(genlight_obj, pop=as.factor(unlist(genlight_obj$pop)), n.pca=12,glPca = pca_obj ,n.da = 3,
                     var.loadings=T, pca.info=T)

# Scatterplot of results
scatter.dapc(location.dapc, grp=as.factor(unlist(genlight_obj$pop)), legend=T)

ploidy.dapc <- dapc(genlight_obj, pop=as.factor(unlist(genlight_obj$ploidy)), n.pca=2,glPca = pca_obj, n.da=3,
                     var.loadings=T, pca.info=T,)
scatter.dapc(ploidy.dapc, grp=as.factor(unlist(genlight_obj$ploidy)), legend=T)


males.ids <- meta_data$BioSample[meta_data$Sex == "male"]
females.ids <- meta_data$BioSample[meta_data$Sex == "female"]

# Extract individuals by their names
male.gl<- genlight_obj[indNames(genlight_obj) %in% males.ids]


# Extract individuals by their names
female.gl <- genlight_obj[indNames(genlight_obj) %in% females.ids]
