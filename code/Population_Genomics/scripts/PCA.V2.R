library(vcfR)
library(adegenet)
library(parallel)
library(ggplot2)
library(ape)
library(ggpubr)

rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
meta.path <- args[2]
vcf.path <- args[3]

meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample" 
setwd(out.dir)
# Convert filtered VCF to genlight object
pruned_vcf <- read.vcfR(vcf.path)
genlight_obj <- vcfR2genlight(pruned_vcf)
# Ensure the order of samples matches between genlight object and metadata
id.order <- genlight_obj@ind.names
meta_data <- meta_data[match(id.order, meta_data$BioSample),]

# Assign ploidy based on metadata
meta_data$ploidy <- ifelse(meta_data$Sex == "female", yes = as.integer(2), no = as.integer(1))
genlight_obj@ploidy <- meta_data$ploidy
# Assign populations based on metadata
genlight_obj@pop <- as.factor(meta_data$Location)
print(genlight_obj)

n <- length(genlight_obj@ind.names)
# Perform PCA
pca_obj <- glPca(genlight_obj, parallel = TRUE, n.cores = 4,nf = n-1, alleleAsUnit = TRUE)

# Prepare data for plotting
pca_scores <- as.data.frame(pca_obj$scores)
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
pca_scores$Sample <- rownames(pca_scores)
pca_scores <- merge(pca_scores, meta_data, by.x = "Sample", by.y = "BioSample")

write.csv(pca_scores,"PCA.csv")
# Create PCA plot for all individuals
pc12_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Location, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Entire Population", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_all)
# Subset males and females
males.ids <- meta_data$BioSample[meta_data$Sex == "male"]
females.ids <- meta_data$BioSample[meta_data$Sex == "female"]
# Extract individuals by their names
male_gl <- genlight_obj[indNames(genlight_obj) %in% males.ids]


# Extract individuals by their names
female_gl <- genlight_obj[indNames(genlight_obj) %in% females.ids]

# Perform PCA for males
n.m  <- length(male_gl@ind.names)
pca_obj_males <- glPca(male_gl, parallel = TRUE, n.cores = 4,nf = n.m -1, alleleAsUnit = TRUE)

pca_scores_males <- as.data.frame(pca_obj_males$scores)
colnames(pca_scores_males) <- paste0("PC", 1:ncol(pca_scores_males))
pca_scores_males$Sample <- rownames(pca_scores_males)
pca_scores_males <- merge(pca_scores_males, meta_data, by.x = "Sample", by.y = "BioSample")

# Create PCA plot for males
pc12_males <- ggplot(pca_scores_males, aes(x = PC1, y = PC2, color = Location, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Males", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_males)
# Perform PCA for females
n.f  <- length(female_gl@ind.names)
pca_obj_females <- glPca(female_gl, parallel = TRUE, n.cores = 4, alleleAsUnit = TRUE,nf= n.f -1)
pca_scores_females <- as.data.frame(pca_obj_females$scores)
colnames(pca_scores_females) <- paste0("PC", 1:ncol(pca_scores_females))
pca_scores_females$Sample <- rownames(pca_scores_females)
pca_scores_females <- merge(pca_scores_females, meta_data, by.x = "Sample", by.y = "BioSample")

# Create PCA plot for females
pc12_females <- ggplot(pca_scores_females, aes(x = PC1, y = PC2, color = Location, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Females", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_females)
# Arrange plots
combined_plot <- ggarrange(
  pc12_all,
  ggarrange(pc12_females, pc12_males,ncol = 2, labels = c("B", "C"),legend = "none"),
  nrow = 2,
  labels = "A",legend="top")

# Print combined plot
print(combined_plot)
setwd(out.dir)

ggsave("PCA.png",combined_plot,width = 11,height=8)


pc12_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Mapped.reads, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Entire Population", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_all)

# Create PCA plot for males
pc12_males <- ggplot(pca_scores_males, aes(x = PC1, y = PC2, color = Mapped.reads, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Males", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))

# Create PCA plot for females
pc12_females <- ggplot(pca_scores_females, aes(x = PC1, y = PC2, color = Mapped.reads, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Females", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_females)
# Arrange plots
combined_plot <- ggarrange(
  pc12_all,
  ggarrange(pc12_females, pc12_males,ncol = 2, labels = c("B", "C"),legend = "none"),
  nrow = 2,
  labels = "A",legend="top")

ggsave("PCAvsMapped.png",combined_plot,width = 11,height=8)

pc12_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Depth.mean, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Entire Population", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_all)
# Create PCA plot for males
pc12_males <- ggplot(pca_scores_males, aes(x = PC1, y = PC2, color = Depth.mean, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Males", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(16), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))

# Create PCA plot for females
pc12_females <- ggplot(pca_scores_females, aes(x = PC1, y = PC2, color =Depth.mean, shape = Sex)) +
  geom_point(size = 3) +
  labs(title = "Females", x = "PC1", y = "PC2") +
  scale_shape_manual(values = c(17), name = "Sex") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face = "bold",size = 20))
print(pc12_females)
# Arrange plots
combined_plot <- ggarrange(
  pc12_all,
  ggarrange(pc12_females, pc12_males,ncol = 2, labels = c("B", "C"),legend = "none"),
  nrow = 2,
  labels = "A",legend="top")

ggsave("PCAvsCoverage.png",combined_plot,width = 11,height=8)
