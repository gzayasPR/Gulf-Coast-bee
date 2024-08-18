library(tidyverse)
rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
name <- args[2]
meta.path <- args[3]
cv.path <- args[4]
chosen.k <- args[5]
q.path <- args[6]
ids.path <- args[7]

setwd(out.dir)
# Load metadata
meta_data <- read.csv(meta.path)
print(head(meta_data ))
names(meta_data)[1] <- "BioSample"
cv_logs <- read.table(cv.path)[3:4]
names(cv_logs) <- c("K_text","CV")

cv_logs$K <- as.numeric(gsub("\\D", "", cv_logs$K_text))
head(cv_logs)
cv.plot.path <- paste(out.dir,"/cv.logs.plot.png",sep="")
png(cv.plot.path)
plot(cv_logs$K,cv_logs$CV)
dev.off()



all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())


ids <- read.table(ids.path)[2]
names(ids) <- "BioSample"
print(head(ids))
IDs <- merge(ids,meta_data,by="BioSample")

for (k in 1:nrow(cv_logs)){
  data <- read_delim(paste0(q.path,".",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")
  data$sample <- IDs$BioSample
  data$Location <- IDs$Location
  data$Sex <- IDs$Sex
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k,-Location,-Sex) -> data
  all_data <- rbind(all_data,data)
}


plots.bar <- all_data %>%
  arrange(desc(Sex)) %>%
  filter(k == chosen.k) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c(1:as.numeric(chosen.k))) + facet_wrap(~Sex)
sex_plot.name <- paste0(out.dir,"/",name,"_Sex_",chosen.k,".png")

ggsave(sex_plot.name,plots.bar)



plots.bar <- all_data %>%
  arrange(desc(Location)) %>%
  filter(k == chosen.k) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c(1:as.numeric(chosen.k))) + facet_wrap(~Location)
Location_plot.name <- paste0(out.dir,"/",name,"_Location_",chosen.k,".png")

ggsave(Location_plot.name,plots.bar)