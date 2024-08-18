
# Load metadata
library(tidyverse)
library(stringr)
rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
name <- args[2]
meta.path <- args[3]
chosen.k <- args[4]
ids.path <- args[5]
reps <- args[6]
ks <- args[7]
print(args)
setwd(out.dir)

# Load metadata
meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample"

all_data <- tibble(Biosample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())


ids <- read.table(ids.path)[1]
names(ids) <- "BioSample"
print(head(ids))
IDs <- merge(ids,meta_data,by="BioSample")

log.files <-list.files(pattern = ".log$", full.names = T)

print(log.files)
bigData<-lapply(1:length(log.files), FUN = function(i) readLines(log.files[i]))
#this will pull out the line that starts with "b" from each file and return it as a list
foundset<-sapply(1:length(log.files), FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])

foundset

#now lets store it in a dataframe
#make a dataframe with an index 1:7, this corresponds to our K values
logs <- data.frame(K = rep(1:as.numeric(ks), each=as.numeric(reps)))


#add to it our likelihood values

logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))


#and now we can calculate our delta K and probability

log_results <- data.frame(K = 1:ks, test = tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x))))
names(log_results) <- c("K","test")

       
print(log_results)
cv.plot.path <- paste(out.dir,"/k.logs.plot.png",sep="")
png(cv.plot.path)
plot(log_results$K,log_results$test)
dev.off()


for (k in 1:as.numeric(ks)){
  data <- read_delim(paste0("myoutfiles",".k",k,".rep1.qopt"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")
  data$BioSample <- ids$BioSample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -BioSample,-k) -> data
  all_data <- rbind(all_data,data)
}

                 
all_data.merged <- merge(all_data,meta_data,by="BioSample")
                                                
 print(head(all_data.merged))                
plots.bar <- all_data.merged  %>%
  arrange(desc(Sex)) %>%
  filter(k == chosen.k) %>%
  ggplot(.,aes(x=BioSample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c(1:as.numeric(chosen.k))) + facet_wrap(~Sex, scales = "free")
sex_plot.name <- paste0(out.dir,"/",name,"_Sex_",chosen.k,".png")

ggsave(sex_plot.name,plots.bar)



plots.bar <- all_data.merged  %>%
  arrange(desc(Location)) %>%
  filter(k == chosen.k) %>%
  ggplot(.,aes(x=BioSample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c(1:as.numeric(chosen.k))) + facet_wrap(~Location, scales = "free")
Location_plot.name <- paste0(out.dir,"/",name,"_Location_",chosen.k,".png")

ggsave(Location_plot.name,plots.bar)
                                                  
