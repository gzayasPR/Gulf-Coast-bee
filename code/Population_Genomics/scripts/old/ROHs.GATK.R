library(detectRUNS)

proj.dir <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/"
setwd(proj.dir)
results.dir <- paste(proj.dir,"results/Population_Genomics/ROH/",sep="")
setwd(results.dir)

genotypeFilePath <- "females2.ped"
mapFilePath <- "females2.map"
slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 10, 
  threshold = 0.05,
  minSNP = 15, 
  ROHet = FALSE, 
  maxOppWindow = 3, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 10000, 
  minDensity = 1/10^3, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 
consecutiveRuns <- consecutiveRUNS.run(
  genotypeFile =genotypeFilePath,
  mapFile = mapFilePath,
  minSNP = 10,
  ROHet = FALSE,
  maxGap = 10^6,
  minLengthBps = 10000,
  maxOppRun = 1,
  maxMissRun = 1
)

plot_Runs(runs = slidingRuns,savePlot=T)
summaryList <- summaryRuns(
  runs = slidingRuns, mapFile = mapFilePath, genotypeFile = genotypeFilePath,  snpInRuns = TRUE)

snpsinruns <- read.table("plink.hom.summary",header=T)
snpsinruns$index <- nrow(snpsinruns)
snpsinruns$PERCENTAGE <- snpsinruns$UNAFF/9

png("manhattan.png")
plot(snpsinruns$index,snpsinruns$PERCENTAGE)
dev.off()
