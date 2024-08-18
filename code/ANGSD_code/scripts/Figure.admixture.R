library(tidyverse)
library(stringr)
library(plyr)

rm(list=ls())
args <- commandArgs(TRUE)

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

all_data <- tibble(BioSample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

ids <- read.table(ids.path)[1]
names(ids) <- "BioSample"
print(head(ids))
IDs <- merge(ids, meta_data, by="BioSample")

for (k in 1:as.numeric(ks)) {
  data <- read_delim(paste0("myoutfiles", ".k", k, ".rep1.qopt"),
                     col_names = paste0("Q", seq(1:k)),
                     delim=" ")
  data$BioSample <- ids$BioSample
  data$k <- k
  
  # This step converts from wide to long.
  data %>% gather(Q, value, -BioSample, -k) -> data
  all_data <- rbind(all_data, data)
}

all_data.merged <- merge(all_data, meta_data, by="BioSample")
print(head(all_data.merged))

plots.bar <- all_data.merged %>%
  arrange(desc(Location)) %>%
  filter(k == chosen.k)

# Rename levels of the Location column
plots.bar <- plots.bar %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Remove NAs
plots.bar <- na.omit(plots.bar)

# Custom sorting within each location and creating an index column
# Sort by Q1, then Q2, then Q3, and so on
if (chosen.k == 2) {
  plots.bar <- plots.bar %>%
    group_by(Location, BioSample) %>%
    mutate(max_Q = max(as.numeric(Q))) %>%
    spread(Q, value) %>%
    arrange(Location, desc(Q1)) %>%
    gather(Q, value, starts_with("Q")) %>%
    arrange(Location, max_Q, desc(Q)) %>%
    mutate(index = as.integer(factor(BioSample, levels = unique(BioSample))))
} else if (chosen.k == 3) {
  plots.bar <- plots.bar %>%
    group_by(Location, BioSample) %>%
    mutate(max_Q = max(as.numeric(Q))) %>%
    spread(Q, value) %>%
    arrange(Location, desc(Q1), desc(Q2)) %>%
    gather(Q, value, starts_with("Q")) %>%
    arrange(Location, max_Q, desc(Q)) %>%
    mutate(index = as.integer(factor(BioSample, levels = unique(BioSample))))
} else if (chosen.k == 4) {
  plots.bar <- plots.bar %>%
    group_by(Location, BioSample) %>%
    mutate(max_Q = max(as.numeric(Q))) %>%
    spread(Q, value) %>%
    arrange(Location, desc(Q1), desc(Q2), desc(Q3)) %>%
    gather(Q, value, starts_with("Q")) %>%
    arrange(Location, max_Q, desc(Q)) %>%
    mutate(index = as.integer(factor(BioSample, levels = unique(BioSample))))
} else if (chosen.k == 5) {
  plots.bar <- plots.bar %>%
    group_by(Location, BioSample) %>%
    mutate(max_Q = max(as.numeric(Q))) %>%
    spread(Q, value) %>%
    arrange(Location, desc(Q1), desc(Q2), desc(Q3), desc(Q4)) %>%
    gather(Q, value, starts_with("Q")) %>%
    arrange(Location, max_Q, desc(Q)) %>%
    mutate(index = as.integer(factor(BioSample, levels = unique(BioSample))))
} else {
  plots.bar <- plots.bar %>%
    group_by(Location, BioSample) %>%
    mutate(max_Q = max(as.numeric(Q))) %>%
    spread(Q, value) %>%
    arrange(Location, desc(Q1), desc(Q2), desc(Q3), desc(Q4), desc(Q5)) %>%
    gather(Q, value, starts_with("Q")) %>%
    arrange(Location, max_Q, desc(Q)) %>%
    mutate(index = as.integer(factor(BioSample, levels = unique(BioSample))))
}

# Print the renamed index
print(plots.bar)

# Plot the bar graph
plot.bar <- ggplot(plots.bar, aes(x=index, y=value, fill=factor(Q))) + 
  geom_bar(stat="identity", position="stack", width=0.95) +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # Remove x-ticks
    axis.ticks.x = element_blank(),  # Remove x-ticks
    strip.background = element_blank(),  # Remove background of facet labels
    strip.placement = "outside",  # Place facet labels outside
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.5, "lines"),  # Remove space between facets
    panel.border = element_blank(),  # Remove facet boxes
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  ) +
  scale_fill_brewer(palette="Set1", name="K",
                    labels=c(1:as.numeric(chosen.k))) + 
  facet_wrap(~Location, scales = "free_x", strip.position = "bottom")  # Display facets at the bottom

Location_plot.name <- paste0(out.dir, "/", name, "_Figure_", chosen.k, ".png")
ggsave(Location_plot.name, plot.bar, width = 11, height= 6)
