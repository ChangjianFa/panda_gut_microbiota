#################### In/Out Degree Analysis
rm(list = ls())
setwd("D:/data/panda/results") 
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

ed="edge_matrix"
edge <- read.csv(paste0(ed,".csv"))

# Prepare node data
nodes <- unique(c(edge$source, edge$target))
data <- data.frame(id = nodes)

s <- table(edge$source)
t <- table(edge$target)

# In-degree positive, out-degree negative
data$In_degree <- as.numeric(t[as.character(data$id)])
data$Out_degree <- -as.numeric(s[as.character(data$id)])
data[is.na(data)] <- 0

# Sort by in-degree
data <- data %>% arrange(desc(In_degree))
data$id <- factor(data$id, levels = data$id)

# Convert to long format for plotting
bardata <- data %>%
  select(id, In_degree, Out_degree) %>%
  pivot_longer(cols = c(In_degree, Out_degree),
               names_to = "variable", values_to = "value")

# Plot
p <- ggplot(bardata, aes(x=id, y=value, fill=variable)) +
  geom_bar(stat='identity', width=0.8) +
  scale_fill_manual(values = c("In_degree"="#FFC60C", "Out_degree"="#00C23F")) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill="white", colour="black", linewidth=0.25),
    axis.line = element_line(colour="black", linewidth=0.25),
    axis.title = element_text(size=13, color="black"),
    axis.text = element_text(size=8, color="black"),
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1)
  ) +
  labs(x="Nodes", y="Edge Number") +
  geom_hline(yintercept=0, color="black") +  # Zero reference line
  scale_y_continuous(expand=c(0.1,0))

p

ggsave(
  file.path("D:/data/panda/results", paste0(ed, "_bar.pdf")),
  width = 8, 
  height = 4
)

