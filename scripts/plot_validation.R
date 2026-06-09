if (!exists("output_dir")) output_dir <- "results"

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# Step 0: load data
pred <- read.csv(file.path(output_dir, "year2_predicted.csv"), row.names=1, check.names=FALSE)
obs  <- read.csv(file.path(output_dir, "year2_observed.csv"),  row.names=1, check.names=FALSE)
metrics <- read.csv(file.path(output_dir, "year2_genus_metrics.csv"))

cat("pred rownames:", rownames(pred), "\n")
cat("pred colnames[1:3]:", colnames(pred)[1:3], "\n")

month_labels <- c("3"="Mar", "8"="Aug", "10"="Oct", "12"="Dec")
new_rownames <- month_labels[rownames(pred)]
cat("new_rownames:", new_rownames, "class:", class(new_rownames), "\n")
cat("anyDuplicated:", anyDuplicated(new_rownames), "\n")

rownames(pred) <- new_rownames
cat("pred renamed OK\n")

rownames(obs) <- new_rownames
cat("obs renamed OK\n")

# Step 1: scatter
pv <- data.frame(
  Predicted = as.vector(as.matrix(pred)),
  Observed   = as.vector(as.matrix(obs)),
  Genus = rep(colnames(obs), each=4),
  Month = rep(rownames(pred), ncol(obs)),
  stringsAsFactors=FALSE
)
cat("pv rows:", nrow(pv), "\n")

p1 <- ggplot(pv, aes(x=Observed, y=Predicted)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray50") +
  geom_point(aes(color=Month), size=1.5, alpha=0.7) +
  scale_color_manual(values=c("Mar"="green4","Aug"="red","Oct"="orange","Dec"="blue")) +
  labs(title="Model Validation: Predicted vs Observed",
       x="Observed Relative Abundance", y="Predicted Relative Abundance") +
  theme_minimal() +
  theme(panel.background=element_rect(fill="white", color="black", linewidth=0.5),
        legend.position="bottom")
ggsave(file.path(output_dir, "validation_scatter.pdf"), p1, width=7, height=7)
cat("Scatter done\n")

# Step 2: time series
all_genera <- names(sort(colMeans(obs), decreasing=TRUE))
n_cols <- min(4, ceiling(sqrt(length(all_genera))))

plot_data <- data.frame(stringsAsFactors=FALSE)
for (g in all_genera) {
  plot_data <- rbind(plot_data, data.frame(
    Month = rep(rownames(pred), 2),
    Abundance = c(obs[, g], pred[, g]),
    Type = rep(c("Observed", "Predicted"), each=4),
    Genus = g, stringsAsFactors=FALSE
  ))
}
plot_data$Month <- factor(plot_data$Month, levels=c("Mar","Aug","Oct","Dec"))
p2 <- ggplot(plot_data, aes(x=Month, y=Abundance, group=Type, color=Type)) +
  geom_line(linewidth=0.8) + geom_point(size=1.5) +
  scale_color_manual(values=c("Observed"="black","Predicted"="red")) +
  facet_wrap(~Genus, scales="free_y", ncol=n_cols) +
  labs(title="Per-Genus Predicted vs Observed Trajectories (Top 12)",
       x="Month", y="Relative Abundance") +
  theme_minimal() +
  theme(panel.background=element_rect(fill="white",color="black",linewidth=0.3),
        strip.background=element_blank(), strip.text=element_text(size=7),
        panel.grid=element_blank(), legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1, size=6))
ggsave(file.path(output_dir, "validation_timeseries.pdf"), p2, width=11, height=10)
cat("Timeseries done\n")

# Step 3: RMSE barplot
metrics <- metrics[order(metrics$RMSE), ]
metrics$Genus <- factor(metrics$Genus, levels=metrics$Genus)
p3 <- ggplot(metrics, aes(x=RMSE, y=Genus)) +
  geom_col(fill="steelblue", width=0.7) +
  geom_text(aes(label=sprintf("%.4f", RMSE)), hjust=-0.1, size=2.5) +
  labs(title="Per-Genus RMSE", x="RMSE", y="") +
  theme_minimal() +
  theme(panel.background=element_rect(fill="white",color="black",linewidth=0.5),
        axis.text.y=element_text(size=7))
ggsave(file.path(output_dir, "validation_rmse.pdf"), p3, width=8, height=10)
cat("RMSE done\n")

# Step 4: heatmap
residuals <- obs - pred
res_long <- melt(as.matrix(residuals))
colnames(res_long) <- c("Month", "Genus", "Residual")
res_long$Month <- factor(res_long$Month, levels=c("Mar","Aug","Oct","Dec"))
p4 <- ggplot(res_long, aes(x=Month, y=Genus, fill=Residual)) +
  geom_tile(color="white", linewidth=0.5) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
  labs(title="Residual Heatmap (Observed - Predicted)", fill="Residual") +
  theme_minimal() +
  theme(axis.text.y=element_text(size=6), axis.text.x=element_text(size=8),
        panel.grid=element_blank())
ggsave(file.path(output_dir, "validation_heatmap.pdf"), p4, width=6, height=10)
cat("Heatmap done\n")

cat("ALL DONE\n")
