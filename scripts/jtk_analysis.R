# jtk_analysis.R — JTK_CYCLE rhythm analysis
rm(list = ls())
setwd("D:/data/panda")

library(MetaCycle)

cat("===== JTK_CYCLE Rhythm Analysis =====\n")

# ============================================================
# 1. Load effect curves, convert to wide-format MetaCycle input
# ============================================================
curves <- read.csv("results/effect_curves.csv", stringsAsFactors = FALSE)
cat("Loaded:", nrow(curves), "rows (", length(unique(curves$Genus)),
    "genera × 4 components × 12 months)\n")

curves_long <- reshape(curves,
  varying = c("Total", "Independent", "Bacterial", "External"),
  v.names = "Value", times = c("Total", "Independent", "Bacterial", "External"),
  timevar = "Component", direction = "long")
curves_long$id <- NULL
rownames(curves_long) <- NULL

curves_long$ID <- paste0(curves_long$Genus, "__", curves_long$Component)
curves_wide <- reshape(curves_long[, c("Month", "ID", "Value")],
  idvar = "Month", timevar = "ID", direction = "wide")
colnames(curves_wide) <- sub("^Value\\.", "", colnames(curves_wide))
rownames(curves_wide) <- curves_wide$Month
curves_wide$Month <- NULL

write.csv(cbind(CycID = colnames(curves_wide), t(curves_wide)),
  "results/jtk_input.csv", row.names = FALSE)

# ============================================================
# 2. Run JTK_CYCLE
# ============================================================
cat("\nRunning JTK_CYCLE (period=12, interval=1)...\n")
meta2d(infile = "results/jtk_input.csv",
  filestyle = "csv", outdir = "results",
  timepoints = 1:12, minper = 12, maxper = 12,
  cycMethod = "JTK", outIntegration = "noIntegration")

# ============================================================
# 3. Read and organize results
# ============================================================
jtk_raw <- read.csv("results/JTKresult_jtk_input.csv")
colnames(jtk_raw)[colnames(jtk_raw) == "ADJ.P"] <- "JTK_P_raw"
colnames(jtk_raw)[colnames(jtk_raw) == "BH.Q"]  <- "JTK_BH_Q"
jtk_raw$Genus     <- sub("__.*", "", jtk_raw$CycID)
jtk_raw$Component <- sub(".*__", "", jtk_raw$CycID)
jtk_raw$Significant <- jtk_raw$JTK_P_raw < 0.05

df <- jtk_raw[, c("CycID", "Genus", "Component", "JTK_P_raw", "JTK_BH_Q",
                  "PER", "LAG", "AMP", "Significant")]
df$PER <- as.numeric(df$PER)
df <- df[order(df$Component, df$JTK_P_raw), ]

# ============================================================
# 4. Wide-format summary: Genus x Component -> P value
# ============================================================
summary_wide <- reshape(df[, c("Genus", "Component", "JTK_P_raw")],
  idvar = "Genus", timevar = "Component", direction = "wide")
colnames(summary_wide) <- sub("^JTK_P_raw\\.", "", colnames(summary_wide))
summary_wide$Total_Sig       <- summary_wide$Total       < 0.05
summary_wide$Independent_Sig <- summary_wide$Independent < 0.05
summary_wide$Bacterial_Sig   <- summary_wide$Bacterial   < 0.05
summary_wide$External_Sig    <- summary_wide$External    < 0.05
summary_wide <- summary_wide[order(summary_wide$Total), ]

write.csv(df, "results/jtk_results.csv", row.names = FALSE)
write.csv(summary_wide, "results/jtk_summary.csv", row.names = FALSE)

# ============================================================
# 5. Report
# ============================================================
cat(sprintf("\n========== JTK Results Summary ==========\n"))
cat(sprintf("Curves tested: %d\n", nrow(df)))
cat(sprintf("P < 0.05:      %d (%.1f%%)\n\n",
  sum(df$Significant), sum(df$Significant) / nrow(df) * 100))

comps <- c("Total", "Independent", "Bacterial", "External")
for (c in comps) {
  sub <- df[df$Component == c, ]
  sig <- sum(sub$Significant)
  cat(sprintf("  %-15s  sig=%2d/%d", c, sig, nrow(sub)))
  if (sig > 0) {
    sig_names <- sub$Genus[sub$Significant]
    cat(sprintf("  (%s)", paste(sig_names, collapse = ", ")))
  }
  cat("\n")
}

# Total significant and non-self-driven pattern
cat("\n--- Total sig, Independent NOT sig, Bacterial sig ---\n")
for (g in unique(df$Genus)) {
  gd <- df[df$Genus == g, ]
  t_s <- gd$Significant[gd$Component == "Total"]
  i_s <- gd$Significant[gd$Component == "Independent"]
  b_s <- gd$Significant[gd$Component == "Bacterial"]
  e_s <- gd$Significant[gd$Component == "External"]
  if (t_s) {
    driver <- c()
    if (!i_s && b_s) driver <- c(driver, "Bacterial")
    if (!i_s && e_s) driver <- c(driver, "External")
    if (i_s) driver <- c(driver, "Self")
    cat(sprintf("  %-30s  T=%.4f  I=%.4f  B=%.4f  E=%.4f  [%s]\n",
      g,
      gd$JTK_P_raw[gd$Component == "Total"],
      gd$JTK_P_raw[gd$Component == "Independent"],
      gd$JTK_P_raw[gd$Component == "Bacterial"],
      gd$JTK_P_raw[gd$Component == "External"],
      paste(driver, collapse = "+")
    ))
  }
}

cat("\n===== Done =====\n")
