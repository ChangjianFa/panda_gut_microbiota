# Preserve tuner/libraries then clean environment
keep_patterns <- c("n_order", "spar", "output_dir", "^res$", "^tag$", "^run_dir$",
                   "^orders$", "^spars$", "^start_time", "^elapsed$", "^df$",
                   "^results$", "^rmse_overall$", "^r2_overall$", "^obs_y2$",
                   "^genus_metrics$", "^n_r2_pos$", "^model_residual$",
                   "^run_id$", "^total$", "^results_list$", "^result$",
                   "^base_results$", "^spar_tag$", "^save_vars$",
                   "^train_residual$")
to_keep <- ls(all.names = TRUE)[grepl(paste(keep_patterns, collapse="|"), ls(all.names = TRUE))]
tryCatch(rm(list = setdiff(ls(all.names = TRUE), to_keep)), error = function(e) NULL)
setwd("D:/data/panda")

if (!exists("n_order")) n_order <- 3
if (!exists("spar"))    spar    <- 0.5

library(deSolve)
library(orthopolynom)

# ============================================================
# 1. Load Year 1 model
# ============================================================
cat("===== Loading Year 1 model =====\n")
load("results/interaction_res.Rdata")
p_v <- res$p_v; p_e <- res$p_e; p <- res$p; n_order <- res$n_order
BETA_est <- res$BETA_est

all_month_y1 <- read.csv("data_processed/all_month.csv", row.names=1, check.names=FALSE)
y1_bac_names <- rownames(all_month_y1)[1:p_v]
colnames(BETA_est) <- y1_bac_names

cat("Model: p=", p, ", bacteria=", p_v, ", env=", p_e, "\n")

# ============================================================
# 2. Load Year 2 (genus set already aligned by intersection preprocessing)
# ============================================================
cat("\n===== Loading Year 2 data =====\n")
all_month_y2 <- read.csv("data_processed/all_month_year2.csv", row.names=1, check.names=FALSE)

# Extract bacteria + env data (must match model exactly)
y2_bac <- all_month_y2[1:p_v, , drop=FALSE]
y2_env <- all_month_y2[(p_v+1):(p_v+p_e), , drop=FALSE]
y2_data <- as.data.frame(t(rbind(y2_bac, y2_env)))

cat("Year 2:", ncol(y2_data), "variables,", nrow(y2_data), "months\n")
stopifnot(all(rownames(y2_bac) == colnames(BETA_est)))

times_y2 <- c(3, 8, 10, 12)
times_full <- 1:12

# ============================================================
# 3. Legendre basis projection
# ============================================================
cat("\n===== Projecting onto Legendre basis =====\n")

get_LOP_M <- function(x, n_order, times, name=NULL) {
  x <- as.numeric(x)
  if (max(abs(x)) < 1e-15) x <- x + 1e-10 * seq_along(x)
  LOP <- legendre.polynomials(n=n_order)
  leg <- polynomial.values(polynomials=LOP, x=scaleX(x, u=-1, v=1))
  fs <- lapply(leg, splinefun, x=times)
  mod <- function(Time, State, Pars, basis) {
    with(as.list(c(State, Pars)), { dy <- basis(Time); return(list(dy)) })
  }
  State <- rep(0, length(leg))
  leg <- sapply(1:length(fs), function(xi)
    ode(func=mod, y=State[[xi]], parms=NULL, basis=fs[[xi]], times=times)[,2])
  if (is.null(name)) {
    colnames(leg) <- paste0("leg", 0:(ncol(leg)-1))
  } else {
    colnames(leg) <- paste0(name, "__leg", 0:(ncol(leg)-1))
  }
  leg
}

x_smooth_y2 <- sapply(1:ncol(y2_data), function(xi) {
  fit <- smooth.spline(times_y2, as.numeric(y2_data[, xi]), spar=spar)
  predict(fit, x=times_full)$y
})
colnames(x_smooth_y2) <- colnames(y2_data)

x_basis_y2 <- Reduce(cbind, lapply(1:ncol(y2_data), function(xi) {
  get_LOP_M(x=x_smooth_y2[, xi], n_order, times_full)
}))

x_basis_est_y2 <- cbind(1, x_basis_y2[,1],
                         x_basis_y2[, !grepl("leg0", colnames(x_basis_y2))])

# ============================================================
# 4. Prediction
# ============================================================
pred_y2 <- x_basis_est_y2 %*% BETA_est
colnames(pred_y2) <- colnames(BETA_est)
rownames(pred_y2) <- times_full

pred_y2_at_obs <- pred_y2[as.character(times_y2), ]
obs_y2 <- as.matrix(y2_data[, 1:p_v])
colnames(obs_y2) <- colnames(BETA_est)

cat("Predicted:", dim(pred_y2_at_obs), "\n")
cat("Observed:", dim(obs_y2), "\n")

# ============================================================
# 5. Validation metrics
# ============================================================
cat("\n===== Validation metrics =====\n")
residuals <- obs_y2 - pred_y2_at_obs
rmse_overall <- sqrt(mean(residuals^2))
ss_res <- sum(residuals^2)
ss_tot <- sum((obs_y2 - mean(obs_y2))^2)
r2_overall <- 1 - ss_res / ss_tot

cat(sprintf("Overall RMSE: %.6f\n", rmse_overall))
cat(sprintf("Overall R²:    %.4f\n", r2_overall))

cat("\n--- Per-genus metrics ---\n")
cat(sprintf("%-35s %10s %10s %12s %12s\n", "Genus", "RMSE", "R²", "Obs_Mean", "Pred_Mean"))
cat(strrep("-", 80), "\n")

genus_metrics <- data.frame(stringsAsFactors=FALSE)
for (i in 1:p_v) {
  rmse_i <- sqrt(mean(residuals[, i]^2))
  ss_res_i <- sum(residuals[, i]^2)
  ss_tot_i <- sum((obs_y2[, i] - mean(obs_y2[, i]))^2)
  r2_i <- if (ss_tot_i > 0) 1 - ss_res_i / ss_tot_i else NA
  obs_m <- mean(obs_y2[, i]); pred_m <- mean(pred_y2_at_obs[, i])
  genus_metrics <- rbind(genus_metrics,
    data.frame(Genus=colnames(obs_y2)[i], RMSE=rmse_i, R2=r2_i,
               Obs_Mean=obs_m, Pred_Mean=pred_m, stringsAsFactors=FALSE))
  cat(sprintf("%-35s %10.6f %10.4f %12.6f %12.6f\n",
              colnames(obs_y2)[i], rmse_i, r2_i, obs_m, pred_m))
}

cat(sprintf("\nMean RMSE: %.6f\n", mean(genus_metrics$RMSE)))
cat(sprintf("Mean R²:   %.4f\n", mean(genus_metrics$R2, na.rm=TRUE)))
cat(sprintf("R² > 0:    %d / %d\n",
            sum(genus_metrics$R2 > 0, na.rm=TRUE), nrow(genus_metrics)))

# Baseline: Year 1 at matching months vs Year 2
cat("\n--- Baseline: Year 1 values at matching months vs Year 2 ---\n")
y1_at_matching <- all_month_y1[1:p_v, c("Month_3","Month_8","Month_10","Month_12")]
baseline_resid <- as.matrix(t(y2_bac)) - as.matrix(t(y1_at_matching))
base_rmse <- sqrt(mean(baseline_resid^2))
ss_base_res <- sum(baseline_resid^2)
base_r2 <- 1 - ss_base_res / ss_tot
cat(sprintf("Baseline RMSE (Y1 direct vs Y2): %.6f\n", base_rmse))
cat(sprintf("Baseline R²:                      %.4f\n", base_r2))

# Save
if (!exists("output_dir")) output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(pred_y2_at_obs, file=file.path(output_dir, "year2_predicted.csv"))
write.csv(obs_y2, file=file.path(output_dir, "year2_observed.csv"))
write.csv(genus_metrics, file=file.path(output_dir, "year2_genus_metrics.csv"), row.names=FALSE)

cat("\n===== Validation complete =====\n")
