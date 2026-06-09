# Preserve tuner/libraries then clean environment
keep_patterns <- c("n_order", "spar", "output_dir", "^res$", "^tag$", "^run_dir$",
                   "^orders$", "^spars$", "^start_time", "^elapsed$", "^df$",
                   "^run_id$", "^total$", "^results_list$", "^result$",
                   "^base_results$", "^spar_tag$", "^save_vars$")
to_keep <- ls(all.names = TRUE)[grepl(paste(keep_patterns, collapse="|"), ls(all.names = TRUE))]
tryCatch(rm(list = setdiff(ls(all.names = TRUE), to_keep)), error = function(e) NULL)
setwd("D:/data/panda")

if (!exists("n_order")) n_order <- 3
if (!exists("spar"))    spar    <- 0.5

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ADSIHT)
library(deSolve)
library(Matrix)
library(orthopolynom)
library(reshape2)
library(patchwork)
library(purrr)

# ============================================================
# Load joint preprocessed data
# ============================================================
cat("===== Loading joint preprocessed data =====\n")
all_month <- read.csv("data_processed/all_month.csv", row.names = 1, check.names = FALSE)
all_month <- as.data.frame(t(all_month))  # Transpose: rows=months, cols=variables

p <- ncol(all_month)
p_v <- p - 5
p_e <- 5
cat("Variables: total=", p, ", bacteria=", p_v, ", env=", p_e, "\n")
cat("Bacteria genera:", colnames(all_month)[1:p_v], "\n")

# ============================================================
# Part 4: ODEsolve - Legendre + ADSIHT interaction modeling
# ============================================================
ODEsolve <- function(data, times, alpha = 5e-3, kappa = 0.99, ic.scale = 0.6, ic.coef = 0.6) {
  orig_times <- times
  data_v <- data[, 1:p_v]
  data_e <- data[, (p_v + 1):p]

  fit <- lapply(1:ncol(data), function(xi)
    smooth.spline(times, as.numeric(data[, xi]), spar = spar))

  x_smooth <- sapply(1:ncol(data), function(xi) predict(fit[[xi]], x = times)$y)
  colnames(x_smooth) <- colnames(data)
  rownames(x_smooth) <- rownames(data)

  select <- c(1:p)
  plt_df1 <- melt(as.matrix(data[, select]))
  plt_df2 <- melt(x_smooth[, select])

  p_smooth <- ggplot() +
    geom_point(plt_df1, mapping = aes(x = Var1, y = value, group = Var2), size = 1) +
    geom_line(plt_df2, mapping = aes(x = Var1, y = value, group = Var2), linewidth = 1) +
    facet_wrap(~Var2, scales = "free_y") +
    labs(title = "Data (points) vs Smoothed (lines)")
  print(p_smooth)

  group <- rep(1:(p_v * p), each = n_order + 1)
  ind_Leg0 <- c(0:(p_v - 1)) * (p * (n_order + 1)) + 1 + rep(0:(p_v - 1)) * (n_order + 1)
  Leg0 <- seq(1, length(group), by = n_order + 1)
  dep_Leg0 <- setdiff(Leg0, ind_Leg0)
  group <- group[-dep_Leg0]

  get_LOP_M <- function(x, n_order, times, name = NULL) {
    x <- as.numeric(x)
    LOP <- legendre.polynomials(n = n_order)
    leg <- polynomial.values(polynomials = LOP, x = scaleX(x, u = -1, v = 1))
    fs <- lapply(leg, splinefun, x = times)
    mod <- function(Time, State, Pars, basis) {
      with(as.list(c(State, Pars)), {
        dy <- basis(Time)
        return(list(dy))
      })
    }
    Pars <- NULL
    State <- rep(0, length(leg))
    leg <- sapply(1:length(fs), function(xi)
      ode(func = mod, y = State[[xi]], parms = Pars,
          basis = fs[[xi]], times = times)[, 2])
    if (is.null(name)) {
      colnames(leg) <- paste0("leg", 0:(ncol(leg) - 1))
    } else {
      colnames(leg) <- paste0(name, "__leg", 0:(ncol(leg) - 1))
    }
    leg
  }

  x_basis <- Reduce(cbind, lapply(1:p, function(xi) {
    get_LOP_M(x = x_smooth[, xi], n_order, times)
  }))

  x_scaled <- scale(x_basis)
  X <- as.matrix(bdiag(rep(list(x_scaled), p_v)))
  X <- X[, -dep_Leg0]
  y_scaled <- scale(data_v)
  Y <- matrix(as.vector(y_scaled), ncol = 1)
  order0_positions <- sapply(as.numeric(names(table(group)[table(group) == max(table(group))])),
                             function(num) which(group == num)[1])

  res_ds <- ADSIHT(x = X, y = Y, group = group, L = n_order,
                   kappa = kappa, ic.scale = ic.scale, ic.coef = ic.coef)

  beta_ds <- res_ds$beta[, which.min(res_ds[["ic"]])]
  beta_ds <- matrix(beta_ds, nrow = p * n_order + 1, ncol = p_v)
  beta0_scaled_ds <- res_ds$intercept[which.min(res_ds[["ic"]])]
  lambda <- res_ds$lambda[which.min(res_ds[["ic"]])]
  cat("Selected lambda:", lambda, "\n")

  ind_0 <- which(sapply(1:p_v, function(i) {
    all(beta_ds[c(n_order * (i - 1) + 1, n_order * i + 1), i] == 0)
  }))
  cat("Equations without independent effect:", length(ind_0), "\n")

  for (equ in ind_0) {
    m <- nrow(x_scaled)
    q <- p * n_order + 1
    X_sub <- X[((equ - 1) * m + 1):(equ * m), ((equ - 1) * q + 1):(equ * q)]
    Y_sub <- Y[((equ - 1) * m + 1):(equ * m), ]
    X_sub <- X_sub[, c(((equ - 1) * n_order + 2):(equ * n_order + 1),
                       which(beta_ds[, equ] != 0))]
    res_sub <- lm(Y_sub ~ X_sub)
    ols_loss <- function(beta) {
      sum((Y_sub - (X_sub %*% beta + beta0_scaled_ds))^2 +
            alpha * (sum(beta^2)))
    }
    init_beta <- res_sub$coefficients[-1]
    init_beta[is.na(init_beta)] <- rep(0, length(init_beta[is.na(init_beta)]))
    ols_result <- optim(par = init_beta, fn = ols_loss, method = "BFGS")
    beta_ds[c(((equ - 1) * n_order + 2):(equ * n_order + 1),
              which(beta_ds[, equ] != 0)), equ] <- ols_result$par
  }

  ind_1 <- setdiff(1:p_v, ind_0)
  for (equ in ind_1) {
    m <- nrow(x_scaled)
    q <- p * n_order + 1
    X_sub <- X[((equ - 1) * m + 1):(equ * m), ((equ - 1) * q + 1):(equ * q)]
    Y_sub <- Y[((equ - 1) * m + 1):(equ * m), ]
    X_sub <- X_sub[, c(which(beta_ds[, equ] != 0))]
    res_sub <- lm(Y_sub ~ X_sub)
    res_sub$coefficients[is.na(res_sub$coefficients)] <-
      rep(0, length(res_sub$coefficients[is.na(res_sub$coefficients)]))
    ols_loss <- function(beta) {
      sum((Y_sub - (X_sub %*% beta + beta0_scaled_ds))^2 +
            alpha * (sum(beta^2)))
    }
    init_beta <- res_sub$coefficients[-1]
    ols_result <- optim(par = init_beta, fn = ols_loss, method = "BFGS")
    beta_ds[c(which(beta_ds[, equ] != 0)), equ] <- ols_result$par
  }
  beta_ds[abs(beta_ds) < lambda] <- 0

  get_B_est <- function(beta, beta0_scaled) {
    beta_est <- matrix(beta[-order0_positions], nrow = p * n_order, ncol = p_v)
    beta_est <- rbind(beta[order0_positions], beta_est)

    sY <- attr(y_scaled, "scaled:scale")
    mY <- attr(y_scaled, "scaled:center")
    sX <- attr(x_scaled, "scaled:scale")
    sX <- c(sX[1], sX[!grepl("leg0", names(sX))])
    mX <- attr(x_scaled, "scaled:center")
    mX <- c(mX[1], mX[!grepl("leg0", names(mX))])

    BETA <- lapply(1:p_v, function(xi) {
      beta_orig <- beta_est[, xi] * (sY[xi] / sX)
      return(beta_orig)
    })
    BETA0 <- lapply(1:p_v, function(xi) {
      orig_int <- beta0_scaled * sY[xi] + mY[xi] - sum(BETA[[xi]] * mX)
      return(orig_int)
    })
    B_est <- t(cbind(unlist(BETA0), t(Reduce(cbind, BETA))))
    B_est
  }
  B_est_ds <- get_B_est(beta = beta_ds, beta0_scaled = beta0_scaled_ds)
  x_basis_est <- cbind(1, x_basis[, 1], x_basis[, !grepl("leg0", colnames(x_basis))])

  plot_est <- function() {
    tmp <- x_basis_est %*% B_est_ds
    colnames(tmp) <- colnames(data)[1:p_v]
    rownames(tmp) <- orig_times
    tmp <- melt(tmp)

    rownames(data) <- orig_times
    tmp2 <- melt(as.matrix(data[, 1:p_v]))
    p <- ggplot() +
      geom_line(data = tmp, mapping = aes(x = Var1, y = value, color = Var2)) +
      geom_point(data = tmp2, mapping = aes(x = Var1, y = value, color = Var2),
                 alpha = 0.3) +
      facet_wrap(~Var2, scales = "free_y") +
      guides(color = "none") +
      labs(title = "Fitted (line) vs Observed (point)")
    p
  }
  print(plot_est())

  residual <- data[, 1:p_v] - x_basis_est %*% B_est_ds
  residual <- sum(residual^2)

  return(list(p = p, p_v = p_v, p_e = p_e, n_order = n_order,
              x_basis = x_basis_est, BETA_est = B_est_ds,
              x_smooth = x_smooth, obs = data, residual = residual,
              x_fit = fit))
}

if (!exists("n_order")) n_order <- 3
times <- 1:12
cat("\n===== Running ODEsolve on Year 1 =====\n")
res <- ODEsolve(data = all_month, times = times)
if (!exists("output_dir")) output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
save(res, file = "results/interaction_res.Rdata")
cat("Model residual:", res$residual, "\n")

# ============================================================
# Part 5: Effect decomposition and visualization
# ============================================================
B_est_ds <- res$BETA_est
colnames(B_est_ds) <- colnames(all_month)[1:p_v]
x_basis_est <- res$x_basis
obs <- res$obs[, 1:p_v]
x_smooth <- res$x_smooth

idx <- seq(3, (p * n_order + 2))
idx_list <- rep(list(split(idx, cut(seq_along(idx), breaks = p, labels = FALSE))),
                times = p_v)
idx_list <- lapply(seq_along(idx_list), function(i) {
  idx_list[[i]][[i]] <- c(1, 2, idx_list[[i]][[i]])
  return(idx_list[[i]])
})

trans_out <- function(out, j) {
  id <- paste0("M", j)
  out <- as.data.frame(out)
  colnames(out) <- c('time', "est", paste0("M", 1:p))
  out <- out[, colSums(out) != 0, drop = FALSE]
  out <- melt(out, id.vars = 'time')
  out$group <- 'Dep'
  out$group[out$variable == id] <- 'Ind'
  out$group[out$variable == 'est'] <- 'est'
  out
}

trans_out1 <- function(out, j) {
  out <- as.data.frame(out)
  out <- melt(out, id.vars = 'time')
  out$group <- 'Dep'
  out$group[out$variable == 'Ind'] <- 'Ind'
  out$group[out$variable == 'est'] <- 'est'
  out
}

get_effect_plot <- function(B_est, x_basis_est, times) {
  trans_effect_est_list <- lapply(1:p_v, function(j) {
    effect_est <- cbind(times,
                        rowSums(sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est[idx_list[[j]][[xi]], j])),
                        sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est[idx_list[[j]][[xi]], j]))
    trans_out(effect_est, j)
  })

  effect_est_all <- do.call(rbind, lapply(seq_along(trans_effect_est_list),
    function(i) {
      df <- trans_effect_est_list[[i]]
      df$equations <- colnames(x_smooth)[i]
      return(df)
    }))

  ori <- as.data.frame(cbind(times, obs))
  ori <- melt(ori, id.vars = 'times')
  ori$equations <- rep(colnames(x_smooth)[1:p_v], each = length(times))

  # Summary effect (Ind, Bacterial, External)
  trans_effect_est_list <- lapply(1:p_v, function(j) {
    effect_est <- cbind(times,
                        rowSums(sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est[idx_list[[j]][[xi]], j])),
                        sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est[idx_list[[j]][[xi]], j]))
    effect_est <- cbind(effect_est[, 1:2],
                        effect_est[, j + 2],
                        rowSums(effect_est[, c(1:p_v + 2)[-j]]),
                        rowSums(effect_est[, (p_v + 3):ncol(effect_est)]))
    colnames(effect_est) <- c('time', 'est', 'Ind', "Bacterial", 'External')
    trans_out1(effect_est, j)
  })

  effect_est_all <- do.call(rbind, lapply(seq_along(trans_effect_est_list),
    function(i) {
      df <- trans_effect_est_list[[i]]
      df$equations <- colnames(x_smooth)[i]
      return(df)
    }))

  dep_data <- effect_est_all[effect_est_all$group == "Dep", ]
  dep_labels <- dep_data %>%
    group_by(variable) %>%
    filter(time == max(time))

  p3 <- ggplot() +
    geom_line(effect_est_all,
              mapping = aes(x = time, y = value, group = variable, color = group),
              linewidth = 0.3, linetype = 1) +
    geom_point(ori, mapping = aes(x = times, y = value),
               color = 'blue', size = 0.5, alpha = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black",
               linewidth = 0.15) +
    geom_text(dep_labels,
              mapping = aes(x = time, y = value, label = variable),
              nudge_x = -0.3, size = 2, show.legend = FALSE, color = "black") +
    scale_color_manual(values = c(est = 'blue', Ind = 'red', Dep = 'green')) +
    theme_minimal() +
    facet_wrap(~equations, scales = "free_y",
               nrow = ceiling(sqrt(p_v)),
               ncol = ceiling(p_v / ceiling(sqrt(p_v)))) +
    labs(x = "Month", y = "Effect Value", color = "") +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size = 10,
                                margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.1, "cm")
    ) +
    scale_x_continuous(breaks = unique(effect_est_all$time))

  return(list(p3 = p3))
}

cat("\n===== Effect decomposition =====\n")
effect_plot <- get_effect_plot(B_est = B_est_ds, x_basis_est = x_basis_est, times = times)
effect_plot$p3
ggsave(file.path(output_dir, "effect_plot.pdf"), effect_plot$p3, width = 14, height = 11)

# ---- Export four effect curves (Total, Ind, Bacterial, External) ----
effect_curves <- data.frame(Month = integer(), Genus = character(),
                            Total = numeric(), Independent = numeric(),
                            Bacterial = numeric(), External = numeric(),
                            stringsAsFactors = FALSE)
for (j in 1:p_v) {
  est_mat <- cbind(times,
                   rowSums(sapply(1:p, function(xi)
                     x_basis_est[, idx_list[[j]][[xi]]] %*%
                       B_est_ds[idx_list[[j]][[xi]], j])),
                   sapply(1:p, function(xi)
                     x_basis_est[, idx_list[[j]][[xi]]] %*%
                       B_est_ds[idx_list[[j]][[xi]], j]))
  est_total <- est_mat[, 2]
  est_ind   <- est_mat[, j + 2]
  est_bac   <- rowSums(est_mat[, c(1:p_v + 2)[-j], drop = FALSE])
  est_env   <- rowSums(est_mat[, (p_v + 3):ncol(est_mat), drop = FALSE])
  for (t in 1:length(times)) {
    effect_curves <- rbind(effect_curves, data.frame(
      Month = times[t],
      Genus = colnames(x_smooth)[j],
      Total = est_total[t],
      Independent = est_ind[t],
      Bacterial = est_bac[t],
      External = est_env[t],
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(effect_curves, file.path(output_dir, "effect_curves.csv"), row.names = FALSE)
cat("Effect curves saved to", file.path(output_dir, "effect_curves.csv"), "\n")

# ============================================================
# Part 6: Network construction
# ============================================================
cat("\n===== Network construction =====\n")

trans_d_effect_est_list <- lapply(1:p_v, function(j) {
  effect_est <- cbind(times,
                      rowSums(sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j])),
                      sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j]))
  fit <- lapply(1:(p + 1), function(xi)
    smooth.spline(effect_est[, 1], as.numeric(effect_est[, -1][, xi])))
  x_deriv <- sapply(1:(p + 1), function(xi)
    predict(fit[[xi]], x = effect_est[, 1], deriv = 1)$y)
  x_deriv <- cbind(effect_est[, 1], x_deriv)
  dep <- which(colSums(x_deriv[, c(-1, -2)]) != 0)
  dep <- dep[which(dep != j)]
  for (di in dep) {
    x <- x_smooth[, di]
    f <- x_deriv[, c(-1, -2)][, di]
    fit_spline <- smooth.spline(x, f, spar = spar)
    x_deriv[, c(-1, -2)][, di] <-
      x_deriv[, c(-1, -2)][, di] - predict(fit_spline, min(x))$y
  }
  if (length(dep) == 1) {
    x_deriv[, -c(1, 2)][, j] <- x_deriv[, 2] - x_deriv[, -c(1, 2)][, dep]
  } else {
    x_deriv[, -c(1, 2)][, j] <- x_deriv[, 2] - rowSums(x_deriv[, -c(1, 2)][, dep])
  }
  x_deriv
})

edge_matrix_to_list <- function(edge_matrix) {
  edge_list <- data.frame()
  for (i in 1:nrow(edge_matrix)) {
    for (j in 1:ncol(edge_matrix)) {
      if (i != j && edge_matrix[i, j] != 0) {
        edge_list <- rbind(edge_list, data.frame(
          source = colnames(edge_matrix)[j],
          target = rownames(edge_matrix)[i],
          cor = edge_matrix[i, j]
        ))
      }
    }
  }
  return(edge_list)
}

edge_matrix <- lapply(1:p_v, function(j) {
  effect_est <- trans_d_effect_est_list[[j]][, -c(1, 2)]
  result <- colMeans(effect_est)
})
edge_matrix <- do.call(rbind, edge_matrix)
colnames(edge_matrix) <- colnames(all_month)
rownames(edge_matrix) <- colnames(all_month)[1:p_v]

edge_matrix <- edge_matrix_to_list(edge_matrix)
edge_matrix <- edge_matrix %>%
  mutate(
    sign = ifelse(cor >= 0, "1", "-1"),
    weight = abs(cor)
  ) %>%
  select(source, target, sign, weight)

cat("Network edges: ", nrow(edge_matrix), "\n")
write.csv(edge_matrix, file = file.path(output_dir, "edge_matrix.csv"), row.names = FALSE, quote = FALSE)

# ============================================================
# Part 7-8: PC1 analysis and contribution quantification
# ============================================================
cat("\n===== PC1 Beta-diversity analysis =====\n")

abundance_mat <- lapply(1:p_v, function(j) {
  effect_est <- rowSums(sapply(1:p, function(xi)
    x_basis_est[, idx_list[[j]][[xi]]] %*%
      B_est_ds[idx_list[[j]][[xi]], j]))
  effect_est
})
abundance_mat <- do.call(cbind, abundance_mat)
colnames(abundance_mat) <- colnames(B_est_ds)

min_val <- min(abundance_mat)
abundance_mat_pos <- if (min_val < 0) abundance_mat + abs(min_val) else abundance_mat
bc_dist <- vegdist(abundance_mat_pos, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
PC1_total <- pcoa_res$points[, 1]

abundance_env <- lapply(1:p_v, function(j) {
  effect_est <- cbind(times,
                      rowSums(sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j])),
                      sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j]))
  effect_est <- effect_est[, (p_v + 3):ncol(effect_est)]
  rowSums(effect_est)
})
abundance_env <- do.call(cbind, abundance_env)
min_val <- min(abundance_env)
abundance_env_pos <- if (min_val < 0) abundance_env + abs(min_val) else abundance_env
bc_dist <- vegdist(abundance_env_pos, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
PC1_env <- pcoa_res$points[, 1]

abundance_bac <- lapply(1:p_v, function(j) {
  effect_est <- cbind(times,
                      rowSums(sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j])),
                      sapply(1:p, function(xi)
                        x_basis_est[, idx_list[[j]][[xi]]] %*%
                          B_est_ds[idx_list[[j]][[xi]], j]))
  effect_est <- effect_est[, 3:(p_v + 3)]
  rowSums(effect_est)
})
abundance_bac <- do.call(cbind, abundance_bac)
min_val <- min(abundance_bac)
abundance_bac_pos <- if (min_val < 0) abundance_bac + abs(min_val) else abundance_bac
bc_dist <- vegdist(abundance_bac_pos, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
PC1_interaction <- pcoa_res$points[, 1]

# Individual environment PC1
PC1_envs <- c()
for (z in 1:p_e) {
  abundance_env_single <- lapply(1:p_v, function(j) {
    effect_est <- cbind(times,
                        rowSums(sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est_ds[idx_list[[j]][[xi]], j])),
                        sapply(1:p, function(xi)
                          x_basis_est[, idx_list[[j]][[xi]]] %*%
                            B_est_ds[idx_list[[j]][[xi]], j]))
    effect_est <- effect_est[, (p_v + 3):ncol(effect_est)]
    effect_est <- effect_est[, z]
    effect_est
  })
  abundance_env_single <- do.call(cbind, abundance_env_single)
  min_val <- min(abundance_env_single)
  abundance_env_single_pos <- if (min_val < 0) abundance_env_single + abs(min_val)
                               else abundance_env_single
  bc_dist <- vegdist(abundance_env_single_pos, method = "bray")
  pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
  pc1 <- pcoa_res$points[, 1]
  PC1_envs <- cbind(PC1_envs, pc1)
}
colnames(PC1_envs) <- colnames(all_month)[(p_v + 1):p]

# ============================================================
# Part 8: Contribution Quantification
# ============================================================
cat("\n===== Contribution quantification =====\n")

var_total <- var(PC1_total)
var_inter <- var(PC1_interaction)
var_env <- var(PC1_env)

contrib_var_inter <- var_inter / var_total
contrib_var_env <- var_env / var_total

cat("Variance method:\n")
cat("  Interaction contribution:", round(contrib_var_inter, 3), "\n")
cat("  Environment contribution:", round(contrib_var_env, 3), "\n\n")

contrib_inter <- cor(PC1_total, PC1_interaction)^2
contrib_env <- cor(PC1_total, PC1_env)^2
contrib_envs <- apply(PC1_envs, 2, function(x) cor(PC1_total, x)^2)

cat("R^2 method:\n")
cat("  Interaction contribution:", round(contrib_inter, 3), "\n")
cat("  Environment contribution:", round(contrib_env, 3), "\n")

all_contrib <- c(Interaction = contrib_inter, contrib_envs)
pct_contrib <- all_contrib / sum(all_contrib) * 100
cat("  Contribution percentages (%):\n")
print(round(pct_contrib, 1))

cat("\nRegression method:\n")
PC1_total_s <- scale(PC1_total)
PC1_inter_s <- scale(PC1_interaction)
PC1_env_s <- scale(PC1_env)

lm_res <- lm(PC1_total_s ~ PC1_inter_s + PC1_env_s)
summary_lm <- summary(lm_res)
coef_table <- summary_lm$coefficients
cat("  Interaction coefficient:", round(coef_table["PC1_inter_s", 1], 3), "\n")
cat("  Environment coefficient:", round(coef_table["PC1_env_s", 1], 3), "\n")
cat("  R^2:", round(summary_lm$r.squared, 3), "\n\n")

# ============================================================
# Part 9: Plot PC1 dynamics
# ============================================================
cat("===== Plotting PC1 dynamics =====\n")

time_vec <- 1:12
plot_data <- data.frame(
  Time = rep(time_vec, 3),
  PC1 = c(PC1_total, PC1_interaction, PC1_env),
  Group = rep(c("Total", "Interaction", "Environment"), each = length(time_vec))
)

p_pc1 <- ggplot(plot_data, aes(x = Time, y = PC1, color = Group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = c("Total" = "black",
                                "Interaction" = "red",
                                "Environment" = "blue")) +
  labs(x = "Month", y = "PC1") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_x_continuous(breaks = unique(plot_data$Time))

print(p_pc1)
ggsave(file.path(output_dir, "PC1_dynamics.pdf"), p_pc1, width = 10, height = 6)

# Save PC1 contribution results
pc1_results <- rbind(
  data.frame(Method = "Variance", Component = "Interaction", Value = contrib_var_inter),
  data.frame(Method = "Variance", Component = "Environment", Value = contrib_var_env),
  data.frame(Method = "R2", Component = "Interaction", Value = contrib_inter),
  data.frame(Method = "R2", Component = "Environment", Value = contrib_env),
  data.frame(Method = "R2", Component = names(pct_contrib), Value = as.numeric(pct_contrib)),
  data.frame(Method = "Regression", Component = "Interaction_coef", Value = coef_table["PC1_inter_s", 1]),
  data.frame(Method = "Regression", Component = "Environment_coef", Value = coef_table["PC1_env_s", 1]),
  data.frame(Method = "Regression", Component = "R2", Value = summary_lm$r.squared)
)
write.csv(pc1_results, file.path(output_dir, "pc1_contributions.csv"), row.names = FALSE)
cat("PC1 contributions saved to", file.path(output_dir, "pc1_contributions.csv"), "\n")

cat("\n===== Year 1 model complete =====\n")
