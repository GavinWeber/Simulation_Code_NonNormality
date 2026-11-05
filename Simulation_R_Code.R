# ==============================================================================
# WHEN IT'S NOT NORMAL — SIMULATION / POWER ANALYSIS
# ==============================================================================

# Sets limit to 4000 MiB
options(future.globals.maxSize = 4000 * 1024^2)

# Force the Linux server to use UTF-8 encoding
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# OUTPUT DIRECTORY

# All generated files will be saved inside:
#
#   output/
#       YYYY-MM-DD/
#
# The date is generated automatically each time the script runs.

output_root <- file.path(getwd(), "output")

output_dir <- file.path(
  output_root,
  format(Sys.Date(), "%Y-%m-%d")
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n============================================================\n")
cat("Output directory:\n")
cat(output_dir, "\n")
cat("============================================================\n\n")

# =================Libraries===========================

library(tidyverse)
library(purrr)
library(writexl)
library(matrixTests)
library(data.table)

library(gt)
library(viridis)

library(future)
library(furrr)
library(gridExtra)
library(cowplot)

library(robustbase)
library(patchwork)
library(webshot2)
library(ragg)
library(listenv)

# Use Slurm's allocation if available; otherwise use four workers
slurm_workers <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
manual_workers <- Sys.getenv("SIM_WORKERS", unset = "4")

requested_workers <-
  if (nzchar(slurm_workers)) slurm_workers else manual_workers

n_workers <- suppressWarnings(as.integer(requested_workers))

if (is.na(n_workers) || n_workers < 1L) {
  n_workers <- 4L
}

# Conservative cap for a shared server
n_workers <- min(n_workers, 4L)

message("Using ", n_workers, " parallel workers")

plan(multisession, workers = n_workers)

# REPRODUCIBLE PARALLEL RNG

RNGkind("L'Ecuyer-CMRG")
set.seed(13)

# MOORS' KURTOSIS

moors_kurtosis <- function(x) {
  
  e <- quantile(
    x,
    probs = c(1/8, 2/8, 3/8, 5/8, 6/8, 7/8),
    names = FALSE
  )
  
  moors_val <- (
    (e[6] - e[4]) +
      (e[3] - e[1])
  ) / (e[5] - e[2])
  
  return(moors_val)
}

# ========================C2 & BRUNNER-MUNZEL TEST FUNCTION==================================

getTestResults <- function(
    X,
    Y,
    alpha = 0.05,
    permutation = TRUE,
    nperm = 10000,
    seed = sample(1:100000, 1)
) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (ncol(X) != ncol(Y)) {
    stop("number of samples should be equal for both groups")
  }
  
  d <- ncol(X)
  
  results <- data.table(
    test = character(),
    lower = numeric(),
    theta = numeric(),
    upper = numeric(),
    variance = numeric(),
    ties = numeric(),
    statistic = numeric(),
    pval = numeric(),
    rejH0 = logical(),
    separation = logical()
  )
  
  for (i in 1:d) {
    
    x <- X[, i]
    y <- Y[, i]
    
    nx <- length(x)
    ny <- length(y)
    Nxy <- nx + ny
    m <- min(nx, ny)
    xy <- c(x, y)
    rxy <- rank(xy)
    rx <- rank(x)
    ry <- rank(y)
    
    rxry <- c(rx, ry)
    rxmax <- rank(x, ties.method = "max")
    rxmin <- rank(x, ties.method = "min")
    
    rxymax <- rank(xy, ties.method = "max")[1:nx]
    rxymin <- rank(xy, ties.method = "min")[1:nx]
    
    pl <- rxy - rxry
    
    mplx <- mean(pl[1:nx])
    mply <- mean(pl[(nx + 1):Nxy])
    
    plx <- (1 / ny) * (rxy[1:nx] - rx)
    ply <- (1 / nx) * (rxy[(nx + 1):Nxy] - ry)
    
    vx <- var(plx)
    vy <- var(ply)
    
    theta <- mean(ply)
    
    # Brunner-Munzel Test

    sigmaBF <-  Nxy * (vx / nx + vy / ny)
    
    df <- (
      vx / nx + vy / ny
      )^2 / (
      vx^2 / (nx^2 * (nx - 1)) +
        vy^2 / (ny^2 * (ny - 1))
    )
    
    df[is.na(df)] <- 1000
    
    thetaBM <- theta
    theta0 <- theta == 0
    theta1 <- theta == 1
    
    thetaBM[theta0] <-
      if (any(x %in% y))
        1 / (2 * nx * ny)
    else
      1 / (nx * ny)
    
    thetaBM[theta1] <-
      if (any(x %in% y))
        (2 * nx * ny - 1) / (2 * nx * ny)
    else
      (nx * ny - 1) / (nx * ny)
    
    sigmaBF[sigmaBF == 0] <-
      if (any(x %in% y))
        Nxy / (2 * nx^2 * ny^2)
    else
      2 * Nxy / (nx^2 * ny^2)
    
    se <- sqrt(sigmaBF / Nxy)
    
    Tn_BM <-
      sqrt(Nxy) *
      (thetaBM - 1 / 2) /
      sqrt(sigmaBF)
    
    pval_BM <-
      2 *
      min(
        pt(Tn_BM, df = df),
        1 - pt(Tn_BM, df = df)
      )
    
    lower_BM <-
      theta -
      se *
      qt(1 - alpha / 2, df = df)
    
    upper_BM <-
      theta +
      se *
      qt(1 - alpha / 2, df = df)
    
    
    # C2 Test
    
    tau <-
      1 / ny *
      (mean(rxymax) - mean(rxymin) - (mean(rxmax) - mean(rxmin))
      )
    
    sigmaN <-
      1 /
      (nx * (nx - 1) * ny * (ny - 1)) *
      (sum((pl -c(rep(mplx, nx),rep(mply, ny)))^2) -
          nx * ny *
          (theta * (1 - theta) - 0.25 * tau)
      )
    
    q <- sigmaN /(theta *(1 - theta))
    
    C2 <- (4 / q) * (theta - 1 / 2)^2
    
    C2[sigmaN == 0] <-
      4 * m * (theta - 1 / 2)^2
    
    pval_C2 <- 1 - pchisq(C2, df = 1)
    
    cAlpha <- qchisq(1 - alpha, df = 1)
    
    ctr <-
      (2 * theta + q * cAlpha) /
      (2 * (1 + q * cAlpha))
    
    margin <-
      (sqrt(cAlpha) *
          sqrt(4 * q * theta * (1 - theta) + q^2 * cAlpha)
      ) / (2 * (1 + q * cAlpha))
    
    lower_r <- ctr - margin
    
    upper_r <- ctr + margin
    
    lower_BK <-
      (2 * m * theta + cAlpha
      ) /
      (2 * (m + cAlpha)) -
      (sqrt(cAlpha) * sqrt(4 * m * theta * (1 - theta) + cAlpha)
      ) /
      (2 *(m + cAlpha))
    
    upper_BK <-
      (2 * m * theta + cAlpha) /
      (2 * (m + cAlpha)) +
      (sqrt(cAlpha) * sqrt(4 *m *theta *(1 - theta) + cAlpha)
      ) / (2 * (m + cAlpha))
    
    
    # Save test results
    
    results_BM <- data.table(
      test = "BM-Test",
      lower = lower_BM,
      theta = thetaBM,
      upper = upper_BM,
      variance = sigmaBF,
      ties = tau,
      statistic = Tn_BM,
      pval = pval_BM,
      rejH0 = pval_BM < alpha,
      separation = theta0 | theta1
    )
    
    results_C2 <- data.table(
      test = "C2-Test",
      lower =
        ifelse(
          theta0 | theta1,
          lower_BK,
          lower_r
        ),
      theta = theta,
      upper =
        ifelse(
          theta0 | theta1,
          upper_BK,
          upper_r
        ),
      variance = sigmaN * Nxy,
      ties = tau,
      statistic = C2,
      pval = pval_C2,
      rejH0 = pval_C2 < alpha,
      separation = theta0 | theta1
    )
    
    results <-
      rbind(
        results,
        results_BM,
        results_C2
      )
    
    if (permutation) {
      
      results_perm <- getPermutationResults(
        x = x,
        y = y,
        Tn_BM = Tn_BM,
        theta = thetaBM,
        se = se,
        nperm = nperm,
        alpha = alpha,
        seed = seed,
        separation = theta0 | theta1
      )
      
      results <-
        rbind(
          results,
          results_perm
        )
    }
  }
  
  results <-
    results[
      order(results$test),
    ]
  
  results$sample <-
    if (permutation) {
      
      rep(
        paste("Sample", 1:d),
        3
      )
      
    } else {
      
      rep(
        paste("Sample", 1:d),
        2
      )
    }
  
  return(results)
}



# ============================SIMULATORS==============================

simulate_floor <- function(
    n,
    Nsim,
    p,
    mu,
    sigma
) {
  
  total_obs <- n * Nsim
  
  # TRUE indicates a non-floor observation
  nonfloor <- rbinom(total_obs, 1, 1 - p) == 1
  
  y <- numeric(total_obs)
  
  if (any(nonfloor)) {
    # Probability that an underlying N(mu, sigma^2) value is <= 0
    lower_prob <- pnorm(0, mean = mu, sd = sigma)
    
    # Sample from N(mu, sigma^2) conditional on X > 0
    u <- lower_prob +
      runif(sum(nonfloor)) * (1 - lower_prob)
    
    u <- pmin(u, 1 - .Machine$double.eps)
    
    y[nonfloor] <- qnorm(
      u,
      mean = mu,
      sd = sigma
    )
  }
  
  matrix(y, nrow = n, ncol = Nsim)
}


simulate_degen <- function(
    n,
    Nsim,
    mean,
    sigma
) {
  
  total_obs <- n * Nsim
  
  matrix(
    rep(mean, total_obs) +
      rnorm(total_obs,0,sigma),
    nrow = n,
    ncol = Nsim
  )
}


simulate_t_heavy <- function(
    n,
    Nsim,
    df,
    location = 0,
    scale = 1
) {
  
  total_obs <- n * Nsim
  
  z <-
    if (is.infinite(df))
      rnorm(total_obs,mean = 0,sd = 1)
  else
    rt(total_obs,df = df)
  
  matrix(
    location +scale * z,
    nrow = n,
    ncol = Nsim
  )
}


simulate_pareto <- function(
    n,
    Nsim,
    xm = 1,
    alpha = 2,
    eps = 1e-12
) {
  
  total_obs <- n * Nsim
  
  u <-pmin(pmax(runif(total_obs),eps),1 - eps)
  
  matrix(
    xm *(u^(-1 / alpha)),
    nrow = n, ncol = Nsim
  )
}

# ================SIMULATOR REGISTRY=======================================

sim_registry <- tribble(
  ~distribution, ~sim_fun,           ~defaults,
  
  "Floor",
  simulate_floor,
  list(
    p = 0.2,
    mu = 1,
    sigma = 0.5
  ),
  
  "Degenerate",
  simulate_degen,
  list(
    mean = 1,
    sigma = 0.2
  ),
  
  "t-heavy",
  simulate_t_heavy,
  list(
    df = 3,
    location = 0,
    scale = 1
  ),
  
  "Pareto",
  simulate_pareto,
  list(
    xm = 1,
    alpha = 2,
    eps = 1e-12
  )
)



# PARAMETER GRID BUILDERS


make_param_grid_floor <- function(
    sim_row,
    p_levels = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8),
    shift_mu_lvls = c(0,0.125,0.45,1),
    sigma = 0.5,
    sample_sizes = c(5,15,30,50,60,80,100)
    ){
  
  sfun <- sim_row$sim_fun[[1]]
  
  expand_grid(
    distribution = "Floor",
    p = p_levels,
    shift_mu = shift_mu_lvls,
    sample_size = sample_sizes
  ) |>
    mutate(
      prop_floor = p,
      mu1_in = 1,
      mu2_in = 1 + shift_mu,
      mix_mean1 = (1 - p) * mu1_in,
      mix_mean2 = (1 - p) * mu2_in,
      shift_mix = mix_mean2 - mix_mean1,
      
      args1 =
        pmap(list(p),
             \(pf)
          list(p = pf, mu = 1,sigma = sigma)
        ),
      
      args2 =
        pmap(list(p, shift_mu),
          \(pf, sh)
          list(p = pf,mu = 1 + sh, sigma = sigma)
        ),
      
      H0_true = shift_mu == 0,
      
      sim_fun = list(sfun)
    )
}


make_param_grid_degen <- function(
    sim_row,
    sigma_levels = c(0.0001,0.0005,0.0008,0.001,0.005,0.01,0.05),
    d_targets = c(0,0.2,0.5,0.8),
    sample_sizes = c(5,15,30,50,60,80,100)
) {
  
  sfun <- sim_row$sim_fun[[1]]
  
  expand_grid(
    distribution = "Degenerate",
    sigma = sigma_levels,
    d_target = d_targets,
    sample_size = sample_sizes
  ) |>
    mutate(
      shift = d_target * sigma,
      mean1 = 1,
      mean2 = 1 + shift,
      args1 =
        pmap(
          list(sigma),
          \(sg)
          list(mean = 1,sigma = sg)
        ),
      args2 =
        pmap(
          list(sigma, shift),
          \(sg, sh)
          list(mean = 1 + sh,sigma = sg)
        ),
      H0_true = d_target == 0,
      sim_fun = list(sfun)
    )
}


make_param_grid_t_heavy <- function(
    sim_row,
    df_levels = c(0.5,1,1.5,2,3,5,Inf),
    d_targets = c(0,0.2,0.5,0.8),
    sample_sizes = c(5,15,30,50,60,80,100),
    location = 0,
    scale = 1
) {
  
  t_sigma_r <- function(
    df,scale = 1
  ) {
    
    q75 <-
      if (is.infinite(df))
        qnorm(0.75)
    else
      qt(0.75, df)
    
    q25 <-
      if (is.infinite(df))
        qnorm(0.25)
    else
      qt(0.25, df)
    
    (scale *(q75 - q25)) / 1.349
  }
  
  
  sfun <- sim_row$sim_fun[[1]]
  
  expand_grid(
    distribution = "t-heavy",
    df = df_levels,
    d_target = d_targets,
    sample_size = sample_sizes
  ) |>
    mutate(
      sigma_r =
        vapply(df,\(v)
          t_sigma_r(v,scale),
          numeric(1)
          ),
      
      shift =
        d_target * sigma_r,
      loc1 = location,
      loc2 = location + shift,
      scale1 = scale,
      scale2 = scale,
      
      args1 =
        pmap(
          list(df),
          \(v)
          list(df = v, location = location, scale = scale)
        ),
      
      args2 =
        pmap(
          list(df, shift),
          \(v, sh)
          list(df = v, location = location + sh, scale = scale)
        ),
      
      H0_true = d_target == 0,
      sim_fun = list(sfun)
    )
}


make_param_grid_pareto <- function(
    sim_row,
    alpha_levels =
      c(0.6,0.75,1,1.10,1.40,1.70,2),
    r_levels =
      c(1,1.5,2,2.5),
    sample_sizes =
      c(5,15,30,50,60,80,100),
    xm_base = 1
) {
  
  sfun <- sim_row$sim_fun[[1]]
  
  expand_grid(
    distribution = "Pareto",
    alpha = alpha_levels,
    r = r_levels,
    sample_size = sample_sizes
  ) |>
    mutate(
      
      xm2 = xm_base,
      xm1 = r * xm2,
      
      H0_true = r == 1,
      
      eps_level =
        dplyr::if_else(
          alpha < 1.2,
          1e-8,
          1e-12
        ),
      
      mean1_real =
        dplyr::if_else(
          alpha > 1,
          alpha * xm1 /
            (alpha - 1),
          NA_real_
        ),
      
      mean2_real =
        dplyr::if_else(
          alpha > 1,
          alpha * xm2 /
            (alpha - 1),
          NA_real_
        ),
      
      mean_ratio = mean1_real / mean2_real,
      
      args1 = pmap(
          list(xm1, alpha, eps_level),
          \(x, a, e)
          list(xm = x, alpha = a, eps = e)
        ),
      
      args2 = pmap(
          list(xm2,alpha,eps_level),
          \(x, a, e)
          list(xm = x,alpha = a,eps = e)
        ),
      
      sim_fun = list(sfun)
    )
}


# ====================CENTRAL DESIGN LEVELS===================================

sample_sizes_ext <-
  c(5,15,30,50,100,500,1500,2500)

prop_levels_floor <-
  c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)

loc_shifts_floor <-
  c(0,0.125,0.45,1)

sigma_levels <-
  c(0.0001,0.0005,0.0008,0.001,0.005,0.01,0.05)

d_targets <-
  c(0,0.2,0.5,0.8)

df_levels <-
  c(0.5,1,1.5,2,3,5,Inf)

d_targets_t <-
  c(0,0.2,0.5,0.8)

alpha_levels <-
  c(0.6,0.75,1,1.10,1.40,1.70,2)

r_levels <-
  c(1,1.5,2,2.5)

xm_base <- 1


# ===========Build Grids===================================


grid_floor <- make_param_grid_floor(
  filter(
    sim_registry,
    distribution == "Floor"
  ),
  p_levels = prop_levels_floor,
  shift_mu_lvls = loc_shifts_floor,
  sample_sizes = sample_sizes_ext
)

grid_degen <- make_param_grid_degen(
  filter(
    sim_registry,
    distribution == "Degenerate"
  ),
  sigma_levels = sigma_levels,
  d_targets = d_targets,
  sample_sizes = sample_sizes_ext
)

grid_t <- make_param_grid_t_heavy(
  filter(
    sim_registry,
    distribution == "t-heavy"
  ),
  df_levels = df_levels,
  d_targets = d_targets_t,
  sample_sizes = sample_sizes_ext
)

grid_pareto <- make_param_grid_pareto(
  filter(
    sim_registry,
    distribution == "Pareto"
  ),
  alpha_levels = alpha_levels,
  r_levels = r_levels,
  sample_sizes = sample_sizes_ext,
  xm_base = xm_base
)


# ===================MONTE CARLO CONTROLS======================


nsim_by_n <- function(n) {
  
  case_when(
    n < 150 ~ 25000L,
    n < 300 ~ 25000L,
    n >= 300 ~ 25000L,
    TRUE ~ 25000L
  )
}



# FAST TEST KERNELS


testing_t_fast <- function(
    X,
    Y,
    alpha = 0.05
) {
  results <-
    matrixTests::col_t_welch(X,Y)
  
  tibble(
    power = mean(results$pvalue < alpha, na.rm = TRUE),
    cohen_d = NA_real_
  )
}


testing_new_tests_fast <- function(
    X,
    Y,
    alpha = 0.05
) {
  
  raw_res <- getTestResults(X, Y, alpha = alpha,permutation = FALSE)
  
  summ_res <-
    raw_res[, .(rate = mean(rejH0, na.rm = TRUE)), by = test]
  
  tibble(
    bm_rate =summ_res$rate[summ_res$test == "BM-Test"],
    
    c2_rate =
      summ_res$rate[summ_res$test == "C2-Test"]
  )
}

# SIMULATION ORCHESTRATOR 


run_simulation <- function(
    sim_fun,
    args1,
    args2,
    sample_size,
    H0_true,
    alpha = 0.05,
    Nsim = NULL
) {
  
  if (is.null(Nsim)) {
    Nsim <- nsim_by_n(sample_size)
  }
  
  X <- do.call(
      sim_fun,
      c(list(n = sample_size, Nsim = Nsim), args1)
    )
  
  Y <- do.call(
      sim_fun,
      c(list(n = sample_size, Nsim = Nsim), args2)
    )
  
 
  if (H0_true) {
    
    pop_flat <- as.numeric(X)
    pop_medcouple <- robustbase::mc(pop_flat)
    pop_moors <-moors_kurtosis(pop_flat)
    
  } else {
    
    pop_medcouple <- NA_real_
    pop_moors <- NA_real_
  }
  
  
  res_t <- testing_t_fast(X,Y,alpha)
  
  res_new <- testing_new_tests_fast(X,Y, alpha)
  
  
  tibble(
    type1_t =
      if (H0_true)
        res_t$power
    else
      NA_real_,
    
    type1_bm =
      if (H0_true)
        res_new$bm_rate
    else
      NA_real_,
    
    type1_c2 =
      if (H0_true)
        res_new$c2_rate
    else
      NA_real_,
    
    power_t =
      if (!H0_true)
        res_t$power
    else
      NA_real_,
    
    power_bm =
      if (!H0_true)
        res_new$bm_rate
    else
      NA_real_,
    
    power_c2 =
      if (!H0_true)
        res_new$c2_rate
    else
      NA_real_,
    
    cohen_d = NA_real_,
    
    pop_medcouple = pop_medcouple,
    
    pop_moors = pop_moors
  )
}


# RESULTS GENERATOR

generate_results <- function(
    grids,
    alpha = 0.05,
    drop_list_cols = TRUE
) {
  
  grid_all <- bind_rows(grids) |>
    mutate(
      Nsim = nsim_by_n(sample_size)
    )
  
  
  rows_out <- furrr::future_map(
      seq_len(nrow(grid_all)),
      \(i) {
        
        row <- grid_all[i, ]
        
        sim_out <-
          run_simulation(
            sim_fun = row$sim_fun[[1]],
            
            args1 = row$args1[[1]],
            
            args2 = row$args2[[1]],
            
            sample_size = row$sample_size,
            
            H0_true = row$H0_true,
            
            alpha = alpha,
            
            Nsim = row$Nsim
          )
        
        
        cond_cols <-
          names(row)
        
        if (drop_list_cols) {
          
          cond_cols <- setdiff(
            cond_cols, c("args1","args2","sim_fun")
            )
        }
        
        
        bind_cols(
          as_tibble(row[, cond_cols, drop = FALSE]), sim_out
          )
      },
      
      .options = furrr_options(seed = TRUE, packages = c("matrixTests", "robustbase")
        )
    )
  
  
  bind_rows(rows_out) |>
    arrange(distribution,sample_size)
}


# =====RUN SIMULATIONS==============

message("Starting floor simulation")
results_floor <-
  generate_results(grid_floor)

message("Starting floor simulation")
results_t <-
  generate_results(grid_t)

message("Starting floor simulation")
results_degen <-
  generate_results(grid_degen)

message("Starting floor simulation")
results_pareto <-
  generate_results(grid_pareto)


results_all <- bind_rows(
    results_floor,
    results_t,
    results_degen,
    results_pareto
  )



# ========= TYPE I ERROR HEATMAPS================

make_type1_heatmap <- function(
    results,
    condition,
    y_lab = NULL,
    fill_lab = "Type I error",
    title = NULL,
    subtitle = NULL,
    palette = "gradient",
    low = "honeydew",
    high = "darkgreen",
    limits = c(0, 0.10),
    percent_lab = TRUE,
    alpha_nom = 0.05,
    band = 0.01
) {
  
  if (!condition %in% names(results)) {
    stop(
      glue::glue(
        "Column '{condition}' not found in results. ", "Available: ",
        paste(names(results), collapse = ", ")
      )
    )
  }
  
  type1_cols <- intersect(
    c("type1_t", "type1_bm", "type1_c2"), names(results)
  )
  
  # 1. Calculate Monte Carlo Standard Error (MCSE) around nominal alpha
  
  # One-sided 95% Monte Carlo reference cutoff under the nominal null
  sim_count <- 25000
  mcse_val <- sqrt((alpha_nom * (1 - alpha_nom)) / sim_count)
  upper_threshold <- alpha_nom + qnorm(0.95) * mcse_val
  
  df_long <- results |>
    dplyr::filter(H0_true) |>
    dplyr::select(
      dplyr::any_of(c("distribution", "sample_size", type1_cols, condition))
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(type1_cols), names_to = "test", values_to = "rate"
    ) |>
    dplyr::mutate(
      test = dplyr::case_when(
        test == "type1_t" ~ "Welch's t-test",
        test == "type1_bm" ~ "Brunner-Munzel test",
        test == "type1_c2" ~ "C2 test",
        TRUE ~ test
      ),
      sample_size = factor(sample_size),
      !!condition := factor(.data[[condition]]),
      in_band = dplyr::between(rate, alpha_nom - band, alpha_nom + band)
    )
  
  lab_fun <- if (percent_lab) scales::percent_format(accuracy = 1) else ggplot2::waiver()
  
  scale_fill_fun <- {
    viridis_opts <- c("viridis","magma","plasma","inferno","cividis","rocket", "mako","A", "B", "C", "D")
    
    if (is.character(palette) && length(palette) == 1 && palette %in% viridis_opts) {
      ggplot2::scale_fill_viridis_c(option = palette, limits = limits, labels = lab_fun, name = fill_lab, oob = scales::squish)
    } else if (is.character(palette) && length(palette) == 1 && palette == "gradient") {
      ggplot2::scale_fill_gradient(low = low, high = high, limits = limits, labels = lab_fun, name = fill_lab, oob = scales::squish)
    } else if (is.character(palette) && length(palette) >= 2) {
      ggplot2::scale_fill_gradientn(colours = palette, limits = limits, labels = lab_fun, name = fill_lab, oob = scales::squish)
    } else {
      ggplot2::scale_fill_viridis_c(option = "viridis", limits = limits, labels = lab_fun, name = fill_lab, oob = scales::squish)
    }
  }
  
  df_long$test <- factor(
    df_long$test,
    levels = c("Welch's t-test", "Brunner-Munzel test", "C2 test")
  )
  
  # 2. Format the explanatory caption
  caption_text <- sprintf(
    paste0(
      "Black outlines flag estimates above the one-sided 95%% ",
      "null reference cutoff: > %.4f\n",
      "(Nominal \u03B1 = %.2f + qnorm(0.95) * MCSE)"
    ),
    upper_threshold,
    alpha_nom
  )
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = sample_size, y = .data[[condition]], fill = rate)) +
    
    # Base heatmap layer
    ggplot2::geom_tile(colour = NA) +
    
    # 3. Highlight layer: draws a box around values > (alpha_nom + mcse_val)
    ggplot2::geom_tile(
      data = function(x) dplyr::filter(x, rate > upper_threshold),
      fill = NA,            # Keep the inside transparent to show the gradient
      colour = "black",     # Outline color
      linewidth = 0.8       # Thickness of the box
    ) +
    
    scale_fill_fun +
    
    ggplot2::facet_wrap(~ test, ncol = 2) +
    
    ggplot2::labs(
      title = title %||% "Type I error vs. condition \u00D7 sample size",
      subtitle = subtitle,
      x = "Sample size (n)",
      y = y_lab %||% condition,
      caption = caption_text  # 4. Add caption to the graphic
    ) +
    
    ggplot2::theme_minimal(base_size = 12) +
    
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 18),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6), size = 18),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 22),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 18),
      plot.caption = ggplot2::element_text(hjust = 0, face = "italic", color = "grey30", margin = ggplot2::margin(t = 15), size = 18),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 16, hjust = 1, margin = margin(b=8)),
      legend.title.position = "top",
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.04),
      legend.justification.inside = c(2, 0),
      legend.direction = "vertical",
      legend.key.height = grid::unit(1, "cm"),
      legend.key.width = grid::unit(.5, "cm")
    )
}

# BUILD HEATMAPS

hm_floor <- make_type1_heatmap(
    results_floor,
    condition = "prop_floor",
    y_lab = "Floor Proportion (p)",
    title = "Floor"
  )

hm_degen <- make_type1_heatmap(
    results_degen,
    condition = "sigma",
    y_lab = "Standard Deviation (σ)",
    title = "Near-Degenerate"
  )

hm_t <- make_type1_heatmap(
    results_t,
    condition = "df",
    y_lab = "Degrees of Freedom (ν)",
    title = "Student's t"
  )

hm_pareto <- make_type1_heatmap(
    results_pareto,
    condition = "alpha",
    y_lab = "Tail Index (α)",
    title = "Pareto",
    subtitle = expression(
        "Mean undefined if " * alpha <= 1 ~ "; variance undefined if 1 < " * alpha * " \u2264 2"
      )
  )

# SAVE TYPE I ERROR FIGURE

ggsave(
  filename = file.path(output_dir,"S9.jpeg"),
  plot = hm_floor,
  width = 12,
  height = 9,
  dpi = 800,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir,"S10.jpeg"),
  plot = hm_degen,
  width = 12,
  height = 9,
  dpi = 800,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir,"S11.jpeg"),
  plot = hm_pareto,
  width = 12,
  height = 9,
  dpi = 800,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir,"S12.jpeg"),
  plot = hm_t,
  width = 12,
  height = 9,
  dpi = 800,
  bg = "white"
)

# ========Worst Case Type I Error======================================

mk_type1_gt <- function(
    tbl,
    title,
    subtitle,
    alpha_nom = 0.05,
    band = 0.01
) {
  
  tbl |>
    gt() |>
    tab_header(title = title,subtitle = subtitle) |>
    cols_align(align = "center",columns = everything()) |>
    fmt_percent(
      columns = any_of(c("Type I (t-test)","Type I (BM)","Type I (C2)")),
      decimals = 2
    )
}


alpha_nom <- 0.05
sim_count <- 25000
band <- sqrt((alpha_nom * (1 - alpha_nom)) / sim_count)
per_dist_limit <- 3


type1_candidates <-
  results_all |>
  filter(H0_true) |>
  mutate(
    liberal_t =
      if ("type1_t" %in% names(results_all))
        type1_t > (alpha_nom + band)
    else
      FALSE,
    
    liberal_bm =
      if ("type1_bm" %in% names(results_all))
        type1_bm > (alpha_nom + band)
    else
      FALSE,
    
    liberal_c2 =
      if ("type1_c2" %in% names(results_all))
        type1_c2 > (alpha_nom + band)
    else
      FALSE,
    
    
    cons_t =
      if ("type1_t" %in% names(results_all))
        type1_t < (alpha_nom - band)
    else
      FALSE,
    
    cons_bm =
      if ("type1_bm" %in% names(results_all))
        type1_bm < (alpha_nom - band)
    else
      FALSE,
    
    cons_c2 =
      if ("type1_c2" %in% names(results_all))
        type1_c2 < (alpha_nom - band)
    else
      FALSE,
    
    
    any_flag =   liberal_t | liberal_bm | liberal_c2 | cons_t | cons_bm | cons_c2,
    
    
    max_dev = pmax(
        if ("type1_t" %in% names(results_all))
          abs(type1_t - alpha_nom)
        else
          0,
        
        if ("type1_bm" %in% names(results_all))
          abs(type1_bm - alpha_nom)
        else
          0,
        
        if ("type1_c2" %in% names(results_all))
          abs(type1_c2 - alpha_nom)
        else
          0,
        
        na.rm = TRUE
      ),
    
    
    Parameter = case_when(
        distribution == "Floor" ~ "p (floor proportion)",
        distribution == "Pareto" ~ "α (tail index)",
        distribution == "t-heavy" ~ "df (Student’s t)",
        distribution == "Degenerate" ~ "σ (SD)"
      ),
    
    
    `Parameter level` = case_when(
        distribution == "Floor" ~ as.numeric(prop_floor),
        distribution == "Pareto" ~ as.numeric(alpha),
        distribution == "t-heavy" ~ as.numeric(df),
        distribution == "Degenerate" ~ as.numeric(sigma)
      )
  ) |>
  filter(any_flag) |>
  group_by(distribution) |>
  arrange(desc(max_dev),.by_group = TRUE)


if (is.finite(per_dist_limit)) {
  
  type1_candidates <- type1_candidates |>
    slice_head(n = per_dist_limit)
}


worst_type1_tbl <- type1_candidates |>
  ungroup() |>
  transmute(
    
    Distribution = distribution,
    
    `Sample size` = sample_size,
    
    Parameter = Parameter,
    
    `Parameter level` = `Parameter level`,
    
    `Type I (t-test)` =
      if ("type1_t" %in% names(type1_candidates))
        type1_t
    else
      NA,
    
    `Type I (BM)` =
      if ("type1_bm" %in% names(type1_candidates))
        type1_bm
    else
      NA,
    
    `Type I (C2)` =
      if ("type1_c2" %in% names(type1_candidates))
        type1_c2
    else
      NA
  ) |>
  arrange(Distribution, desc(`Type I (t-test)`)) |>
  select(
    where(~ !all(is.na(.x)))
  )


type1_table <- mk_type1_gt(
    worst_type1_tbl,
    title = "Type I error — Worst Cases",
    subtitle = "Rows with liberal or conservative deviations from α across all tested procedures"
  )


gtsave(
  data = type1_table, 
  filename = file.path(output_dir, "Worst_Case_Type1_Error.html")
)


# =================POWER CURVES================================

plot_power_facet_grid <- function(
    results,
    family,
    shape_param,
    effect_param,
    levels = NULL,
    calibrated_only = FALSE, # Changed default to FALSE
    alpha = 0.05
) {
  
  # Ensure Nsim exists for the standard error calculation
  if (!"Nsim" %in% names(results)) {
    results$Nsim <- 25000
  }
  
  # 0. Map display names to actual distribution names
  distribution_name <- switch(
    family,
    "Floor" = "Floor",
    "Near-Degenerate Normal" = "Degenerate",
    "Student's t" = "t-heavy",
    "Pareto" = "Pareto",
    family
  )
  
  # 1. Alternative hypothesis data
  df_h1 <- results |>
    dplyr::filter(distribution == distribution_name, !H0_true)
  
  if (!nrow(df_h1)) {
    stop(sprintf("\nNo H1 rows found for '%s'.\n", family))
  }
  
  # 2. Null hypothesis calibration using MCSE threshold
  type1_cols <- intersect(c("type1_t", "type1_bm", "type1_c2"), names(results))
  
  sim_count <- 25000
  mcse_alpha <- sqrt((alpha * (1 - alpha)) / sim_count)
  strict_threshold <- alpha + 1.645 * mcse_alpha
  
  df_h0 <- results |>
    dplyr::filter(distribution == distribution_name, H0_true) |>
    dplyr::select(sample_size, dplyr::all_of(shape_param), dplyr::all_of(type1_cols)) |>
    dplyr::mutate(
      ok_t = if ("type1_t" %in% names(results)) type1_t <= strict_threshold else NA,
      ok_bm = if ("type1_bm" %in% names(results)) type1_bm <= strict_threshold else NA,
      ok_c2 = if ("type1_c2" %in% names(results)) type1_c2 <= strict_threshold else NA
    ) |>
    dplyr::select(sample_size, dplyr::all_of(shape_param), ok_t, ok_bm, ok_c2)
  
  # 3. Join H1 and H0 data
  df_joined <- df_h1 |>
    dplyr::left_join(df_h0, by = c("sample_size", shape_param))
  
  # 4. Convert power columns to long format
  powers_to_pivot <- intersect(c("power_t", "power_bm", "power_c2"), names(df_joined))
  
  if (!length(powers_to_pivot)) {
    stop(sprintf("\nNo power columns found for '%s'.\n", family))
  }
  
  df_long <- df_joined |>
    dplyr::select(
      distribution, sample_size, Nsim, 
      dplyr::all_of(shape_param), dplyr::all_of(effect_param),
      dplyr::all_of(powers_to_pivot), ok_t, ok_bm, ok_c2
    ) |>
    tidyr::pivot_longer(dplyr::all_of(powers_to_pivot), names_to = "which_power", values_to = "power") |>
    dplyr::mutate(
      test = dplyr::case_when(
        which_power == "power_t" ~ "Welch's t-test",
        which_power == "power_bm" ~ "Brunner-Munzel test",
        which_power == "power_c2" ~ "C2 test",
        TRUE ~ NA_character_
      ),
      ok = dplyr::case_when(
        test == "Welch's t-test" ~ ok_t,
        test == "Brunner-Munzel test" ~ ok_bm,
        test == "C2 test" ~ ok_c2,
        TRUE ~ FALSE
      ),
      mcse = sqrt((power * (1 - power)) / sim_count)
    )
  
  # FORCE EQUIDISTANT X-AXIS SPACING
  # Convert numeric sample_size to an ordered factor
  df_long$sample_size <- factor(df_long$sample_size, levels = sort(unique(df_long$sample_size)))
  
  # 5. Filter to requested parameter levels
  if (!is.null(levels)) {
    df_long <- df_long |> dplyr::filter(.data[[shape_param]] %in% levels)
  }
  
  # 6. Calibration filter
  if (calibrated_only) {
    df_long <- df_long |> dplyr::filter(ok %in% TRUE)
  }
  
  # 7. Useful empty-data error
  if (!nrow(df_long)) {
    stop(sprintf("No rows left to plot for '%s' after filtering.", family))
  }
  
  # 8. Lock requested parameter ordering
  if (!is.null(levels)) {
    df_long[[shape_param]] <- factor(df_long[[shape_param]], levels = levels)
  }
  
  # Extract max SE for the caption
  max_mcse <- max(df_long$mcse, na.rm = TRUE)
  caption_text <- sprintf("Maximum Monte Carlo standard error \u2248 %.4f (Nsim = %d)", max_mcse, sim_count)
  
  if (calibrated_only) {
    caption_text <- paste0(caption_text, sprintf("\nTests exceeding Type I error of %.4f are excluded.", strict_threshold))
  }
  
  # 9. TEST COLORS
  test_cols <- c("Welch's t-test" = "#0072B2", "Brunner-Munzel test" = "#D55E00", "C2 test" = "#009E73")
  shape_vals <- c("Welch's t-test" = 21, "Brunner-Munzel test" = 22, "C2 test" = 23)
  line_vals <- c("Welch's t-test" ="solid", "Brunner-Munzel test" ="dotted", "C2 test" ="twodash")
  
  # 10. Facet labels
  my_labeller <- ggplot2::labeller(
    prop_floor = function(x) paste("p =", x),
    alpha = function(x) paste("α =", x),
    df = function(x) paste("ν =", x),
    sigma = function(x) paste("σ =", x),
    shift_mu = function(x) paste("Effect Size =", x),
    d_target = function(x) paste("Effect Size =", x),
    r = function(x) paste("Effect Size =", x),
    .default = ggplot2::label_value
  )
  
  
  # 11. Build plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample_size, y = power, color = test, group = interaction(test, .data[[shape_param]]))) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkred", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(linetype = test), linewidth = 0.8, alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(shape = test), size = 2.5, fill = "white", stroke = 1.2) +
    ggplot2::scale_color_manual(name = "Test", values = test_cols, breaks = names(test_cols)) +
    ggplot2::scale_shape_manual(name = "Test", values = shape_vals, breaks = names(shape_vals)) +
    ggplot2::scale_linetype_manual(name = "Test", values = line_vals, breaks = names(line_vals)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    
    # NEW: Discrete x-axis for perfectly even spacing
    ggplot2::scale_x_discrete() +
    
    ggplot2::facet_grid(rows = ggplot2::vars(.data[[shape_param]]), cols = ggplot2::vars(.data[[effect_param]]), labeller = my_labeller) +
    ggplot2::labs(
      title = family,
      subtitle = "Red line indicates 80% power",
      x = "Sample Size (n)",
      y = "Power",
      caption = caption_text
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 18),
      strip.background = ggplot2::element_rect(fill = "gray90"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 22),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 18, color = "darkred"),
      plot.caption = ggplot2::element_text(hjust = 0, face = "italic", color = "grey30", margin = ggplot2::margin(t = 15), size = 18),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 18),
      legend.text = ggplot2::element_text(size = 18)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1)
    )
  
  return(p)
}

# BUILD FACET GRID POWER PLOTS

grid_plot_floor <-
  plot_power_facet_grid(
    results = results_floor,
    family = "Floor",
    shape_param = "prop_floor",
    effect_param = "shift_mu",
    levels =
      c(0.8,0.5,0.2)
  )


grid_plot_degen <-
  plot_power_facet_grid(
    results = results_degen,
    family = "Near-Degenerate Normal",
    shape_param = "sigma",
    effect_param = "d_target",
    levels = c(0.0001,0.001,0.05)
  )


grid_plot_t <-
  plot_power_facet_grid(
    results = results_t,
    family = "Student's t",
    shape_param = "df",
    effect_param = "d_target",
    levels =
      c(0.5,2,5)
  )


grid_plot_pareto <-
  plot_power_facet_grid(
    results = results_pareto,
    family = "Pareto",
    shape_param = "alpha",
    effect_param = "r",
    levels =
      c(0.6,1.10,2)
  )


grid_plot_floor
grid_plot_degen
grid_plot_t
grid_plot_pareto


# ===== Save Power Curves================================

ggsave(
  file.path(output_dir, "S1.jpeg"),
  plot = grid_plot_floor,
  width = 12,
  height = 8,
  dpi = 1200,
  bg = "white"
)

ggsave(
  file.path(output_dir, "S2.jpeg"),
  plot = grid_plot_degen,
  width = 12,
  height = 8,
  dpi = 1200,
  bg = "white"
)

ggsave(
  file.path(output_dir, "S4.jpeg"),
  plot = grid_plot_t,
  width = 12,
  height = 8,
  dpi = 1200,
  bg = "white"
)

ggsave(
  file.path(output_dir, "S3.jpeg"),
  plot = grid_plot_pareto,
  width = 12,
  height = 8,
  dpi = 1200,
  bg = "white"
)


# =======CALIBRATED POWER & MIN N TABLES=======================================================

alpha_nom <- 0.05
sim_count <- 25000
mcse_null <- sqrt((alpha_nom * (1 - alpha_nom)) / sim_count)
one_sided95_band <- qnorm(0.95) * mcse_null
band <- one_sided95_band

calibration_mode <- "liberal_only"

ok_fun <- function(x) {
  !is.na(x) & x <= (alpha_nom + one_sided95_band)
}


pp <- rlang::sym(pareto_param)


h0_flags <- results_all |>
  filter(H0_true) |>
  transmute(
    distribution,
    sample_size,
    
    prop_floor = coalesce( prop_floor, NA_real_),
    
    sigma = coalesce(sigma, NA_real_),
    
    df = coalesce(df, NA_real_),
    
    alpha = coalesce(alpha, NA_real_),
    
    ok_t = ok_fun(type1_t),
    
    ok_bm =
      if (
        "type1_bm" %in% names(results_all)
      )
        ok_fun(type1_bm)
    else
      NA,
    
    ok_c2 =
      if (
        "type1_c2" %in% names(results_all)
      )
        ok_fun(type1_c2)
    else
      NA
  )


worst_power_tbl <- results_all |>
  filter(!H0_true) |>
  left_join(
    h0_flags,
    by = c("distribution", "sample_size", "prop_floor", "sigma", "df", "alpha")
  ) |>
  mutate(
    ok_t = coalesce(ok_t, FALSE),
    ok_bm = coalesce(ok_bm, FALSE),
    ok_c2 = coalesce(ok_c2, FALSE),

    `Parameter level` = case_when(
        distribution == "Floor" ~ prop_floor,
        distribution == "Degenerate" ~ sigma,
        distribution == "t-heavy" ~ df,
        distribution == "Pareto" ~ !!pp
      ),
    
    
    `Parameter` = case_when(
        distribution == "Floor" ~ "p (floor proportion)",
        distribution == "Degenerate" ~ "σ (SD)",
        distribution == "t-heavy" ~ "df (Student’s t)",
        distribution == "Pareto" ~
          ifelse(pareto_param == "alpha", "α (tail index)", "r (scale ratio)")
      ),
    
    
    Effect =case_when(
        distribution == "Floor" ~ shift_mu,
        distribution == "Degenerate" ~ d_target,
        distribution == "t-heavy" ~ d_target,
        distribution == "Pareto" ~ r
      ),
    
    
    power_t_cal = if_else(ok_t, power_t, NA_real_),
    
    power_bm_cal =
      if ("power_bm" %in% names(results_all))
        if_else(ok_bm, power_bm, NA_real_)
    else
      NA_real_,
    
    power_c2_cal =
      if ("power_c2" %in% names(results_all))
        if_else(ok_c2, power_c2, NA_real_)
    else
      NA_real_
  ) |>
  group_by(distribution, `Parameter level`, sample_size) |>
  mutate(
    has_nonzero = any(Effect != 0 & is.finite(Effect)
      )
  ) |>
  filter(
    if_else(has_nonzero, Effect != 0, TRUE)
  ) |>
  group_modify(
    ~ {
      
      if (
        .y$distribution == "Pareto"
      ) {
        
        dplyr::slice_min(
          .x,
          order_by = Effect,
          n = 1,
          with_ties = FALSE
        )
        
      } else {
        
        dplyr::slice_max(
          .x,
          order_by = Effect,
          n = 1,
          with_ties = FALSE
        )
      }
    }
  ) |>
  ungroup() |>
  filter(
    pmin(power_t_cal, power_bm_cal, power_c2_cal, na.rm = TRUE) < 0.30
  ) |>
  group_by(distribution) |>
  arrange(pmin(power_t_cal, power_bm_cal, power_c2_cal, na.rm = TRUE),
    .by_group = TRUE
  ) |>
  slice_head(n = per_dist_limit) |>
  ungroup() |>
  transmute(
    
    Distribution = distribution,
    
    `Sample size` = sample_size,
    
    `Parameter` = `Parameter`,
    
    `Parameter level` = `Parameter level`,
    
    `Effect / contrast` = Effect,
    
    `Power (t-test)` = power_t_cal,
    
    `Power (BM)` = power_bm_cal,
    
    `Power (C2)` = power_c2_cal
  ) |>
  select(
    where(~ !all(is.na(.x)))
  )


worst_power_tbl2 <- worst_power_tbl |>
  mutate(
    min_pwr =
      pmin( `Power (t-test)`,
        if (
          "Power (BM)" %in% names(worst_power_tbl)
        )
          `Power (BM)`
        else
          NA,
        
        if (
          "Power (C2)" %in% names(worst_power_tbl)
        )
          `Power (C2)`
        else
          NA,
        
        na.rm = TRUE
      )
  )


power_extremes <- worst_power_tbl2 |>
  select(-min_pwr) |>
  gt() |>
  tab_header(
    title = "Power — Worst Cases"
  ) |>
  fmt_number(
    columns = any_of(c("Power (t-test)", "Power (BM)", "Power (C2)")),
    decimals = 4
  ) |>
  tab_style(style = cell_text(weight = "bold"),
    
    locations = cells_body(
        columns = "Power (t-test)",
        rows = abs(`Power (t-test)` - worst_power_tbl2$min_pwr) < 1e-9
      )
  )



if (
  "Power (BM)" %in% names(worst_power_tbl2)
) {
  
  power_extremes <- power_extremes |>
    tab_style(
      style = cell_text(weight = "bold"),
      
      locations = cells_body(columns = "Power (BM)",
          
          rows = abs(`Power (BM)` - worst_power_tbl2$min_pwr) < 1e-9
        )
    )
}


if (
  "Power (C2)" %in% names(worst_power_tbl2)
) {
  
  power_extremes <- power_extremes |>
    tab_style(style = cell_text(weight = "bold"),
      
      locations = cells_body(columns = "Power (C2)",
          
          rows = abs(`Power (C2)` - worst_power_tbl2$min_pwr) < 1e-9
        )
    )
}


gtsave(
  data = power_extremes, 
  filename = file.path(output_dir, "Worst_Case_Power.html")
)


# ===============Power Tables=========================

render_dist_power_table <- function(
    data,
    dist_name
) {
  
  sub_data <- data |>  
    filter(
      Distribution == dist_name
    )
  
  
  tbl <- sub_data |>
    select(
      -Distribution,
      -min_pwr
    ) |>
    gt() |>
    tab_header(
      title = paste( "Power — Worst Cases:", dist_name),
      
      subtitle = paste( "Filtered poor performance scenarios for", dist_name)
    ) |>
    cols_align(align = "center", columns = everything()) |>
    fmt_number(
      columns =
        any_of(c("Power (t-test)", "Power (BM)", "Power (C2)")),
      decimals = 4
    ) |>
    tab_style(style = cell_text(weight = "bold"),
      
      locations =
        cells_body(columns = "Power (t-test)",
          
          rows = abs(sub_data$`Power (t-test)` -  sub_data$min_pwr) < 1e-9
        )
    )
  
  
  if (
    "Power (BM)" %in% names(sub_data)
  ) {
    
    tbl <- tbl |>
      tab_style(style = cell_text(weight = "bold"),
        
        locations = cells_body(columns = "Power (BM)",
            
            rows = abs(sub_data$`Power (BM)` - sub_data$min_pwr) < 1e-9
          )
      )
  }
  
  
  if (
    "Power (C2)" %in%
    names(sub_data)
  ) {
    
    tbl <- tbl |>
      tab_style(
        style = cell_text(weight = "bold"),
        
        locations = cells_body(columns = "Power (C2)",
            
            rows = abs(sub_data$`Power (C2)` - sub_data$min_pwr) < 1e-9
          )
      )
  }
  
  
  return(tbl)
}


floor_power_tbl <- render_dist_power_table(worst_power_tbl2, "Floor")

pareto_power_tbl <- render_dist_power_table(worst_power_tbl2, "Pareto")

t_heavy_power_tbl <- render_dist_power_table(worst_power_tbl2, "t-heavy")

degen_power_tbl <- render_dist_power_table(worst_power_tbl2,"Degenerate")

gtsave(
  data = floor_power_tbl, filename = file.path(output_dir, "Floor_Power.html")
)
gtsave(
  data = pareto_power_tbl, filename = file.path(output_dir, "Pareto_Power.html")
)
gtsave(
  data = t_heavy_power_tbl, filename = file.path(output_dir, "Students_t_Power.html")
)
gtsave(
  data = degen_power_tbl, filename = file.path(output_dir, "Near_Degen_Power.html")
)





# ===================== MINIMUM n TABLES =====================

pretty_param <- function(param) {
  labels <- c(
    prop_floor = "Floor Proportion (p)",
    sigma      = "Standard Deviation (σ)",
    df         = "Degrees of Freedom (ν)",
    alpha      = "Tail Index (α)"
  )
  
  if (param %in% names(labels)) {
    unname(labels[[param]])
  } else {
    param
  }
}


make_min_n_table_wide <- function(
    results_list,
    families = c("Floor", "Degenerate", "t-heavy", "Pareto"),
    params = c("prop_floor", "sigma", "df", "alpha"),
    effect_params = c("shift_mu", "d_target", "d_target", "r"),
    target = 0.80,
    alpha_nom = 0.05,
    sim_count = 25000,
    calibration = c("liberal_only", "strict", "none"),
    numeric_param_sort = FALSE
) {
  
  calibration <- match.arg(calibration)
  
  mcse_null <- sqrt(alpha_nom * (1 - alpha_nom) / sim_count)
  one_sided95_band <- qnorm(0.95) * mcse_null
  two_sided95_band <- qnorm(0.975) * mcse_null
  
  ok_fun <- switch(
    calibration,
    
    liberal_only = function(x) {
      !is.na(x) & x <= alpha_nom + one_sided95_band
    },
    
    strict = function(x) {
      !is.na(x) &
        dplyr::between(
          x,
          alpha_nom - two_sided95_band,
          alpha_nom + two_sided95_band
        )
    },
    
    none = function(x) {
      rep(TRUE, length(x))
    }
  )
  
  one_family <- function(res, fam, param, eff_param) {
    
    grp_label <- paste0(
      fam, " — ", pretty_param(param)
    )
    
    # Start with the alternative-hypothesis power results.
    power_data <- res |>
      dplyr::filter(
        distribution == fam,
        !H0_true
      ) |>
      dplyr::select(
        sample_size,
        dplyr::all_of(c(param, eff_param)),
        dplyr::any_of(
          c("power_t", "power_bm", "power_c2")
        )
      )
    
    if (!"power_t" %in% names(power_data)) {
      stop("Missing power_t for ", fam)
    }
    
    # Check column availability BEFORE entering summarise().
    # This avoids confusing the data frame with the Student's t
    # degrees-of-freedom column named df.
    has_bm <- "power_bm" %in% names(power_data)
    has_c2 <- "power_c2" %in% names(power_data)
    
    if (calibration == "none") {
      
      # Uncalibrated: minimum n depends only on power.
      power_data <- power_data |>
        dplyr::mutate(
          ok_t = TRUE,
          ok_bm = TRUE,
          ok_c2 = TRUE
        )
      
    } else {
      
      # Calibrated: match each power condition to its H0
      # Type I error at the same parameter level and n.
      has_type1_t <- "type1_t" %in% names(res)
      has_type1_bm <- "type1_bm" %in% names(res)
      has_type1_c2 <- "type1_c2" %in% names(res)
      
      h0 <- res |>
        dplyr::filter(
          distribution == fam,
          H0_true
        ) |>
        dplyr::select(
          sample_size,
          dplyr::all_of(param),
          dplyr::any_of(
            c("type1_t", "type1_bm", "type1_c2")
          )
        ) |>
        dplyr::mutate(
          ok_t = if (has_type1_t) {
            ok_fun(type1_t)
          } else {
            rep(FALSE, dplyr::n())
          },
          
          ok_bm = if (has_type1_bm) {
            ok_fun(type1_bm)
          } else {
            rep(FALSE, dplyr::n())
          },
          
          ok_c2 = if (has_type1_c2) {
            ok_fun(type1_c2)
          } else {
            rep(FALSE, dplyr::n())
          }
        ) |>
        dplyr::select(
          sample_size,
          dplyr::all_of(param),
          ok_t,
          ok_bm,
          ok_c2
        )
      
      power_data <- power_data |>
        dplyr::left_join(
          h0,
          by = c("sample_size", param)
        )
    }
    
    grab_min <- function(ok, pwr, n) {
      
      eligible <-
        !is.na(ok) &
        ok &
        !is.na(pwr) &
        pwr >= target &
        !is.na(n)
      
      if (any(eligible)) {
        min(n[eligible])
      } else {
        NA_real_
      }
    }
    
    power_data |>
      dplyr::group_by(
        .data[[param]],
        .data[[eff_param]]
      ) |>
      dplyr::summarise(
        `Welch t-test min n` =
          grab_min(ok_t, power_t, sample_size),
        
        `BM min n` =
          if (has_bm) {
            grab_min(ok_bm, power_bm, sample_size)
          } else {
            NA_real_
          },
        
        `C2 min n` =
          if (has_c2) {
            grab_min(ok_c2, power_c2, sample_size)
          } else {
            NA_real_
          },
        
        .groups = "drop"
      ) |>
      dplyr::mutate(
        Group = grp_label,
        `Parameter level` = .data[[param]],
        `Effect Size` = .data[[eff_param]]
      ) |>
      dplyr::select(
        Group,
        `Parameter level`,
        `Effect Size`,
        `Welch t-test min n`,
        `BM min n`,
        `C2 min n`
      )
  }
  
  tbl <- purrr::pmap_dfr(
    list(
      res = results_list,
      fam = families,
      param = params,
      eff_param = effect_params
    ),
    one_family
  ) |>
    dplyr::select(
      where(~ !all(is.na(.x)))
    ) |>
    dplyr::mutate(
      dplyr::across(
        ends_with("min n"),
        ~ ifelse(is.na(.x), "—", as.character(.x))
      )
    )
  
  if (numeric_param_sort) {
    tbl <- tbl |>
      dplyr::arrange(
        Group,
        suppressWarnings(
          as.numeric(`Parameter level`)
        ),
        `Effect Size`
      )
  } else {
    tbl <- tbl |>
      dplyr::arrange(
        Group,
        `Parameter level`,
        `Effect Size`
      )
  }
  
  calibration_label <- switch(
    calibration,
    
    liberal_only = sprintf(
      "type I ≤ %.4f",
      alpha_nom + one_sided95_band
    ),
    
    strict = sprintf(
      "type I ∈ [%.4f, %.4f]",
      alpha_nom - two_sided95_band,
      alpha_nom + two_sided95band
    ),
    
    none = "none; power only"
  )
  
  gt::gt(
    tbl,
    groupname_col = "Group"
  ) |>
    gt::cols_label(
      `Parameter level` = "Parameter level",
      `Effect Size` = "Effect Size",
      `Welch t-test min n` = "Welch t-test",
      `BM min n` = "Brunner-Munzel",
      `C2 min n` = "C2 test"
    ) |>
    gt::tab_header(
      title = sprintf(
        "Minimum n for %.0f%% power (%s)",
        100 * target,
        calibration_label
      )
    ) |>
    gt::tab_options(
      table.font.size = gt::px(13)
    )
}


render_min_n_dist_table_1 <- function(
    full_tbl,
    group_pattern,
    subtitle
) {
  
  sub_data <- full_tbl |>
    dplyr::filter(
      grepl(group_pattern, Group, fixed = TRUE)
    ) |>
    dplyr::mutate(
      parameter_label = sub(
        "^.* — ",
        "",
        Group
      ),
      `Parameter level` = paste0(
        parameter_label,
        " = ",
        `Parameter level`
      )
    )
  
  gt::gt(
    sub_data,
    groupname_col = "Parameter level"
  ) |>
    gt::cols_align(
      align = "center",
      columns = everything()
    ) |>
    gt::cols_label(
      `Effect Size` = "Effect Size",
      `Welch t-test min n` = "Welch t-test",
      `BM min n` = "Brunner-Munzel",
      `C2 min n` = "C2 test"
    ) |>
    gt::tab_header(
      title = paste(
        "Minimum n for 80% Power:",
        group_pattern
      ),
      subtitle = subtitle
    ) |>
    gt::tab_options(
      table.align = "center"
    ) |>
    gt::cols_hide(
      columns = c(Group, parameter_label)
    ) |>
    gt::tab_style(
      style = gt::cell_text(
        weight = "bold",
        align = "center"
      ),
      locations = gt::cells_row_groups()
    ) |>
    gt::text_transform(
      locations = gt::cells_body(),
      fn = function(x) {
        ifelse(
          x == "—",
          "<span style='color:lightgrey;'>—</span>",
          x
        )
      }
    )
}


# Keep this order aligned with families, params,
# and effect_params in make_min_n_table_wide().
results_list <- list(
  results_floor,
  results_degen,
  results_t,
  results_pareto
)

# Change this to "liberal_only" or "strict" when wanted.
calibration_mode <- "liberal_only"

min_n <- make_min_n_table_wide(
  results_list,
  calibration = calibration_mode
)

min_n_data <- min_n$`_data`

subtitle_text <- switch(
  calibration_mode,
  
  none =
    "Uncalibrated: power only; Type I error not screened",
  
  liberal_only =
    sprintf(
      "Calibrated: matched null Type I error \u2264 %.4f (one-sided 95%% Monte Carlo cutoff)",
      0.05 + qnorm(0.95) * sqrt(0.05 * 0.95 / 25000)
    ),
  
  strict =
    "Calibrated: Type I error within specified band"
)

file_tag <- switch(
  calibration_mode,
  none = "uncalibrated",
  liberal_only = "liberal_only",
  strict = "strict"
)

degen_min_n_tbl <- render_min_n_dist_table_1(
  min_n_data,
  "Degenerate",
  subtitle_text
)

floor_min_n_tbl <- render_min_n_dist_table_1(
  min_n_data,
  "Floor",
  subtitle_text
)

pareto_min_n_tbl <- render_min_n_dist_table_1(
  min_n_data,
  "Pareto",
  subtitle_text
)

t_heavy_min_n_tbl <- render_min_n_dist_table_1(
  min_n_data,
  "t-heavy",
  subtitle_text
)

gt::gtsave(
  data = degen_min_n_tbl,
  filename = file.path(
    output_dir,
    paste0("S13", ".pdf")
  )
)

gt::gtsave(
  data = floor_min_n_tbl,
  filename = file.path(
    output_dir,
    paste0("S12", ".pdf")
  )
)

gt::gtsave(
  data = pareto_min_n_tbl,
  filename = file.path(
    output_dir,
    paste0("S15",  ".pdf")
  )
)

gt::gtsave(
  data = t_heavy_min_n_tbl,
  filename = file.path(
    output_dir,
    paste0("S16",".pdf")
  )
)

# =====================Supplementary Results=====================

appendix_prepare <-
  function(
    results_all,
    alpha_nom = 0.05,
    band = 0.01
  ) {
    
    df <- results_all |>
      mutate(
        Parameter = case_when(
            
            distribution == "Floor" ~ "p (floor proportion)",
            
            distribution == "Degenerate" ~ "σ (SD)",
            
            distribution == "t-heavy" ~ "df (Student’s t)",
            
            distribution == "Pareto" ~ "α (tail index)"
          ),
        
        Parameter_level = case_when(
            
            distribution == "Floor" ~ as.character(prop_floor),
            
            distribution == "Degenerate" ~ as.character(sigma),
            
            distribution == "t-heavy" ~ as.character(df),
            
            distribution == "Pareto" ~ as.character(alpha)
          ),
        
        Effect_or_contrast = case_when(
            
            distribution == "Floor" ~ shift_mu,
            
            distribution == "Degenerate" ~ d_target,
            
            distribution == "t-heavy" ~ d_target,
            
            distribution == "Pareto" ~ r
          )
      )
    
    
    h0_flags <- df |>
      filter(H0_true) |>
      transmute(
        
        distribution,
        sample_size,
        
        prop_floor = coalesce(prop_floor, NA_real_),
        
        sigma = coalesce(sigma,NA_real_),
        
        df =coalesce(df,NA_real_),
        
        alpha =coalesce(alpha,NA_real_),
        
        ok_t = is.na(type1_t) | type1_t <= (alpha_nom + band),
        
        ok_bm = if ("type1_bm" %in% names(df))
            is.na(type1_bm) | type1_bm <= (alpha_nom + band)
        else
          NA,
        
        ok_c2 = if ("type1_c2" %in% names(df))
            is.na(type1_c2) | type1_c2 <= (alpha_nom + band)
        else
          NA
      )
    
    
    df2 <- df |>
      left_join(h0_flags,
        by = c(
            "distribution",
            "sample_size",
            "prop_floor",
            "sigma",
            "df",
            "alpha"
          )
      ) |>
      mutate(
        
        ok_t = coalesce(ok_t, FALSE),
        
        ok_bm = coalesce(ok_bm, FALSE),
        
        ok_c2 = coalesce(ok_c2, FALSE),
        
        power_t_cal = if_else( !H0_true & ok_t, power_t, NA_real_),
        
        power_bm_cal =
          if ("power_bm" %in% names(df))
            
            if_else(
              !H0_true & ok_bm,
              power_bm,
              NA_real_
            )
        else
          NA_real_,
        
        power_c2_cal =
          if ("power_c2" %in% names(df))
            
            if_else(
              !H0_true & ok_c2,
              power_c2,
              NA_real_
            )
        else
          NA_real_,
        
        
        min_cal_pow =
          pmin(
            power_t_cal,
            power_bm_cal,
            power_c2_cal,
            na.rm = TRUE
          ),
        
        max_cal_pow =
          pmax(
            power_t_cal,
            power_bm_cal,
            power_c2_cal,
            na.rm = TRUE
          ),
        
        
        Smaller_calib = case_when(
            
            power_t_cal == min_cal_pow ~ "t-test",
            
            power_bm_cal == min_cal_pow ~ "BM",
            
            power_c2_cal == min_cal_pow ~ "C2",
            
            TRUE ~ NA_character_
          ),
        
        
        Larger_calib = case_when(
            
            power_t_cal == max_cal_pow ~ "t-test",
            
            power_bm_cal == max_cal_pow ~ "BM",
            
            power_c2_cal == max_cal_pow ~ "C2",
            
            TRUE ~ NA_character_
          )
      )
    
    
    df_out <- df2 |>
      transmute(
        
        Distribution = distribution,
        
        Sample_size = sample_size,
        
        H0_true = H0_true,
        
        Parameter,
        
        Parameter_level,
        
        Effect_or_contrast,
        
        Type1_t = type1_t,
        
        Type1_bm =
          if ( "type1_bm" %in% names(df2))
            type1_bm
        else
          NA,
        
        Type1_c2 =
          if ("type1_c2" %in% names(df2))
            type1_c2
        else
          NA,
        
        Power_t = power_t,
        
        Power_bm =
          if ("power_bm" %in% names(df2))
            power_bm
        else
          NA,
        
        Power_c2 =
          if ("power_c2" %in% names(df2))
            power_c2
        else
          NA,
        
        Power_t_cal = power_t_cal,
        
        Power_bm_cal = power_bm_cal,
        
        Power_c2_cal = power_c2_cal,
        
        Cohens_d_abs = cohen_d,
        
        Calibrated_t = ok_t,
        
        Calibrated_bm = ok_bm,
        
        Calibrated_c2 = ok_c2,
        
        Smaller_calib = Smaller_calib,
        
        Larger_calib = Larger_calib,
        
        pop_medcouple = pop_medcouple,
        
        pop_moors = pop_moors,
        
        prop_floor = prop_floor,
        
        sigma = sigma,
        
        df = df,
        
        alpha = alpha,
        
        r = r
      ) |>
      arrange(
        Distribution,
        Parameter,
        suppressWarnings(
          as.numeric(
            Parameter_level
          )
        ),
        Sample_size
      ) |>
      select(
        where(~ !all(is.na(.x))
        )
      )
    
    
    return(df_out)
  }


appendix_all <- appendix_prepare(
    results_all,
    alpha_nom = alpha_nom,
    band = band
  )


appendix_floor <- appendix_all |>
  filter(Distribution == "Floor")

appendix_degen <- appendix_all |>
  filter(Distribution == "Degenerate")

appendix_t <- appendix_all |>
  filter(Distribution == "t-heavy")

appendix_pareto <- appendix_all |>
  filter(Distribution == "Pareto")

# WORST SIZE

worst_size <-
  results_all |>
  filter(H0_true) |>
  mutate(
    
    dev_t = if_else(
        is.finite(type1_t),
        abs(type1_t - alpha_nom),
        NA_real_
      ),
    
    dev_bm =
      if ("type1_bm" %in% names(results_all))
        if_else(
          is.finite(type1_bm),
          abs(type1_bm - alpha_nom),
          NA_real_
        )
    else
      NA_real_,
    
    dev_c2 =
      if ( "type1_c2" %in% names(results_all))
        if_else(
          is.finite(type1_c2),
          abs(type1_c2 - alpha_nom),
          NA_real_
        )
    else
      NA_real_
  ) |>
  group_by(distribution) |>
  slice_max(
    order_by = pmax(
        coalesce(dev_t, -Inf),
        
        coalesce(dev_bm, -Inf),
        
        coalesce(dev_c2, -Inf)
      ),
    
    n = 5,
    
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    Parameter = case_when(
        
        distribution == "Floor" ~ "p (floor proportion)",
        
        distribution == "Degenerate" ~ "σ (SD)",
        
        distribution == "t-heavy" ~ "df (Student’s t)",
        
        distribution == "Pareto" ~ "α (tail index)"
      ),
    
    Parameter_level = case_when(
      
        distribution == "Floor" ~ as.character(prop_floor),
        
        distribution == "Degenerate" ~ as.character(sigma),
        
        distribution == "t-heavy" ~ as.character(df),
        
        distribution == "Pareto" ~ as.character(alpha)
      )
  ) |>
  transmute(
    
    Distribution = distribution,
    
    Sample_size = sample_size,
    
    Parameter,
    
    Parameter_level,
    
    `Type I (t)` = type1_t,
    
    `Type I (BM)` =
      if ("type1_bm" %in% names(results_all))
        type1_bm
    else
      NA,
    
    `Type I (C2)` =
      if ("type1_c2" %in% names(results_all))
        type1_c2
    else
      NA,
    
    `|Type I - α| (t)` = dev_t,
    
    `|Type I - α| (BM)` = dev_bm,
    
    `|Type I - α| (C2)` = dev_c2
  ) |>
  arrange( Distribution, desc(`|Type I - α| (t)`)) |>
  select(
    where( ~ !all(is.na(.x)))
  )

# WORST POWER

worst_power <- appendix_all |>
  filter(!H0_true) |>
  mutate(
    
    min_power_cal = pmin(
        
        Power_t_cal,
        
        if (
          "Power_bm_cal" %in% names(appendix_all)
        )
          Power_bm_cal
        else
          NA,
        
        if (
          "Power_c2_cal" %in% names(appendix_all)
        )
          Power_c2_cal
        else
          NA,
        
        na.rm = TRUE
      )
  ) |>
  group_by(Distribution) |>
  filter(
    is.finite(min_power_cal)
  ) |>
  slice_min(
    order_by = min_power_cal,
    n = 5,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(Distribution,min_power_cal) |>
  select(-min_power_cal)


# ==========Codebook==========================

codebook <- tribble(
  
  ~Column, ~Meaning,
  
  "Distribution",
  "Family of data-generating mechanism (Floor, Degenerate, t-heavy, Pareto).",
  
  "Sample_size",
  "Per-group sample size n used in the Monte Carlo experiment.",
  
  "H0_true",
  "TRUE if null hypothesis (no effect/contrast) is true for the row’s condition.",
  
  "Parameter",
  "Primary condition parameter for the family (e.g., p, σ, df, α).",
  
  "Parameter_level",
  "Value of the primary condition parameter as displayed.",
  
  "Effect_or_contrast",
  "Magnitude defining H1: shift_mu (Floor), d_target (Degenerate/t-heavy), r (Pareto).",
  
  "Type1_t",
  "Empirical Type I error (Welch t) estimated on H0 rows.",
  
  "Type1_bm",
  "Empirical Type I error (Brunner-Munzel) estimated on H0 rows.",
  
  "Type1_c2",
  "Empirical Type I error (C2 test) estimated on H0 rows.",
  
  "Power_t",
  "Empirical power (Welch t) estimated on H1 rows (uncalibrated).",
  
  "Power_bm",
  "Empirical power (Brunner-Munzel) estimated on H1 rows (uncalibrated).",
  
  "Power_c2",
  "Empirical power (C2 test) estimated on H1 rows (uncalibrated).",
  
  "Power_*_cal",
  "Power set to NA unless H0 size ≤ α+band for that (n, param).",
  
  "Cohens_d_abs",
  "Average absolute Cohen’s d from Welch t kernel across replicates.",
  
  "Calibrated_*",
  "Whether test size passed calibration (≤ α+band) at matching (n, param).",
  
  "Smaller_calib",
  "Which test has the smallest calibrated power for the row (t-test/BM/C2/NA).",
  
  "Larger_calib",
  "Which test has the largest calibrated power for the row (t-test/BM/C2/NA).",
  
  "prop_floor",
  "Raw family parameter: floor proportion (Floor only).",
  
  "sigma",
  "Raw family parameter: SD (Degenerate only).",
  
  "df",
  "Raw family parameter: degrees of freedom (t-heavy only).",
  
  "alpha",
  "Raw family parameter: Pareto tail index α (Pareto only).",
  
  "r",
  "Raw family parameter: Pareto scale ratio xm1/xm2 (Pareto only)."
)

# SAVE ALL EXCEL OUTPUTS

writexl::write_xlsx(
  
  list(
    
    "All_Results" = appendix_all,
    
    "Floor" = appendix_floor,
    
    "Degenerate" = appendix_degen,
    
    "t-heavy" = appendix_t,
    
    "Pareto" = appendix_pareto,
    
    "Worst_Size" = worst_size,
    
    "Worst_Power" = worst_power
  ),
  
  path = file.path(output_dir, "Supplementary_AllResults.xlsx")
)


writexl::write_xlsx(
  
  list(
    "Codebook" = codebook
  ),
  
  path = file.path(output_dir, "Supplementary_Codebook.xlsx")
)

# OPTIONAL: SAVE R OBJECTS

# This gives you a complete R-readable copy of the simulation results
# in the same dated output folder.

saveRDS(
  list(
    results_floor = results_floor,
    results_degen = results_degen,
    results_t = results_t,
    results_pareto = results_pareto,
    results_all = results_all,
    appendix_all = appendix_all,
    worst_size = worst_size,
    worst_power = worst_power,
    min_n_data = min_n_data
  ),
  
  file =
    file.path(
      output_dir,
      "Simulation_Results.rds"
    )
)


# ========================FINISHED===============================

cat("\n")
cat("============================================================\n")
cat("SIMULATION COMPLETE\n")
cat("============================================================\n")
cat("All output files saved to:\n")
cat(output_dir, "\n")
cat("============================================================\n")

