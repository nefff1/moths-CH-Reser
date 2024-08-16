# SETUP ------------------------------- ########################################
################################################################################.

# ... packages #################################################################
################################################################################.

library(tidyverse); theme_set(theme_classic())
library(rstan); options(mc.cores = 4)
library(brms)
library(bayestestR)
library(mgcv)
library(cowplot)
library(ggnewscale)
library(ggh4x)
library(openxlsx)

# ... set global parameters ####################################################
################################################################################.

min_mass <- min(d_mod_z$mass_tot[d_mod_z$mass_tot > 0])

log_plus_trans <- scales::trans_new(name = "log_plus",
                                    transform = \(x) log(x + min_mass, 10),
                                    inverse = \(x) out <- 10 ^ (x) - min_mass)

v_textsize <- c(axis.title = 8, axis.text = 7, legend.title = 8, additional.text = 6) 

n_iter <- 2000 # number of MCMC iterations

# ... define functions #########################################################
################################################################################.

f_analysis_A_height <- function(formula, data_z, scalings, family, iter, seed, 
                                hours_sel = NULL, run_loo = FALSE) {
  
  data_z <- droplevels(data_z)
  
  set.seed(seed)
  
  out <- list()
  
  response <- as.character(formula[2])
  
  
  l_data <- make_standata(update(formula, ". ~ A * height + ."), 
                          data = data_z)
  
  if (length(hours_sel) == 0){
    hours_sel <- NULL
  }
  
  if (!is.null(hours_sel)){
    l_data_sub <- make_standata(update(formula, ". ~ s(active_hours)"),
                                data = data_z[hours_sel, ])
    
    l_data$N_SEL <- length(hours_sel)
    l_data$SEL <-  hours_sel
    l_data$Xs_add <- l_data_sub$Xs
    l_data$nb_add <- 1
    l_data$knots_add <- l_data$knots_1
    dimnames(l_data$knots_add)[[1]] <- "Zs_add_1"
    l_data$Zs_add_1 <- l_data_sub$Zs_1_1
  }
  
  
  if (family == "zero_inflated_negbinomial"){
    if (!is.null(hours_sel)){
      mod <- stan(file = "Stan_splines/Stan_nb_spline_s1p1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    } else {
      mod <- stan(file = "Stan_splines/Stan_nb_spline_s1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    }
  } else if (family == "hurdle_gamma"){
    if (!is.null(hours_sel)){
      mod <- stan(file = "Stan_splines/Stan_hg_spline_s1p1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    } else {
      mod <- stan(file = "Stan_splines/Stan_hg_spline_s1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    }
  }
  
  # extract parameter data
  par_names <- names(mod)
  pos <- regexpr("\\[", par_names) 
  pos <- ifelse(pos == -1, 200, pos - 1)
  pars <- unique(substr(par_names, 1, pos))
  pars <- pars[pars %in% c("b", "Intercept", "bs", "bs_add", "zs_1_1",
                           "sds_1_1", "zs_2_1", "sds_2_1", "zs_add_1",
                           "sds_add_1", "shape", "zi", "sd_1", "sd_2",
                           "sd_3", "sd_4", "s_1_1", "s_2_1", "s_add_1")]
  fit <- rstan::extract(mod, pars = pars)
  out$fit <- fit
  if (run_loo) out$loo_lm <- loo(fit)
  
  terms <- brmsterms(update(formula, ". ~ A * height + ."))
  
  # predictions for fixed effects
  out$l_pred_fe <- f_pred_fixed(fit = fit, data = data_z, 
                                terms = terms, scalings = scalings)
  
  
  # predictions for smooths
  out$l_pred_sm <- f_pred_smooths(fit = fit, data = data_z, 
                                  l_data = l_data,
                                  terms = terms, scalings = scalings,
                                  hours_sel = hours_sel)
  
  
  out
}

f_pred_fixed <- function(fit, data, terms, scalings) {
  
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  
  # effects without interactions -----------------------------------------------.
  
  l_out <- list()
  for (var_i in all.vars(terms$dpars$mu$fe)){
    if (is.numeric(deframe(data[, var_i]))) {
      
      d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                      length.out = 150)) %>% 
        mutate(expl_c = expl - mean(data[, var_i, drop = T]),
               expl_or = expl * scalings$sd[scalings$var == var_i] +
                 scalings$mean[scalings$var == var_i])
      
      
      m_pred <- matrix(rep(fit$Intercept, each = 150), nrow = 150) + 
        matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                  !grepl(":", colnames(mm_full)))]
      
      d_pred <- cbind(d_pred, m_pred)
      
      d_pred <- d_pred %>% 
        select(-c(expl, expl_c)) %>% 
        pivot_longer(cols = -expl_or,
                     names_to = "iter", values_to = "pred") %>%
        rename(!! sym(var_i) := expl_or) %>% 
        mutate(pred_exp = exp(pred)) %>%
        group_by(!! sym(var_i)) %>% 
        summarise(log_lower = ci(pred, .95)$CI_low,
                  log_upper = ci(pred, .95)$CI_high,
                  log_estimate = mean(pred),
                  lower = ci(pred_exp, .95)$CI_low,
                  upper = ci(pred_exp, .95)$CI_high,
                  estimate = mean(pred_exp),
                  .groups = "drop")
      
    } else {
      
      d_pred <- data.frame(expl = unique(data[, var_i]))
      
      comps <- strsplit(paste(deparse(terms$dpars$mu$fe), collapse = ""), split = " \\+ ")[[1]]
      
      comps <- comps[grepl(var_i, comps)]
      
      mm_data <- model.matrix(as.formula(paste("~", comps)), data = data)[, -1, drop = F]
      mm_pred <- model.matrix(as.formula(paste("~", comps)), data = d_pred)[, -1, drop = F]
      
      for (i in seq_len(ncol(mm_pred))){
        mm_pred[, i] <- mm_pred[, i] - mean(mm_data[, i])
      }
      
      m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred)) + 
        as.matrix(mm_pred) %*% t(as.matrix(fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                           !grepl(":", colnames(mm_full)))]))
      
      d_pred <- cbind(d_pred, m_pred)
      
      d_pred <- d_pred %>% 
        pivot_longer(cols = - !! enquo(var_i),
                     names_to = "iter", values_to = "pred") %>%
        mutate(pred_exp = exp(pred)) %>%
        group_by(!! sym(var_i)) %>% 
        summarise(log_lower = ci(pred, .95)$CI_low,
                  log_upper = ci(pred, .95)$CI_high,
                  log_estimate = mean(pred),
                  lower = ci(pred_exp, .95)$CI_low,
                  upper = ci(pred_exp, .95)$CI_high,
                  estimate = mean(pred_exp),
                  .groups = "drop")
    }
    
    l_out[[var_i]] <- d_pred
  }
  
  # effects with interactions --------------------------------------------------.
  for (int_i in colnames(mm_full)[grepl(":", colnames(mm_full))]) {
    vars_i <- strsplit(int_i, ":")[[1]]
    
    d_pred <- data.frame()
    for (var_i in vars_i){
      
      if (nrow(d_pred) == 0){
        d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                        length.out = 150)) %>%
          mutate(expl_c = expl - mean(data[, var_i, drop = T]),
                 expl_or = expl * scalings$sd[scalings$var == var_i] +
                   scalings$mean[scalings$var == var_i]) %>% 
          rename(!! sym(var_i) := expl,
                 !! sym(paste0(var_i, "_C")) := expl_c,
                 !! sym(paste0(var_i, "_OR")) := expl_or)
      } else {
        d_pred <- expand_grid(d_pred, 
                              data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                                    length.out = 150)) %>% 
                                mutate(expl_c = expl - mean(data[, var_i, drop = T]),
                                       expl_or = expl * scalings$sd[scalings$var == var_i] +
                                         scalings$mean[scalings$var == var_i]) %>% 
                                rename(!! sym(var_i) := expl,
                                       !! sym(paste0(var_i, "_C")) := expl_c,
                                       !! sym(paste0(var_i, "_OR")) := expl_or))
      }
      
    }
    
    
    m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred))
    
    # add main effects
    for (var_i in vars_i){
      m_pred <- m_pred +
        matrix(d_pred[, paste0(var_i, "_C"), drop = T]) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                                            !grepl(":", colnames(mm_full)))]
    }
    
    # add interactive effects
    m_pred <- m_pred + 
      matrix(d_pred %>% 
               select(ends_with("_C")) %>% 
               apply(., 1, prod)) %*% fit$b[, which(grepl(int_i, colnames(mm_full)))]
    
    
    d_pred <- cbind(d_pred, m_pred)
    
    d_pred <-
      d_pred %>% 
      select(- c(!!! syms(vars_i)),
             - c(!!! syms(paste0(vars_i, "_C")))) %>% 
      pivot_longer(cols = -ends_with("_OR"),
                   names_to = "iter", values_to = "pred") %>%
      rename_all(~ gsub("_OR$", "", .)) %>% 
      mutate(pred_exp = exp(pred)) %>%
      group_by(!!! syms(vars_i)) %>% 
      summarise(log_lower = ci(pred, .95)$CI_low,
                log_upper = ci(pred, .95)$CI_high,
                log_estimate = mean(pred),
                lower = ci(pred_exp, .95)$CI_low,
                upper = ci(pred_exp, .95)$CI_high,
                estimate = mean(pred_exp),
                .groups = "drop")
    
    l_out[[int_i]] <- d_pred
  }
  
  l_out
}

f_pred_smooths <-  function(fit, data, l_data, terms, scalings, hours_sel){
  l_out <- list()
  
  if (!is.null(hours_sel)){
    vars <- c(all.vars(terms$dpars$mu$sm), "active_hours")
  } else {
    vars <- all.vars(terms$dpars$mu$sm)
  }
  
  
  for (var_i in vars){ 
    
    if (var_i == "active_hours"){
      d_target_unscaled <- data[hours_sel, var_i] %>% 
        rename(var_target = !! sym(var_i))
      
      d_target_scaled <- l_data$Xs_add  %>% 
        as.data.frame()
      
      v_bs <- fit$bs_add[, 1]
      v_s <- fit$s_add_1
      
    } else {
      d_target_unscaled <- data[, var_i] %>% 
        rename(var_target = !! sym(var_i))
      
      d_target_scaled <- l_data$Xs[, which(colnames(l_data$Xs) == paste0("s", var_i, "_1"))]  %>% 
        as.data.frame()
      
      count_i <- which(all.vars(terms$dpars$mu$sm) == var_i)
      
      v_bs <- fit$bs[, count_i]
      v_s <- fit[[paste0("s_", count_i, "_1")]]
      
    }
    
    
    names(d_target_scaled) <- "var_target"
    
    sm_unscaled <- smoothCon(s(var_target),
                             data  = d_target_unscaled,
                             knots = NULL,
                             absorb.cons = TRUE, 
                             modCon = 3,
                             diagonal.penalty = TRUE)
    
    sm_scaled <- smoothCon(s(var_target),
                           data = d_target_scaled,
                           knots = NULL,
                           absorb.cons = TRUE, 
                           modCon = 3,
                           diagonal.penalty = TRUE)
    
    d_newdata_unscaled <- data.frame(var_target = seq(min(d_target_unscaled$var_target), 
                                                      max(d_target_unscaled$var_target), 
                                                      length.out = 365))
    pred_sm <- PredictMat(sm_unscaled[[1]], d_newdata_unscaled)
    
    d_newdata_scaled <- data.frame(var_target = pred_sm[, 9])
    
    pred_sm <- pred_sm[, -9][, c(1,8:2)]
    
    pred <- matrix(rep(fit$Intercept, each = nrow(d_newdata_scaled)),
                   nrow = nrow(d_newdata_scaled)) + 
      as.matrix(d_newdata_scaled$var_target) %*%  v_bs + 
      pred_sm %*% t(v_s)
    
    d_pred <- data.frame(d_newdata_unscaled, pred) %>%
      pivot_longer(-var_target, names_to = "iter", values_to = "pred")
    
    d_pred <- d_pred %>%
      mutate(var_target = var_target * scalings$sd[scalings$var == var_i] +
               scalings$mean[scalings$var == var_i]) %>%
      rename(!! sym(var_i) := var_target) %>% 
      mutate(pred_exp = exp(pred)) %>%
      group_by(!! sym(var_i)) %>% 
      summarise(log_lower = ci(pred, .95)$CI_low,
                log_upper = ci(pred, .95)$CI_high,
                log_estimate = mean(pred),
                lower = ci(pred_exp, .95)$CI_low,
                upper = ci(pred_exp, .95)$CI_high,
                estimate = mean(pred_exp),
                .groups = "drop")
    
    l_out[[var_i]] <- d_pred
  }
  l_out
}

f_A_height_change_numbers <- function(fit, data, terms) {
  
  int_i <- "A:height"
  
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  
  # also considers main effects, not just interactive effects!
  
  vars_i <- strsplit(int_i, ":")[[1]]
  
  d_pred <- data.frame()
  for (var_i in vars_i){
    
    if (nrow(d_pred) == 0){
      d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                      length.out = 150)) %>%
        mutate(expl_c = expl - mean(data[, var_i, drop = T])) %>%
        rename(!! sym(var_i) := expl,
               !! sym(paste0(var_i, "_C")) := expl_c)
    } else {
      d_pred <- expand_grid(d_pred, 
                            data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                                  length.out = 150)) %>% 
                              mutate(expl_c = expl - mean(data[, var_i, drop = T])) %>%
                              rename(!! sym(var_i) := expl,
                                     !! sym(paste0(var_i, "_C")) := expl_c))
    }
    
  }
  
  
  m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred))
  
  # add main effects
  for (var_i in vars_i){
    m_pred <- m_pred +
      matrix(d_pred[, paste0(var_i, "_C"), drop = T]) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                                          !grepl(":", colnames(mm_full)))]
  }
  
  # add interactive effects
  m_pred <- m_pred + 
    matrix(d_pred %>% 
             select(ends_with("_C")) %>% 
             apply(., 1, prod)) %*% fit$b[, which(grepl(int_i, colnames(mm_full)))]
  
  
  d_pred <- cbind(d_pred, m_pred)
  
  
  d_pred %>% 
    filter(A %in% range(A),
           height %in% range(height)) %>% 
    mutate(height = ifelse(height == min(height), "lowest", "highest")) %>% 
    select(- c(!!! syms(paste0(vars_i, "_C")))) %>% 
    arrange(height, A) %>% 
    select(-A) %>% 
    pivot_longer(-height, names_to = "iter", values_to = "pred") %>% 
    group_by(height, iter) %>% 
    summarise(diff_mean = mean(pred[2] - pred[1]),
              .groups = "drop") %>% 
    group_by(height) %>% 
    summarise(diff_log_estimate = mean(diff_mean),
              diff_log_lower = ci(diff_mean)$CI_low,
              diff_log_upper = ci(diff_mean)$CI_high,
              factor_estimate = mean(exp(diff_mean)),
              factor_lower = ci(exp(diff_mean))$CI_low,
              factor_upper = ci(exp(diff_mean))$CI_high,
              .groups = "drop")
  
}

f_apply_Rhat <- function(fit){
  # Rhat function basend on rstan package
  f_Rhat <- \(sims) {
    
    f_tomatrix <- \(obj_draws){
      matrix(as.numeric(obj_draws), ncol = 8, byrow = F) # 4 chains split in two
    }
    
    f_rhat_rfun <- \(sims) {
      chains <- ncol(sims)
      n_samples <- nrow(sims)
      chain_mean <- numeric(chains)
      chain_var <- numeric(chains)
      for (i in seq_len(chains)) {
        chain_mean[i] <- mean(sims[, i])
        chain_var[i] <- var(sims[, i])
      }
      var_between <- n_samples * var(chain_mean)
      var_within <- mean(chain_var)
      sqrt((var_between/var_within + n_samples - 1)/n_samples)
    }
    
    f_z_scale <- \(x){
      S <- length(x)
      r <- rank(x, ties.method = 'average')
      z <- qnorm((r - 1/2)/S)
      z[is.na(x)] <- NA
      if (!is.null(dim(x))) {
        z <- array(z, dim = dim(x), dimnames = dimnames(x))
      }
      z
    }
    
    bulk_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims)))
    sims_folded <- abs(sims - median(sims))
    tail_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims_folded)))
    max(bulk_rhat, tail_rhat)
  }
  
  d_Rhat <- data.frame()
  for (var_i in names(fit)[grepl("Intercept|^b|^bs|^s_", names(fit))]){
    if (length(dim(fit[[var_i]])) > 1){ 
      for (i in seq_len(ncol(fit[[var_i]]))){
        d_Rhat <- d_Rhat |> 
          bind_rows(data.frame(var = paste(var_i, i, sep = "_"),
                               rhat = f_Rhat(fit[[var_i]][, i])))
      }
    } else {
      d_Rhat <- d_Rhat |> 
        bind_rows(data.frame(var = var_i,
                             rhat = f_Rhat(fit[[var_i]])))
    }
  }
  
  d_Rhat
}

f_A_height_plot <- function(pred,
                            data_raw = NULL, mean_yday = NULL, response = NULL,
                            quantiles = c(0, .25, .5, .75, 1)){
  
  # simplified version to take quantiles
  sel <- round(length(unique(pred$height)) * quantiles)
  sel[sel == 0] <- 1
  sel_height <- unique(pred$height)[sel]
  
  pred <- pred %>% 
    filter(height %in% sel_height) %>% 
    rowwise() %>%
    mutate(height_cat = paste0("Q-", quantiles[which(sel_height == height)] * 100, "%")) %>% 
    ungroup() |> 
    mutate(height_cat = factor(height_cat, levels = unique(height_cat),
                               labels = paste0(unique(height_cat), "\n(   ", round(sel_height), "m)")))
  
  
  if (!is.null(mean_yday)){
    pred <- pred %>%
      mutate(Date = as.Date(format(date_decimal(A)), "%Y-%m-%d") + 
               mean_yday) 
  } else {
    pred <- pred %>%
      mutate(Date = A) 
  }
  
  p <- pred %>%
    ggplot(aes(x = Date, col = height_cat))
  
  if (!is.null(data_raw)){
    p <- p +
      geom_point(data = data_raw, aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4")
  }
  
  p <- p +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "grey20", alpha = .3, lty = 2) +
    geom_line(aes(y = exp(log_estimate)), size = 1) + # more stable mean estimates
    scale_colour_manual(values = RColorBrewer::brewer.pal(length(quantiles),
                                                          "RdYlBu"))
  
  
  if (!is.null(response)){
    if (response == "mass_tot"){
      p <-  p +
        scale_y_continuous(trans = log_plus_trans)
    } else {
      p <-  p +
        scale_y_continuous(trans = "log1p")
    }
  }
  
  p +
    labs(colour = "Elevation", x = "Year")
}

f_A_height_incr_decr <- function(fit, data, scalings, formula = NULL,
                                    c_factors = c(1.1, 1.5, 2, 4)){
  
  
  f_P_sum <- \(x, resolution = 1000){
    
    seq_neg <- seq(-0.9, -(1 - 1 / c_factors[1]), length.out = resolution)
    seq_pos <-  (1 / (1+seq_neg)) - 1
    
    data.frame(threshold = seq_neg,
               test = "decrease") |> 
      rowwise() |> 
      mutate(P_sum = mean(x <= threshold)) |> 
      ungroup() |> 
      bind_rows(data.frame(threshold = seq_pos,
                           test = "increase") |> 
                  rowwise() |> 
                  mutate(P_sum = mean(x >= threshold)) |> 
                  ungroup()) 
    
  }
  
  if (is.null(formula)){
    formula <- response ~ A * height + s(yday) + P_2day + T_2day +
      C(traptype, "contr.sum") +
      C(bulbtype, "contr.sum") +
      n_trap +
      C(sample_previous, "contr.sum") +
      (1 | gr) +
      (1 | trap_ID) +
      (1 | night_ID) +
      (1 | trap_ID_A)
  }
  
  terms <- brmsterms(formula)
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  vars_i <- c("A", "height")
  
  d_pred <- data.frame()
  for (var_i in vars_i){
    
    if (nrow(d_pred) == 0){
      d_pred <-
        data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                              length.out = 3)) %>%
        mutate(expl_c = expl - mean(data[, var_i, drop = T]),
               expl_or = expl * scalings$sd[scalings$var == var_i] +
                 scalings$mean[scalings$var == var_i]) %>% 
        rename(!! sym(var_i) := expl,
               !! sym(paste0(var_i, "_C")) := expl_c,
               !! sym(paste0(var_i, "_OR")) := expl_or)
    } else {
      d_pred <- expand_grid(d_pred, 
                            data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                                  length.out = 3)) %>% 
                              mutate(expl_c = expl - mean(data[, var_i, drop = T]),
                                     expl_or = expl * scalings$sd[scalings$var == var_i] +
                                       scalings$mean[scalings$var == var_i]) %>% 
                              rename(!! sym(var_i) := expl,
                                     !! sym(paste0(var_i, "_C")) := expl_c,
                                     !! sym(paste0(var_i, "_OR")) := expl_or))
    }
    
  }
  
  
  m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred))
  
  if (ncol(m_pred) < 4000) warning("Attention! Not 4 chains!")
  
  # add main effects
  for (var_i in vars_i){
    m_pred <- m_pred +
      matrix(d_pred[, paste0(var_i, "_C"), drop = T]) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                                          !grepl(":", colnames(mm_full)))]
  }
  
  # add interactive effects
  m_pred <- m_pred + 
    matrix(d_pred %>% 
             select(ends_with("_C")) %>% 
             apply(., 1, prod)) %*% fit$b[, colnames(mm_full) == "A:height"]
  
  
  d_pred <- cbind(d_pred, m_pred)
  
  d_pred |> 
    filter(A_OR %in% c(min(A_OR), max(A_OR))) |> 
    mutate_at(vars(!any_of(c("A", "A_C", "A_OR", "height", "height_C", "height_OR"))),
              ~ exp(.)) |>
    group_by(height_OR) |>
    summarise_at(.vars = vars(!any_of(c("A", "A_C", "A_OR", "height", "height_C", "height_OR"))),
                 ~ (.[A_OR == max(A_OR)] - .[A_OR == min(A_OR)]) / .[A_OR == min(A_OR)]) |> 
    pivot_longer(-height_OR, names_to = "iter", values_to = "prop_change") |> 
    mutate(height = ifelse(height_OR == max(height_OR),
                           "max", ifelse(height_OR == min(height_OR),
                                         "min", "median")),
           height = factor(height, levels = c("min", "median", "max"),
                           labels = c("Min. elevation", "Med. elevation", "Max. elevation"))) |> 
    group_by(height) |> 
    group_map(~ data.frame(f_P_sum(.$prop_change)) |> 
                mutate(height = unique(.$height)), .keep = T) |> 
    bind_rows() |> 
    mutate(threshold_cut = cut(threshold, breaks = c(-Inf, rev(1/c_factors - 1), c_factors - 1, Inf)),
           P_sum = ifelse(threshold_cut == "(-0.0909,0.1]", NA, P_sum)) |> 
    # to have a line closing the curve on the right and left:
    group_by(height) |> 
    group_modify(~add_row(., threshold = -(1 - 1 / 1.1), test = "decrease", P_sum = 0)) |> 
    group_modify(~add_row(., threshold = 0.1, test = "increase", P_sum = 0)) |>
    ungroup() |> 
    group_by(height, test) |> 
    mutate(xwidth = (c(diff(threshold)[1], diff(threshold)) +
                       c(diff(threshold), diff(threshold)[length(threshold)])) / 2) |> 
    ungroup()
}

f_A_height_plot_comb <- function(l_pred,
                                 response = NULL,
                                 quantiles = c(0, .25, .5, .75, 1),
                                 name = NULL){
  
  # simplified version to take quantiles
  sel <- round(length(unique(l_pred[[1]]$height)) * quantiles)
  sel[sel == 0] <- 1
  sel_height <- unique(l_pred[[1]]$height)[sel]
  
  
  d_plot <- lapply(l_pred,
                   function(x) x %>% 
                     filter(height %in% sel_height) %>% 
                     rowwise() %>%
                     mutate(height_cat = paste0("Q", quantiles[which(sel_height == height)])) %>% 
                     ungroup() |> 
                     mutate(height_cat = factor(height_cat, levels = unique(height_cat),
                                                labels = paste0(unique(height_cat), 
                                                                "\n(   ", round(sel_height), "m)")))) %>% 
    bind_rows(.id = "Trait") %>% 
    mutate(Trait = factor(Trait, levels = names(l_pred)),
           name = name)
  
  
  p <-
    d_plot %>%
    ggplot(aes(x = A, col = height_cat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "grey20", alpha = .3, lty = 2) +
    geom_line(aes(y = exp(log_estimate)), size = 1) +
    scale_colour_manual(values = RColorBrewer::brewer.pal(length(quantiles),
                                                          "RdYlBu"))
  
  if (!is.null(name)) {
    p <- p + facet_nested_wrap(~ name + Trait, scales = "free_y", nrow = 1)
  } else {
    p <- p + facet_wrap(~ Trait, scales = "free_y", nrow = 1)
  }
  
  if (!is.null(response)){
    if (response == "mass_tot"){
      p <-  p +
        scale_y_continuous(trans = log_plus_trans)
    } else {
      p <-  p +
        scale_y_continuous(trans = "log1p")
    }
  }
  
  p +
    labs(colour = "Elevation", x = "Year")
}

# RUN MODELS -------------------- ##############################################
################################################################################.

# ... overall models ###########################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A <- f_analysis_A_height(formula = abu_tot ~
                                 s(yday) + P_2day + T_2day +
                                 C(traptype, "contr.sum") +
                                 C(bulbtype, "contr.sum") +
                                 n_trap +
                                 C(sample_previous, "contr.sum") +
                                 (1 | gr) +
                                 (1 | trap_ID) +
                                 (1 | night_ID) +
                                 (1 | trap_ID_A),
                               data_z = d_mod_z,
                               scalings = filter(d_scalings, data == "full"),
                               family = "zero_inflated_negbinomial",
                               hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                               iter = n_iter, seed = 126)

# richness ---------------------------------------------------------------------.
l_ric_A <- f_analysis_A_height(formula = sric ~
                                 s(yday) + P_2day + T_2day +
                                 C(traptype, "contr.sum") +
                                 C(bulbtype, "contr.sum") +
                                 n_trap +
                                 C(sample_previous, "contr.sum") +
                                 (1 | gr) + (1 | trap_ID) +
                                 (1 | night_ID) +
                                 (1 | trap_ID_A),
                               data_z = d_mod_z,
                               scalings = filter(d_scalings, data == "full"),
                               family = "zero_inflated_negbinomial",
                               hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                               iter = n_iter, seed = 8989)

# biomass ----------------------------------------------------------------------.
l_mass_A <- f_analysis_A_height(formula = mass_tot ~
                                  s(yday) + P_2day + T_2day +
                                  C(traptype, "contr.sum") +
                                  C(bulbtype, "contr.sum") +
                                  n_trap +
                                  C(sample_previous, "contr.sum") +
                                  (1 | gr) +
                                  (1 | trap_ID) +
                                  (1 | night_ID) +
                                  (1 | trap_ID_A),
                                data_z = d_mod_z,
                                scalings = filter(d_scalings, data == "full"),
                                family = "hurdle_gamma",
                                hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                                iter = n_iter, seed = 811)

# ... per body size group ######################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_small <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                         !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "small"]),
                                     iter = n_iter, seed = 323)
l_abu_A_medium <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                          !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "medium"]),
                                      iter = n_iter, seed = 963)
l_abu_A_large <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                         !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "large"]),
                                     iter = n_iter, seed = 97)

# richness ---------------------------------------------------------------------.
l_ric_A_small <- f_analysis_A_height(formula = sric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                         !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "small"]),
                                     iter = n_iter, seed = 8460)
l_ric_A_medium <- f_analysis_A_height(formula = sric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                          !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "medium"]),
                                      iter = n_iter, seed = 510)
l_ric_A_large <- f_analysis_A_height(formula = sric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                         !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "large"]),
                                     iter = n_iter, seed = 879)

# biomass ----------------------------------------------------------------------.
l_mass_A_small <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                          !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "small"]),
                                      iter = n_iter, seed = 2512)
l_mass_A_medium <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | gr) +
                                         (1 | trap_ID) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                           !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "medium"]),
                                       iter = n_iter, seed = 1717)
l_mass_A_large <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                          !d_mod_mass_z$estimate[d_mod_mass_z$mass_cat == "large"]),
                                      iter = n_iter, seed = 112)

# ... per temperature niche group ##############################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_cold <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                        !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                    iter = n_iter, seed = 1688)
l_abu_A_intermediate <- f_analysis_A_height(formula = abu_tot ~ height +
                                              s(yday) + P_2day + T_2day +
                                              C(traptype, "contr.sum") +
                                              C(bulbtype, "contr.sum") +
                                              n_trap +
                                              C(sample_previous, "contr.sum") +
                                              (1 | gr) +
                                              (1 | trap_ID) +
                                              (1 | night_ID) +
                                              (1 | trap_ID_A),
                                            data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                            scalings = filter(d_scalings, data == "full"),
                                            family = "zero_inflated_negbinomial",
                                            hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                            iter = n_iter, seed = 190)
l_abu_A_warm <- f_analysis_A_height(formula = abu_tot ~ height +
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                        !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                    iter = n_iter, seed = 1881)

# richness ---------------------------------------------------------------------.
l_ric_A_cold <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                        !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                    iter = n_iter, seed = 871)
l_ric_A_intermediate <- f_analysis_A_height(formula = sric ~
                                              s(yday) + P_2day + T_2day +
                                              C(traptype, "contr.sum") +
                                              C(bulbtype, "contr.sum") +
                                              n_trap +
                                              C(sample_previous, "contr.sum") +
                                              (1 | gr) +
                                              (1 | trap_ID) +
                                              (1 | night_ID) +
                                              (1 | trap_ID_A),
                                            data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                            scalings = filter(d_scalings, data == "full"),
                                            family = "zero_inflated_negbinomial",
                                            hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                            iter = n_iter, seed = 111)
l_ric_A_warm <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                        !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                    iter = n_iter, seed = 901)

# richness ---------------------------------------------------------------------.
l_mass_A_cold <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                         !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                     iter = n_iter, seed = 964)
l_mass_A_intermediate <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                               C(traptype, "contr.sum") +
                                               C(bulbtype, "contr.sum") +
                                               n_trap +
                                               C(sample_previous, "contr.sum") +
                                               (1 | gr) +
                                               (1 | trap_ID) +
                                               (1 | night_ID) +
                                               (1 | trap_ID_A),
                                             data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                             scalings = filter(d_scalings, data == "full"),
                                             family = "hurdle_gamma",
                                             hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                 !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                             iter = n_iter, seed = 111)
l_mass_A_warm <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | gr) +
                                       (1 | trap_ID) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                         !d_mod_Tavg_z$estimate[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                     iter = n_iter, seed = 255)

# ... per specialisation group #################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_spec_m <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Monophagous"]),
                                      iter = n_iter, seed = 346)
# l_abu_A_spec_o <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Oligophagous"]),
                                      iter = n_iter, seed = 665)
# l_abu_A_spec_p <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Polyphagous"]),
                                      iter = n_iter, seed = 863)

# richness ---------------------------------------------------------------------.
l_ric_A_spec_m <- f_analysis_A_height(formula = sric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Monophagous"]),
                                      iter = n_iter, seed = 8132)
l_ric_A_spec_o <- f_analysis_A_height(formula = sric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Oligophagous"]),
                                      iter = n_iter, seed = 551)
l_ric_A_spec_p <- f_analysis_A_height(formula = sric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | gr) +
                                        (1 | trap_ID) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                          !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Polyphagous"]),
                                      iter = n_iter, seed = 918)

# biomass ----------------------------------------------------------------------.
l_mass_A_spec_m <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | gr) +
                                         (1 | trap_ID) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                           !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Monophagous"]),
                                       iter = n_iter, seed = 810)
l_mass_A_spec_o <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | gr) +
                                         (1 | trap_ID) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                           !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Oligophagous"]),
                                       iter = n_iter, seed = 743)
l_mass_A_spec_p <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | gr) +
                                         (1 | trap_ID) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                           !d_mod_spec_z$estimate[d_mod_spec_z$Spec == "Polyphagous"]),
                                       iter = n_iter, seed = 822)

# ... per overwintering stage ##################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_egg <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                                    iter = n_iter, seed = 121)
l_abu_A_larva <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                                    iter = n_iter, seed = 442)
l_abu_A_pupa <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                                    iter = n_iter, seed = 9877)
l_abu_A_adult <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                                    iter = n_iter, seed = 2111)

# richness ---------------------------------------------------------------------.
l_ric_A_egg <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                                    iter = n_iter, seed = 12)
l_ric_A_larva <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                                    iter = n_iter, seed = 8797)
l_ric_A_pupa <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                                    iter = n_iter, seed = 194)
l_ric_A_adult <- f_analysis_A_height(formula = sric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                                    iter = n_iter, seed = 654)

# biomass ----------------------------------------------------------------------.
l_mass_A_egg <- f_analysis_A_height(formula = mass_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                                    iter = n_iter, seed = 537)
l_mass_A_larva <- f_analysis_A_height(formula = mass_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                                    iter = n_iter, seed = 634)
l_mass_A_pupa <- f_analysis_A_height(formula = mass_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                                    iter = n_iter, seed = 19)
l_mass_A_adult <- f_analysis_A_height(formula = mass_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | gr) +
                                      (1 | trap_ID) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                        !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                                    iter = n_iter, seed = 944)

# ... sensitivity analyses #####################################################
################################################################################.

# abundance --------------------------------------------------------------------.

# full observation hours data:
l_abu_A_ne <- f_analysis_A_height(formula = abu_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | gr) +
                                    (1 | trap_ID) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_ne_z,
                                  scalings = filter(d_scalings, data == "no_estimates"),
                                  family = "zero_inflated_negbinomial",
                                  hours_sel = which(d_mod_ne_z$traptype == "p"),
                                  iter = n_iter, seed = 54987)

# fixed traps only:
l_abu_A_lf <- f_analysis_A_height(formula = abu_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    C(sample_previous, "contr.sum") +
                                    (1 | gr) +
                                    (1 | trap_ID) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_lf_z,
                                  scalings = filter(d_scalings, data == "LF"),
                                  family = "zero_inflated_negbinomial",
                                  iter = n_iter, seed = 1213)


# manual traps only:
l_abu_A_p <- f_analysis_A_height(formula = abu_tot ~
                                   s(yday) + P_2day + T_2day +
                                   C(bulbtype, "contr.sum") +
                                   n_trap +
                                   C(sample_previous, "contr.sum") +
                                   (1 | gr) +
                                   (1 | trap_ID) +
                                   (1 | night_ID) +
                                   (1 | trap_ID_A),
                                 data_z = d_mod_p_z,
                                 scalings = filter(d_scalings, data == "p"),
                                 family = "zero_inflated_negbinomial",
                                 hours_sel = which(!d_mod_p_z$estimate),
                                 iter = n_iter, seed = 6894)

# richness ---------------------------------------------------------------------.

# full observation hours data:
l_ric_A_ne <- f_analysis_A_height(formula = sric ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | gr) +
                                    (1 | trap_ID) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_ne_z,
                                  scalings = filter(d_scalings, data == "no_estimates"),
                                  family = "zero_inflated_negbinomial",
                                  hours_sel = which(d_mod_ne_z$traptype == "p"),
                                  iter = n_iter, seed = 1211)
# fixed traps only:
l_ric_A_lf <- f_analysis_A_height(formula = sric ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    C(sample_previous, "contr.sum") +
                                    (1 | gr) + (1 | trap_ID) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_lf_z,
                                  scalings = filter(d_scalings, data == "LF"),
                                  family = "zero_inflated_negbinomial",
                                  iter = n_iter, seed = 45461)
# manual traps only:
l_ric_A_p <- f_analysis_A_height(formula = sric ~
                                   s(yday) + P_2day + T_2day +
                                   C(bulbtype, "contr.sum") +
                                   n_trap +
                                   C(sample_previous, "contr.sum") +
                                   (1 | gr) + (1 | trap_ID) +
                                   (1 | night_ID) +
                                   (1 | trap_ID_A),
                                 data_z = d_mod_p_z,
                                 scalings = filter(d_scalings, data == "p"),
                                 family = "zero_inflated_negbinomial",
                                 hours_sel = which(!d_mod_p_z$estimate),
                                 iter = n_iter, seed = 6894)

# biomass ----------------------------------------------------------------------.

# full observation hours data:
l_mass_A_ne <- f_analysis_A_height(formula = mass_tot ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     n_trap +
                                     C(sample_previous, "contr.sum") +
                                     (1 | gr) +
                                     (1 | trap_ID) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_ne_z,
                                   scalings = filter(d_scalings, data == "no_estimates"),
                                   family = "hurdle_gamma",
                                   hours_sel = which(d_mod_ne_z$traptype == "p"),
                                   iter = n_iter, seed = 9101)

# fixed traps only:
l_mass_A_lf <- f_analysis_A_height(formula = mass_tot ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     C(sample_previous, "contr.sum") +
                                     (1 | gr) +
                                     (1 | trap_ID) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_lf_z,
                                   scalings = filter(d_scalings, data == "LF"),
                                   family = "hurdle_gamma",
                                   iter = n_iter, seed = 565)

# manual traps only:
l_mass_A_p <- f_analysis_A_height(formula = mass_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | gr) +
                                    (1 | trap_ID) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_p_z,
                                  scalings = filter(d_scalings, data == "p"),
                                  family = "hurdle_gamma",
                                  hours_sel = which(!d_mod_p_z$estimate),
                                  iter = n_iter, seed = 6894)


# CREATE OUTPUTS -------------- ################################################
################################################################################.

# ... Text / numbers ###########################################################
################################################################################.

f_A_height_change_numbers(l_abu_A$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_ric_A$fit,
                          d_mod_z,
                          brmsterms(update(sric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_mass_A$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_abu_A_small$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_small$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_large$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_abu_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(sric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(sric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_spec_o$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(sric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | gr) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

# ... Rhat #####################################################################
################################################################################.

d_Rhat_abu <- f_apply_Rhat(l_abu_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_abu_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_ric <- f_apply_Rhat(l_ric_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_ric_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_mass <- f_apply_Rhat(l_mass_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_mass_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_abu |> 
  mutate(response = "Abundance") |> 
  bind_rows(d_Rhat_ric |> 
              mutate(response = "Richness"))|> 
  bind_rows(d_Rhat_mass |> 
              mutate(response = "Biomass")) |> 
  group_by(response, trait, data) |> 
  summarise(prop_thresh = mean(rhat < 1.1)) |> 
  pivot_wider(values_from = prop_thresh, names_from = response)


# ... Figure 1 #################################################################
################################################################################.

plot_grid(
  f_A_height_plot(pred = l_abu_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    labs(x = "Year", y = "Abundance") +
    theme(legend.position = "none",
          axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"])),
  f_A_height_plot(pred = l_ric_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "sric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    labs(x = "Year", y = "Richness") +
    theme(legend.position = "none",
          axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"])),
  f_A_height_plot(pred = l_mass_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       trans = log_plus_trans) +
    labs(x = "Year", y = "Biomass (g)", colour = "Elevation") +
    guides(colour = guide_legend(reverse = T)) +
    theme(axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"]),
          legend.title = element_text(size = v_textsize["legend.title"]),
          legend.text = element_text(size = v_textsize["axis.text"])),
  nrow = 1, rel_widths = c(1.05, 1, 1.6)
) 

ggsave("Output/Figures/Year_trends.pdf", width = 180, height = 65,
       units = "mm", dpi = 600)

# ... Figure 2 / Figure S3 #####################################################
################################################################################.

trans_spec <- scales::trans_new(name = "trans_spec",
                                transform = \(x) ifelse(x > 0 , -(1  / (x + 1) - 1),
                                                        x),
                                inverse = \(x) ifelse(x > 0, 1 / (1 - x) - 1, x))
c_factors <- c(1.1, 1.5, 2, 4)

# abundance --------------------------------------------------------------------.

d_incr_decr_abu <- f_A_height_incr_decr(l_abu_A$fit, d_mod_z, 
                                           filter(d_scalings, data == "full")) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_small$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "small"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_medium$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "medium"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_large$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "large"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_cold$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_intermediate$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_warm$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_spec_m$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Monophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_spec_o$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Oligophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_spec_p$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Polyphagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_egg$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "egg"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_larva$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "larva"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_pupa$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "pupa"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_abu_A_adult$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "adult"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

# richness ---------------------------------------------------------------------.

d_incr_decr_ric <- f_A_height_incr_decr(l_ric_A$fit, d_mod_z, 
                                           filter(d_scalings, data == "full")) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_small$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "small"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_medium$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "medium"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_large$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "large"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_cold$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_intermediate$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_warm$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_spec_m$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Monophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_spec_o$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Oligophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_spec_p$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Polyphagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_egg$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "egg"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_larva$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "larva"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_pupa$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "pupa"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_ric_A_adult$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "adult"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

# combined graph (abundance/richness) ------------------------------------------.

d_incr_decr_abu |> 
  mutate(response = "Abundance") |> 
  bind_rows(d_incr_decr_ric |> 
              mutate(response = "Richness")) |> 
  mutate(P_cut = cut(P_sum, breaks = c(0, .8, .9, .95, 1),
                     labels = c("I1", "I2", "I3", "I4"),
                     include.lowest = T, right = F),
         test = factor(test, levels = c("decrease", "increase"),
                       labels = c("Decr.", "Incr."))) |> 
  group_by(response, data, test) |> 
  arrange(P_sum) |> 
  mutate(first = !duplicated(P_cut) & P_cut != "I1") |> 
  ungroup() |> 
  ggplot(aes(x = threshold, y = P_sum)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            aes(fill = trait %in% c("Body size", "Specialisation")),
            colour = NA) +
  scale_fill_manual(values = c("white", "grey80"), guide = "none") +
  geom_hline(yintercept = c(0, 1), colour = "grey20") +
  new_scale("fill") +
  geom_tile(aes(fill = threshold_cut, height = P_sum, y = P_sum/2, width = xwidth)) +
  geom_line(size = .25) +
  geom_point(data = function(x) filter(x, first), aes(size = P_cut)) +
  facet_nested(trait + data ~ response + height + test, scales = "free_x") +
  scale_x_continuous(trans = trans_spec,
                     breaks = c(rev(1/c_factors - 1), c_factors - 1),
                     labels = c(paste0("1/", rev(c_factors)),
                                c_factors)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "PiYG")[-c(5,7)], guide = "none") +
  scale_size_manual(values = c(I2 = .6, I3 = .8, I4 = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"])) +
  labs(x = "50-year change", y = "Cumulative propability (%)")

ggsave("Output/Figures/Incr_decr_abu_ric.pdf", width = 180, height = 205,
       units = "mm", dpi = 600)

# biomass ----------------------------------------------------------------------.

d_incr_decr_mass <- f_A_height_incr_decr(l_mass_A$fit, d_mod_z, 
                                            filter(d_scalings, data == "full")) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_small$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "small"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_medium$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "medium"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_large$fit, 
                                    d_mod_mass_z %>% filter(mass_cat == "large"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_cold$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_intermediate$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_warm$fit, 
                                    d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_spec_m$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Monophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_spec_o$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Oligophagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_spec_p$fit, 
                                    d_mod_spec_z %>% filter(Spec == "Polyphagous"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_egg$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "egg"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_larva$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "larva"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_pupa$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "pupa"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_A_height_incr_decr(l_mass_A_adult$fit, 
                                    d_mod_hib_z %>% filter(overwintering_stage == "adult"), 
                                    filter(d_scalings, data == "full")) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_incr_decr_mass |> 
  mutate(P_cut = cut(P_sum, breaks = c(0, .8, .9, .95, 1),
                     labels = c("I1", "I2", "I3", "I4"),
                     include.lowest = T, right = F),
         test = factor(test, levels = c("decrease", "increase"),
                       labels = c("Decr.", "Incr."))) |> 
  group_by(data, test) |> 
  arrange(P_sum) |> 
  mutate(first = !duplicated(P_cut) & P_cut != "I1") |> 
  ungroup() |> 
  ggplot(aes(x = threshold, y = P_sum)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            aes(fill = trait %in% c("Body size", "Specialisation")),
            colour = NA) +
  scale_fill_manual(values = c("white", "grey80"), guide = "none") +
  geom_hline(yintercept = c(0, 1), colour = "grey20") +
  new_scale("fill") +
  geom_tile(aes(fill = threshold_cut, height = P_sum, y = P_sum/2, width = xwidth)) +
  geom_line(size = .25) +
  geom_point(data = function(x) filter(x, first), aes(size = P_cut)) +
  facet_nested(trait + data ~ height + test, scales = "free_x") +
  scale_x_continuous(trans = trans_spec,
                     breaks = c(rev(1/c_factors - 1), c_factors - 1),
                     labels = c(paste0("1/", rev(c_factors)),
                                c_factors)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "PiYG")[-c(5,7)], guide = "none") +
  scale_size_manual(values = c(I2 = .6, I3 = .8, I4 = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"])) +
  labs(x = "50-year change", y = "Cumulative propability (%)")

ggsave("Output/Figures/Incr_decr_mass.pdf", width = 90, height = 195,
       units = "mm", dpi = 600)

# ... Figure 3 / Figures S5 & S6 ###############################################
################################################################################.

p_legend <- cowplot::get_legend(f_A_height_plot_comb(list(small = l_abu_A_small$l_pred_fe$`A:height`),
                                                     response = "abu_tot") +
                                  theme(legend.title = element_text(size = v_textsize["legend.title"]),
                                        legend.text = element_text(size = v_textsize["axis.text"])))
p_empty <- ggplot() + theme_nothing()

# abundance --------------------------------------------------------------------.
plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_abu_A_small$l_pred_fe$`A:height`,
                                medium = l_abu_A_medium$l_pred_fe$`A:height`,
                                large = l_abu_A_large$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Body size") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_abu_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_abu_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_abu_A_warm$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Temperature niche") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_abu_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_abu_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_abu_A_spec_p$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Food specialisation") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .305)),
  f_A_height_plot_comb(list(egg = l_abu_A_egg$l_pred_fe$`A:height`,
                            larva = l_abu_A_larva$l_pred_fe$`A:height`,
                            pupa = l_abu_A_pupa$l_pred_fe$`A:height`,
                            adult = l_abu_A_adult$l_pred_fe$`A:height`),
                       response = "abu_tot",
                       name = "Overwintering stage") +
    labs(y = "Abundance") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, .4, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_abu.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# richness ---------------------------------------------------------------------.

plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_ric_A_small$l_pred_fe$`A:height`,
                                medium = l_ric_A_medium$l_pred_fe$`A:height`,
                                large = l_ric_A_large$l_pred_fe$`A:height`),
                           response = "sric",
                           name = "Body size") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_ric_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_ric_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_ric_A_warm$l_pred_fe$`A:height`),
                           response = "sric",
                           name = "Temperature niche") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_ric_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_ric_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_ric_A_spec_p$l_pred_fe$`A:height`),
                           response = "sric",
                           name = "Food specialisation") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .32)),
  f_A_height_plot_comb(list(egg = l_ric_A_egg$l_pred_fe$`A:height`,
                            larva = l_ric_A_larva$l_pred_fe$`A:height`,
                            pupa = l_ric_A_pupa$l_pred_fe$`A:height`,
                            adult = l_ric_A_adult$l_pred_fe$`A:height`),
                       response = "sric",
                       name = "Overwintering stage") +
    labs(y = "Richness") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, 0, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_ric.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# biomass ----------------------------------------------------------------------.

plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_mass_A_small$l_pred_fe$`A:height`,
                                medium = l_mass_A_medium$l_pred_fe$`A:height`,
                                large = l_mass_A_large$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Body size") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_mass_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_mass_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_mass_A_warm$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Temperature niche") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_mass_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_mass_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_mass_A_spec_p$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Food specialisation") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .305)),
  f_A_height_plot_comb(list(egg = l_mass_A_egg$l_pred_fe$`A:height`,
                            larva = l_mass_A_larva$l_pred_fe$`A:height`,
                            pupa = l_mass_A_pupa$l_pred_fe$`A:height`,
                            adult = l_mass_A_adult$l_pred_fe$`A:height`),
                       response = "mass_tot",
                       name = "Overwintering stage") +
    labs(y = "Biomass (g)") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, .4, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_mass.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# ... Figure S4 ################################################################
################################################################################.

plot_grid(
  # abundance ------------------------------------------------------------------.
  f_A_height_plot(pred = l_abu_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Abundance") +
    theme(legend.position = "none",
          axis.text.x = element_blank()) +
    facet_grid(~ "Full model"),
  f_A_height_plot(pred = l_abu_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Fixed only"),
  f_A_height_plot(pred = l_abu_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Manual only"),
  f_A_height_plot(pred = l_abu_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Full data only"),
  
  # richness -------------------------------------------------------------------.
  f_A_height_plot(pred = l_ric_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "sric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p",
                       labels = c("0", "10", "    100",  "1000")) +
    coord_cartesian(ylim = c(NA, 300),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Richness") +
    theme(legend.position = "none",
          axis.text.x = element_blank()),
  f_A_height_plot(pred = l_ric_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "sric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 300),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  f_A_height_plot(pred = l_ric_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "sric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 300),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  f_A_height_plot(pred = l_ric_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "sric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 300),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  
  # biomass --------------------------------------------------------------------.
  f_A_height_plot(pred = l_mass_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c(0, 0.0001, " 0.001", 0.01, 0.1, 1),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, .5),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Biomass (g)") +
    theme(legend.position = "none"),
  f_A_height_plot(pred = l_mass_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, .5),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  f_A_height_plot(pred = l_mass_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, .5),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  f_A_height_plot(pred = l_mass_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, .5),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  
  nrow = 3, rel_widths = c(1, .75, .75, .75))


ggsave("Output/Figures/Year_trends_sensibility.jpeg", width = 180, height = 180,
       units = "mm", dpi = 400)

# ... Tables S2-S17 ############################################################
################################################################################.

# full models ------------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") +
    (1 | gr) +
    (1 | trap_ID) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 9, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height.xlsx", overwrite = T)

# full observation hours data --------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_ne")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") +
    (1 | gr) +
    (1 | trap_ID) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_ne_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_ne.xlsx", overwrite = T)

# fixed traps only -------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_lf")))
  
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    C(sample_previous, "contr.sum") +
    (1 | gr) +
    (1 | trap_ID) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_lf_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_LF.xlsx", overwrite = T)


# manual traps only ------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_p")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day + 
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") + 
    (1 | gr) + 
    (1 | trap_ID) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_p_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_p.xlsx", overwrite = T)

# ------------------------------------------------------------------------------.
# trait subsets ----------------------------------------------------------------.
# ------------------------------------------------------------------------------.

d_traits <- data.frame(trait = "mass", traitvalue = "small") |> 
  add_row(trait = "mass", traitvalue = "medium") |> 
  add_row(trait = "mass", traitvalue = "large") |> 
  add_row(trait = "Tavg", traitvalue = "cold") |> 
  add_row(trait = "Tavg", traitvalue = "intermediate") |> 
  add_row(trait = "Tavg", traitvalue = "warm") |> 
  add_row(trait = "spec", traitvalue = "m") |> 
  add_row(trait = "spec", traitvalue = "o") |> 
  add_row(trait = "spec", traitvalue = "p") |> 
  add_row(trait = "hib", traitvalue = "egg") |> 
  add_row(trait = "hib", traitvalue = "larva") |> 
  add_row(trait = "hib", traitvalue = "pupa") |> 
  add_row(trait = "hib", traitvalue = "adult")

c_traits <- c(mass = "mass_cat",
              Tavg = "Tavg_mean_cat",
              spec = "Spec",
              hib = "overwintering_stage")
c_traitvalues <- c(m = "Monophagous",
                   o = "Oligophagous",
                   p = "Polyphagous")

for (trait_i in unique(d_traits$trait)){
  wb_out <- createWorkbook()
  for (resp_i in c("abu", "ric", "mass")){
    for (traitvalue_i in d_traits$traitvalue[d_traits$trait == trait_i]){
      
      if (trait_i == "spec"){
        mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_spec_", traitvalue_i)))
      } else {
        mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_", traitvalue_i)))
      }
      
      formula_i <- response ~ s(yday) + P_2day + T_2day +
        C(traptype, "contr.sum") +
        C(bulbtype, "contr.sum") +
        n_trap +
        C(sample_previous, "contr.sum") +
        (1 | gr) +
        (1 | trap_ID) +
        (1 | night_ID) +
        (1 | trap_ID_A)
      
      data <- eval(parse(text = paste0("d_mod_", trait_i, "_z")))
      
      if (trait_i == "spec"){
        data <- data |> 
          filter(!! sym(c_traits[trait_i]) == c_traitvalues[traitvalue_i])
      } else {
        data <- data |> 
          filter(!! sym(c_traits[trait_i]) == traitvalue_i)
      }
      
      
      table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                                data = data, model = "year")
      
      addWorksheet(wb_out, paste0(resp_i, "_", traitvalue_i))
      
      writeDataTable(wb_out, paste0(resp_i, "_", traitvalue_i), table_i)
      
      for (c_i in c(1, 2)){
        goon <- T
        start <- 1
        while(goon){
          
          if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
            stop <- nrow(table_i)
          } else {
            stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
          }
          
          
          if (is.na(stop)) break
          
          if (stop < start) stop <- start
          
          mergeCells(wb_out, paste0(resp_i, "_", traitvalue_i), cols = c_i, rows = seq(start + 1, stop + 1))
          
          if (stop == nrow(table_i)) goon <- F
          
          start <- stop + 1
        }
      }
      
      
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = seq_len(nrow(table_i) + 1), 
               cols = seq_len(ncol(table_i)), gridExpand = T,
               style = createStyle(fontSize = 8, fontName = "Helvetica", 
                                   valign = 'center', border = "TopBottomLeftRight"))
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 1, cols = 1:100, gridExpand = T,
               style = createStyle(textDecoration = 'bold'), stack = T)
      rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                            (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                               as.numeric(table_i$`Upper 95%-CI`) < 0))
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
               style = createStyle(textDecoration = 'bold'), stack = T)
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 2:100, cols = 4:6, gridExpand = T,
               style = createStyle(halign = 'right'), stack = T)
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 1, cols = 4:6, gridExpand = T,
               style = createStyle(halign = 'center'), stack = T)
      
      setColWidths(wb_out, paste0(resp_i, "_", traitvalue_i), seq_len(ncol(table_i)), widths = "auto")
    }
  }
  
  saveWorkbook(wb_out, file = paste0("Output/Tables/Modsummaries_A_height_", trait_i, ".xlsx"), overwrite = T)
}
