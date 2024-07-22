# SETUP ------------------------------- ########################################
################################################################################.

# ... packages #################################################################
################################################################################.

library(tidyverse); theme_set(theme_classic())
library(rstan); options(mc.cores = 4)
library(brms)
library(bayestestR)
library(mgcv)
library(lubridate)
library(cowplot)

# ... set global parameters ####################################################
################################################################################.

min_mass <- min(d_mod_z$mass_tot[d_mod_z$mass_tot > 0])

log_plus_trans <- scales::trans_new(name = "log_plus",
                                    transform = \(x) log(x + min_mass, 10),
                                    inverse = \(x) out <- 10 ^ (x) - min_mass)

n_iter <- 2000 # number of MCMC iterations

v_covlabels <- c(height = "Elevation (m)",yday = "Day of the year", 
                 P_2day = "Precipitation (mm)", T_2day = "Temperature (°C)",
                 traptype = "Trap type", bulbtype = "Lamp type", 
                 n_trap = "Number of traps", 
                 sample_previous = "Sampling in previous night",
                 active_hours = "Hours active",
                 prop_forest_500 = "Prop. forests",
                 prop_grassland_500 = "Prop. grasslands",
                 prop_crop_500 = "Prop. croplands", 
                 prop_sealed_500 = "Prop. sealed area")
v_covlabels_short <- c(height = "Elevation", 
                       yday = "Day of year", 
                       P_2day = "Prec.", T_2day = "Temp.",
                       traptype = "Trap type", bulbtype = "Lamp type", 
                       n_trap = "Nr. of traps", 
                       sample_previous = "Sampl. prev. night",
                       active_hours = "Hrs. active",
                       gr = "Spat-temp. cluster", trap_ID = "Site ID",
                       prop_forest_500 = "Prop. forests",
                       prop_grassland_500 = "Prop. grasslands",
                       prop_crop_500 =  "Prop. croplands",
                       prop_sealed_500 = "Prop. sealed area",
                       trap_ID_A = "Site ID × Year",
                       night_ID = "Sampling night × location",
                       Intercept = "Intercept")

v_traplabels <- c(LF = "Fixed (t1)", "LF-Changins" = "Fixed (t2)", p = "Manual")
v_bulblabels <-  c("80MLL" = "80W HWL", "150-160MLL" = "150–160W HWL", 
                   "125HQL/150-160MLL" = "150–160W HWL/125W HQL",
                   "125HQL" = "125W HQL")

v_textsize <- c(axis.title = 8, axis.text = 7, legend.title = 8, additional.text = 6) 

# ... define functions #########################################################
################################################################################.

f_analysis_LU <- function(formula, data_z, data, scalings, 
                          family, hours_sel, iter, seed,
                          LU_vars, run_loo = FALSE, extract_random = F){
  
  set.seed(seed)
  
  out <- list()
  
  response <- as.character(formula[2])
  
  formula_full <- update(formula,
                         paste0(". ~ s(height) + ", 
                                paste(LU_vars, collapse = " + "),
                                "+ ."))
  
  # prepare data
  
  l_data <- make_standata(formula_full, 
                          data = data_z)
  
  
  l_data_sub <- make_standata(update(formula_full, ". ~ s(active_hours)"),
                              data = data_z[hours_sel, ])
  
  l_data$N_SEL <- length(hours_sel)
  l_data$SEL <-  hours_sel
  l_data$Xs_add <- l_data_sub$Xs
  l_data$nb_add <- 1
  l_data$knots_add <- l_data$knots_1
  dimnames(l_data$knots_add)[[1]] <- "Zs_add_1"
  l_data$Zs_add_1 <- l_data_sub$Zs_1_1
  
  
  # run model
  
  if (family == "zero_inflated_negbinomial"){
    mod <- stan(file = "Stan_splines/Stan_nb_spline_s2p1_r4.stan",
                data = l_data,
                chains = 4, cores = 4,
                iter = iter)
  } else if (family == "hurdle_gamma"){
    mod <- stan(file = "Stan_splines/Stan_hg_spline_s2p1_r4.stan",
                data = l_data,
                chains = 4, cores = 4,
                iter = iter)
  }
  
  # extract parameter data
  par_names <- names(mod)
  pos <- regexpr("\\[", par_names) 
  pos <- ifelse(pos == -1, 200, pos - 1)
  pars <- unique(substr(par_names, 1, pos))
  pars <- pars[pars %in% c("b", "Intercept", "bs", "bs_add", "zs_1_1",
                           "sds_1_1", "zs_2_1", "sds_2_1", "zs_add_1",
                           "sds_add_1", "shape", "zi", "hu", "sd_1", "sd_2",
                           "sd_3", "sd_4", "s_1_1", "s_2_1", "s_add_1")]
  if (extract_random){
    pars_plus <- unique(substr(par_names, 1, pos))
    pars_plus <- pars_plus[grepl("^r_", pars_plus)]
    pars <- c(pars, pars_plus)
  }
  
  fit <- rstan::extract(mod, pars = pars)
  out$fit <- fit
  
  if (run_loo) out$loo <- loo(mod)
  
  terms <- brmsterms(formula_full)
  
  # predictions for smooths
  out$l_pred_sm <- f_pred_smooths(fit = fit, data = data_z, 
                                  l_data = l_data,
                                  terms = terms, scalings = scalings,
                                  hours_sel = hours_sel)
  
  
  # predictions for fixed effects
  out$l_pred_fe <-   f_pred_fixed_covariate(fit = fit, data = data_z, 
                                            terms = terms, scalings = scalings,
                                            vars = c("P_2day", "T_2day", "traptype",
                                                     "bulbtype", "n_trap", "sample_previous",
                                                     LU_vars))

  out
}

f_pred_fixed_covariate <- function(fit, data, terms, scalings, vars) {
  
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  
  l_out <- list()
  
  vars_loop <- all.vars(terms$dpars$mu$fe)
  vars_loop <- vars_loop[vars_loop %in% vars]
  
  for (var_i in vars_loop){
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
      
      # centre variables
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

f_pred_fixed_diff <- function(fit, data, scalings, formula = NULL) {
  
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
  
  # effects without interactions -----------------------------------------------.
  
  d_out <- data.frame()
  vars <- all.vars(terms$dpars$mu$fe)
  
  # special case: height_cat only in interactions
  if (!grepl(" height_cat ", formula[3])){
    vars <- vars[vars != "height_cat"]
  }
  
  for (var_i in vars){
    if (is.numeric(deframe(data[, var_i]))) {
      
      d_pred <- data.frame(expl = min(data[, var_i]) + 
                             diff(range(data[, var_i])) * c(.25, .75)) %>% 
        mutate(expl_c = expl - mean(data[, var_i, drop = T]),
               expl_or = expl * scalings$sd[scalings$var == var_i] +
                 scalings$mean[scalings$var == var_i])
      
      
      m_pred <- matrix(rep(fit$Intercept, each = 2), nrow = 2) + 
        matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                  !grepl(":", colnames(mm_full)))]
      
      
      c_diff <- apply(m_pred, 2, diff) 
      
      d_out_tmp <- data.frame(diff = c_diff) |> 
        summarise(estimate_mean = mean(diff),
                  lower_95 = ci(diff, .95)$CI_low,
                  upper_95 = ci(diff, .95)$CI_high,
                  lower_80 = ci(diff, .80)$CI_low,
                  upper_80 = ci(diff, .80)$CI_high,
                  lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
                  upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
                  
                  exp_estimate_mean = mean(exp(diff)),
                  exp_lower_95 = ci(exp(diff), .95)$CI_low,
                  exp_upper_95 = ci(exp(diff), .95)$CI_high,
                  exp_lower_80 = ci(exp(diff), .80)$CI_low,
                  exp_upper_80 = ci(exp(diff), .80)$CI_high,
                  exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
                  exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high) |> 
        mutate(var = var_i,
               diff_or = diff(d_pred$expl_or))
      
      if (grepl(paste0(var_i, ":height_cat"), formula[3])){
        d_out_tmp <- d_out_tmp |> 
          mutate(height_cat = "low")
      }
      d_out <- d_out_tmp |> 
        select(var, diff_or, any_of("height_cat"), everything()) |> 
        (\(x) bind_rows(d_out, x))()
      
      # height cat interaction effect
      if (grepl(paste0(var_i, ":height_cat"), formula[3])){
        m_pred <- matrix(rep(fit$Intercept, each = 2), nrow = 2) + 
          matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                    !grepl(":", colnames(mm_full)))] +
          matrix(d_pred$expl_c) %*% fit$b[, colnames(mm_full) == paste0(var_i, ":height_cathigh")]
        
        
        c_diff <- apply(m_pred, 2, diff) 
        
        d_out <- data.frame(diff = c_diff) |> 
          summarise(estimate_mean = mean(diff),
                    lower_95 = ci(diff, .95)$CI_low,
                    upper_95 = ci(diff, .95)$CI_high,
                    lower_80 = ci(diff, .80)$CI_low,
                    upper_80 = ci(diff, .80)$CI_high,
                    lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
                    upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
                    
                    exp_estimate_mean = mean(exp(diff)),
                    exp_lower_95 = ci(exp(diff), .95)$CI_low,
                    exp_upper_95 = ci(exp(diff), .95)$CI_high,
                    exp_lower_80 = ci(exp(diff), .80)$CI_low,
                    exp_upper_80 = ci(exp(diff), .80)$CI_high,
                    exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
                    exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high) |> 
          mutate(var = var_i,
                 height_cat = "high",
                 diff_or = diff(d_pred$expl_or)) |> 
          select(var, diff_or, height_cat, everything()) |> 
          (\(x) bind_rows(d_out, x))()
      }
      
    } else {
      
      d_pred <- data.frame(expl = unique(data[, var_i]))
      
      comps <- strsplit(paste(deparse(terms$dpars$mu$fe), collapse = ""), split = " \\+ ")[[1]]
      
      comps <- comps[grepl(var_i, comps)]
      
      mm_data <- model.matrix(as.formula(paste("~", comps)), data = data)[, -1, drop = F]
      mm_pred <- model.matrix(as.formula(paste("~", comps)), data = d_pred)[, -1, drop = F]
      
      # centre variables
      for (i in seq_len(ncol(mm_pred))){
        mm_pred[, i] <- mm_pred[, i] - mean(mm_data[, i])
      }
      
      m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred)) + 
        as.matrix(mm_pred) %*% t(as.matrix(fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                           !grepl(":", colnames(mm_full)))]))
      
      d_diff_raw <- data.frame()
      for (i in seq_len(nrow(m_pred))[-1]){
        d_diff_raw <- data.frame(diff = apply(m_pred[c(1, i), ], 2, diff),
                                 varlevel = as.character(d_pred[, var_i][i]),
                                 stringsAsFactors = F) |> 
          (\(x) bind_rows(d_diff_raw, x))()
        
      }
      
      d_out <- d_diff_raw |> 
        group_by(varlevel) |> 
        summarise(estimate_mean = mean(diff),
                  lower_95 = ci(diff, .95)$CI_low,
                  upper_95 = ci(diff, .95)$CI_high,
                  lower_80 = ci(diff, .80)$CI_low,
                  upper_80 = ci(diff, .80)$CI_high,
                  lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
                  upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
                  
                  exp_estimate_mean = mean(exp(diff)),
                  exp_lower_95 = ci(exp(diff), .95)$CI_low,
                  exp_upper_95 = ci(exp(diff), .95)$CI_high,
                  exp_lower_80 = ci(exp(diff), .80)$CI_low,
                  exp_upper_80 = ci(exp(diff), .80)$CI_high,
                  exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
                  exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high,
                  .groups = "drop") |> 
        mutate(var = var_i) |> 
        select(var, varlevel, everything()) |> 
        (\(x) bind_rows(d_out, x))()
    }
  }
  
  d_out <- d_out |> 
    select(var, diff_or, varlevel, any_of("height_cat"), everything())
  
  d_out
}

f_pred_fixed_diff_specific <- function(fit, data, scalings, formula = NULL,
                                       var_i, diff_or) {
  
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
  
  # effects without interactions -----------------------------------------------.
  
  d_out <- data.frame()
  vars <- all.vars(terms$dpars$mu$fe)
  
  # special case: height_cat only in interactions
  if (!grepl(" height_cat ", formula[3])){
    vars <- vars[vars != "height_cat"]
  }
  
  d_pred <- data.frame(expl = diff_or / scalings$sd[scalings$var == var_i] * c(0, 1)) |> 
    mutate(expl_c = expl - mean(data[, var_i, drop = T]),
           expl_or = expl * scalings$sd[scalings$var == var_i] +
             scalings$mean[scalings$var == var_i])
  
  
  m_pred <- matrix(rep(fit$Intercept, each = 2), nrow = 2) + 
    matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                              !grepl(":", colnames(mm_full)))]
  
  
  c_diff <- apply(m_pred, 2, diff) 
  
  d_out_tmp <- data.frame(diff = c_diff) |> 
    summarise(estimate_mean = mean(diff),
              lower_95 = ci(diff, .95)$CI_low,
              upper_95 = ci(diff, .95)$CI_high,
              lower_80 = ci(diff, .80)$CI_low,
              upper_80 = ci(diff, .80)$CI_high,
              lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
              upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
              
              exp_estimate_mean = mean(exp(diff)),
              exp_lower_95 = ci(exp(diff), .95)$CI_low,
              exp_upper_95 = ci(exp(diff), .95)$CI_high,
              exp_lower_80 = ci(exp(diff), .80)$CI_low,
              exp_upper_80 = ci(exp(diff), .80)$CI_high,
              exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
              exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high) |> 
    mutate(var = var_i,
           diff_or = diff(d_pred$expl_or))
  
  if (grepl(paste0(var_i, ":height_cat"), formula[3])){
    d_out_tmp <- d_out_tmp |> 
      mutate(height_cat = "low")
  }
  d_out <- d_out_tmp |> 
    select(var, diff_or, any_of("height_cat"), everything()) |> 
    (\(x) bind_rows(d_out, x))()
  
  # height cat interaction effect
  if (grepl(paste0(var_i, ":height_cat"), formula[3])){
    m_pred <- matrix(rep(fit$Intercept, each = 2), nrow = 2) + 
      matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                !grepl(":", colnames(mm_full)))] +
      matrix(d_pred$expl_c) %*% fit$b[, colnames(mm_full) == paste0(var_i, ":height_cathigh")]
    
    
    c_diff <- apply(m_pred, 2, diff) 
    
    d_out <- data.frame(diff = c_diff) |> 
      summarise(estimate_mean = mean(diff),
                lower_95 = ci(diff, .95)$CI_low,
                upper_95 = ci(diff, .95)$CI_high,
                lower_80 = ci(diff, .80)$CI_low,
                upper_80 = ci(diff, .80)$CI_high,
                lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
                upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
                
                exp_estimate_mean = mean(exp(diff)),
                exp_lower_95 = ci(exp(diff), .95)$CI_low,
                exp_upper_95 = ci(exp(diff), .95)$CI_high,
                exp_lower_80 = ci(exp(diff), .80)$CI_low,
                exp_upper_80 = ci(exp(diff), .80)$CI_high,
                exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
                exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high) |> 
      mutate(var = var_i,
             height_cat = "high",
             diff_or = diff(d_pred$expl_or)) |> 
      select(var, diff_or, height_cat, everything()) |> 
      (\(x) bind_rows(d_out, x))()
  }
  
  d_out   
  
}

f_pred_smooths_diff <- function(fit, data, scalings, formula = NULL) {
  
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
  
  # effects without interactions -----------------------------------------------.
  
  d_out <- data.frame()
  vars <- all.vars(terms$dpars$mu$sm)
  
  l_data <- make_standata(formula, 
                          data = data |> 
                            mutate(response = 1)) # dummy variable
  
  for (var_i in vars){
    
    d_target_unscaled <- data[, var_i] %>% 
      rename(var_target = !! sym(var_i))
    
    d_target_scaled <- l_data$Xs[, which(colnames(l_data$Xs) == paste0("s", var_i, "_1"))]  %>% 
      as.data.frame()
    
    count_i <- which(all.vars(terms$dpars$mu$sm) == var_i)
    
    v_bs <- fit$bs[, count_i]
    v_s <- fit[[paste0("s_", count_i, "_1")]]
    
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
    
    v_50perc <- diff(range(d_target_unscaled$var_target)) * .5
    # slowly move the 50% window through the range (resolution 100, i.e. 100 steps)
    v_plus <- v_50perc * seq(0, 1, length.out = 100)
    
    m_diff <- c()
    for (plus_i in v_plus){
      d_newdata_unscaled <- data.frame(var_target = c(min(d_target_unscaled$var_target) + plus_i, 
                                                      min(d_target_unscaled$var_target) + v_50perc + plus_i))
      
      pred_sm <- PredictMat(sm_unscaled[[1]], d_newdata_unscaled)
      
      d_newdata_scaled <- data.frame(var_target = pred_sm[, 9])
      
      pred_sm <- pred_sm[, -9][, c(1,8:2)]
      
      m_pred <- matrix(rep(fit$Intercept, each = nrow(d_newdata_scaled)),
                       nrow = nrow(d_newdata_scaled)) + 
        as.matrix(d_newdata_scaled$var_target) %*%  v_bs + 
        pred_sm %*% t(v_s)
      
      c_diff <- apply(m_pred, 2, diff) 
      
      m_diff <- cbind(m_diff, c_diff)
    }
    
    c_diff_mean <- rowMeans(m_diff)
    
    d_out <- data.frame(diff = c_diff_mean) |> 
      summarise(estimate_mean = mean(diff),
                lower_95 = ci(diff, .95)$CI_low,
                upper_95 = ci(diff, .95)$CI_high,
                lower_80 = ci(diff, .80)$CI_low,
                upper_80 = ci(diff, .80)$CI_high,
                lower_776 = ci(diff, 1 - sqrt(.05))$CI_low,
                upper_776 = ci(diff, 1 - sqrt(.05))$CI_high,
                
                exp_estimate_mean = mean(exp(diff)),
                exp_lower_95 = ci(exp(diff), .95)$CI_low,
                exp_upper_95 = ci(exp(diff), .95)$CI_high,
                exp_lower_80 = ci(exp(diff), .80)$CI_low,
                exp_upper_80 = ci(exp(diff), .80)$CI_high,
                exp_lower_776 = ci(exp(diff), 1 - sqrt(.05))$CI_low,
                exp_upper_776 = ci(exp(diff), 1 - sqrt(.05))$CI_high) |> 
      mutate(var = var_i,
             diff_or = v_50perc * scalings$sd[scalings$var == var_i]) |> 
      select(var, diff_or, everything()) |> 
      (\(x) bind_rows(d_out, x))()
    
  }
  
  
  d_out <- d_out |> 
    select(var, diff_or, everything())
  
  d_out
}

f_local_maxima <- function(x){
  out <- which(diff(sign(diff(x))) == -2) + 1
  
  if (length(out) == 0){
    NA
  } else {
    out
  }
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

f_plot_pred <- function(x, data, response, hours_sel = NULL, line.size = 1){
  var_i <- names(x)[1]
  
  label <- ifelse(var_i %in% names(v_covlabels),
                  v_covlabels[var_i],
                  var_i)
  
  if (var_i == "traptype"){
    x <- x %>% 
      mutate(traptype = factor(traptype, levels = names(v_traplabels),
                               labels = v_traplabels))
    data <- data %>% 
      mutate(traptype = factor(traptype, levels = names(v_traplabels),
                               labels = v_traplabels))
  } else if (var_i == "bulbtype"){
    x <- x %>% 
      mutate(bulbtype = factor(bulbtype, levels = names(v_bulblabels),
                               labels = v_bulblabels))
    data <- data %>% 
      mutate(bulbtype = factor(bulbtype, levels = names(v_bulblabels),
                               labels = v_bulblabels))
  }
  
  
  if (var_i == "active_hours"){
    p <- ggplot(x, aes_string(x = var_i))  +
      geom_point(data = data[hours_sel, ], aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4") + 
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = .5) +
      geom_line(aes(y = estimate), size = line.size, colour = "darkgreen") 
  } else if (is.numeric(x[, var_i, drop = T])){
    p <- ggplot(x, aes_string(x = var_i)) +
      geom_point(data = data, aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4") + 
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = .5) +
      geom_line(aes(y = estimate), size = line.size, colour = "darkgreen")
  } else {
    p <- ggplot(x, aes_string(x = var_i)) +
      geom_jitter(data = data, aes_string(y = response), 
                  alpha = .05, size = .4, col = "salmon4", width = .25) + 
      geom_segment(aes_string(y = "lower", yend = "upper", xend = var_i), 
                   col = "grey20", 
                   arrow = arrow(angle = 90, ends = "both", length = unit(.1, "inches"))) +
      geom_point(aes(y = estimate), size = 2, colour = "darkgreen")
  }
  
  p <- p +
    ylab("Estimate") +
    xlab(label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, NA))
  
  if (response == "mass_tot"){
    p <-  p +
      scale_y_continuous(trans = log_plus_trans)
  } else {
    p <-  p +
      scale_y_continuous(trans = "log1p")
  }
  
  p
}

f_summarytable <- function(mod, formula, covlabels = v_covlabels_short, data, 
                           LU_vars, ...){

  formula <- update(formula,
                    paste0(". ~ s(height) + ", 
                           paste(LU_vars, collapse = " + "),
                           "+ ."))
  
  
  terms <- brmsterms(formula)
  mm_full <- model.matrix(terms$dpars$mu$fe, data = droplevels(data))[, -1]
  
  # Fixed components -----------------------------------------------------------.
  
  out <- data.frame(parameter = "Fixed effect",
                    type = "Fixed", 
                    var = colnames(mm_full),
                    estimate = apply(mod$b, 2, mean),
                    lower = apply(mod$b, 2, \(x) ci(x)$CI_low),
                    upper = apply(mod$b, 2, \(x) ci(x)$CI_high))
  
  # Spline components ----------------------------------------------------------.
  
  out <- data.frame(parameter = "Fixed effect",
                    type = "Spline", 
                    var = all.vars(terms$dpars$mu$sm),
                    estimate = apply(mod$bs, 2, mean),
                    lower = apply(mod$bs, 2, \(x) ci(x)$CI_low),
                    upper = apply(mod$bs, 2, \(x) ci(x)$CI_high)) %>% 
    bind_rows(out, .)
  
  for (i in seq_len(length(all.vars(terms$dpars$mu$sm)))){
    out <- data.frame(parameter = "SD",
                      type = "Spline", 
                      var = all.vars(terms$dpars$mu$sm)[i],
                      estimate = mean(mod[[paste0("sds_", i, "_1")]]),
                      lower = ci(mod[[paste0("sds_", i, "_1")]])$CI_low,
                      upper = ci(mod[[paste0("sds_", i, "_1")]])$CI_high) %>% 
      bind_rows(out, .)
    
    
  }
  
  # additional spline components (active_hours) --------------------------------.
  
  if ("bs_add" %in% names(mod)){
    out <- data.frame(parameter = "Fixed effect",
                      type = "Spline", 
                      var = "active_hours",
                      estimate = mean(mod$bs_add),
                      lower = ci(mod$bs_add)$CI_low,
                      upper = ci(mod$bs_add)$CI_high) %>% 
      bind_rows(out, .)
    
    out <- data.frame(parameter = "SD",
                      type = "Spline", 
                      var = "active_hours",
                      estimate = mean(mod$sds_add_1),
                      lower = ci(mod$sds_add_1)$CI_low,
                      upper = ci(mod$sds_add_1)$CI_high) %>% 
      bind_rows(out, .)
  }
  
  # random factors -------------------------------------------------------------.
  
  v_r <- terms$dpars$mu$re$group
  
  for (i in seq_len(length(v_r))){
    out <- data.frame(parameter = "SD",
                      type = "Random", 
                      var = v_r[i],
                      estimate = mean(mod[[paste0("sd_", i)]]),
                      lower = ci(mod[[paste0("sd_", i)]])$CI_low,
                      upper = ci(mod[[paste0("sd_", i)]])$CI_high) %>% 
      bind_rows(out, .)
    
  }
  
  
  out <-
    out %>% 
    mutate_at(vars(estimate, lower, upper),
              ~ as.character(formatC(., digits = 3, format = "fg", flag = "#"))) %>% 
    mutate(varplus = ifelse(grepl('^C.*"contr.sum")[0-9]$', var), 
                            paste0("(contr. sum ", substr(var, nchar(var), 
                                                          nchar(var)), ")"), ""),
           var = ifelse(grepl('^C.*"contr.sum")[0-9]$', var), 
                        substr(var, 3, nchar(var) - 15), var),
           varplus = ifelse(grepl("\\.L$", var), "(linear)", varplus),
           varplus = ifelse(grepl("\\.Q$", var), "(quadratic)", varplus),
           varplus = ifelse(grepl("\\.C$", var), "(cubic)", varplus),
           var = ifelse(grepl("\\.L$|\\.Q$|\\.C$", var), 
                        substr(var, 1, nchar(var) - 2), var)) %>% 
    rowwise() |> 
    mutate(var = ifelse(grepl(":", var),
                        paste(covlabels[strsplit(var, ":")[[1]]],
                              collapse = " × "),
                        covlabels[var]),
           var = ifelse(varplus != "", paste(var, varplus), var)) %>% 
    ungroup() |> 
    arrange(parameter, type) %>% 
    select(Parameter = parameter, 
           Type = type, 
           Variable = var, 
           Estimate = estimate, 
           "Lower 95%-CI" = lower, 
           "Upper 95%-CI" = upper) |> 
    as.data.frame()
  
}

# RUN MODELS -------------------- ##############################################
################################################################################.

# ... overall models ###########################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_LU_500 <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 923,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              run_loo = T)

# richness ---------------------------------------------------------------------.
l_ric_LU_500 <- f_analysis_LU(formula = sric ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 87127,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              run_loo = T)

# biomass ----------------------------------------------------------------------.
l_mass_LU_500 <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 1411,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              run_loo = T)


# ... per overwintering stage ##################################################
################################################################################.

# abundance --------------------------------------------------------------------. 
l_abu_LU_egg <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "egg"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                              iter = n_iter, seed = 897444,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_abu_LU_larva <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "larva"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                              iter = n_iter, seed = 913,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_abu_LU_pupa <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "pupa"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                              iter = n_iter, seed = 817,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_abu_LU_adult <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "adult"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                              iter = n_iter, seed = 1253,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))


# richness ---------------------------------------------------------------------.
l_ric_LU_egg <- f_analysis_LU(formula = sric ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "egg"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                              iter = n_iter, seed = 183456,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_ric_LU_larva <- f_analysis_LU(formula = sric ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "larva"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                              iter = n_iter, seed = 512344,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_ric_LU_pupa <- f_analysis_LU(formula = sric ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "pupa"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                              iter = n_iter, seed = 51,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_ric_LU_adult <- f_analysis_LU(formula = sric ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "adult"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                              iter = n_iter, seed = 662,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

# biomass ----------------------------------------------------------------------.
l_mass_LU_egg <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "egg"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "egg"]),
                              iter = n_iter, seed = 918,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_mass_LU_larva <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "larva"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "larva"]),
                              iter = n_iter, seed = 234,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_mass_LU_pupa <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "pupa"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "pupa"]),
                              iter = n_iter, seed = 15144,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

l_mass_LU_adult <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod_hib %>% filter(overwintering_stage == "adult"),
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                  !d_mod_hib_z$estimate[d_mod_hib_z$overwintering_stage == "adult"]),
                              iter = n_iter, seed = 621,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"))

# CREATE OUTPUTS -------------- ################################################
################################################################################.

# define model formula
LU_vars <- c("Forest" = "prop_forest_500",
             "Grassland" = "prop_grassland_500",
             "Crop" = "prop_crop_500",
             "Sealed" = "prop_sealed_500")
formula_full <- response ~ 
  s(yday) + P_2day + T_2day +
  C(traptype, "contr.sum") +
  C(bulbtype, "contr.sum") +
  n_trap +
  C(sample_previous, "contr.sum") +
  (1 | gr) +
  (1 | trap_ID) +
  (1 | night_ID) +
  (1 | trap_ID_A)
formula_full <- update(formula_full,
                       paste0(". ~ s(height) + ", 
                              paste(LU_vars, collapse = " + "),
                              "+ ."))

# ... Text #####################################################################
################################################################################.

# factor estimates for fixed variables:
f_pred_fixed_diff(l_abu_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                  formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, varlevel, smry)
f_pred_fixed_diff(l_ric_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                  formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, varlevel, smry)
f_pred_fixed_diff(l_mass_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                  formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, varlevel, smry)

# factor estimates for smoothing terms:
f_pred_smooths_diff(l_abu_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                    formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, smry)
f_pred_smooths_diff(l_ric_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                    formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, smry)
f_pred_smooths_diff(l_mass_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                    formula = formula_full) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(var, diff_or, smry)


# where are the local maxima of the smoothing terms:
l_abu_LU_500$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))
l_ric_LU_500$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))
l_mass_LU_500$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))

l_abu_LU_egg$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))
l_ric_LU_egg$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))
l_mass_LU_egg$l_pred_sm$yday |> 
  slice(f_local_maxima(estimate)) |> 
  mutate(Date = parse_date_time(as.character(round(yday)), orders = "j"))

l_abu_LU_500$l_pred_sm$active_hours |> 
  slice(f_local_maxima(estimate))
l_ric_LU_500$l_pred_sm$active_hours |> 
  slice(f_local_maxima(estimate))
l_mass_LU_500$l_pred_sm$active_hours |> 
  slice(f_local_maxima(estimate))

l_ric_LU_500$l_pred_sm$height |> 
  slice(f_local_maxima(estimate))

# factor estimates for a 5 degree Celsius difference:
f_pred_fixed_diff_specific(l_abu_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                           formula = formula_full,
                           var_i = "T_2day", diff_or = 5) |> 
  mutate(response = "abu_tot") |> 
  bind_rows(f_pred_fixed_diff_specific(l_ric_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                                       formula = formula_full,
                                       var_i = "T_2day", diff_or = 5) |> 
              mutate(response = "sric")) |> 
  bind_rows(f_pred_fixed_diff_specific(l_mass_LU_500$fit, d_mod_z, d_scalings |> filter(data == "full"),
                                       formula = formula_full,
                                       var_i = "T_2day", diff_or = 5) |> 
              mutate(response = "mass_tot")) |> 
  mutate(smry = paste0(round(exp(estimate_mean), 2), " (95%-CI: ",
                       round(exp_lower_95, 2), "–",
                       round(exp_upper_95, 2), ")")) |> 
  select(response, var, diff_or, smry)

# ... Rhat #####################################################################
################################################################################.

f_apply_Rhat(l_abu_LU_500$fit) |> 
  mutate(response = "Abundance") |> 
  bind_rows(f_apply_Rhat(l_ric_LU_500$fit) |> 
              mutate(response = "Richness"))|> 
  bind_rows(f_apply_Rhat(l_mass_LU_500$fit) |> 
              mutate(response = "Biomass")) |> 
  mutate(model = "full") |> 
  bind_rows(f_apply_Rhat(l_abu_LU_egg$fit) |> 
              mutate(response = "Abundance") |> 
              bind_rows(f_apply_Rhat(l_abu_LU_larva$fit) |> 
                          mutate(response = "Abundance")) |> 
              bind_rows(f_apply_Rhat(l_abu_LU_pupa$fit) |> 
                          mutate(response = "Abundance")) |> 
              bind_rows(f_apply_Rhat(l_abu_LU_adult$fit) |> 
                          mutate(response = "Abundance")) |> 
              bind_rows(f_apply_Rhat(l_ric_LU_egg$fit) |> 
                          mutate(response = "Richness")) |> 
              bind_rows(f_apply_Rhat(l_ric_LU_larva$fit) |> 
                          mutate(response = "Richness")) |> 
              bind_rows(f_apply_Rhat(l_ric_LU_pupa$fit) |> 
                          mutate(response = "Richness")) |> 
              bind_rows(f_apply_Rhat(l_ric_LU_adult$fit) |> 
                          mutate(response = "Richness")) |> 
              bind_rows(f_apply_Rhat(l_mass_LU_egg$fit) |> 
                          mutate(response = "Biomass")) |> 
              bind_rows(f_apply_Rhat(l_mass_LU_larva$fit) |> 
                          mutate(response = "Biomass")) |> 
              bind_rows(f_apply_Rhat(l_mass_LU_pupa$fit) |> 
                          mutate(response = "Biomass")) |> 
              bind_rows(f_apply_Rhat(l_mass_LU_adult$fit) |> 
                          mutate(response = "Biomass")) |> 
              mutate(model = "overwintering_stages")) |> 
  group_by(response, model) |> 
  summarise(prop_thresh = mean(rhat < 1.1)) 

# ... Table 1 ##################################################################
################################################################################.

f_pred_fixed_diff(l_abu_LU_500$fit, d_mod_z, d_scalings |>
                    filter(data == "full"),
                  formula = formula_full) |> 
  bind_rows(f_pred_smooths_diff(l_abu_LU_500$fit, d_mod_z, d_scalings |> 
                                  filter(data == "full"),
                                formula = formula_full)) |> 
  mutate(response = "Abundance") |> 
  bind_rows(f_pred_fixed_diff(l_ric_LU_500$fit, d_mod_z, d_scalings |>
                                filter(data == "full"),
                              formula = formula_full) |> 
              bind_rows(f_pred_smooths_diff(l_ric_LU_500$fit, d_mod_z, d_scalings |> 
                                              filter(data == "full"),
                                            formula = formula_full)) |> 
              mutate(response = "Richness")) |> 
  bind_rows(f_pred_fixed_diff(l_mass_LU_500$fit, d_mod_z, d_scalings |>
                                filter(data == "full"),
                              formula = formula_full) |> 
              bind_rows(f_pred_smooths_diff(l_mass_LU_500$fit, d_mod_z, d_scalings |> 
                                              filter(data == "full"),
                                            formula = formula_full)) |> 
              mutate(response = "Biomass (g)")) |>
  filter(var %in% c('prop_forest_500', 'prop_grassland_500', 'prop_crop_500',
                    'prop_sealed_500', 'P_2day', 'T_2day', 'height')) |> 
  mutate(exp_estimate_mean = exp(estimate_mean),
         diff_or = round(diff_or, 3),
         var = recode(var, !!! v_covlabels),
         var = factor(var, levels = v_covlabels)) |> 
  arrange(response, var) |> 
  mutate(var = as.character(var)) |> 
  select(response, var, diff_or, exp_estimate_mean, exp_lower_95, exp_upper_95) |> 
  mutate_at(vars(exp_estimate_mean, exp_lower_95, exp_upper_95),
            ~ round(., 2)) |> 
  write_csv("Output/Tables/factors_summary.csv")

# ... Figure 1 #################################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_plots_cov1_abu <- c(l_abu_LU_500$l_pred_fe, l_abu_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("height", "P_2day", "T_2day", "yday")) %>%
  map(f_plot_pred, data = d_mod, response = "abu_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

p <-
  l_plots_cov1_abu$yday +
  geom_line(
    data = l_abu_LU_egg$l_pred_sm$yday %>% 
      mutate(hib = "e") %>% 
      bind_rows(l_abu_LU_larva$l_pred_sm$yday %>% 
                  mutate(hib = "l")) %>% 
      bind_rows(l_abu_LU_pupa$l_pred_sm$yday %>% 
                  mutate(hib = "p")) %>% 
      bind_rows(l_abu_LU_adult$l_pred_sm$yday %>% 
                  mutate(hib = "a")) %>% 
      mutate(hib = factor(hib, levels = c("e", "l", "p", "a"))),
    aes(colour = hib, y = estimate),
    size = .5) +
  scale_colour_manual(values = c("#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"),
                      name = "Hib. stage") +
  theme(legend.position = c(0.5, 1),
        legend.justification = c(-.1, .9),
        legend.background = element_rect(colour = NA, 
                                         fill = alpha("white", .5)),
        legend.key.height = unit(.25, "cm"),
        legend.title = element_text(size = v_textsize["axis.text"]),
        legend.text = element_text(size = v_textsize["additional.text"]))
p$layers <- p$layers[c(1,
                       length(p$layers), 
                       seq(2, length(p$layers) - 1))]
l_plots_cov1_abu$yday <- p


l_plots_cov1_abu <- lapply(l_plots_cov1_abu,
                           \(x) x +
                             theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                   axis.text = element_text(size = v_textsize["axis.text"])) +
                             labs(y = "Abundance")+ 
                             scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p"))

l_plots_cov1_abu[c("P_2day", "T_2day", "yday")] <- 
  lapply(l_plots_cov1_abu[c("P_2day", "T_2day", "yday")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_cov1_abu <- 
  lapply(l_plots_cov1_abu[c("height", "P_2day", "T_2day", "yday")],
         \(x) x + 
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank()))

grid_abu <- plot_grid(plotlist = l_plots_cov1_abu,
                      nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                      align = "h")

# richness ---------------------------------------------------------------------.
l_plots_cov1_ric <- c(l_ric_LU_500$l_pred_fe, l_ric_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("height", "P_2day", "T_2day", "yday")) %>%
  map(f_plot_pred, data = d_mod, response = "sric", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

p <-
  l_plots_cov1_ric$yday +
  geom_line(
    data = l_ric_LU_egg$l_pred_sm$yday %>% 
      mutate(hib = "e") %>% 
      bind_rows(l_ric_LU_larva$l_pred_sm$yday %>% 
                  mutate(hib = "l")) %>% 
      bind_rows(l_ric_LU_pupa$l_pred_sm$yday %>% 
                  mutate(hib = "p")) %>% 
      bind_rows(l_ric_LU_adult$l_pred_sm$yday %>% 
                  mutate(hib = "a")) %>% 
      mutate(hib = factor(hib, levels = c("e", "l", "p", "a"))),
    aes(colour = hib, y = estimate),
    size = .5) +
  scale_colour_manual(values = c("#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"),
                      name = "Hib. stage", guide = "none")
p$layers <- p$layers[c(1,
                       length(p$layers), 
                       seq(2, length(p$layers) - 1))]
l_plots_cov1_ric$yday <- p


l_plots_cov1_ric <- lapply(l_plots_cov1_ric,
                           \(x) x +
                             theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                   axis.text = element_text(size = v_textsize["axis.text"])) +
                             labs(y = "Richness")+ 
                             scale_y_continuous(breaks = c(0, 10, 100), 
                                                labels = c("0", "10", "    100"), # for layout of the final plot
                                                trans = "log1p"))

l_plots_cov1_ric[c("P_2day", "T_2day", "yday")] <- 
  lapply(l_plots_cov1_ric[c("P_2day", "T_2day", "yday")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_cov1_ric <- 
  lapply(l_plots_cov1_ric[c("height", "P_2day", "T_2day", "yday")],
         \(x) x + 
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank()))

grid_ric <- plot_grid(plotlist = l_plots_cov1_ric,
                      nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                      align = "h")

# biomass ----------------------------------------------------------------------.
l_plots_cov1_mass <- c(l_mass_LU_500$l_pred_fe, l_mass_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("height", "P_2day", "T_2day", "yday")) %>%
  map(f_plot_pred, data = d_mod, response = "mass_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

p <-
  l_plots_cov1_mass$yday +
  geom_line(
    data = l_mass_LU_egg$l_pred_sm$yday %>% 
      mutate(hib = "e") %>% 
      bind_rows(l_mass_LU_larva$l_pred_sm$yday %>% 
                  mutate(hib = "l")) %>% 
      bind_rows(l_mass_LU_pupa$l_pred_sm$yday %>% 
                  mutate(hib = "p")) %>%
      bind_rows(l_mass_LU_adult$l_pred_sm$yday %>% 
                  mutate(hib = "a")) %>% 
      mutate(hib = factor(hib, levels = c("e", "l", "p", "a"))),
    aes(colour = hib, y = estimate),
    size = .5) +
  scale_colour_manual(values = c("#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"),
                      name = "Hib. stage", guide = F) 
p$layers <- p$layers[c(1,
                       length(p$layers), 
                       seq(2, length(p$layers) - 1))]
l_plots_cov1_mass$yday <- p


l_plots_cov1_mass <- lapply(l_plots_cov1_mass,
                            \(x) x +
                              theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                    axis.text = element_text(size = v_textsize["axis.text"])) +
                              labs(y = "Biomass (g)")+ 
                              scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                                 labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                                 trans = log_plus_trans))

l_plots_cov1_mass[c("P_2day", "T_2day", "yday")] <- 
  lapply(l_plots_cov1_mass[c("P_2day", "T_2day", "yday")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

# change order
l_plots_cov1_mass <- l_plots_cov1_mass[c("height", "P_2day", "T_2day", "yday")]

grid_mass <- plot_grid(plotlist = l_plots_cov1_mass,
                       nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                       align = "h")



# All in one plot --------------------------------------------------------------.
p <- plot_grid(grid_abu,
               grid_ric,
               grid_mass,
               ncol = 1, rel_heights = c(1, 1, 1.2))

ggsave(p, file = "Output/Figures/Covariates_main_comb.pdf",
       width = 170, height = 180, units = "mm", dpi = 600)

# ... Figure 2 #################################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_plots_LU_abu <- c(l_abu_LU_500$l_pred_fe, l_abu_LU_500$l_pred_sm) %>% 
  keep(grepl("prop_", names(.))) %>%
  map(f_plot_pred, data = d_mod, response = "abu_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

l_plots_LU_abu <- lapply(l_plots_LU_abu,
                         \(x) x +
                           theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                 axis.text = element_text(size = v_textsize["axis.text"])) +
                           labs(y = "Abundance")+ 
                           scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p"))

l_plots_LU_abu[2:4] <- 
  lapply(l_plots_LU_abu[2:4],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_LU_abu <- 
  lapply(l_plots_LU_abu,
         \(x) x + 
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank()))

grid_abu <- plot_grid(plotlist = l_plots_LU_abu,
                      nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                      align = "h")

# richness ---------------------------------------------------------------------.
l_plots_LU_ric <- c(l_ric_LU_500$l_pred_fe, l_ric_LU_500$l_pred_sm) %>% 
  keep(grepl("prop_", names(.))) %>%
  map(f_plot_pred, data = d_mod, response = "sric", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))


l_plots_LU_ric <- lapply(l_plots_LU_ric,
                         \(x) x +
                           theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                 axis.text = element_text(size = v_textsize["axis.text"])) +
                           labs(y = "Richness")+ 
                           scale_y_continuous(breaks = c(0, 10, 100), 
                                              labels = c("0", "10", "    100"), # for layout of the final plot
                                              trans = "log1p"))

l_plots_LU_ric[2:4] <- 
  lapply(l_plots_LU_ric[2:4],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_LU_ric <- 
  lapply(l_plots_LU_ric,
         \(x) x + 
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank()))

grid_ric <- plot_grid(plotlist = l_plots_LU_ric,
                      nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                      align = "h")

# biomass ----------------------------------------------------------------------.
l_plots_LU_mass <- c(l_mass_LU_500$l_pred_fe, l_mass_LU_500$l_pred_sm) %>% 
  keep(grepl("prop_", names(.))) %>%
  map(f_plot_pred, data = d_mod, response = "mass_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

l_plots_LU_mass <- lapply(l_plots_LU_mass,
                          \(x) x +
                            theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                  axis.text = element_text(size = v_textsize["axis.text"])) +
                            labs(y = "Biomass (g)")+ 
                            scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                               labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                               trans = log_plus_trans))

l_plots_LU_mass[2:4] <- 
  lapply(l_plots_LU_mass[2:4],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

grid_mass <- plot_grid(plotlist = l_plots_LU_mass,
                       nrow = 1, rel_widths = c(1.3, 1, 1, 1),
                       align = "h")


# All in one plot
p <- plot_grid(grid_abu,
               grid_ric,
               grid_mass,
               ncol = 1, rel_heights = c(1, 1, 1.2))

ggsave(p, file = "Output/Figures/LU_comb.pdf",
       width = 170, height = 180, units = "mm", dpi = 600)

# ... Figure S7 ################################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_plots_cov2_abu <- c(l_abu_LU_500$l_pred_fe, l_abu_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("traptype", "bulbtype", "n_trap", "active_hours", "sample_previous")) %>%
  map(f_plot_pred, data = d_mod, response = "abu_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

l_plots_cov2_abu <- lapply(l_plots_cov2_abu,
                           \(x) x +
                             theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                   axis.text = element_text(size = v_textsize["axis.text"])) +
                             labs(y = "Abundance") + 
                             scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p"))

l_plots_cov2_abu[c("bulbtype", "n_trap", "active_hours", "sample_previous")] <- 
  lapply(l_plots_cov2_abu[c("bulbtype", "n_trap", "active_hours", "sample_previous")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_cov2_abu <- lapply(l_plots_cov2_abu,
                           \(x) x + 
                             theme(axis.title.x = element_blank(),
                                   axis.text.x = element_blank()))

grid_abu_sup <- plot_grid(plotlist = l_plots_cov2_abu[c("traptype", "bulbtype",
                                                        "n_trap", "active_hours", 
                                                        "sample_previous")],
                          nrow = 1, rel_widths = c(1.2, 1, 1, 1, 1),
                          align = "h")

# richness ---------------------------------------------------------------------.
l_plots_cov2_ric <- c(l_ric_LU_500$l_pred_fe, l_ric_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("traptype", "bulbtype", "n_trap", "active_hours", "sample_previous")) %>%
  map(f_plot_pred, data = d_mod, response = "sric", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

l_plots_cov2_ric <- lapply(l_plots_cov2_ric,
                           \(x) x +
                             theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                   axis.text = element_text(size = v_textsize["axis.text"])) +
                             labs(y = "Richness") + 
                             scale_y_continuous(breaks = c(0, 10, 100), 
                                                labels = c("0", "10", "    100"), # for layout of the final plot
                                                trans = "log1p"))

l_plots_cov2_ric[c("bulbtype", "n_trap", "active_hours", "sample_previous")] <- 
  lapply(l_plots_cov2_ric[c("bulbtype", "n_trap", "active_hours", "sample_previous")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

l_plots_cov2_ric <- lapply(l_plots_cov2_ric,
                           \(x) x + 
                             theme(axis.title.x = element_blank(),
                                   axis.text.x = element_blank()))


grid_ric_sup <- plot_grid(plotlist = l_plots_cov2_ric[c("traptype", "bulbtype",
                                                        "n_trap", "active_hours", 
                                                        "sample_previous")],
                          nrow = 1, rel_widths = c(1.2, 1, 1, 1, 1),
                          align = "h")

# biomass ----------------------------------------------------------------------.
l_plots_cov2_mass <- c(l_mass_LU_500$l_pred_fe, l_mass_LU_500$l_pred_sm) %>% 
  keep(names(.) %in% c("traptype", "bulbtype", "n_trap", "active_hours", "sample_previous")) %>%
  map(f_plot_pred, data = d_mod, response = "mass_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))

l_plots_cov2_mass <- lapply(l_plots_cov2_mass,
                            \(x) x +
                              theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                    axis.text = element_text(size = v_textsize["axis.text"])) +
                              labs(y = "Biomass (g)") + 
                              scale_y_continuous(breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                                 labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
                                                 trans = log_plus_trans))

l_plots_cov2_mass[c("bulbtype", "n_trap", "active_hours", "sample_previous")] <- 
  lapply(l_plots_cov2_mass[c("bulbtype", "n_trap", "active_hours", "sample_previous")],
         \(x) x + 
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))


grid_mass_sup <- plot_grid(plotlist = l_plots_cov2_mass[c("traptype", "bulbtype",
                                                          "n_trap", "active_hours", 
                                                          "sample_previous")],
                           nrow = 1, rel_widths = c(1.2, 1, 1, 1, 1),
                           align = "h")


# All in one plot --------------------------------------------------------------.
p <- plot_grid(grid_abu_sup,
               grid_ric_sup,
               grid_mass_sup,
               ncol = 1, rel_heights = c(1, 1, 1.5))

ggsave(p, file = "Output/Figures/Covariates_sup_comb.pdf",
       width = 180, height = 210, units = "mm", dpi = 600)

# ... Figure S6 ################################################################
################################################################################.

# run models and also export random effects

l_abu_LU_500_R <- f_analysis_LU(formula = abu_tot ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 923,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              extract_random = T)

l_ric_LU_500_R <- f_analysis_LU(formula = sric ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "zero_inflated_negbinomial",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 87127,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              extract_random = T)

l_mass_LU_500_R <- f_analysis_LU(formula = mass_tot ~
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
                              data = d_mod,
                              scalings = filter(d_scalings, data == "full"),
                              family = "hurdle_gamma",
                              hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate),
                              iter = n_iter, seed = 1411,
                              LU_vars = c("Forest" = "prop_forest_500",
                                          "Grassland" = "prop_grassland_500",
                                          "Crop" = "prop_crop_500",
                                          "Sealed" = "prop_sealed_500"),
                              extract_random = T)

# caculate mus and ysims
LU_vars <- c("Forest" = "prop_forest_500",
             "Grassland" = "prop_grassland_500",
             "Crop" = "prop_crop_500",
             "Sealed" = "prop_sealed_500")
formula_full <- response ~
  s(yday) + P_2day + T_2day +
  C(traptype, "contr.sum") +
  C(bulbtype, "contr.sum") +
  n_trap +
  C(sample_previous, "contr.sum") +
  (1 | gr) +
  (1 | trap_ID) +
  (1 | night_ID) +
  (1 | trap_ID_A)
formula_full <- update(formula_full,
                       paste0(". ~ s(height) + ",
                              paste(LU_vars, collapse = " + "),
                              "+ ."))

set.seed(23)

# abundance
m_mu <- f_mu_s2p1_r4(formula = formula_full,
                     fit = l_abu_LU_500_R$fit,
                     data = d_mod_z,
                     hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))


m_ysim <- c()
sel <- sample(seq_len(ncol(m_mu)), 100)
for (sel_i in sel){
  out <- c()
  for (y_i in seq_len(nrow(m_mu))){
    out <- c(out, rnbinom(n = 1,
                          size = l_abu_LU_500_R$fit$shape[sel_i],
                          mu = exp(m_mu[y_i, sel_i])) *
               rbinom(n = 1, size = 1, prob = 1 - l_abu_LU_500_R$fit$zi[sel_i]))
  }
  m_ysim <- rbind(m_ysim, out)
}


d_check_abu <- data.frame(y = d_mod_z$abu_tot,
                          rep_id = 99999,
                          type = "Empirical") |>
  bind_rows(data.frame(y = c(t(m_ysim)),
                       rep_id = rep(sel, each = ncol(m_ysim)),
                       type = "Predictive"))

# richness
m_mu <- f_mu_s2p1_r4(formula = update(formula_full, sric ~ .),
                     fit = l_ric_LU_500_R$fit,
                     data = d_mod_z,
                     hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))


m_ysim <- c()
sel <- sample(seq_len(ncol(m_mu)), 100)
for (sel_i in sel){
  out <- c()
  for (y_i in seq_len(nrow(m_mu))){
    out <- c(out, rnbinom(n = 1,
                          size = l_ric_LU_500_R$fit$shape[sel_i],
                          mu = exp(m_mu[y_i, sel_i])) *
               rbinom(n = 1, size = 1, prob = 1 - l_ric_LU_500_R$fit$zi[sel_i]))
  }
  m_ysim <- rbind(m_ysim, out)
}


d_check_ric <- data.frame(y = d_mod_z$sric,
                          rep_id = 99999,
                          type = "Empirical") |>
  bind_rows(data.frame(y = c(t(m_ysim)),
                       rep_id = rep(sel, each = ncol(m_ysim)),
                       type = "Predictive"))

# biomass
m_mu <- f_mu_s2p1_r4(formula = update(formula_full, mass_tot ~ .),
                     fit = l_mass_LU_500_R$fit,
                     data = d_mod_z,
                     hours_sel = which(d_mod_z$traptype == "p" & !d_mod_z$estimate))


m_ysim <- c()
sel <- sample(seq_len(ncol(m_mu)), 100)
for (sel_i in sel){
  out <- c()
  for (y_i in seq_len(nrow(m_mu))){
    out <- c(out, rgamma(n = 1,
                         shape = l_mass_LU_500_R$fit$shape[sel_i],
                          rate = l_mass_LU_500_R$fit$shape[sel_i] * exp(-m_mu[y_i, sel_i])) * 
               rbinom(n = 1, size = 1, prob = 1 - l_mass_LU_500_R$fit$hu[sel_i]))
  }
  m_ysim <- rbind(m_ysim, out)
}


d_check_mass <- data.frame(y = d_mod_z$mass_tot,
                          rep_id = 99999,
                          type = "Empirical") |>
  bind_rows(data.frame(y = c(t(m_ysim)),
                       rep_id = rep(sel, each = ncol(m_ysim)),
                       type = "Predictive"))

# plotting ---------------------------------------------------------------------.

# abundance:
p_check_abu <-  d_check_abu %>% 
  ggplot(aes(x = y, y = ..density..)) +
  geom_histogram(data = \(x) filter(x, type == "Empirical"), fill = "grey50", bins = 30) +
  stat_density(aes(group = rep_id, col = type, 
                   size = type, alpha = type), 
               geom = "line", position = "identity") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = c("turquoise3", "darkorchid4"), 
                      labels = c("Empirical data", "Posterior predictive"),
                      name = "Distribution") +
  scale_alpha_manual(values = c(1, .2), guide = F) +
  scale_size_manual(values = c(1, .5), guide = F) +
  xlab("Abundance") +
  ylab("Density") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# richness:
p_check_ric <- d_check_ric %>% 
  ggplot(aes(x = y, y = ..density..)) +
  geom_histogram(data = \(x) filter(x, type == "Empirical"), fill = "grey50", bins = 30) +
  stat_density(aes(group = rep_id, col = type, 
                   size = type, alpha = type), 
               geom = "line", position = "identity") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = c("turquoise3", "darkorchid4"), 
                      labels = c("Empirical data", "Posterior predictive"),
                      name = "Distribution") +
  scale_alpha_manual(values = c(1, .2), guide = F) +
  scale_size_manual(values = c(1, .5), guide = F) +
  xlab("Richness") +
  ylab("Density") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# biomass:
p_check_mass <- d_check_mass %>% 
  ggplot(aes(x = y, y = ..density..)) +
  geom_histogram(data = \(x) filter(x, type == "Empirical"), fill = "grey50", bins = 30) +
  stat_density(aes(group = rep_id, col = type, 
                   size = type, alpha = type), 
               geom = "line", position = "identity") +
  scale_x_continuous(trans = log_plus_trans, 
                     breaks = c(0, 0.00001, 0.0001, 0.001, .01, .1, 1, 10, 100, 1000),
                     labels = c(0, 0.00001, 0.0001, 0.001, .01, .1, 1, 10, 100, 1000)) +
  scale_colour_manual(values = c("turquoise3", "darkorchid4"), 
                      labels = c("Empirical data", "Posterior predictive"),
                      name = "Distribution") +
  scale_alpha_manual(values = c(1, .2), guide = "none") +
  scale_size_manual(values = c(1, .5), guide = "none") +
  xlab("Biomass [g]") +
  ylab("Density") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# all in one plot:
plot_grid(p_check_abu, p_check_ric, p_check_mass, nrow = 1, 
          rel_widths = c(1, 1, 1.7), labels = letters[1:3], align = "h")

ggsave("Output/Figures/Dist_emp_pred_LU.jpeg", width = 240, height = 90,
       units = "mm", dpi = 400)

# ... Table S3 #################################################################
################################################################################.

formula <- response ~ 
  s(yday) + P_2day + T_2day +
  C(traptype, "contr.sum") +
  C(bulbtype, "contr.sum") +
  n_trap +
  C(sample_previous, "contr.sum") +
  (1 | gr) +
  (1 | trap_ID) +
  (1 | night_ID) +
  (1 | trap_ID_A)

f_summarytable(l_abu_LU_500$fit, formula, v_covlabels_short,
               data = d_mod_z,
               LU_vars = c("Forest" = "prop_forest_500",
                           "Grassland" = "prop_grassland_500",
                           "Crop" = "prop_crop_500",
                           "Sealed" = "prop_sealed_500")) |> 
  print()
f_summarytable(l_ric_LU_500$fit, formula, v_covlabels_short,
               data = d_mod_z,
               LU_vars = c("Forest" = "prop_forest_500",
                           "Grassland" = "prop_grassland_500",
                           "Crop" = "prop_crop_500",
                           "Sealed" = "prop_sealed_500")) |> 
  print()
f_summarytable(l_mass_LU_500$fit, formula, v_covlabels_short,
               data = d_mod_z,
               LU_vars = c("Forest" = "prop_forest_500",
                           "Grassland" = "prop_grassland_500",
                           "Crop" = "prop_crop_500",
                           "Sealed" = "prop_sealed_500")) |> 
  print()
