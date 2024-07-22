functions {
  /* hurdle gamma log-PDF of a single response
  * Args:
  *   y: the response value
  *   alpha: shape parameter of the gamma distribution
  *   beta: rate parameter of the gamma distribution
  *   hu: hurdle probability
  * Returns:
  *   a scalar to be added to the log posterior
  */
  real hurdle_gamma_lpdf(real y, real alpha, real beta, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
      gamma_lpdf(y | alpha, beta);
    }
  }
  /* hurdle gamma log-PDF of a single response
  * logit parameterization of the hurdle part
  * Args:
  *   y: the response value
  *   alpha: shape parameter of the gamma distribution
  *   beta: rate parameter of the gamma distribution
  *   hu: linear predictor for the hurdle part
  * Returns:
  *   a scalar to be added to the log posterior
  */
  real hurdle_gamma_logit_lpdf(real y, real alpha, real beta, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
      gamma_lpdf(y | alpha, beta);
    }
  }
  // hurdle gamma log-CCDF and log-CDF functions
  real hurdle_gamma_lccdf(real y, real alpha, real beta, real hu) {
    return bernoulli_lpmf(0 | hu) + gamma_lccdf(y | alpha, beta);
  }
  real hurdle_gamma_lcdf(real y, real alpha, real beta, real hu) {
    return log1m_exp(hurdle_gamma_lccdf(y | alpha, beta, hu));
  }
  
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_SEL;
  vector[N] Y;  // response variable
  int SEL[N_SEL]; 
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  matrix[N_SEL, 1] Xs_add;
  
  // data for spline 1
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  
  // data for spline 2
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  
  // data for additional spline
  int nb_add;  // number of bases
  int knots_add[nb_add];  // number of knots
  // basis function matrices
  matrix[N_SEL, knots_add[1]] Zs_add_1;
  
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  vector[1] bs_add;  // spline coefficients
  // parameters for spline 1
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline 2
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  
  // parameters for additional spline
  // standarized spline coefficients
  vector[knots_add[1]] zs_add_1;
  real<lower=0> sds_add_1;  // standard deviations of spline coefficients
  
  real<lower=0> shape;  // shape parameter
  real<lower=0,upper=1> hu;  // hurdle probability
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // actual spline coefficients
  vector[knots_add[1]] s_add_1;
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  vector[N_3] r_3_1;  // actual group-level effects
  vector[N_4] r_4_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual spline coefficients
  s_add_1 = sds_add_1 * zs_add_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_3_1 = (sd_3[1] * (z_3[1]));
  r_4_1 = (sd_4[1] * (z_4[1]));
  lprior += normal_lpdf(b | 0, 5);
  lprior += student_t_lpdf(Intercept | 3, 3, 2.5);
  lprior += normal_lpdf(bs | 0, 5);
  lprior += normal_lpdf(bs_add | 0, 5);
  lprior += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_2_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_add_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(shape | 0.01, 0.01);
  lprior += beta_lpdf(hu | 1, 1);
  lprior += cauchy_lpdf(sd_1 | 0, 1)
  - 1 * cauchy_lccdf(0 | 0, 1);
  lprior += cauchy_lpdf(sd_2 | 0, 1)
  - 1 * cauchy_lccdf(0 | 0, 1);
  lprior += cauchy_lpdf(sd_3 | 0, 1)
  - 1 * cauchy_lccdf(0 | 0, 1);
  lprior += cauchy_lpdf(sd_4 | 0, 1)
  - 1 * cauchy_lccdf(0 | 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
    for (m in 1:N_SEL){
      mu[SEL[m]] += Xs_add[m, 1] * bs_add[1] + Zs_add_1[m, ] * s_add_1;
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = shape * exp(-(mu[n]));
    }
    for (n in 1:N) {
      target += hurdle_gamma_lpdf(Y[n] | shape, mu[n], hu);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_add_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  vector[N] log_lik;
  
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  
  vector[N] mu = Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
  for (m in 1:N_SEL){
    mu[SEL[m]] += Xs_add[m, 1] * bs_add[1] + Zs_add_1[m, ] * s_add_1;
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
  
  for (n in 1:N) log_lik[n] = hurdle_gamma_lpdf(Y[n] | shape, mu[n], hu);
}