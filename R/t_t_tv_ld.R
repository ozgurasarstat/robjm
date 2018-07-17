t_t_tv_ld = "

data{
  int<lower = 1> ntot;        // total number observations
  int id[ntot];               // id matrix, takes values like 1, 2, 3, ...
  vector[ntot] y;             // response matrix
  int<lower = 1> p;           // number of covariates in the fixed effects design matrix
  int<lower = 1> q;           // number of covariates in the random effects design matrix
  int<lower = 1> ngroup;      // number of subjects/clusters/groups
  matrix[ntot, p] x;          // fixed effects design matrix
  matrix[ntot, q * ngroup] d; // random effects design matrix, block diagonal
  int<lower = 1> s;           // number of columns of a matrix
  matrix[ntot, s] a;          // spline matrix
  vector[5] priors;           // prior hyperparameters, order: theta, omega, sigma_bstar, sigmaZ, beta
}

transformed data{
  //QR decomposition
  matrix[ntot, p] Q_star;
  matrix[p, p] R_star;
  matrix[p, p] R_star_inv;

  Q_star = qr_Q(x)[, 1:p] * sqrt(ntot - 1);
  R_star = qr_R(x)[1:p, ] / sqrt(ntot - 1);
  R_star_inv = inverse(R_star);
}

parameters{
  vector[p] theta;              // fixed effects coefficients
  matrix[ngroup, q] Bstar;
  corr_matrix[q] Omega;
  vector<lower = 0>[q] sigma_Bstar;
  vector<lower = 0>[ngroup] V;
  real<lower = 0> sigma_Zstar;      // scale parameter of measurement error
  vector<lower = 0>[ntot] W;
  vector[s] beta;
  real<lower = 0, upper = 100> delta0;
  real<lower = 2, upper = 100> phi;
}

transformed parameters{
  cov_matrix[q] Sigma;
  vector[ntot] linpred;
  matrix[ngroup, q] B;
  matrix[ngroup * q, 1] B_mat;
  vector<lower = 0>[ntot] delta;
  vector[q] zero_Bstar = rep_vector(0, q);
  delta = to_vector((2 + delta0) * exp(a * beta));

  for(i in 1:ngroup){
    B[i, ] = Bstar[i, ] * sqrt(V[i]);
  }
  B_mat = to_matrix(B', ngroup * q, 1);
  linpred = Q_star * theta + to_vector(d * B_mat);
  Sigma = quad_form_diag(Omega, sigma_Bstar);
}

model{

  theta ~ cauchy(0, priors[1]);

  for(i in 1:ngroup)
    Bstar[i] ~ multi_normal(zero_Bstar, Sigma);


  Omega ~ lkj_corr(priors[2]);
  sigma_Bstar ~ cauchy(0, priors[3]);
  sigma_Zstar ~ cauchy(0, priors[4]);

  V ~ inv_gamma(phi/2, phi/2 -1);
  beta ~ cauchy(0, priors[5]);


  for(i in 1:ntot){
    W[i] ~ inv_gamma(delta[i]/2, delta[i]/2 - 1);
    y[i] ~ normal(linpred[i], sigma_Zstar * sqrt(W[i]));
  }
}

generated quantities{
  vector[p] alpha;
  real sigmasq;
  
  alpha = R_star_inv * theta;
  sigmasq = sigma_Zstar^2;
}

"
