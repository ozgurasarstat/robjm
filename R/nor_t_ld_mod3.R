nor_t_ld_mod3 = "

data{
int<lower = 1> ntot;        // total number observations
int<lower = 1> id[ntot];    // id matrix, takes values like 1, 2, 3, ...
vector[ntot] y;             // response matrix
int<lower = 1> p;           // number of covariates in the fixed effects design matrix
int<lower = 1> q;           // number of covariates in the random effects design matrix
int<lower = 1> ngroup;      // number of subjects/clusters/groups
matrix[ntot, p] x;          // fixed effects design matrix
matrix[ntot, q * ngroup] d; // random effects design matrix, block diagonal
vector[4] priors;           // prior hyperparameters, order: theta, omega, sigma_b, sigma_Z
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
vector[p] theta;                           // fixed effects coefficients with QR decomposition
matrix[ngroup, q] B;                       // B matrix
corr_matrix[q] Omega;                      // correlation matrix for B
vector<lower = 0>[q] sigma_B;              // scale parameters for B
real<lower = 0> sigma_Zstar;               // scale parameter of measurement error
vector<lower = 0>[ntot] W;                 // scaling r.v. for Z
real<lower = 0.01, upper = 0.5> delta_inv; // inverse of the parameter for W
}

transformed parameters{
cov_matrix[q] Sigma;
vector[ntot] linpred;
matrix[ngroup * q, 1] B_mat;
vector[q] zero_B = rep_vector(0, q);
real<lower = 2, upper = 100> delta;

delta = 1/delta_inv;

B_mat = to_matrix(B', ngroup * q, 1);

linpred = Q_star * theta + to_vector(d * B_mat);

Sigma = quad_form_diag(Omega, sigma_B);

}

model{

theta ~ cauchy(0, priors[1]);

for(i in 1:ngroup){
B[i] ~ multi_normal(zero_B, Sigma);
}

Omega ~ lkj_corr(priors[2]);
sigma_B ~ cauchy(0, priors[3]);
sigma_Zstar ~ cauchy(0, priors[4]);

W ~ inv_gamma(delta/2, delta/2);
//delta_inv ~ uniform(0.01, 0.5);//the prior is uniform with -infty and infty, constrained above

for(i in 1:ntot) y[i] ~ normal(linpred[i], sigma_Zstar * sqrt(W[i]));

}

generated quantities{
vector[p] alpha;
real sigmasq;

alpha = R_star_inv * theta; // convert estimates based on QR to original formulation
sigmasq = sigma_Zstar^2;        // report sigma^2 for Z
}
"
