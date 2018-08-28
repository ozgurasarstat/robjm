nor_nor_ld = "

data{
int<lower = 1> ntot;        // total number observations
int id[ntot];               // id matrix, takes values like 1, 2, 3, ...
vector[ntot] y;             // response matrix
int<lower = 1> p;           // number of covariates in the fixed effects design matrix
int<lower = 1> q;           // number of covariates in the random effects design matrix
int<lower = 1> ngroup;      // number of subjects/clusters/groups
matrix[ntot, p] x;          // fixed effects design matrix
matrix[ntot, q * ngroup] d; // random effects design matrix, block diagonal
vector[4] priors; // prior hyperparameters, order: theta, Omega, sigma_B, sigma_Z
}

transformed data{
matrix[ntot, p] Q_star;
matrix[p, p] R_star;
matrix[p, p] R_star_inv;

Q_star = qr_Q(x)[, 1:p] * sqrt(ntot - 1.0);
R_star = qr_R(x)[1:p, ] / sqrt(ntot - 1.0);
R_star_inv = inverse(R_star);
}

parameters{
vector[p] theta;              // fixed effects coefficients - qr
matrix[ngroup, q] B;          // random effects coefficients
corr_matrix[q] Omega;             // correlation matrix for random effects
vector<lower = 0>[q] sigma_B; // scale parameters for random effects
real<lower = 0> sigma_Z;      // scale parameter of measurement error
}

transformed parameters{
cov_matrix[q] Sigma;
vector[ntot] linpred;
matrix[ngroup * q, 1] Bmat;
vector[q] zero_B = rep_vector(0, q);

Bmat = to_matrix(B', ngroup * q, 1);

linpred = Q_star * theta + to_vector(d * Bmat);
Sigma = quad_form_diag(Omega, sigma_B);
}

model{

theta ~ cauchy(0, priors[1]);

for(i in 1:ngroup){
B[i] ~ multi_normal(zero_B, Sigma);
}

Omega ~ lkj_corr(priors[2]);
sigma_B ~ cauchy(0, priors[3]);
sigma_Z ~ cauchy(0, priors[4]);

y ~ normal(linpred, sigma_Z);

}

generated quantities{
vector[p] alpha;
real sigmasq;
alpha = R_star_inv * theta;
sigmasq = sigma_Z^2;
}

"
