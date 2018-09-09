t_t_tv_jm = "

data{
//longitudinal sub-model data
int<lower = 1> ntot;        // total number observations
int id[ntot];               // id matrix, takes values like 1, 2, 3, ...
vector[ntot] y;             // response matrix
int<lower = 1> p;           // number of covariates in the fixed effects design matrix
int<lower = 1> q;           // number of covariates in the random effects design matrix
int<lower = 1> ngroup;      // number of subjects/clusters/groups
matrix[ntot, p] x;          // fixed effects design matrix
matrix[ntot, q * ngroup] d; // random effects design matrix, block diagonal
vector[5] priors_long; // prior hyperparameters, order: alpha, Omega, sigma_B, sigma_Z, beta
int<lower = 1> ncol_a;           // number of columns of a matrix
matrix[ntot, ncol_a] a;          // spline matrix for tv dof

//quadratures
int<lower = 1> Q; //number of Gauss-Legendre quadratures

//survival sub-model data
int<lower = 1> ntot_quad;         // total number of observations for quadratures
vector<lower = 0.000001>[ngroup] S;      // survival times
vector<lower = 0, upper = 1>[ngroup] E;  // event indicators
int<lower = 1> ncol_e;            // number of columns in the e matrix, spline matrix for baseline hazard
matrix[ngroup, ncol_e] e;         // spline matrix for baseline hazard
matrix[ntot_quad, ncol_e] e_quad; // extended spline matrix for baseline hazard for quadrature approx.
int<lower = 1> ncol_c;            // number of columns in the c matrix; survival submodel fixed effects
matrix[ngroup, ncol_c] c;         // survival sub-model fixed effects matrix
matrix[ntot_quad, ncol_c] c_quad; // extended survival sub-model fixed effects matrix for quadrature approx.

matrix[ngroup, p] x_T;                // x matrix at survival times
matrix[ntot_quad, p] x_quad;          // x matrix for quadrature approx 
matrix[ngroup, q * ngroup] d_T;       // d matrix at survival times
matrix[ntot_quad, q * ngroup] d_quad; // d matrix for qaudrature approx

vector[3] priors_surv;  //prior hyperparameters, order: zeta, omega, eta

vector[ntot_quad] wt_quad; // extended quadrature weights to be used during Gauss-Legendre approx.

}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
//longitudinal sub-model
vector[p] alpha;              // fixed effects coefficients
matrix[ngroup, q] Bstar;          // random effects coefficients
corr_matrix[q] Omega;         // correlation matrix for random effects
vector<lower = 0>[q] sigma_B; // scale parameters for random effects
real<lower = 0> sigma_Z;      // scale parameter of measurement error
real<lower = 0.01, upper = 0.5> phi_inv; // inverse d.o.f. for B
real<lower = 0.01, upper = 0.5> delta0_inv;
vector[ncol_a] beta;
vector<lower = 0>[ngroup] V;
vector<lower = 0>[ntot] W;

//survival sub-model
vector[ncol_e] zeta;  // spline coefficents for baseline hazard
vector[ncol_c] omega; // fixed effects parameters
real eta;             // association parameter
}

transformed parameters{
cov_matrix[q] Sigma; 
vector[ntot] linpred;
matrix[ngroup, q] B;         
matrix[ngroup * q, 1] Bmat;
real<lower = 2, upper = 100> phi;
real<lower = 2, upper = 100> delta0;
vector<lower = 0>[ntot] delta;

vector[ngroup] lsd_expr1;
vector[ngroup] lsd_expr1_spl;
vector[ngroup] lsd_expr1_fix;
vector[ngroup] lsd_expr1_ystar;
vector[ngroup] lsd_expr2;
vector[ntot_quad] lsd_expr2_quad;
vector[ntot_quad] lsd_expr2_quad_spl;
vector[ntot_quad]lsd_expr2_quad_fix;
vector[ntot_quad] lsd_expr2_quad_ystar;
vector[ngroup] lsd;

//longitudinal sub-model

phi = 1/phi_inv;
delta0 = 1/delta0_inv;

delta = to_vector(delta0 * exp(a * beta));

for(i in 1:ngroup){
B[i, ] = Bstar[i, ] * sqrt(V[i]);
}

Bmat = to_matrix(B', ngroup * q, 1);
linpred = x * alpha + to_vector(d * Bmat);

Sigma = quad_form_diag(Omega, sigma_B);

//survival sub-model, lsd: log-survival density
lsd_expr1_spl = e * zeta ; 
lsd_expr1_fix = c * omega; 
lsd_expr1_ystar = x_T * alpha + to_vector(d_T * Bmat);

lsd_expr1 = E .* (lsd_expr1_spl +lsd_expr1_fix + rep_vector(eta, ngroup) .* lsd_expr1_ystar);

lsd_expr2_quad_spl = e_quad * zeta; 
lsd_expr2_quad_fix = c_quad * omega; 
lsd_expr2_quad_ystar = x_quad * alpha + to_vector(d_quad * Bmat); 

lsd_expr2_quad = wt_quad .* exp(lsd_expr2_quad_spl + lsd_expr2_quad_fix + rep_vector(eta, ntot_quad) .* lsd_expr2_quad_ystar);

for(i in 1:ngroup){
lsd_expr2[i] = 0.5 * S[i] * sum(lsd_expr2_quad[((i-1)*Q+1):(i*Q)]);
}

lsd = lsd_expr1 - lsd_expr2;

}

model{

alpha ~ cauchy(0, priors_long[1]);

for(i in 1:ngroup){
Bstar[i] ~ multi_normal(zero_B, Sigma);
}

Omega ~ lkj_corr(priors_long[2]);
sigma_B ~ cauchy(0, priors_long[3]);
sigma_Z ~ cauchy(0, priors_long[4]);

V ~ inv_gamma(phi/2, phi/2);
//phi_inv ~ uniform(0.01, 0.5);//the prior is uniform with -infty and infty, constrained above
for(i in 1:ntot){
W[i] ~ inv_gamma(delta[i]/2, delta[i]/2);
}
beta ~ cauchy(0, priors_long[5]);
//prior for delta0_inv is uniform -infty, infty; constrained above

for(i in 1:ntot) y[i] ~ normal(linpred[i], sigma_Z * sqrt(W[i]));

zeta ~ cauchy(0, priors_surv[1]);
omega ~ cauchy(0, priors_surv[2]);
eta ~ cauchy(0, priors_surv[3]);

target += lsd;

}

generated quantities{
real sigmasq;
sigmasq = sigma_Z^2;
}

"