t_t_jm_mod1_weibull = "

data{
//longitudinal sub-model data
int<lower = 1> ntot;        // total number observations
int id[ntot];               // id matrix, takes values like 1, 2, 3, ...
vector[ntot] y;             // response matrix
int<lower = 1> p;           // number of covariates in the fixed effects design matrix
int<lower = 1> q;           // number of covariates in the random effects design matrix
int<lower = 1> ngroup;      // number of subjects/clusters/groups
matrix[ntot, p] x;          // fixed effects design matrix
//matrix[ntot, q * ngroup] d; // random effects design matrix, block diagonal
matrix[ntot, q] d;
vector[4] priors_long; // prior hyperparameters, order: alpha, Omega, sigma_B, sigma_Z

int d_ind[ngroup, 2];
int Q_ind[ngroup, 2];

int nrepeat[ngroup];

//quadratures
int<lower = 1> Q; //number of Gauss-Legendre quadratures
//survival sub-model data
int<lower = 1> ntot_quad;         // total number of observations for quadratures
vector<lower = 0.000001>[ngroup] S;      // survival times
vector<lower = 0, upper = 1>[ngroup] E;  // event indicators
//int<lower = 1> ncol_e;            // number of columns in the e matrix, spline matrix for baseline hazard
//matrix[ngroup, ncol_e] e;         // spline matrix for baseline hazard
//matrix[ntot_quad, ncol_e] e_quad; // extended spline matrix for baseline hazard for quadrature approx.
int<lower = 1> ncol_c;            // number of columns in the c matrix; survival submodel fixed effects
matrix[ngroup, ncol_c] c;         // survival sub-model fixed effects matrix
matrix[ntot_quad, ncol_c] c_quad; // extended survival sub-model fixed effects matrix for quadrature approx.

matrix[ngroup, p] x_T;                // x matrix at survival times
matrix[ntot_quad, p] x_quad;          // x matrix for quadrature approx 
//matrix[ngroup, q * ngroup] d_T;       // d matrix at survival times
matrix[ngroup, q] d_T;
//matrix[ntot_quad, q * ngroup] d_quad; // d matrix for qaudrature approx
matrix[ntot_quad, q] d_quad;

vector[4] priors_surv;  //prior hyperparameters, order: log_lambda, log_nu, omega, eta
vector[ntot_quad] wt_quad; // extended quadrature weights to be used during Gauss-Legendre approx.
vector[ntot_quad] t_quad;
}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{

//longitudinal sub-model
vector[p] alpha;              // fixed effects coefficients
//matrix[ngroup, q] Bstar;          // random effects coefficients
vector[q] Bstar[ngroup];
corr_matrix[q] Omega;         // correlation matrix for random effects
vector<lower = 0>[q] sigma_B; // scale parameters for random effects
real<lower = 0> sigma_Z;      // scale parameter of measurement error
real<lower = 0.01, upper = 0.5> phi_inv; // inverse d.o.f. for B and Z
vector<lower = 0>[ngroup] V;

//survival sub-model
//vector[ncol_e] zeta;  // spline coefficents for baseline hazard
real log_lambda;
real log_nu;
vector[ncol_c] omega; // fixed effects parameters
real eta;             // association parameter
}

transformed parameters{
cov_matrix[q] Sigma; 
//vector[ntot] linpred;
matrix[ngroup, q] B;         
//matrix[ngroup * q, 1] Bmat;
vector[ntot] d_B;
vector[ngroup] d_T_B;
vector[ntot_quad] d_quad_B;
//real<lower = 2, upper = 100> phi;

vector<lower = 0>[ntot] V_ext;

vector[ngroup] lsd_expr1;
vector[ngroup] lsd_expr1_bh;
vector[ngroup] lsd_expr1_fix;
vector[ngroup] lsd_expr1_ystar;
vector[ngroup] lsd_expr2;
vector[ntot_quad] lsd_expr2_quad;
vector[ntot_quad] lsd_expr2_quad_bh;
vector[ntot_quad]lsd_expr2_quad_fix;
vector[ntot_quad] lsd_expr2_quad_ystar;
vector[ngroup] lsd;

//longitudinal sub-model
//phi = 1/phi_inv;
//for(i in 1:ngroup){
//B[i, ] = Bstar[i, ] * sqrt(V[i]);
//}

for(i in 1:q){
B[, i] = to_vector(Bstar[, i]) .* sqrt(V);
}

//for(i in 1:ntot) V_ext[i] = V[id[i]];

//Bmat = to_matrix(B', ngroup * q, 1);

for(i in 1:ngroup){
V_ext[d_ind[i, 1]:d_ind[i, 2]] = rep_vector(V[i], nrepeat[i]);
d_B[d_ind[i, 1]:d_ind[i, 2]] = d[d_ind[i, 1]:d_ind[i, 2], ] * to_vector(B[i]);
d_T_B[i] = sum(d_T[i] .* B[i]);
d_quad_B[Q_ind[i, 1]:Q_ind[i, 2]] = d_quad[Q_ind[i, 1]:Q_ind[i, 2]] * to_vector(B[i]);
}

//linpred = x * alpha + to_vector(d * Bmat);
Sigma = quad_form_diag(Omega, sigma_B);

//survival sub-model, lsd: log-survival density
lsd_expr1_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(S);  
lsd_expr1_fix = c * omega; 
//lsd_expr1_ystar = x_T * alpha + to_vector(d_T * Bmat);
lsd_expr1_ystar = x_T * alpha + d_T_B;

lsd_expr1 = E .* (lsd_expr1_bh + 
                  lsd_expr1_fix + 
                  eta * lsd_expr1_ystar);

lsd_expr2_quad_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad);  
lsd_expr2_quad_fix = c_quad * omega; 
//lsd_expr2_quad_ystar = x_quad * alpha + to_vector(d_quad * Bmat); 
lsd_expr2_quad_ystar = x_quad * alpha + d_quad_B; 

lsd_expr2_quad = wt_quad .* exp(lsd_expr2_quad_bh + 
                                lsd_expr2_quad_fix + 
                                eta * lsd_expr2_quad_ystar);

for(i in 1:ngroup){
lsd_expr2[i] = 0.5 * S[i] * sum(lsd_expr2_quad[((i-1)*Q+1):(i*Q)]);
}

lsd = lsd_expr1 - lsd_expr2;
}

model{

alpha[1] ~ cauchy(0, priors_long[1] * 4);
for(i in 2:p) alpha[i] ~ cauchy(0, priors_long[1]);

Bstar ~ multi_normal(zero_B, Sigma);

Omega ~ lkj_corr(priors_long[2]);
sigma_B ~ cauchy(0, priors_long[3]);
sigma_Z ~ cauchy(0, priors_long[4]);

V ~ inv_gamma((1/phi_inv)/2, (1/phi_inv)/2);
//phi_inv ~ uniform(0.01, 0.5);//the prior is uniform with -infty and infty, constrained above

y ~ normal(x * alpha + d_B, sigma_Z * sqrt(V_ext));

//zeta ~ cauchy(0, priors_surv[1]);
log_lambda ~ cauchy(0, priors_surv[1]);
log_nu ~ cauchy(0, priors_surv[2]);
omega ~ cauchy(0, priors_surv[3]);
eta ~ cauchy(0, priors_surv[4]);

target += lsd;
}

generated quantities{
real sigmasq;
real lambda;
real nu;
real phi;

lambda = exp(log_lambda);
nu = exp(log_nu);
sigmasq = sigma_Z^2;
phi = 1/phi_inv;
}

"