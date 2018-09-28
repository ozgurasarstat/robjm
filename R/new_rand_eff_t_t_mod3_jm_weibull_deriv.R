new_rand_eff_t_t_mod3_jm_weibull_deriv = "

data{
int<lower = 1> ntot;
int<lower = 1> ngroup;
int<lower = 1> p;
int<lower = 1> q;

int id[ntot];
vector[ntot] y;
matrix[ntot, p] x;
matrix[ntot, q * ngroup] d;

int<lower = 1> Q;
int<lower = 1> ntot_quad;

vector<lower = 0.0>[ngroup] S; 
int<lower = 1> ncol_c;           
matrix[ngroup, ncol_c] c;         
matrix[ntot_quad, ncol_c] c_quad; 
matrix[ngroup, p] x_T;               
matrix[ntot_quad, p] x_quad;          
matrix[ngroup, q * ngroup] d_T;       
matrix[ntot_quad, q * ngroup] d_quad; 

vector[ntot_quad] t_quad;
vector[ntot_quad] wt_quad; 

int<lower = 1> p_deriv;
int<lower = 1> q_deriv;

matrix[ngroup, p_deriv] x_deriv_T;
matrix[ntot_quad, p_deriv] x_deriv_quad;
matrix[ngroup, q_deriv * ngroup] d_deriv_T;
matrix[ntot_quad, q_deriv * ngroup] d_deriv_quad;

int deriv_alpha_ind[p_deriv];
int deriv_B_ind[q_deriv];

vector[p] alpha;
matrix[q, q] Sigma;
real sigma_Z;
real log_lambda;
real log_nu;
vector[ncol_c] omega;
vector[2] eta;
real phi;
real delta;

}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
matrix[ngroup, q] Bstar;
vector<lower = 0>[ngroup] V;
vector<lower = 0>[ntot] W;
}

transformed parameters{

vector[ntot] linpred;
matrix[ngroup, q] B;         
matrix[ngroup * q, 1] Bmat;
vector[p_deriv] alpha_deriv;
matrix[ngroup, q_deriv] B_deriv;
matrix[ngroup * q_deriv, 1] Bmat_deriv;

//vector[ngroup] lsd_expr1;
//vector[ngroup] lsd_expr1_bh;
//vector[ngroup] lsd_expr1_fix;
//vector[ngroup] lsd_expr1_ystar;
vector[ngroup] lsd_expr2;
vector[ntot_quad] lsd_expr2_quad;
vector[ntot_quad] lsd_expr2_quad_bh;
vector[ntot_quad] lsd_expr2_quad_fix;
vector[ntot_quad] lsd_expr2_quad_ystar;
vector[ntot_quad] lsd_expr2_quad_ystar_deriv;
//vector[ngroup] lsd;

//longitudinal sub-model
for(i in 1:ngroup){
B[i, ] = Bstar[i, ] * sqrt(V[i]);
}
Bmat = to_matrix(B', ngroup * q, 1);
linpred = x * alpha + to_vector(d * Bmat);

for(i in 1:p_deriv) alpha_deriv[i] = alpha[deriv_alpha_ind[i]];
for(i in 1:q_deriv) B_deriv[, i] = B[, deriv_B_ind[i]];
Bmat_deriv = to_matrix(B_deriv', ngroup * q_deriv, 1);

//survival sub-model, lsd: log-survival density
//lsd_expr1_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(S); 
//lsd_expr1_fix = c * omega; 
//lsd_expr1_ystar = x_T * alpha + to_vector(d_T * Bmat);

//lsd_expr1 = E .* (lsd_expr1_bh + lsd_expr1_fix + rep_vector(eta, ngroup) .* lsd_expr1_ystar);

lsd_expr2_quad_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad); 
lsd_expr2_quad_fix = c_quad * omega; 
lsd_expr2_quad_ystar = x_quad * alpha + to_vector(d_quad * Bmat);
lsd_expr2_quad_ystar_deriv = x_deriv_quad * alpha_deriv + to_vector(d_deriv_quad * Bmat_deriv);

lsd_expr2_quad = wt_quad .* exp(lsd_expr2_quad_bh + lsd_expr2_quad_fix + 
                                rep_vector(eta[1], ntot_quad) .* lsd_expr2_quad_ystar + 
                                rep_vector(eta[2], ntot_quad) .* lsd_expr2_quad_ystar_deriv);

for(i in 1:ngroup){
lsd_expr2[i] = 0.5 * S[i] * sum(lsd_expr2_quad[((i-1)*Q+1):(i*Q)]);
}

lsd_expr2 = -1.0 * lsd_expr2; 
//lsd = lsd_expr1 - lsd_expr2;
}

model{

for(i in 1:ngroup){
Bstar[i] ~ multi_normal(zero_B, Sigma);
}

V ~ inv_gamma(phi/2, phi/2);
W ~ inv_gamma(delta/2, delta/2);

for(i in 1:ntot){
y[i] ~ normal(linpred[i], sigma_Z * sqrt(W[i]));
}

target += lsd_expr2;

}

"