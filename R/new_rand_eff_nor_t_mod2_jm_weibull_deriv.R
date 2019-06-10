new_rand_eff_nor_t_mod2_jm_weibull_deriv = "

data{
int<lower = 1> ntot;
int<lower = 1> ngroup;
int<lower = 1> p;
int<lower = 1> q;

int id[ntot];
vector[ntot] y;
matrix[ntot, p] x;
//matrix[ntot, q * ngroup] d;
matrix[ntot, q] d;

int<lower = 1> Q;
int<lower = 1> ntot_quad;

int d_ind[ngroup, 2];
int Q_ind[ngroup, 2];

int nrepeat[ngroup];

vector<lower = 0.0>[ngroup] S; 
int<lower = 1> ncol_c;           
matrix[ngroup, ncol_c] c;         
matrix[ntot_quad, ncol_c] c_quad; 
matrix[ngroup, p] x_T;               
matrix[ntot_quad, p] x_quad;          
//matrix[ngroup, q * ngroup] d_T;       
matrix[ngroup, q] d_T;       
//matrix[ntot_quad, q * ngroup] d_quad; 
matrix[ntot_quad, q] d_quad; 

vector[ntot_quad] t_quad;
vector[ntot_quad] wt_quad; 

int<lower = 1> p_deriv;
int<lower = 1> q_deriv;

matrix[ngroup, p_deriv] x_deriv_T;
matrix[ntot_quad, p_deriv] x_deriv_quad;
//matrix[ngroup, q_deriv * ngroup] d_deriv_T;
matrix[ngroup, q_deriv] d_deriv_T;
//matrix[ntot_quad, q_deriv * ngroup] d_deriv_quad;
matrix[ntot_quad, q_deriv] d_deriv_quad;

int deriv_alpha_ind[p_deriv];
int deriv_B_ind[q_deriv];

vector[p] alpha;
matrix[q, q] Sigma;
real sigma_Z;
real log_lambda;
real log_nu;
vector[ncol_c] omega;
real eta1;
real eta2;
//real phi;
real delta;

}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
//matrix[ngroup, q] Bstar;
vector[q] Bstar[ngroup];
//vector<lower = 0>[ngroup] V;
vector<lower = 0>[ngroup] W;
}

transformed parameters{

matrix[ngroup, q] B;         
//matrix[ngroup * q, 1] Bmat;
vector[p_deriv] alpha_deriv;
matrix[ngroup, q_deriv] B_deriv;
//matrix[ngroup * q_deriv, 1] Bmat_deriv;
vector[ntot] d_B;
//vector[ngroup] d_T_B;
vector[ntot_quad] d_quad_B;
//vector[ngroup] d_deriv_T_B;
vector[ntot_quad] d_deriv_quad_B;

vector<lower = 0>[ntot] W_ext;

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
for(i in 1:q){
B[, i] = to_vector(Bstar[, i]);// .* sqrt(V);
}

//Bmat = to_matrix(B', ngroup * q, 1);

for(i in 1:p_deriv) alpha_deriv[i] = alpha[deriv_alpha_ind[i]];
for(i in 1:q_deriv) B_deriv[, i] = B[, deriv_B_ind[i]];
//Bmat_deriv = to_matrix(B_deriv', ngroup * q_deriv, 1);

for(i in 1:ngroup){
W_ext[d_ind[i, 1]:d_ind[i, 2]] = rep_vector(W[i], nrepeat[i]);

d_B[d_ind[i, 1]:d_ind[i, 2]] = d[d_ind[i, 1]:d_ind[i, 2], ] * to_vector(B[i]);
d_quad_B[Q_ind[i, 1]:Q_ind[i, 2]] = d_quad[Q_ind[i, 1]:Q_ind[i, 2], ] * to_vector(B[i]);
d_deriv_quad_B[Q_ind[i, 1]:Q_ind[i, 2]] = d_deriv_quad[Q_ind[i, 1]:Q_ind[i, 2]] * to_vector(B_deriv[i]);
}

//linpred = x * alpha + to_vector(d * Bmat);

//survival sub-model, lsd: log-survival density
//lsd_expr1_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(S); 
//lsd_expr1_fix = c * omega; 
//lsd_expr1_ystar = x_T * alpha + to_vector(d_T * Bmat);

//lsd_expr1 = E .* (lsd_expr1_bh + lsd_expr1_fix + rep_vector(eta, ngroup) .* lsd_expr1_ystar);

lsd_expr2_quad_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad); 
lsd_expr2_quad_fix = c_quad * omega; 
//lsd_expr2_quad_ystar = x_quad * alpha + to_vector(d_quad * Bmat);
lsd_expr2_quad_ystar = x_quad * alpha + d_quad_B;
//lsd_expr2_quad_ystar_deriv = x_deriv_quad * alpha_deriv + to_vector(d_deriv_quad * Bmat_deriv);
lsd_expr2_quad_ystar_deriv = x_deriv_quad * alpha_deriv + d_deriv_quad_B;

lsd_expr2_quad = wt_quad .* exp(lsd_expr2_quad_bh + lsd_expr2_quad_fix + 
                                eta1 * lsd_expr2_quad_ystar + 
                                eta2 * lsd_expr2_quad_ystar_deriv);

for(i in 1:ngroup){
lsd_expr2[i] = 0.5 * S[i] * sum(lsd_expr2_quad[Q_ind[i, 1]:Q_ind[i, 2]]);
}

//lsd_expr2 = -1.0 * lsd_expr2; 
//lsd = lsd_expr1 - lsd_expr2;
}

model{

Bstar ~ multi_normal(zero_B, Sigma);

//V ~ inv_gamma(phi/2, phi/2);
W ~ inv_gamma(delta/2, delta/2);

y ~ normal(x * alpha + d_B, sigma_Z * sqrt(W_ext));

target += lsd_expr2;

}

"