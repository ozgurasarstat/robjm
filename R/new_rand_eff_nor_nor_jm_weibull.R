new_rand_eff_nor_nor_jm_weibull = "

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

vector<lower = 0.000001>[ngroup] S; 
int<lower = 1> ncol_c;           
matrix[ngroup, ncol_c] c;         
matrix[ntot_quad, ncol_c] c_quad; 
matrix[ngroup, p] x_T;               
matrix[ntot_quad, p] x_quad;          
matrix[ngroup, q * ngroup] d_T;       
matrix[ntot_quad, q * ngroup] d_quad; 

vector[ntot_quad] t_quad;
vector[ntot_quad] wt_quad; 

vector[p] alpha;
matrix[q, q] Sigma;
real sigma_Z;
real log_lambda;
real log_nu;
vector[ncol_c] omega;
real eta;

}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
matrix[ngroup, q] B;
}

transformed parameters{
vector[ntot] linpred;
matrix[ngroup * q, 1] Bmat;

//vector[ngroup] lsd_expr1;
//vector[ngroup] lsd_expr1_bh;
//vector[ngroup] lsd_expr1_fix;
//vector[ngroup] lsd_expr1_ystar;
vector[ngroup] lsd_expr2;
vector[ntot_quad] lsd_expr2_quad;
vector[ntot_quad] lsd_expr2_quad_bh;
vector[ntot_quad] lsd_expr2_quad_fix;
vector[ntot_quad] lsd_expr2_quad_ystar;
//vector[ngroup] lsd;

//longitudinal sub-model
Bmat = to_matrix(B', ngroup * q, 1);
linpred = x * alpha + to_vector(d * Bmat);

//survival sub-model, lsd: log-survival density
//lsd_expr1_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(S); 
//lsd_expr1_fix = c * omega; 
//lsd_expr1_ystar = x_T * alpha + to_vector(d_T * Bmat);

//lsd_expr1 = E .* (lsd_expr1_bh +lsd_expr1_fix + rep_vector(eta, ngroup) .* lsd_expr1_ystar);

lsd_expr2_quad_bh = log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad); 
lsd_expr2_quad_fix = c_quad * omega; 
lsd_expr2_quad_ystar = x_quad * alpha + to_vector(d_quad * Bmat); 

lsd_expr2_quad = wt_quad .* exp(lsd_expr2_quad_bh + lsd_expr2_quad_fix + rep_vector(eta, ntot_quad) .* lsd_expr2_quad_ystar);

for(i in 1:ngroup){
lsd_expr2[i] = 0.5 * S[i] * sum(lsd_expr2_quad[((i-1)*Q+1):(i*Q)]);
}

lsd_expr2 = -1.0 * lsd_expr2; 
//lsd = lsd_expr1 - lsd_expr2;
}

model{

y ~ normal(linpred, sigma_Z);

for(i in 1:ngroup){
B[i, ] ~ multi_normal(zero_B, Sigma);
}

target += lsd_expr2;

}

"