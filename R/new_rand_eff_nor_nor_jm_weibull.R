new_rand_eff_nor_nor_jm_weibull = "

data{
int<lower = 1> ntot;
int<lower = 1> ngroup;
int<lower = 1> p;
int<lower = 1> q;
int<lower = 1> ncol_c;

matrix[ntot, 1] Y;
matrix[ntot, p] x;
matrix[ntot, q * ngroup] d;
matrix[ngroup, 1] S;
matrix[ngroup, ncol_c] c;

matrix[p, 1] alpha;
matrix[q, q] Sigma;
matrix[p, 1] sigma;
matrix[1, 1] lambda;
matrix[1, 1] nu;
matrix[ncol_c, 1] omega;
matrix[1, 1] gamma;

int<lower = 1> Q;


}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
matrix[ngroup, q] B;
}

transformed parameters{

}

model{

Y ~ normal(linpred, sigma);

for(i in 1:ngroup) B[i, ] ~ normal(zero_B, Sigma);

}

"