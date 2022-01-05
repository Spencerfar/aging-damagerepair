functions {
	real partial_sum(int[] y_slice,
					 int start, int end,
					 matrix X, real[,,] spline_2d, int[] index, int[] damage, matrix offsets,
					 vector beta_int, vector beta_r, vector beta_d, matrix b,
					 matrix S_repair_raw_0, matrix S_repair_raw_s, matrix S_repair_raw_t, matrix S_repair_raw_st,
					 matrix S_damage_raw_0, matrix S_damage_raw_s, matrix S_damage_raw_t, matrix S_damage_raw_st,
					 int num_basis, vector tau_repair_we, vector tau_damage_we, vector tau_repair_ba, vector tau_damage_ba,
					 vector p_spline) {
		
		vector[end-start+1] lambda_r;
		vector[end-start+1] lambda_d;
		matrix[end-start+1, 8] spline = rep_matrix(0.0, end-start+1, 8);
		
		matrix[num_basis, num_basis] S_repair_0 = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_s = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_t = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_st = rep_matrix(0.0, num_basis, num_basis);
		
		matrix[num_basis, num_basis] S_damage_0 = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_s = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_t = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_st = rep_matrix(0.0, num_basis, num_basis);
		
		real result = 0;
		
		S_repair_0[1, 1] = S_repair_raw_0[1, 1];
		S_repair_0[1, 2:num_basis] = S_repair_raw_0[1, 1:num_basis-1] + S_repair_raw_0[1,2:num_basis]*tau_repair_ba[1];
		S_repair_0[2:num_basis, 1] = S_repair_raw_0[1:num_basis-1, 1] + S_repair_raw_0[2:num_basis,1]*tau_repair_we[1];
		
		S_repair_s[1, 1] = S_repair_raw_s[1, 1];
		S_repair_s[1, 2:num_basis] = S_repair_raw_s[1, 1:num_basis-1] + S_repair_raw_s[1,2:num_basis]*tau_repair_ba[2];
		S_repair_s[2:num_basis, 1] = S_repair_raw_s[1:num_basis-1, 1] + S_repair_raw_s[2:num_basis,1]*tau_repair_we[2];
		
		S_repair_t[1, 1] = S_repair_raw_t[1, 1];
		S_repair_t[1, 2:num_basis] = S_repair_raw_t[1, 1:num_basis-1] + S_repair_raw_t[1,2:num_basis]*tau_repair_ba[3];
		S_repair_t[2:num_basis, 1] = S_repair_raw_t[1:num_basis-1, 1] + S_repair_raw_t[2:num_basis,1]*tau_repair_we[3];
		
		S_repair_st[1, 1] = S_repair_raw_st[1, 1];
		S_repair_st[1, 2:num_basis] = S_repair_raw_st[1, 1:num_basis-1] + S_repair_raw_st[1,2:num_basis]*tau_repair_ba[4];
		S_repair_st[2:num_basis, 1] = S_repair_raw_st[1:num_basis-1, 1] + S_repair_raw_st[2:num_basis,1]*tau_repair_we[4];
		
		
		S_damage_0[1, 1] = S_damage_raw_0[1, 1];
		S_damage_0[1, 2:num_basis] = S_damage_raw_0[1, 1:num_basis-1] + S_damage_raw_0[1,2:num_basis]*tau_damage_ba[1];
		S_damage_0[2:num_basis, 1] = S_damage_raw_0[1:num_basis-1, 1] + S_damage_raw_0[2:num_basis,1]*tau_damage_we[1];
		
		S_damage_s[1, 1] = S_damage_raw_s[1, 1];
		S_damage_s[1, 2:num_basis] = S_damage_raw_s[1, 1:num_basis-1] + S_damage_raw_s[1,2:num_basis]*tau_damage_ba[2];
		S_damage_s[2:num_basis, 1] = S_damage_raw_s[1:num_basis-1, 1] + S_damage_raw_s[2:num_basis,1]*tau_damage_we[2];
		
		S_damage_t[1, 1] = S_damage_raw_t[1, 1];
		S_damage_t[1, 2:num_basis] = S_damage_raw_t[1, 1:num_basis-1] + S_damage_raw_t[1,2:num_basis]*tau_damage_ba[3];
		S_damage_t[2:num_basis, 1] = S_damage_raw_t[1:num_basis-1, 1] + S_damage_raw_t[2:num_basis,1]*tau_damage_we[3];
		
		S_damage_st[1, 1] = S_damage_raw_st[1, 1];
		S_damage_st[1, 2:num_basis] = S_damage_raw_st[1, 1:num_basis-1] + S_damage_raw_st[1,2:num_basis]*tau_damage_ba[4];
		S_damage_st[2:num_basis, 1] = S_damage_raw_st[1:num_basis-1, 1] + S_damage_raw_st[2:num_basis,1]*tau_damage_we[4];
		
		
		for(i in 1:num_basis) {
			for(j in 1:num_basis) {
				if((i>1) && (j>1)) {
					S_repair_0[i,j] = p_spline[1]*(S_repair_0[i-1, j] + S_repair_raw_0[i,j]*tau_repair_we[1]) +
						p_spline[2]*(S_repair_0[i, j-1] + S_repair_raw_0[i,j]*tau_repair_ba[1]);
					
					S_repair_s[i,j] = p_spline[1]*(S_repair_s[i-1, j] + S_repair_raw_s[i,j]*tau_repair_we[2]) +
						p_spline[2]*(S_repair_s[i, j-1] + S_repair_raw_s[i,j]*tau_repair_ba[2]);
					
					S_repair_t[i,j] = p_spline[1]*(S_repair_t[i-1, j] + S_repair_raw_t[i,j]*tau_repair_we[3]) +
						p_spline[2]*(S_repair_t[i, j-1] + S_repair_raw_t[i,j]*tau_repair_ba[3]);
					
					S_repair_st[i,j] = p_spline[1]*(S_repair_st[i-1, j] + S_repair_raw_st[i,j]*tau_repair_we[4]) +
						p_spline[2]*(S_repair_st[i, j-1] + S_repair_raw_st[i,j]*tau_repair_ba[4]);

					S_damage_0[i,j] = p_spline[1]*(S_damage_0[i-1, j] + S_damage_raw_0[i,j]*tau_damage_we[1]) +
						p_spline[2]*(S_damage_0[i, j-1] + S_damage_raw_0[i,j]*tau_damage_ba[1]);
					
					S_damage_s[i,j] = p_spline[1]*(S_damage_s[i-1, j] + S_damage_raw_s[i,j]*tau_damage_we[2]) +
						p_spline[2]*(S_damage_s[i, j-1] + S_damage_raw_s[i,j]*tau_damage_ba[2]);
					
					S_damage_t[i,j] = p_spline[1]*(S_damage_t[i-1, j] + S_damage_raw_t[i,j]*tau_damage_we[3]) +
						p_spline[2]*(S_damage_t[i, j-1] + S_damage_raw_t[i,j]*tau_damage_ba[3]);
					
					S_damage_st[i,j] = p_spline[1]*(S_damage_st[i-1, j] + S_damage_raw_st[i,j]*tau_damage_we[4]) +
						p_spline[2]*(S_damage_st[i, j-1] + S_damage_raw_st[i,j]*tau_damage_ba[4]);
				}
				
				spline[,1] = spline[,1] + to_vector(spline_2d[start:end,i,j]) * S_repair_0[i,j];
				spline[,2] = spline[,2] + to_vector(spline_2d[start:end,i,j]) * S_repair_s[i,j];
				spline[,3] = spline[,3] + to_vector(spline_2d[start:end,i,j]) * S_repair_t[i,j];
				spline[,4] = spline[,4] + to_vector(spline_2d[start:end,i,j]) * S_repair_st[i,j];
				
				spline[,5] = spline[,5] + to_vector(spline_2d[start:end,i,j]) * S_damage_0[i,j];
				spline[,6] = spline[,6] + to_vector(spline_2d[start:end,i,j]) * S_damage_s[i,j];
				spline[,7] = spline[,7] + to_vector(spline_2d[start:end,i,j]) * S_damage_t[i,j];
				spline[,8] = spline[,8] + to_vector(spline_2d[start:end,i,j]) * S_damage_st[i,j];
				
			}
		}
		
		lambda_r = log(1+exp(beta_int[1] + X[start:end]*beta_r + spline[,1] + spline[,2] .* X[start:end,1] + spline[,3] .* X[start:end,2] + spline[,4].*X[start:end,3] + b[index[start:end],1])) .* offsets[start:end, 1] .* offsets[start:end, 3];
		
		lambda_d = log(1+exp(beta_int[2] + X[start:end]*beta_d + spline[,5] + spline[,6] .* X[start:end,1] + spline[,7].*X[start:end,2] + spline[,8].*X[start:end,3] + b[index[start:end],2])) .* (offsets[start:end, 2] - offsets[start:end, 1]) .* offsets[start:end, 3];
		
		
		for (i in 1:(end-start+1)) {
	    	
			if (offsets[start:end, 1][i] > 0)
				result += poisson_lpmf(y_slice[i] | lambda_r[i]);
			
			if ((offsets[start:end, 2][i] - offsets[start:end, 1][i]) > 0)
				result += poisson_lpmf(damage[start:end][i] | lambda_d[i]);
			
		}
		return result;
		
	}
	
}
data {
    //dimensions
    int<lower=1> m;     // number of individuals
    int<lower=1> n;     // number of observations
    int<lower=1> p;     // number of predictors (X)
	
    //input variables
    matrix[n, p] X;     // predictors for damage/repair
    int<lower=0, upper=m> index[n];
    matrix[n, 3] offsets;
	
    //response
    int<lower=0> repair[n];  // count response
    int<lower=0> damage[n];  // count response
	
    //spline input
    int num_basis;
	
	//2d spline
	real spline_2d[n, num_basis, num_basis];
	real spline_2d_deriv[n, num_basis, num_basis];
	
}
transformed data {
    vector[2] zero_mu;
	vector[2] p_alpha = rep_vector(1.5, 2);
    zero_mu = rep_vector(0.0, 2);
	
}
parameters {
	
    //fixed effect longitudinal parameters
    vector[p] beta_r;         // repair
    vector[p] beta_d;         // damage
    vector[2] beta_int;       // intercepts
	
	vector<lower=0>[4] tau_repair_we;
	vector<lower=0>[4] tau_damage_we;
	vector<lower=0>[4] tau_repair_ba;
	vector<lower=0>[4] tau_damage_ba;
	
	simplex[2] p_spline;
    
	
	//random effects parameters
	cholesky_factor_corr[2] Lcorr;  
	vector<lower=0>[2] tau_rand;
    matrix[m,2] b;
	
	//2D splines
	matrix[num_basis, num_basis] S_repair_raw_0;
	matrix[num_basis, num_basis] S_repair_raw_s;
	matrix[num_basis, num_basis] S_repair_raw_t;
	matrix[num_basis, num_basis] S_repair_raw_st;
	
	matrix[num_basis, num_basis] S_damage_raw_0;
	matrix[num_basis, num_basis] S_damage_raw_s;
	matrix[num_basis, num_basis] S_damage_raw_t;
	matrix[num_basis, num_basis] S_damage_raw_st;
	
}
model {
	
	//spline parameters
	to_vector(S_repair_raw_0) ~ std_normal();
	to_vector(S_repair_raw_s) ~ std_normal();
	to_vector(S_repair_raw_t) ~ std_normal();
	to_vector(S_repair_raw_st) ~ std_normal();
	
	to_vector(S_damage_raw_0) ~ std_normal();
	to_vector(S_damage_raw_s) ~ std_normal();
	to_vector(S_damage_raw_t) ~ std_normal();
	to_vector(S_damage_raw_st) ~ std_normal();
	
	//RW prior
	tau_repair_we ~ std_normal();
	tau_damage_we ~ std_normal();
	tau_repair_ba ~ std_normal();
	tau_damage_ba ~ std_normal();
	
	p_spline ~ dirichlet(p_alpha);
	
	//longitudinal priors
	beta_int ~ normal(0.0, 5.0);
	beta_r ~ std_normal();
	beta_d ~ std_normal();
	
	
	//correlated random effect priors
	for(i in 1:m) {
		b[i] ~ multi_normal_cholesky(zero_mu, diag_pre_multiply(tau_rand, Lcorr));
	}
	tau_rand ~ cauchy(0.0, 1.0);
	Lcorr ~ lkj_corr_cholesky(2.0);
	
	target += reduce_sum(partial_sum, repair, 400, X, spline_2d, index, damage, offsets, beta_int, beta_r, beta_d, b, 
						 S_repair_raw_0, S_repair_raw_s, S_repair_raw_t, S_repair_raw_st,
						 S_damage_raw_0, S_damage_raw_s, S_damage_raw_t, S_damage_raw_st,
						 num_basis, tau_repair_we, tau_damage_we, tau_repair_ba, tau_damage_ba, p_spline);
}
generated quantities {
    
    vector[n] sampled_repair;
    vector[n] sampled_damage;
    vector[n+m] sampled_n;
	
	vector[n] lambda_r;
    vector[n] lambda_d;
	
	vector[n] deriv_r;
	vector[n] deriv_d;
	
	vector[n] deriv_r_f;
	vector[n] deriv_d_f;

	vector[n] deriv_r_w;
	vector[n] deriv_d_w;
	{
		int jj;

		matrix[n, 8] spline = rep_matrix(0.0, n, 8);
		matrix[n, 8] spline_deriv = rep_matrix(0.0, n, 8);
		
		matrix[num_basis, num_basis] S_repair_0 = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_s = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_t = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_st = rep_matrix(0.0, num_basis, num_basis);

		matrix[num_basis, num_basis] S_damage_0 = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_s = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_t = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_st = rep_matrix(0.0, num_basis, num_basis);

		matrix[num_basis, num_basis] S_repair_0_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_s_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_t_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_repair_st_deriv = rep_matrix(0.0, num_basis, num_basis);

		matrix[num_basis, num_basis] S_damage_0_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_s_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_t_deriv = rep_matrix(0.0, num_basis, num_basis);
		matrix[num_basis, num_basis] S_damage_st_deriv = rep_matrix(0.0, num_basis, num_basis);
		
		S_repair_0[1, 1] = S_repair_raw_0[1, 1];
		S_repair_0[1, 2:num_basis] = S_repair_raw_0[1, 1:num_basis-1] + S_repair_raw_0[1,2:num_basis]*tau_repair_ba[1];
		S_repair_0[2:num_basis, 1] = S_repair_raw_0[1:num_basis-1, 1] + S_repair_raw_0[2:num_basis,1]*tau_repair_we[1];
		
		S_repair_s[1, 1] = S_repair_raw_s[1, 1];
		S_repair_s[1, 2:num_basis] = S_repair_raw_s[1, 1:num_basis-1] + S_repair_raw_s[1,2:num_basis]*tau_repair_ba[2];
		S_repair_s[2:num_basis, 1] = S_repair_raw_s[1:num_basis-1, 1] + S_repair_raw_s[2:num_basis,1]*tau_repair_we[2];
		
		S_repair_t[1, 1] = S_repair_raw_t[1, 1];
		S_repair_t[1, 2:num_basis] = S_repair_raw_t[1, 1:num_basis-1] + S_repair_raw_t[1,2:num_basis]*tau_repair_ba[3];
		S_repair_t[2:num_basis, 1] = S_repair_raw_t[1:num_basis-1, 1] + S_repair_raw_t[2:num_basis,1]*tau_repair_we[3];
		
		S_repair_st[1, 1] = S_repair_raw_st[1, 1];
		S_repair_st[1, 2:num_basis] = S_repair_raw_st[1, 1:num_basis-1] + S_repair_raw_st[1,2:num_basis]*tau_repair_ba[4];
		S_repair_st[2:num_basis, 1] = S_repair_raw_st[1:num_basis-1, 1] + S_repair_raw_st[2:num_basis,1]*tau_repair_we[4];
		
		
		S_damage_0[1, 1] = S_damage_raw_0[1, 1];
		S_damage_0[1, 2:num_basis] = S_damage_raw_0[1, 1:num_basis-1] + S_damage_raw_0[1,2:num_basis]*tau_damage_ba[1];
		S_damage_0[2:num_basis, 1] = S_damage_raw_0[1:num_basis-1, 1] + S_damage_raw_0[2:num_basis,1]*tau_damage_we[1];
		
		S_damage_s[1, 1] = S_damage_raw_s[1, 1];
		S_damage_s[1, 2:num_basis] = S_damage_raw_s[1, 1:num_basis-1] + S_damage_raw_s[1,2:num_basis]*tau_damage_ba[2];
		S_damage_s[2:num_basis, 1] = S_damage_raw_s[1:num_basis-1, 1] + S_damage_raw_s[2:num_basis,1]*tau_damage_we[2];
		
		S_damage_t[1, 1] = S_damage_raw_t[1, 1];
		S_damage_t[1, 2:num_basis] = S_damage_raw_t[1, 1:num_basis-1] + S_damage_raw_t[1,2:num_basis]*tau_damage_ba[3];
		S_damage_t[2:num_basis, 1] = S_damage_raw_t[1:num_basis-1, 1] + S_damage_raw_t[2:num_basis,1]*tau_damage_we[3];
		
		S_damage_st[1, 1] = S_damage_raw_st[1, 1];
		S_damage_st[1, 2:num_basis] = S_damage_raw_st[1, 1:num_basis-1] + S_damage_raw_st[1,2:num_basis]*tau_damage_ba[4];
		S_damage_st[2:num_basis, 1] = S_damage_raw_st[1:num_basis-1, 1] + S_damage_raw_st[2:num_basis,1]*tau_damage_we[4];
		
		
		for(i in 1:num_basis) {
			for(j in 1:num_basis) {
				if((i>1) && (j>1)) {
					S_repair_0[i,j] = p_spline[1]*(S_repair_0[i-1, j] + S_repair_raw_0[i,j]*tau_repair_we[1]) +
						p_spline[2]*(S_repair_0[i, j-1] + S_repair_raw_0[i,j]*tau_repair_ba[1]);
					
					S_repair_s[i,j] = p_spline[1]*(S_repair_s[i-1, j] + S_repair_raw_s[i,j]*tau_repair_we[2]) +
						p_spline[2]*(S_repair_s[i, j-1] + S_repair_raw_s[i,j]*tau_repair_ba[2]);
					
					S_repair_t[i,j] = p_spline[1]*(S_repair_t[i-1, j] + S_repair_raw_t[i,j]*tau_repair_we[3]) +
						p_spline[2]*(S_repair_t[i, j-1] + S_repair_raw_t[i,j]*tau_repair_ba[3]);
					
					S_repair_st[i,j] = p_spline[1]*(S_repair_st[i-1, j] + S_repair_raw_st[i,j]*tau_repair_we[4]) +
						p_spline[2]*(S_repair_st[i, j-1] + S_repair_raw_st[i,j]*tau_repair_ba[4]);
					
					S_damage_0[i,j] = p_spline[1]*(S_damage_0[i-1, j] + S_damage_raw_0[i,j]*tau_damage_we[1]) +
	                      p_spline[2]*(S_damage_0[i, j-1] + S_damage_raw_0[i,j]*tau_damage_ba[1]);
					
					S_damage_s[i,j] = p_spline[1]*(S_damage_s[i-1, j] + S_damage_raw_s[i,j]*tau_damage_we[2]) +
						p_spline[2]*(S_damage_s[i, j-1] + S_damage_raw_s[i,j]*tau_damage_ba[2]);
					
					S_damage_t[i,j] = p_spline[1]*(S_damage_t[i-1, j] + S_damage_raw_t[i,j]*tau_damage_we[3]) +
						p_spline[2]*(S_damage_t[i, j-1] + S_damage_raw_t[i,j]*tau_damage_ba[3]);
					
					S_damage_st[i,j] = p_spline[1]*(S_damage_st[i-1, j] + S_damage_raw_st[i,j]*tau_damage_we[4]) +
						p_spline[2]*(S_damage_st[i, j-1] + S_damage_raw_st[i,j]*tau_damage_ba[4]);
				}
				
				spline[,1] = spline[,1] + to_vector(spline_2d[,i,j]) * S_repair_0[i,j];
				spline[,2] = spline[,2] + to_vector(spline_2d[,i,j]) * S_repair_s[i,j];
				spline[,3] = spline[,3] + to_vector(spline_2d[,i,j]) * S_repair_t[i,j];
				spline[,4] = spline[,4] + to_vector(spline_2d[,i,j]) * S_repair_st[i,j];
				
				spline[,5] = spline[,5] + to_vector(spline_2d[,i,j]) * S_damage_0[i,j];
				spline[,6] = spline[,6] + to_vector(spline_2d[,i,j]) * S_damage_s[i,j];
				spline[,7] = spline[,7] + to_vector(spline_2d[,i,j]) * S_damage_t[i,j];
				spline[,8] = spline[,8] + to_vector(spline_2d[,i,j]) * S_damage_st[i,j];

				
				spline_deriv[,1] = spline_deriv[,1] + to_vector(spline_2d_deriv[,i,j]) * S_repair_0[i,j];
				spline_deriv[,2] = spline_deriv[,2] + to_vector(spline_2d_deriv[,i,j]) * S_repair_s[i,j];
				spline_deriv[,3] = spline_deriv[,3] + to_vector(spline_2d_deriv[,i,j]) * S_repair_t[i,j];
				spline_deriv[,4] = spline_deriv[,4] + to_vector(spline_2d_deriv[,i,j]) * S_repair_st[i,j];
				
				spline_deriv[,5] = spline_deriv[,5] + to_vector(spline_2d_deriv[,i,j]) * S_damage_0[i,j];
				spline_deriv[,6] = spline_deriv[,6] + to_vector(spline_2d_deriv[,i,j]) * S_damage_s[i,j];
				spline_deriv[,7] = spline_deriv[,7] + to_vector(spline_2d_deriv[,i,j]) * S_damage_t[i,j];
				spline_deriv[,8] = spline_deriv[,8] + to_vector(spline_2d_deriv[,i,j]) * S_damage_st[i,j];
				
			}
		}
		
		lambda_r = log(1+exp(beta_int[1] + X*beta_r + spline[,1] + spline[,2] .* X[,1] + spline[,3] .* X[,2] + spline[,4] .* X[,3] + b[index,1]));
		lambda_d = log(1+exp(beta_int[2] + X*beta_d + spline[,5] + spline[,6] .* X[,1] + spline[,7] .* X[,2] + spline[,8] .* X[,3] + b[index,2]));
	
	
		jj = 1;
    for (i in 1:n) {
		
		real temp_rate;
		
		if (lambda_r[i] * offsets[i, 1] * offsets[i, 3] > 0)
			sampled_repair[i] = poisson_rng(lambda_r[i] * offsets[i, 1] * offsets[i, 3]);
		else
			sampled_repair[i] = 0.0;
		
		if (lambda_d[i] * (offsets[i, 2] - offsets[i, 1]) * offsets[i, 3] > 0)
			sampled_damage[i] = poisson_rng(lambda_d[i] * (offsets[i, 2] - offsets[i, 1]) * offsets[i, 3]);
		else
			sampled_damage[i] = 0.0;
		
		deriv_r[i] = (beta_r[2] + beta_r[3]*X[i,1] + beta_r[8]*X[i,5] + beta_r[9]*X[i,4] + beta_r[10]*X[i,7] + beta_r[11]*X[i,6] + spline[i,3] + spline[i,4]*X[i,1])*exp(lambda_r[i])/(exp(lambda_r[i]) + 1);
	
		deriv_d[i] = (beta_d[2] + beta_r[3]*X[i,1] + beta_d[8]*X[i,5] + beta_d[9]*X[i,4] + beta_d[10]*X[i,7] + beta_d[11]*X[i,6] + spline[i,7] + spline[i,8]*X[i,1])*exp(lambda_d[i])/(exp(lambda_d[i]) + 1);
		
		deriv_r_f[i] = (beta_r[12] + beta_r[13]*X[i,1] + beta_r[14]*X[i, 4] + beta_r[15]*X[i, 5])*exp(lambda_r[i])/(exp(lambda_r[i]) + 1);
		
		deriv_d_f[i] = (beta_d[12] + beta_d[13]*X[i,1] + beta_d[14]*X[i, 4] + beta_d[15]*X[i, 5])*exp(lambda_d[i])/(exp(lambda_d[i]) + 1);
		
		
		deriv_r_w[i] = (beta_r[4] + beta_r[6]*X[i,1] + beta_r[9]*X[i, 2] + beta_r[11]*X[i, 3] + beta_r[14]*X[i,12] + spline_deriv[i,1] + spline_deriv[i,2]*X[i,1] + spline_deriv[i,3]*X[i,2] + spline_deriv[i,4]*X[i,3])*exp(lambda_r[i])/(exp(lambda_r[i]) + 1);
		
		deriv_d_w[i] = (beta_d[4] + beta_d[6]*X[i,1] + beta_d[9]*X[i, 2] + beta_d[11]*X[i, 3] + beta_d[14]*X[i,12] + spline_deriv[i,5] + spline_deriv[6,2]*X[i,1] + spline_deriv[7,3]*X[i,2] + spline_deriv[8,4]*X[i,3])*exp(lambda_d[i])/(exp(lambda_d[i]) + 1);
		
		if (i == 1) {
			
			sampled_n[jj] =  offsets[i,1];
			
		}
		else {
			
			if(index[i] != index[i-1]) {
				
				sampled_n[jj] = offsets[i-1,1] + fmin(sampled_damage[i-1], (23.0-offsets[i-1,1])) - fmin(offsets[i-1,1], sampled_repair[i-1]);
				
				jj = jj+1;
				sampled_n[jj] =  offsets[i,1];
			}
			else {
				
				sampled_n[jj] = offsets[i-1,1] + fmin(sampled_damage[i-1], (23.0-offsets[i-1,1])) - fmin(offsets[i-1,1], sampled_repair[i-1]);
				
			}
			
			if (i == n) {
				jj = jj + 1;
				
				sampled_n[jj] = offsets[i,1] + fmin(sampled_damage[i], (23.0-offsets[i,1]))
					- fmin(offsets[i,1], sampled_repair[i]);
				
				break;
			}
			
		}
		
	jj = jj+1;
	
    }

	}
	
}

