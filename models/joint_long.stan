data {
    //dimensions
    int<lower=1> m;     // number of individuals
    int<lower=1> n;     // number of observations
    int<lower=1> p;     // number of predictors (X)
    int<lower=1> ps;     // number of survival predictors (surv_X)
    int<lower=1> q;     // number of predictors r/d random effect
	int<lower=1> N;     // number of deficits

    //input variables
    matrix[n, p] X;     // predictors for damage/repair
    matrix[n, ps] surv_X;     // predictors for survival
    matrix[n, q] Z;     // predictors for damage/repair random effects
	
    int<lower=0,upper=m> index[n];
	int index_size[m];
    real<lower=0> offsets[n, 3];

    //response
    int<lower=0> repair[n];  // count response
    int<lower=0> damage[n];  // count response
    vector[n] Tstart;  // Time to event response
    vector[n] T;  // Time to event response
    vector[n] status;  // Time to event response

    //spline input
    int num_basis; //number of spline basis
    int num_quad; //number of points for gaussian quadrature
	
    real quad_splines[n, num_quad, num_basis]; //msplines for integrating hazard with gaussian quadrature
    matrix[n, num_basis] msplines; //msplines for hazard
    vector[num_quad] quad_nodes; //gaussian quadrature points
    vector[num_quad] quad_weights; // gaussian quadrature weights
    
}
transformed data {
	vector[4] zero_mu = rep_vector(0, 4);
}
parameters {
	
	//fixed effect longitudinal parameters
	vector[p] beta_r;         // repair coefficients
	vector[p] beta_d;         // damage coefficients
	vector[2] beta_int;       // intercepts
	
	//random effects parameters
	matrix[m, 2*q] b;
	cholesky_factor_corr[4] Lcorr;  
	vector<lower=0>[4] tau_rand;
    
	//survival parameters
	real gamma_int; // intercept
	vector[ps] gamma; //coeficients
	vector[2] gamma_rand;// gamma random effect
	
	//spline parameters
	simplex[num_basis] spline_gamma; 
    
}
transformed parameters {
  
    vector[n] log_lambda_r;
    vector[n] log_lambda_d;
    
    vector[n] hazard;
    vector[n] cum_hazard;
    
    for (i in 1:n) {
	  vector[num_quad] quad_x;
	
	  //repair rate
	  log_lambda_r[i] = log(log(1 + exp(beta_int[1] + dot_product(X[i], beta_r) + dot_product(Z[i], b[index[i],1:2]))));
	
	  //damage rate
	  log_lambda_d[i] = log(log(1 + exp(beta_int[2] + dot_product(X[i], beta_d) + dot_product(Z[i], b[index[i],3:4]))));
	  
	  hazard[i] = dot_product(spline_gamma, msplines[i]) .* exp(gamma_int
			      + gamma_rand[1]*(beta_int[1] + beta_r[1]*T[i] + beta_r[2]*X[i,2] + beta_r[3]*X[i,3] + b[index[i], 1] + b[index[i], 2]*T[i])
			      + gamma_rand[2]*(beta_int[2] + beta_d[1]*T[i] + beta_d[2]*X[i,2] + beta_d[3]*X[i,3] + b[index[i], 3] + b[index[i], 4]*T[i])
																+ dot_product(surv_X[i], gamma));
	  
	
	  //need to convert T and Tstart to age scale
	  quad_x = (T[i] - Tstart[i])/2*quad_nodes + (T[i] + Tstart[i])/2;
	  cum_hazard[i] = (T[i] - Tstart[i])/2*dot_product(to_vector(to_row_vector(spline_gamma) * to_matrix(quad_splines[i])') .* exp(gamma_int
	                      + gamma_rand[1]*(beta_int[1] + beta_r[1]*quad_x + beta_r[2]*X[i,2] + beta_r[3]*X[i,3] + b[index[i], 1] + b[index[i], 2]*quad_x)
						  + gamma_rand[2]*(beta_int[2] + beta_d[1]*quad_x + beta_d[2]*X[i,2] + beta_d[3]*X[i,3] + b[index[i], 3] + b[index[i], 4]*quad_x)
														 + dot_product(surv_X[i], gamma)), quad_weights);
    }
    
}
model {
    
	//longitudinal priors
	beta_int ~ normal(0.0, 3.0);
	beta_r ~ std_normal();
	beta_d ~ std_normal();
	
	//survival priors
	gamma_int ~ normal(0.0, 3.0);
	gamma ~ std_normal();
	gamma_rand ~ std_normal();
	
	//correlated random effect priors
	for(i in 1:m) {
		b[i] ~ multi_normal_cholesky(zero_mu, diag_pre_multiply(tau_rand, Lcorr));
	}
	tau_rand ~ cauchy(0.0, 1.0);
	Lcorr ~ lkj_corr_cholesky(2.0);
	
	//spline prior
	spline_gamma ~ dirichlet(rep_vector(1.0, num_basis));
    
	//likelihood 
	target += status .* log(hazard) - cum_hazard;
	target += poisson_log_lpmf(repair | log_lambda_r + to_vector(log(offsets[, 1])) + to_vector(log(offsets[, 3]))); 
	target += poisson_log_lpmf(damage | log_lambda_d + log(to_vector(offsets[, 2]) - to_vector(offsets[, 1])) + to_vector(log(offsets[, 3])));
}
generated quantities {

	// sampled repair, damage, and deficit counts
    vector[n] sampled_repair;
    vector[n] sampled_damage;
    vector[n+m] sampled_n;
	
	// rates
    vector[n] lambda_r;
    vector[n] lambda_d;

	// time derivative of rates
	vector[n] deriv_r;
	vector[n] deriv_d;

	// Frailty Index derivative of rates
	vector[n] deriv_r_f;
	vector[n] deriv_d_f;

	// survival function
    vector[n] survival;
	vector[m] output_age;
	vector[m] output_status;
	
	//
	int j;
	j = 1;

	// compute rates and survival function
    lambda_r = exp(log_lambda_r);
    lambda_d = exp(log_lambda_d);
    survival = exp(-cum_hazard);
  
    for (i in 1:n) {

		// sample repair and damage counts
		sampled_repair[i] = poisson_rng(exp(log_lambda_r[i] + log(offsets[i, 1]) + log(offsets[i, 3])));
		sampled_damage[i] = poisson_rng(exp(log_lambda_d[i] + log(offsets[i, 2] - offsets[i, 1]) + log(offsets[i, 3])));

		// time derivatives
		deriv_r[i] = (beta_r[1] + b[index[i],2])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));
		deriv_d[i] = (beta_d[1] + b[index[i],4])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));

		// f derivatives
		deriv_r_f[i] = (beta_r[2])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));
		deriv_d_f[i] = (beta_d[2])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));
		
		if (i == 1) {
			//set first deficit count as observed count at baseline
			sampled_n[j] = offsets[i,1];
		}
		else {

			// if this is a new individual
			if(i < n && index[i] != index[i-1]) {
				
				sampled_n[j] = offsets[i-1,1] + fmin(sampled_damage[i-1], (116.0-offsets[i-1,1]))
					- fmin(offsets[i-1,1], sampled_repair[i-1]);
				
				j = j+1;
				sampled_n[j] =  offsets[i,1];
			}
			// if this is the same indiviudal
			else {
				
				sampled_n[j] = offsets[i-1,1] + fmin(sampled_damage[i-1], (116.0-offsets[i-1,1]))
					- fmin(offsets[i-1,1], sampled_repair[i-1]);
				
			}
			// if this is the last point
			if (i == n) {
				j += 1;
				
				sampled_n[j] = offsets[i,1] + fmin(sampled_damage[i], (116.0-offsets[i,1]))
					- fmin(offsets[i,1], sampled_repair[i]);
				break;
			}
			
			
			
		}
		j = j+1;
		
    }

	// posterior predictive check
	{
		
		int indiv = 1;
		
		for(k in 1:m) {
			
			int alive[index_size[k]];
			int cumprod_alive[index_size[k]];
			
			int ii = 1;
			int flag = 0;
			for(i in indiv:indiv+index_size[k]-1) {
				
				real prob_death = exp(-cum_hazard[i]);
				int dead = bernoulli_rng(prob_death);
				
				if(dead == 0) {
					output_age[k] = T[i];
					output_status[k] = 1;//status[T > T[i]]
					flag = 1;
					break;
				}
				
				
			}
			//int ii = 1;
			while(flag == 0) {
				real prob_death = exp(-(cum_hazard[indiv + index_size[k]-1] + ii*hazard[indiv + index_size[k]-1]*0.1));// + ii*(cum_hazard[indiv + index_size[k]-3] - cum_hazard[indiv + index_size[k]-4])));
				int dead = bernoulli_rng(prob_death);
				if(dead==0) {
					output_age[k] = T[indiv + index_size[k]-1] + ii*0.1;//(T[indiv + index_size[k]-1] - T[indiv + index_size[k]-2]);
					output_status[k] = 1;
					break;
				}
				ii = ii + 1;
			}
			
			indiv = indiv + index_size[k];
		}
		
	}
	
}

