data {
    //dimensions
    int<lower=1> m;     // number of individuals
    int<lower=1> n;     // number of observations
    int<lower=1> p;     // number of predictors (X)
    int<lower=1> ps;     // number of survival predictors (surv_X)
    int<lower=1> q;     // number of random effects 
	int<lower=1> N; // number of deficits
	
    //input variables
    matrix[n, p] X;     // predictors for damage/repair
    matrix[n, ps] surv_X; // predictors for survival
    matrix[n, q] Z;     // predictors for damage/repair random effects
	
    int<lower=0,upper=m> index[n];
	int index_size[m];
    real<lower=0> offsets[n, 3];
    matrix[n,p] X_T;// predictors with event time for time

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

	//mean and std of predictors
	vector[p] means;
	vector[p] stds;
	
}
transformed data {
	vector[2*q] zero_mu = rep_vector(0, 2*q);
}
parameters {
    
	//fixed effect longitudinal parameters
	vector[p] beta_r;         // repair coefficients
	vector[p] beta_d;         // damage coefficients
	vector[2] beta_int;       // intercepts
	
	//random effects parameters
	matrix[m, 2*q] b;
	cholesky_factor_corr[2*q] Lcorr;  
	vector<lower=0>[2*q] tau_rand;
    
	//survival parameters
	real gamma_int; // intercept
	vector[ps] gamma; // coefficients
	vector[4] gamma_repair;
	vector[4] gamma_damage;
	
	//sex-specific spline parameters
	simplex[num_basis] spline_gamma;
	simplex[num_basis] spline_gamma_sex; 
    
}
transformed parameters {
  
    vector[n] log_lambda_r;
    vector[n] log_lambda_d;
    
    vector[n] hazard;
    vector[n] cum_hazard;
    
    for (i in 1:n) {
		vector[num_quad] quad_x;
		matrix[num_quad, p] quad_X;
		vector[num_quad] unnormed_x;
		vector[num_quad] unnormed_s;
		vector[num_quad] unnormed_t;
		vector[num_basis] spline_coef;
		real repair_hazard;
		vector[num_quad] repair_cum_hazard;
		real damage_hazard;
		vector[num_quad] damage_cum_hazard;
		
		// sex spline for male/female
		if (X[i,1] > 0)
			spline_coef = spline_gamma;
		else
			spline_coef = spline_gamma_sex;
		
		//repair rate
		log_lambda_r[i] = log(log(1 + exp(beta_int[1] + dot_product(X[i], beta_r) + dot_product(Z[i], b[index[i],1:q]))));
		
		//damage rate
		log_lambda_d[i] = log(log(1 + exp(beta_int[2] + dot_product(X[i], beta_d) + dot_product(Z[i], b[index[i],q+1:2*q]))));
	
		//need to convert T and Tstart to age scale
		quad_x = (T[i] - Tstart[i])/2*quad_nodes + (T[i] + Tstart[i])/2;
		
		
		// need to create new predictors over the gaussian quadrature prediction range
		quad_X[,1] = rep_vector(X[i, 1], num_quad);
		quad_X[,2] = rep_vector(X[i, 2], num_quad);
		quad_X[,4] = rep_vector(X[i, 4], num_quad);
		quad_X[,5] = rep_vector(X[i, 5], num_quad);
		quad_X[,8] = rep_vector(X[i, 8], num_quad);
		quad_X[,10] = rep_vector(X[i, 10], num_quad);
		quad_X[,11] = rep_vector(X[i, 11], num_quad);
		quad_X[,12] = rep_vector(X[i, 12], num_quad);
		quad_X[,13] = rep_vector(X[i, 13], num_quad);
		
		//need to normalize predictors involving time 
		unnormed_x = quad_x*stds[3]+means[3];
		unnormed_s = rep_vector(X[i,1]*stds[1]+means[1],num_quad);
		unnormed_t = rep_vector(X[i,2]*stds[2]+means[2],num_quad);
		quad_X[,3] = quad_x;
		quad_X[,6] = (unnormed_x .*unnormed_s - means[6])/stds[6];
		quad_X[,7] = (unnormed_x .*unnormed_t - means[7])/stds[7];
		quad_X[,9] = (unnormed_x .*unnormed_s .*unnormed_t - means[9])/stds[9];
		
		
		//female
		if (X[i,1] < 0) {
			//control
			if (X[i,2] < 0) {
				
				repair_hazard = gamma_repair[1]*(beta_int[1] + dot_product(X_T[i], beta_r) + b[index[i], 1] + b[index[i], 2]*T[i]);
				repair_cum_hazard = gamma_repair[1]*(beta_int[1] + quad_X*beta_r + b[index[i], 1] + b[index[i], 2]*quad_x);
				
				damage_hazard = gamma_damage[1]*(beta_int[2] + dot_product(X_T[i], beta_d) + b[index[i], 3] + b[index[i], 4]*T[i]);
				damage_cum_hazard = gamma_damage[1]*(beta_int[2] + quad_X*beta_d + b[index[i], 3] + b[index[i], 4]*quad_x);
				
			}
			//drug
			else{
				
				repair_hazard = (gamma_repair[1]+gamma_repair[3])*(beta_int[1] + dot_product(X_T[i], beta_r) + b[index[i], 1] + b[index[i], 2]*T[i]);
				repair_cum_hazard = (gamma_repair[1]+gamma_repair[3])*(beta_int[1] + quad_X*beta_r + b[index[i], 1] + b[index[i], 2]*quad_x);
				
				damage_hazard = (gamma_damage[1]+gamma_damage[3])*(beta_int[2] + dot_product(X_T[i], beta_d) + b[index[i], 3] + b[index[i], 4]*T[i]);
				damage_cum_hazard = (gamma_damage[1]+gamma_damage[3])*(beta_int[2] + quad_X*beta_d + b[index[i], 3] + b[index[i], 4]*quad_x);
				
		}
		}
		//male
		else {
			//control
			if (X[i,2] < 0) {
				repair_hazard = (gamma_repair[1]+gamma_repair[2])*(beta_int[1] + dot_product(X_T[i], beta_r) + b[index[i], 1] + b[index[i], 2]*T[i]);
				repair_cum_hazard = (gamma_repair[1]+gamma_repair[2])*(beta_int[1] + quad_X*beta_r + b[index[i], 1] + b[index[i], 2]*quad_x);
				
				damage_hazard = (gamma_damage[1] + gamma_damage[2])*(beta_int[2] + dot_product(X_T[i], beta_d) + b[index[i], 3] + b[index[i], 4]*T[i]);
				damage_cum_hazard =  (gamma_damage[1] + gamma_damage[2])*(beta_int[2] + quad_X*beta_d + b[index[i], 3] + b[index[i], 4]*quad_x);
			}
			//drug
			else{
				repair_hazard = (gamma_repair[1]+gamma_repair[2]+gamma_repair[3]+gamma_repair[4])*(beta_int[1] + dot_product(X_T[i], beta_r)  + b[index[i], 1] + b[index[i], 2]*T[i]);
				repair_cum_hazard = (gamma_repair[1]+gamma_repair[2]+gamma_repair[3]+gamma_repair[4])*(beta_int[1] + quad_X*beta_r + b[index[i], 1] + b[index[i], 2]*quad_x);
				
				damage_hazard = (gamma_damage[1]+gamma_damage[2]+gamma_damage[3]+gamma_damage[4])*(beta_int[2] + dot_product(X_T[i], beta_d) + b[index[i], 3] + b[index[i], 4]*T[i]);
				damage_cum_hazard =  (gamma_damage[1]+gamma_damage[2]+gamma_damage[3]+gamma_damage[4])*(beta_int[2] + quad_X*beta_d + b[index[i], 3] + b[index[i], 4]*quad_x);
			}
		}
		
		
		hazard[i] = dot_product(spline_coef, msplines[i]) .* exp(gamma_int 
																 + repair_hazard + damage_hazard + dot_product(surv_X[i], gamma));
		
		cum_hazard[i] = (T[i] - Tstart[i])/2*dot_product(to_vector(to_row_vector(spline_coef) * to_matrix(quad_splines[i])') .*
														 exp(gamma_int + repair_cum_hazard + damage_cum_hazard
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
	gamma_repair ~ std_normal();
	gamma_damage ~ std_normal();
	
	
	//correlated random effect priors
	for(i in 1:m) {
		b[i] ~ multi_normal_cholesky(zero_mu, diag_pre_multiply(tau_rand, Lcorr));
	}
	tau_rand ~ cauchy(0.0, 1.0);
	Lcorr ~ lkj_corr_cholesky(2.0);
	
	//spline priors
	spline_gamma ~ dirichlet(rep_vector(1.0, num_basis));
	spline_gamma_sex ~ dirichlet(rep_vector(1.0, num_basis));
    
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

	// hazard rates for repair and damage rates for sex and treatment groups
	real female_control_repair;
	real male_control_repair;
	real female_drug_repair;
	real male_drug_repair;
	
	real female_control_damage;
	real male_control_damage;
	real female_drug_damage;
	real male_drug_damage;

	// 
	int j;
	j = 1;


	//compute hazard rates for groups
	female_control_repair = gamma_repair[1];
	male_control_repair = gamma_repair[1] + gamma_repair[2];
	female_drug_repair = gamma_repair[1] + gamma_repair[3];
	male_drug_repair = gamma_repair[1] + gamma_repair[2] + gamma_repair[3] + gamma_repair[4];
	
	female_control_damage = gamma_damage[1];
	male_control_damage = gamma_damage[1] + gamma_damage[2];
	female_drug_damage = gamma_damage[1] + gamma_damage[3];
	male_drug_damage = gamma_damage[1] + gamma_damage[2] + gamma_damage[3] + gamma_damage[4];

	// compute rates and survival function
    lambda_r = exp(log_lambda_r);
    lambda_d = exp(log_lambda_d);
    survival = exp(-cum_hazard);
  
    for (i in 1:n) {

		// sample repair and damage counts
		sampled_repair[i] = poisson_rng(exp(log_lambda_r[i] + log(offsets[i, 1]) + log(offsets[i, 3])));
		sampled_damage[i] = poisson_rng(exp(log_lambda_d[i] + log(offsets[i, 2] - offsets[i, 1]) + log(offsets[i, 3])));

		// time derivatives
		deriv_r[i] = (beta_r[3] + X[i,1]*beta_r[6] + X[i,2]*beta_r[7] + X[i,8]*beta_r[9] + b[index[i],2])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));
		deriv_d[i] = (beta_d[3] + X[i,1]*beta_d[6] + X[i,2]*beta_d[7] + X[i,8]*beta_d[9] + b[index[i],4])*exp(lambda_d[i])/(1 + exp(lambda_d[i]));

		// f derivatives
		deriv_r_f[i] = (beta_r[4] + beta_r[11]*X[i,1] + beta_r[12]*X[i,2] + beta_r[13]*X[i,8])*exp(lambda_r[i])/(1 + exp(lambda_r[i]));
		deriv_d_f[i] = (beta_d[4] + beta_d[11]*X[i,1] + beta_d[12]*X[i,2] + beta_d[13]*X[i,8])*exp(lambda_d[i])/(1 + exp(lambda_d[i]));
		
	
		if (i == 1) {
			//set first deficit count as observed count at baseline
			sampled_n[j] = offsets[i,1];
		}
		else {

			// if this is a new individual
			if(i < n && index[i] != index[i-1]) {
				
				sampled_n[j] = offsets[i-1,1] + fmin(sampled_damage[i-1], (N-offsets[i-1,1]))
					- fmin(offsets[i-1,1], sampled_repair[i-1]);
				
				j = j+1;
				sampled_n[j] =  offsets[i,1];
			}
			// if this is the same indiviudal
			else {
				
				sampled_n[j] = offsets[i-1,1] + fmin(sampled_damage[i-1], (N-offsets[i-1,1]))
					- fmin(offsets[i-1,1], sampled_repair[i-1]);
				
			}
			// if this is the last point
			if (i == n) {
				j = j+1;
				
				sampled_n[j] = offsets[i,1] + fmin(sampled_damage[i], (N-offsets[i,1]))
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
				int dead = bernoulli_rng(fmax(0.0, fmin(1.0, prob_death)));
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


