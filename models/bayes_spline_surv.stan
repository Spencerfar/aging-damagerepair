data {
    //dimensions
    int<lower=1> n_censored;     // number of observations
    int<lower=1> n_dead;     // number of observations
    
    //testing data
    int<lower=1> Ntest;
    
    //spline input
    int num_basis;
	real<lower = 0> spline_alpha;
	matrix[n_dead, num_basis] isplines_lower;
	matrix[n_dead, num_basis] isplines_upper;
	matrix[n_censored, num_basis] isplines_right;
	matrix[Ntest, num_basis] isplines_test;
    
}
transformed data{
    vector[num_basis] alpha_hyperparam;
    alpha_hyperparam = rep_vector(spline_alpha, num_basis);
}
parameters {
	// intercept
    real gamma_int;  
    
    //spline parameters
	simplex[num_basis] spline_gamma;
}
transformed parameters {
	
	vector[n_censored] cum_hazard;
	vector[n_dead] cum_hazard_lower;
	vector[n_dead] cum_hazard_upper;
	
	cum_hazard = isplines_right*spline_gamma * exp(gamma_int);
	
	cum_hazard_lower = isplines_lower*spline_gamma * exp(gamma_int);
	cum_hazard_upper = isplines_upper*spline_gamma * exp(gamma_int);
    
}
model {
	
	//intercept prior
	target += normal_lpdf(gamma_int | 0, 10);
	
	//spline priors
	target += dirichlet_lpdf(spline_gamma | alpha_hyperparam);
	
	
	//likelihood 
	target += -cum_hazard;
	target += log(exp(-cum_hazard_lower) - exp(-cum_hazard_upper));
    
}
generated quantities {
	
	vector[Ntest + 1] survival;
	survival[1] = 1;

	survival[2:Ntest+1] = exp(-isplines_test*spline_gamma * exp(gamma_int));

}

