data {
    int<lower=1> K; // number of mixture components
    int<lower=0> Ni; // number of planets w/ incl measurements
    real xi[Ni]; // x values for Ni planets w/ incl measurements
}
parameters {
    simplex[K] theta; // mixing proportions
    real gamma[K]; // power-law indices
    real xl[K]; // xlower
    real xu; // xupper
}
model {
    real ps[K];
    real xlower[K];
    real xupper[K];
    real integral[K];
    real sini_onethird;
    
    xlower[1] <- xl[1]; // xl for first component
    xlower[2] <- xl[2]; // xl for second component
    xupper[1] <- xu; // xu for first component
    xupper[2] <- xu; // xu for second component
    
    gamma ~ uniform(-10.,10.); // prior on gamma
    xl ~ uniform(0., 25.);
    xu ~ uniform(xl[2], 25.); // prior on xu

    for (n in 1:Ni) { //loop over planets  
    	for (k in 1:K){            
	    if (n == 1){
	        integral[k] <- (xupper[k]^(gamma[k]-1.) - xlower[k]^(gamma[k]-1.))/(gamma[k]-1.);
	    }
            if ((xi[n] > xlower[k]) && (xi[n] < xupper[k]))
            {
                ps[k] <- theta[k] * (xi[n])^(gamma[k]-2.)/integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
}
