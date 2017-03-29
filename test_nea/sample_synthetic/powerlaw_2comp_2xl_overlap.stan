data {
    int<lower=1> K; // number of mixture components
    int<lower=0> N; // number of planets w/o incl measurements
    int<lower=0> Ni; // number of planets w/ incl measurements
    real x[N]; // x values for N planets w/o incl measurements
    real xi[Ni]; // x values for Ni planets w/ incl measurements
    real ind[N];
    real indi[Ni];
}
parameters {
    simplex[K] theta; // mixing proportions
    real gamma[K]; // power-law indices
    real xl[K]; // xlower
    real xu; // xupper
    real cosi[N]; // cos inclinations of all planets
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
    cosi ~ uniform(-1.,1.); // prior in cosi

    // for systems without inclination measurements
    for (n in 1:N) { //loop over planets
        sini_onethird <- (1 - cosi[n]*cosi[n])^(0.16667); // converting cosi to sini^(1/3)
        
        for (k in 1:K) { //loop over mixture components
	    if (n==1)
            {
                integral[k] <- (xupper[k]^(gamma[k]-ind[n]) - xlower[k]^(gamma[k]-ind[n]))/(gamma[k]-ind[n]);
            }
            if ((x[n]/sini_onethird > xlower[k]) && (x[n]/sini_onethird < xupper[k]))
            {
                ps[k] <- theta[k] * (x[n]/sini_onethird)^(gamma[k]-1.-ind[n])/integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
    
    for (n in 1:Ni) { //loop over planets
                
	for (k in 1:K){
	    if (n==1)
            {
                integral[k] <- (xupper[k]^(gamma[k]-indi[n]) - xlower[k]^(gamma[k]-indi[n]))/(gamma[k]-indi[n]);
            }
            if ((xi[n] > xlower[k]) && (xi[n] < xupper[k]))
            {
                ps[k] <- theta[k] * (xi[n])^(gamma[k]-1.-indi[n])/integral[k];
            }
            else
            {
                ps[k] <- 0.;
            }
        }     
        increment_log_prob( log(sum(ps)) );
    }
}
