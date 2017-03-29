data {
    int<lower=0> N; // number of planets w/o incl measurements
    int<lower=0> Ni; // number of planets w/ incl measurements
    real x[N]; // x values of planets w/o incl measurements
    real xi[Ni]; // x values of planets w/ incl mesaurements
    real ind[N];
    real indi[Ni];
}
parameters {
    real gamma; // power law index
    real<lower=0.> xl; // lower truncation value of x
    real<upper=30.> xu; // upper truncation value of x
    real cosi[N]; // cos inclinations of all planets
}
model {
    real sum_log_x;
    real sini_onethird;

    sum_log_x <- 0.0;
    
    cosi ~ uniform(-1.,1.); // prior in cosi
    
    // for systems without inclination measurements
    for (i in 1:N){
        sini_onethird <- (1 - cosi[i]*cosi[i])^(0.16667);
        if ((x[i]/sini_onethird > xl) && (x[i]/sini_onethird < xu))
        {
            sum_log_x <- sum_log_x + log((gamma-ind[i]) / (xu^(gamma-ind[i]) - xl^(gamma-ind[i]))) + (gamma-1.-ind[i]) * log(x[i]/sini_onethird);
        }
        else
        {
            sum_log_x <- -1/0.;
        } 
    }
    
    // for systems with inclination measurements
    for (i in 1:Ni){
                
        if ((xi[i] > xl) && (xi[i] < xu))
        {
            sum_log_x <- sum_log_x + log((gamma-indi[i]) / (xu^(gamma-indi[i]) - xl^(gamma-indi[i]))) + (gamma-1.-indi[i]) * log(xi[i]);
        }
        else
        {
            sum_log_x <- -1/0.;
        } 
        
    }
    
    // specify priors in population papers
    gamma ~ uniform(-10.,10.);
    xl ~ uniform(0.,30.);
    xu ~ uniform(xl,30.);
    
    increment_log_prob( sum_log_x );
}
