data {
    int<lower=0> Ni; // number of planets w/ incl measurements
    real xi[Ni]; // x values of planets w/ incl mesaurements
}
parameters {
    real gamma; // power law index
    real<lower=0.> xl; // lower truncation value of x
    real<upper=999.> xu; // upper truncation value of x
}
model {
    real sum_log_x;
    real sini_onethird;

    sum_log_x <- 0.0;
    
    for (i in 1:Ni){
                
        if ((xi[i] > xl) && (xi[i] < xu))
        {
            sum_log_x <- sum_log_x + log((gamma-1.) / (xu^(gamma-1.) - xl^(gamma-1.))) + (gamma-2.) * log(xi[i]);
        }
        else
        {
            sum_log_x <- -1/0.;
        } 
        
    }
    
    // specify priors in population papers
    gamma ~ uniform(-10.,10.);
    xl ~ uniform(0.,xu);
    xu ~ uniform(xl,25.);
    
    increment_log_prob( sum_log_x );
}
