import csv
import generate
import getData
import pystan

filename = "planets_RV.csv"
modelfile = "powerlaw_2comp_2xl_overlap.stan"
nChains = 10
nIterations = 500000
nThin = 250
nJobs = -1

data = getData.NEA_to_dict(filename)
Stan_data = getData.create_Stan_input_mixture(data)
Stan_ICs = generate.ICs_2comp_2xl_overlap([ [0.5,0.5], [0.,-1.], [1.4, 2.5], 27. ], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

print(postSamp)

with open("postSamp_rvall_2comp_2xl_overlap.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
