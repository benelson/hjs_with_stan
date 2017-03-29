import csv
import generate
import getData
import pystan

filename = "TableKeplerGiants.tex"
modelfile = "powerlaw_2comp_2xl_overlap.stan"
nChains = 10
nIterations = 8000000
nThin = 4000
nJobs = -1

data = getData.Santerne_to_dict(filename)
Stan_data = getData.create_Stan_input_mixture(data)
Stan_ICs = generate.ICs_2comp_2xl_overlap([ [0.5,0.5], [1.,1.], [1.5, 2.0], 500. ], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

print(postSamp)

with open("postSamp_Kepler_2comp_2xl_overlap.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
