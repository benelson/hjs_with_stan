import csv
import generate
import getData
import pystan

filename = "TableKeplerGiants.tex"
modelfile = "powerlaw_1comp_sigmas.stan"
nChains = 10
nIterations = 500000
nThin = 250
nJobs = -1

data = getData.Santerne_to_dict(filename)
Stan_data = getData.create_Stan_input_sigmas(data, transit=1)
Stan_ICs = generate.ICs_1comp_sigmas([-1.,1.7,27.], Stan_data, nChains)

print(Stan_data)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)
postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

with open("postSamp_Kepler_1comp_sigmas.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
