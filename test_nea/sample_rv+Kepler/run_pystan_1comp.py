import numpy as np
import csv
import generate
import getData
import pystan

def combine_two_dicts(x, y):
    z = {}
    for k in set(x).intersection(y):
        if isinstance(x[k], int):
            z[k] = x[k] + y[k]
        else:
            z[k] = np.array(list(x[k]) + list(y[k]))
    return z

filename_Kepler = "TableKeplerGiants.tex"
filename_rv = "planets_RV.csv"
modelfile = "powerlaw_1comp.stan"
nChains = 10
nIterations = 400000
nThin = 200
nJobs = -1

data_Kepler = getData.Santerne_to_dict(filename_Kepler)
data_rv = getData.NEA_to_dict(filename_rv)

Stan_data_Kepler = getData.create_Stan_input(data_Kepler, transit=1)
Stan_data_rv = getData.create_Stan_input(data_rv, transit=0)
Stan_data = combine_two_dicts(Stan_data_Kepler, Stan_data_rv)

Stan_ICs = generate.ICs_1comp([-1.,1.5,27.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_1comp(fit.extract())

with open("postSamp_rv+Kepler_1comp.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
