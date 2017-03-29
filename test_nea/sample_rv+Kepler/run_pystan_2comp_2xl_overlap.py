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
modelfile = "powerlaw_2comp_2xl_overlap.stan"
nChains = 10
nIterations = 500000
nThin = 250
nJobs = -1

data_Kepler = getData.Santerne_to_dict(filename_Kepler)
data_rv = getData.NEA_to_dict(filename_rv)

Stan_data_Kepler = getData.create_Stan_input_mixture(data_Kepler, transit=1)
Stan_data_rv = getData.create_Stan_input_mixture(data_rv, transit=0)
Stan_data = combine_two_dicts(Stan_data_Kepler, Stan_data_rv)
Stan_data['K'] = 2 # fix to combine_to_dicts bug

Stan_ICs = generate.ICs_2comp_2xl_overlap( [[0.5,0.5], [0.,-1.], [1.4, 2.5], 25.], Stan_data["x"], nChains)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

with open("postSamp_rv+Kepler_2comp_2xl_overlap.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
