import numpy as np
import csv
import generate
import getData
import pystan

def combine_two_dicts(x, y):
    z = {}
    for k in set(x).intersection(y):
        if k is "K":
            z[k] = x[k]
        elif isinstance(x[k], int):
            z[k] = x[k] + y[k]
        else:
            z[k] = np.array(list(x[k]) + list(y[k]))
    return z

filename_Kepler = "TableKeplerGiants.tex"
filename_rv = "planets_rvall.in"
filename_HAT = "planets_HAT.in"
filename_WASP = "planets_WASP.in"
modelfile = "powerlaw_2comp_2xl_overlap.stan"
nChains = 20
nIterations = 300000
nThin = 300
nJobs = 1

data_Kepler = getData.Santerne_to_dict(filename_Kepler)
data_rv = getData.exoplanetsORG_to_dict(filename_rv)
data_HAT = getData.exoplanetsORG_to_dict(filename_HAT)
data_WASP = getData.exoplanetsORG_to_dict(filename_WASP)

Stan_data_Kepler = getData.create_Stan_input_mixture(data_Kepler, transit=1)
Stan_data_rv = getData.create_Stan_input_mixture(data_rv, transit=0)
Stan_data_HAT = getData.create_Stan_input_mixture(data_HAT, transit=1)
Stan_data_WASP = getData.create_Stan_input_mixture(data_WASP, transit=1)

data_tmp = combine_two_dicts(Stan_data_Kepler, Stan_data_rv)
data_tmp = combine_two_dicts(data_tmp, Stan_data_HAT)
Stan_data = combine_two_dicts(data_tmp, Stan_data_WASP)

Stan_ICs = generate.ICs_2comp_2xl_overlap([[0.25,0.75], [-1.,-2.], [1.05, 2.25], 25.], Stan_data["x"], nChains)
print(Stan_data)
print(Stan_ICs)

fit = pystan.stan(file=modelfile, data=Stan_data, iter=nIterations, chains=nChains, thin=nThin, init=Stan_ICs, n_jobs=nJobs)

postSamp = getData.stan_output_to_posterior_samples_2comp_2xl_overlap(fit.extract())

with open("postSamp_all_2comp_2xl_overlap.txt", "w") as f:
    writer = csv.writer(f, delimiter=' ')
    for item in postSamp:
        writer.writerow(list(item))
