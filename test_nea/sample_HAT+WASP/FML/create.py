import numpy as np
from collections import OrderedDict
from itertools import chain

def input_format(dict):
    params = [dict['period'], dict['half-amp'], dict['eccentricity'], dict['argperi'], dict['meanomaly']]
    params_by_planet = list(map(list, zip(*params)))
    params_chained = list(chain.from_iterable(params_by_planet))
    return params_by_planet, params_chained


def dict_keys(nPlanets=1,nOff=0,slope=False):
    assert nPlanets<10
    assert nOff<10
    param = ['p','k','e','w','M']
    param_tf = ['p', 'k', 'esinw', 'ecosw', 'wpM']
    num = ['1','2','3','4','5','6','7','8','9']
    off = 'c'
    if slope:
        acc = 's'
        
    param_keys = [] # parameters from file read-in/integration
    param_tf_IS_keys = [] # parameters for importance sampling
        
    for i in range(nPlanets):
        for j in param:
            temp = j + num[i]
            param_keys.append(temp)
        for j in param_tf:
            temp = j + num[i]
            param_tf_IS_keys.append(temp)

    if nPlanets>=2:
        param_tf_IS_keys = []
        param_tf_IS_keys.append('p1')
        param_tf_IS_keys.append('k1')
        param_tf_IS_keys.append('wpM1')
        param_tf_IS_keys.append('k2')
        param_tf_IS_keys.append('wpM2')

    if slope:
        param_keys.append(acc)
        param_tf_IS_keys.append(acc)

    for i in range(nOff):
        temp = off + num[i]
        param_keys.append(temp)
        param_tf_IS_keys.append(temp)
        
    return param_keys, param_tf_IS_keys

def posterior_samples_from_file(filename, param_keys):
    dic = OrderedDict()
    for i in range(len(param_keys)):
        dic[param_keys[i]] = []
        
    num = 0

    with open(filename) as infile:
        for line in infile:
            num += 1
            line = line.split()
            for i in range(nPlanets):
                dic[param_keys[i*5+0]].append(float(line[i*9+3]))
                dic[param_keys[i*5+1]].append(float(line[i*9+4]))
                dic[param_keys[i*5+2]].append(float(line[i*9+7]))
                dic[param_keys[i*5+3]].append(float(line[i*9+9]))
                dic[param_keys[i*5+4]].append(float(line[i*9+11]))
            for i in range(nOff):
                dic[param_keys[nPlanets*5+i]].append(float(line[nPlanets*9+3+i]))  
              
    return dic, num


def posterior_samples_from_emcee(samples, param_keys):
    num = len(np.transpose(samples))
    assert len(samples)==len(param_keys)
    dic = OrderedDict(zip(param_keys, samples))
    return dic, num


'''
def rebound_structure_from_amewM(postSamp, nPlanets, nOffs=0):
    planets = [ [] for i in range(len(postSamp['a1'])) ]
    offs = [ [] for i in range(len(postSamp['a1'])) ]
    for i in range(len(planets)):
        for j in range(nPlanets):
            n = str(j+1)
            planets[i].append({'a': postSamp['a'+n][i],\
                                'm': postSamp['m'+n][i],\
                                'e': postSamp['e'+n][i],\
                                'w': postSamp['w'+n][i] * np.pi/180.,\
                                'M': postSamp['M'+n][i] * np.pi/180.})
    for i in range(len(planets)):
        for j in range(nOffs):
            offs[i].append(postSamp['c'+str(j+1)][i])
        
    return planets, offs
'''

def obs_data_from_file(filename):
    data = []
    with open(filename) as infile:
        for line in infile:
            line = list(map(float, line.split()))
            data.append(line)
    data = np.array(data).T
    data[0] = data[0] - data[0][0]
    return {'times': data[0], 'rvs': data[1], 'errs': data[2], 'obs_indices': list(map(int, data[3])) }


def post_samp_from_file(filename):
    with open(filename) as f:
        f.readline()
        lines = f.read().splitlines()

    for i, line in enumerate(lines):
        lines[i] = list(map(float, line.split(' ')))

    return np.transpose(lines)[:-3]

