import numpy as np
import re
from math import isnan

MJtoMsun = 0.000954265748
RJtoAU = 0.000477894503
RstarToAU = 0.00465047

def has_numbers(inputString):
    '''Checks if inputString contains any numerical value'''
    return any(char.isdigit() for char in inputString)

def convert(inputList):
    '''Attempts to convert elements of list to floats. If not, keep as string.'''
    outputList = []
    for x in inputList:
        try:
            x = float(x)
        except ValueError:
            x = np.nan
        outputList.append(x)
    return outputList
    
def clean_array(array):
    '''Cleans the array from a file generated from exoplanets.org into a more readable format'''
    tmp = [ i for i in array if has_numbers(i) ]
    tmp = np.transpose( [ [x for x in i.split(',')] for i in tmp] )
    return [ np.array(convert(x)) for x in tmp ]

def create_array_from_file(filename):
    '''Creates a cleaned array from input file'''
    inds = [1, 3,4,5,6, 7,8,9,10,  11,12,13,14, 15,16,17,18, 19,20,21,22, 23,24,25,26, 27,28,29,30, 32,33,34,35, \
             37,38,39,40, 42,43,44,45 ] 
    array = []
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"):
                array.append(line.split(","))
    cleaned = np.transpose(array[1:])
    cleaned = [ np.array(convert(cleaned[i])) for i in inds ]
    return cleaned
#    return clean_array(array)

def NEA_to_dict(filename):
    tmp = create_array_from_file(filename)
    for i in range(42)[4::4]:
        tmp[i] = (tmp[i-2] - tmp[i-1])/2.
    tmp[20] = np.array([ 0.1 if isnan(i) else i for i in tmp[20] ]) #rad_sigma   
    tmp[21] = np.array([ 1.2 if isnan(i) else i for i in tmp[21] ]) #rad                                                                                                                                  
    tmp[24] = np.array([ 0.1 if isnan(i) else i for i in tmp[24] ]) #rad_sigma                                                                                                                                  
    keys = [ "names",
             "per", "per_upper", "per_lower", "per_sigma", \
             "sma", "sma_upper", "sma_lower", "sma_sigma", \
             "ecc", "ecc_upper", "ecc_lower", "ecc_sigma", \
             "inc", "inc_upper", "inc_lower", "inc_sigma", \
             "Mpl", "Mpl_upper", "Mpl_lower", "Mpl_sigma",\
             "rad", "rad_upper", "rad_lower", "rad_sigma", \
             "tem", "tem_upper", "tem_lower", "tem_sigma",\
             "Mstar", "Mstar_upper", "Mstar_lower", "Mstar_sigma", \
             "amp", "amp_upper", "amp_lower", "amp_sigma", \
             "feh", "feh_upper", "feh_lower", "feh_sigma", \
             "vsi", "vsi_upper", "vsi_lower", "vsi_sigma"]
    
    return dict(zip(keys, tmp))


def clean_Santerne_array(array):
    '''Cleans the array from a Santerne et al. 2016 tex file into a more readable format'''
    mean = []
    err = []
    for i,x in enumerate(array):
        mean.append([])
        err.append([])
        for y in x.rsplit("&"):
            z = y.rsplit("$")
            if z[0].isspace():
                mean[i].append(z[1].rsplit("<")[1])
                err[i].append(np.nan)
            else: 
                mean[i].append(z[0])
                if len(z)>1:
                    if not z[2].isspace():
                        err[i].append(float(z[2]))
                    elif not z[1].isspace():
                        errors = list(map(float, re.findall("\W\d+\.?\d+", z[1])))
                        errors = list(map(abs,errors))
                        powers = [float(x.replace("{","").replace("}","").replace("^","E").replace("0","")) for x in re.findall("10\^\{\-\d+\}", z[1])]
                        if powers:
                            errors = [ errors[j] * powers[j] for j in range(len(powers)) ]
                        err[i].append(sum(errors)/len(errors))
                else:
                    err[i].append(np.nan)


    mean = np.transpose(mean)
    err = np.transpose(err)
    tmp = np.vstack((mean, err))
    
    return [ np.array(convert(x)) if i!=13 else x for i,x in enumerate(tmp) ]

    
def create_array_from_Santerne_file(filename):
    '''Creates a cleaned array from input file'''

    with open(filename) as f:
        array = f.read().splitlines()
        return clean_Santerne_array(array)

def Santerne_to_dict(filename):
    tmp = create_array_from_Santerne_file(filename)
    keys = [ "KICID", "KOIID", "KeplerID",\
             "per", "aRstar", "RpRstar", "rad", "Mpl", \
             "Teff","FeH","Rstar","Mstar","Kp", \
             "status", "references",
             "fill1", "fill2", "fill3",\
             "per_sigma", "aRstar_sigma", "RpRstar_sigma", "rad_sigma", "Mpl_sigma",\
             "Teff_sigma", "FeH_sigma", "Rstar_sigma", "Mstar_sigma", "Kp_sigma",\
             "fill4","fill5"]
    tmp = dict(zip(keys, tmp))
    per_inds = [ i for i,x in enumerate(tmp["per"]) if x<30. ]
    
    incs = np.array([90. for i in range(len(tmp["per"]))])
#    smas = RstarToAU*tmp["aRstar"]*tmp["Rstar"]
    pers = tmp["per"]
    rads = tmp["rad"]
    Mpls = tmp["Mpl"]
    Mstars = tmp["Mstar"]
    statuses = tmp["status"]

    persigma = np.array([0.000001 for i in pers])
    radsigma = tmp["rad_sigma"]
    Mplsigma = tmp["Mpl_sigma"]
    Mstarsigma = tmp["Mstar_sigma"]
    Mstarsigma = np.array([ 0.1 if isnan(i) or i==0. else i for i in Mstarsigma ])

    teff = tmp["Teff"]
    feh = tmp["FeH"]

    return {"inc": incs[per_inds],\
            "per": pers[per_inds],\
            "rad": rads[per_inds],\
            "Mpl": Mpls[per_inds],\
            "status": statuses[per_inds],\
            "per_sigma": persigma[per_inds],\
            "rad_sigma": radsigma[per_inds],\
            "Mpl_sigma": Mplsigma[per_inds],\
            "teff": teff[per_inds],\
            "feh": feh[per_inds]}

def inc_to_sinionethird(li):
    tmp = []
    for i in li:
        tmp.append(np.sin(i*np.pi/180.)**(0.33333))
    return np.array(tmp)

def create_Stan_input(dic, transit=0):
    ilist = inc_to_sinionethird(dic["inc"])
    xlist = (dic["per"]/365.242)**0.66667 * ( 0.462 * (dic["Mpl"] * MJtoMsun)**0.3333 ) / (dic["rad"] * RJtoAU)
    
    x = []
    xi = []
    
    for j,i in enumerate(ilist):
        if isnan(i):
            x.append(xlist[j])
        else:
            xi.append(xlist[j])
            
    return {'N': len(x),
            'Ni': len(xi),
            'x': x,
            'xi': xi,
            'ind': [ transit for i in range(len(x)) ],
            'indi': [ transit for i in range(len(xi)) ] }


def create_Stan_input_alltransit(dic):
    ilist = inc_to_sinionethird(dic["inc"])
    xlist = (dic["per"]/365.242)**0.66667 * ( 0.462 * (dic["Mpl"] * MJtoMsun)**0.3333 ) / (dic["rad"] * RJtoAU)
    
    x = []
    xi = []

    for j,i in enumerate(ilist):
        if isnan(i):
            x.append(xlist[j])
        else:
            xi.append(xlist[j])

    return {'N': len(x),
            'Ni': len(xi),
            'x': x,
            'xi': xi }

def create_Stan_input_mixture(dic, transit=0):
    ilist = inc_to_sinionethird(dic["inc"])
    xlist = (dic["per"]/365.242)**0.66667 * ( 0.462 * (dic["Mpl"] * MJtoMsun)**0.3333 ) / (dic["rad"] * RJtoAU)

    x = []
    xi = []

    for j,i in enumerate(ilist):
        if isnan(i):
            x.append(xlist[j])
        else:
            xi.append(xlist[j])
    return {'K': 2,
            'N': len(x),
            'Ni': len(xi),
            'x': x,
            'xi': xi,
            'ind': [ transit for i in range(len(x)) ],
            'indi': [ transit for i in range(len(xi)) ] }


def create_Stan_input_mixture_alltransit(dic):
    ilist = inc_to_sinionethird(dic["inc"])
    xlist = (dic["per"]/365.242)**0.66667 * ( 0.462 * (dic["Mpl"] * MJtoMsun)**0.3333 ) / (dic["rad"] * RJtoAU)

    xi = []

    for j,i in enumerate(ilist):
        xi.append(xlist[j])
    return {'K': 2,
            'Ni': len(xi),
            'xi': xi }


def create_Stan_input_sigmas(dic, transit=0):
    ilist = inc_to_sinionethird(dic["inc"])

    per, rad, Mpl = [ [] for i in range(3) ]
    per_sigma, rad_sigma, Mpl_sigma = [ [] for i in range(3) ]
    peri, radi, Mpli = [ [] for i in range(3) ] 
    per_sigmai, rad_sigmai, Mpl_sigmai = [ [] for i in range(3) ]

    for j,i in enumerate(ilist):
        if isnan(i):
            per.append(dic["per"][j])
            rad.append(dic["rad"][j])
            Mpl.append(dic["Mpl"][j])
#            Mstar.append(dic["Mstar"][j])
            per_sigma.append(dic["per_sigma"][j])
            rad_sigma.append(dic["rad_sigma"][j])
            Mpl_sigma.append(dic["Mpl_sigma"][j])
#            Mstar_sigma.append(dic["Mstar_sigma"][j])
        else:
            peri.append(dic["per"][j])
            radi.append(dic["rad"][j])
            Mpli.append(dic["Mpl"][j])
#            Mstari.append(dic["Mstar"][j])
            per_sigmai.append(dic["per_sigma"][j])
            rad_sigmai.append(dic["rad_sigma"][j])
            Mpl_sigmai.append(dic["Mpl_sigma"][j])
#            Mstar_sigmai.append(dic["Mstar_sigma"][j])

    return {'N': len(per),
            'Ni': len(peri),
            'per': per, 'per_sigma': per_sigma,
            'peri': peri, 'per_sigmai': per_sigmai,
            'rad': rad, 'rad_sigma': rad_sigma,
            'radi': radi, 'rad_sigmai': rad_sigmai,
            'Mpl': Mpl, 'Mpl_sigma': Mpl_sigma,
            'Mpli': Mpli,'Mpl_sigmai': Mpl_sigmai,
            'ind': [ transit for i in range(len(per)) ],
            'indi': [ transit for i in range(len(peri)) ]}
#            'Mstar': Mstar, 'Mstar_sigma': Mstar_sigma,
#            'Mstari': Mpli,'Mstar_sigmai': Mstar_sigmai }


def create_Stan_input_mixture_sigmas(dic, transit=0):
    ilist = inc_to_sinionethird(dic["inc"])

    per, rad, Mpl = [ [] for i in range(3) ]
    per_sigma, rad_sigma, Mpl_sigma = [ [] for i in range(3) ]
    peri, radi, Mpli = [ [] for i in range(3) ]
    per_sigmai, rad_sigmai, Mpl_sigmai = [ [] for i in range(3) ]

    for j,i in enumerate(ilist):
        if isnan(i):
            per.append(dic["per"][j])
            rad.append(dic["rad"][j])
            Mpl.append(dic["Mpl"][j])
#            Mstar.append(dic["Mstar"][j])                                                                                                                                         
            per_sigma.append(dic["per_sigma"][j])
            rad_sigma.append(dic["rad_sigma"][j])
            Mpl_sigma.append(dic["Mpl_sigma"][j])
#            Mstar_sigma.append(dic["Mstar_sigma"][j])                                                                                                                             
        else:
            peri.append(dic["per"][j])
            radi.append(dic["rad"][j])
            Mpli.append(dic["Mpl"][j])
#            Mstari.append(dic["Mstar"][j])                                                                                                                                        
            per_sigmai.append(dic["per_sigma"][j])
            rad_sigmai.append(dic["rad_sigma"][j])
            Mpl_sigmai.append(dic["Mpl_sigma"][j])
#            Mstar_sigmai.append(dic["Mstar_sigma"][j])                                                                                                                            

    return {'K': 2,
            'N': len(per),
            'Ni': len(peri),
            'per': per, 'per_sigma': per_sigma,
            'peri': peri, 'per_sigmai': per_sigmai,
            'rad': rad, 'rad_sigma': rad_sigma,
            'radi': radi, 'rad_sigmai': rad_sigmai,
            'Mpl': Mpl, 'Mpl_sigma': Mpl_sigma,
            'Mpli': Mpli,'Mpl_sigmai': Mpl_sigmai,
            'ind': [ transit for i in range(len(per)) ],
            'indi': [ transit for i in range(len(peri)) ] }
#            'Mstar': Mstar, 'Mstar_sigma': Mstar_sigma,                                                                                                                           
#            'Mstari': Mpli,'Mstar_sigmai': Mstar_sigmai }    


def create_Stan_input_mixture_sigmas_alltransit(dic):
    ilist = inc_to_sinionethird(dic["inc"])

    peri, radi, Mpli = [ [] for i in range(3) ]
    per_sigmai, rad_sigmai, Mpl_sigmai = [ [] for i in range(3) ]

    for j,i in enumerate(ilist):
        peri.append(dic["per"][j])
        radi.append(dic["rad"][j])
        Mpli.append(dic["Mpl"][j])
#       Mstari.append(dic["Mstar"][j])                                                                                                                                        
        per_sigmai.append(dic["per_sigma"][j])
        rad_sigmai.append(dic["rad_sigma"][j])
        Mpl_sigmai.append(dic["Mpl_sigma"][j])
#       Mstar_sigmai.append(dic["Mstar_sigma"][j])                                                                                                                            

    return {'K': 2,
            'Ni': len(peri),
            'peri': peri, 'per_sigmai': per_sigmai,
            'radi': radi, 'rad_sigmai': rad_sigmai,
            'Mpli': Mpli,'Mpl_sigmai': Mpl_sigmai }
#            'Mstar': Mstar, 'Mstar_sigma': Mstar_sigma,                                                                                                                           
#            'Mstari': Mpli,'Mstar_sigmai': Mstar_sigmai }     


def stan_output_to_posterior_samples_1comp(dic):
    return np.array([ dic['gamma'], dic['xl'], dic['xu'] ]).T

def stan_output_to_posterior_samples_2comp_1xl(dic):
    return np.array([ dic['theta'].T[0], dic['theta'].T[1], dic['gamma'].T[0], dic['gamma'].T[1], dic['xl'], dic['xu'] ]).T

def stan_output_to_posterior_samples_2comp_2xl(dic):
    return np.array([ dic['theta'].T[0], dic['theta'].T[1], dic['gamma'].T[0], dic['gamma'].T[1], dic['xl'], dic['xt'], dic['xu'] ]).T

def stan_output_to_posterior_samples_2comp_overlap(dic):
    return np.array([ dic['theta'].T[0], dic['theta'].T[1], dic['gamma'].T[0], dic['gamma'].T[1], dic['xu'] ]).T

def stan_output_to_posterior_samples_2comp_2xl_overlap(dic):
    return np.array([ dic['theta'].T[0], dic['theta'].T[1], dic['gamma'].T[0], dic['gamma'].T[1], dic['xl'].T[0], dic['xl'].T[1], dic['xu'] ]).T
