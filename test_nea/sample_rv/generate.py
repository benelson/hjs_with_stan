import numpy as np

def ICs_1comp(IC, xlist, nChains):
    gamma = IC[0]
    xl = IC[1]
    xu = IC[2]

    ICs = []
    for i in range(nChains):
        ICs.append({ "gamma": gamma, "xl": xl, "xu": xu, "cosi": [0. for i in range(len(xlist))] })
    return ICs

def ICs_1comp_alltransit(IC, nChains):
    gamma = IC[0]
    xl = IC[1]
    xu = IC[2]

    ICs = []
    for i in range(nChains):
        ICs.append({ "gamma": gamma, "xl": xl, "xu": xu })
    return ICs

def ICs_2comp_2xl_overlap(IC, xlist, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xl = [ IC[2] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xu = [ IC[3] for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(len(xlist))] })
    return ICs



def ICs_1comp_sigmas(IC, data, nChains):
    gamma = IC[0]
    xl = IC[1]
    xu = IC[2]

    per = data["per"]
    log_rad = list(np.log10(data["rad"]))
    log_Mpl = list(np.log10(data["Mpl"]))
#    Mstar = data["Mstar"]

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))
#    Mstari = data["Mstari"]

    N = len(per)

    ICs = []
    for i in range(nChains):
        ICs.append({ "gamma": gamma, "xl": xl, "xu": xu, "cosi": [0. for i in range(N)],\
                     "per_true": per, "log_rad_true": log_rad, "log_Mpl_true": log_Mpl ,\
                     "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs


def ICs_2comp_2xl_overlap_sigmas(IC, data, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xl = [ IC[2] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xu = [ IC[3]*1.05 for i in range(nChains) ]

    per = data["per"]
    log_rad = list(np.log10(data["rad"]))
    log_Mpl = list(np.log10(data["Mpl"]))

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))

    N = len(per)

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i], "cosi": [0. for i in range(N)],\
                         "per_true": per, "log_rad_true": log_rad, "log_Mpl_true": log_Mpl ,\
                         "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs


def ICs_2comp_2xl_overlap_alltransit(IC, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] for i in range(nChains) ]
    xl = [ IC[2] for i in range(nChains) ]
    xu = [ IC[3] for i in range(nChains) ]

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i] })
    return ICs


def ICs_2comp_2xl_overlap_sigmas_alltransit(IC, data, nChains):
    K = len(IC[0])
    theta = [ IC[0] for i in range(nChains) ]
    gamma = [ IC[1] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xl = [ IC[2] + np.random.randn(2)*0.01 for i in range(nChains) ]
    xu = [ IC[3]*1.05 for i in range(nChains) ]

    peri = data["peri"]
    log_radi = list(np.log10(data["radi"]))
    log_Mpli = list(np.log10(data["Mpli"]))

    ICs = []
    for i in range(nChains):
        ICs.append({ "theta": theta[i], "gamma": gamma[i], "xl": xl[i], "xu": xu[i],\
                         "peri_true": peri, "log_radi_true": log_radi, "log_Mpli_true": log_Mpli })
    return ICs
