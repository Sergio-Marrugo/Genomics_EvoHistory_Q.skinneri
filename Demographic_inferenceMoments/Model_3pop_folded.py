	"""
Three-population demographic models.
"""
import numpy as np
import moments


def out_of_Africa(params, ns, pop_ids=["Chagos", "Ningaloo", "Indo"]):
    """
    The Gutenkunst et al (2009) out-of-Africa that has been reinferred a
    number of times.

    :param params: List of parameters, in the order (nuA, TA, nuB, TB, nuEu0,
        nuEuF, nuAs0, nuAsF, TF, mAfB, mAfEu, mAfAs, mEuAs).
    :type params: list of floats
    :param ns: List of population sizes in each population, in order given
        by `pop_ids`.
    :type ns: list of ints
    :param pop_ids: List of population IDs, defaults to ["Chagos", "Ningaloo", "Indo"].
    :type pop_ids: list of strings, optional
    """
    if pop_ids is not None and len(pop_ids) != 3:
        raise ValueError("pop_ids must be a list of three population IDs")
    if len(ns) != 3:
        raise ValueError("ns must have length 3")
    (
        nuA,
        TA,
        nuB,
        TB,
        nuEu0,
        nuEuF,
        nuAs0,
        nuAsF,
        TF,
        mAfB,
        mAfEu,
        mAfAs,
        mEuAs,
    ) = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    # integrate modern human branch with expansion
    fs.integrate([nuA], TA)
    # split into African and Eurasian branches and integrate
    fs = fs.split(0, ns[0], ns[1] + ns[2])
    mig_mat = [[0, mAfB], [mAfB, 0]]
    fs.integrate([nuA, nuB], TB, m=mig_mat)
    # split further into Indo and Ningaloo
    fs = fs.split(1, ns[1], ns[2])
    nu_func = lambda t: [
        nuA,
        nuEu0 * np.exp(np.log(nuEuF / nuEu0) * t / TF),
        nuAs0 * np.exp(np.log(nuAsF / nuAs0) * t / TF),
    ]
    mig_mat = [[0, mAfEu, mAfAs], [mAfEu, 0, mEuAs], [mAfAs, mEuAs, 0]]
    fs.integrate(nu_func, TF, m=mig_mat)
    fs.pop_ids = pop_ids
    return fs

def split_symmig_adjacent(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
    in between populations 1 and 3, which do not come in to contact with one another.
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    nu1: Size of population 1 after split (Chagos)
    nuA: Size of population (2,3) after split from Chagos.
    nu2: Size of population 2 after split (Ningaloo)
    nu3: Size of population 3 after split (Indo)
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m): Chagos-Ningaloo
    m2: Migration rate between populations 2 and 3: Ningaloo-Indo
    
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 (Ningaloo) and 3 (Indo) (in units of 2*Na generations).
   """
    #9 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs.integrate([nu1, nuA], T1, dt_fac=0.01, m=np.array([[0, mA], [mA, 0]]))
    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    mig_mat = [[0, m1, 0], [m1, 0, m2], [0, m2, 0]]
    fs.integrate([nu1, nu2, nu3], T2, dt_fac=0.01, m=mig_mat)
   
    fs.pop_ids=["1", "2", "3"]
    return fs


def sim_split_sym_mig_all(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
    m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
	#6 parameters
    nu1, nu2, nu3, m1, m2, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    mig_mat = [[0, m1, 0], [m1, 0, m2], [0, m2, 0]]
    fs.integrate([nu1, nu2, nu3], T1, dt_fac=0.01, m=mig_mat)
    fs.pop_ids=["1", "2", "3"]
    return fs

