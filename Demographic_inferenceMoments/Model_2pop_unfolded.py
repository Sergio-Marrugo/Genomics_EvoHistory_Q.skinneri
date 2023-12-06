import numpy
import moments
from moments import Numerics, Manips, Integration
from moments.Spectrum_mod import Spectrum


######################
####  SI  and SI+ ####
######################

def SI(params, ns):
    """
    nu1= pop size for Pop 1
    nu2=pop size for Pop 2
    T1= time of split
    O= oriented SNPs
    """
    nu1,nu2,T1,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
    
   
def SI_b(params, ns):
    """
    nu1= pop size for Pop 1
    s=proportion of Pop 1 which invaded Pop 2 (i.e. original bottleneck)
    nu2= final size of Pop 2
    T1= time of split
    """
    nu1,nu2,s,T1,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs 


def SI_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of split
    """
    nu_ae,nu1,nu2,Tae,T1,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsMv
    return fs
    
    
def SI_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh Pop 1 population after the split)
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split and start of population growth in Pop 2
    """
    nu_ae,nu1,nu2,s,Tae,T1,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]

# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def SI_NeC(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change
    
    """
    nu1,nu2,nu1b,nu2b,T1,T2,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
    fs = O*fsO+(1-O)*fsM
    return fs
    
def SI_ae_NeC(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	Tae= timing of ancestral population expansion
	T1= time of split
        
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
  

def SI_2N(params, ns):
    """
    nu1= pop size for Pop 1
	nu2=pop size for Pop 2
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under linked selection 
    """
    nu1,nu2,T1,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
 
 
def SI_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu1,nu2,s,T1,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
    

def SI_ae_2N(params, ns):
    """
    nu1= pop size for Pop 1
	Tae= time of ancestral pop expansion 
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under linked selection 
    """
    nu_ae,nu1,nu2,Tae,T1,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu1*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def SI_ae_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu_ae,nu1,nu2,s,Tae,T1,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
    
def SI_NeC_2N(params, ns):
    """
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
        T2= time of NeC
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    
    """
    nu1,nu2,nu1b,nu2b,T1,T2,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O) 
     
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
  
def SI_ae_NeC_2N(params, ns):
    """
	nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	Tae= timing of ancestral population expansion
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
     
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
	
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
 


######################
####    IM  and IM+   ####
######################

def IM(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2
	T1= time of split
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
    
def IM_b(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of split
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,s,T1,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
	


def IM_ae(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	Tae= timing of ancestral population expansion
	T1= time of split
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,Tae,T1,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def IM_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh Pop 1 population after the split)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split and start of population growth in the Baltic Sea
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu_ae,nu1,nu2,s,Tae,T1,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def IM_NeC(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def IM_ae_NeC(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	Tae= timing of ancestral population expansion
	T1= time of split
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
    
def IM_2N(params, ns):
    """
    nu1= pop size for Pop 1
    nu2=pop size for Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu1,nu2,m12,m21,T1,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def IM_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
    nu2= final size of Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under selection
    """
    nu1,nu2m12,m21,s,T1,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0,m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def IM_ae_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    Tae= time of ancestral pop expansion 
    T1= time of split
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu_ae,m12,m21,nu1,nu2,Tae,T1,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu1*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def IM_ae_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu_ae,nu1,nu2,m12,m21,s,Tae,T1,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2            
    
def IM_NeC_2N(params, ns):
    """
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
        T2= time of NeC
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)

    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
  
def IM_ae_NeC_2N(params, ns):
    """
	nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	Tae= timing of ancestral population expansion
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
     
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2


def IM_2M(params, ns):
    """
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2
    
def IM_ae_2M(params, ns):
    """
	nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	T1= time of split
	Tae= timing of ancestral population expansion
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,Tae,T1,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
	
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def IM_b_2M(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	= time of population split
	T1= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO) 
    
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def IM_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of Pop 1)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2
    
    
def IM_NeC_2M(params, ns):
    """
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
    T2= time of NeC
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))

    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))

    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2
         

def IM_ae_NeC_2M(params, ns):
    """
	nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	nu1b= pop size for Pop1 at time T2
	nu2b=pop size for Pop2 at time T2
	T1= time of split
	T2= Time at which Ne are allowed to change	
	Tae= timing of ancestral population expansion
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

    
def IM_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,T1,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3   
 	
	
def IM_ae_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    Tae= pop size of ancestral pop after size change
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae, nu1,nu2,T1,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3    

def IM_b_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu1,nu2,T1,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate(nuhrf_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate(nuhrf_func, T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def IM_ae_b_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu_ae,nu1,nu2,T1,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate(nuhrf_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate(nuhrf_func, T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def IM_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    nu1b= pop size for Pop1 at time T2
    nu2b=pop size for Pop2 at time T2
    T2= time of NeC
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def IM_ae_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3
	
######################
####   SC and SC +  ####
######################
def SC(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2 
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,T2,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def SC_ae(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact

	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polarized SNPs  
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
    
def SC_b(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of secondary contact
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,s,T1,T2,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs    

def SC_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh Pop 1 population after the split)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact and start of population growth in the Baltic Sea
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def SC_NeC(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2 
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs


	
def SC_ae_NeC(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Pop2)
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12_0= migration rate from Pop1 to Pop2
	m21_0= migration rate from Pop2 to Pop1
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

    
def SC_2N(params, ns):
    """
    nu1= pop size for Pop 1
    nu2=pop size for Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu1,nu2,m12,m21,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def SC_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
    nu2= final size of Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under selection
    """
    nu1,nu2,m12,m21,s,T1,T2,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
	
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
    

def SC_ae_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    Tae= time of ancestral pop expansion 
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu_ae,m12,m21,nu1,nu2,Tae,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu1*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1*hrf, nu1*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def SC_ae_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu_ae,nu1,nu2,m12,m21,s,Tae,T1,T2,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
    
def SC_NeC_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    m12,m21,nu1,nu2,nu1b,nu2b,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2    

def SC_ae_NeC_2N(params, ns):
    """

    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    Tae= time of ancestral pop expansion 
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu_ae,m12,m21,nu1,nu2,nu1b,nu2b,Tae,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def SC_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2 


def SC_ae_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
   """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def SC_b_2M(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,T2,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2


def SC_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of Pop 1)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

  

def SC_NeC_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2


def SC_ae_NeC_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
   """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def SC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0,0], [0, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3
	
	
def SC_ae_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    Tae= pop size of ancestral pop after size change
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae, nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3    

def SC_b_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu1,nu2,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def SC_ae_b_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def SC_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    nu1b= pop size for Pop1 at time T2
    nu2b=pop size for Pop2 at time T2
    T2= time of NeC
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def SC_ae_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3
    
######################
####   AM and AM +  ####
######################
def AM(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2 
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,T2,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def AM_ae(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
	nu1= pop size for Pop1
	nu2=pop size for Pop2
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
    
def AM_b(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of secondary contact
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,s,T1,T2,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs    

def AM_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh Pop 1 population after the split)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact and start of population growth in the Baltic Sea
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01,m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
	
def AM_NeC(params, ns):
    """
    nu1= pop size for Pop1
	nu2=pop size for Pop2 
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
	
def AM_ae_NeC(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Pop2)
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12_0= migration rate from Pop1 to Pop2
	m21_0= migration rate from Pop2 to Pop1
	m12= migration rate from Pop1 to Pop2

	m21= migration rate from Pop2 to Pop1
    O= Proportion of miss-polirized SNPs  
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs	
    
def AM_2N(params, ns):
    """
    nu1= pop size for Pop 1
    nu2=pop size for Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu1,nu2,m12,m21,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def AM_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
    nu2= final size of Pop 2
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under selection
    """
    nu1,nu2,m12,m21,s,T1,T2,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
	
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
    

def AM_ae_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    Tae= time of ancestral pop expansion 
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu_ae,m12,m21,nu1,nu2,Tae,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def AM_ae_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu_ae,nu1,nu2,m12,m21,s,Tae,T1,T2,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2
	
def AM_NeC_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    m12,m21,nu1,nu2,nu1b,nu2b,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2    

def SC_ae_NeC_2N(params, ns):
    """
    nu1= pop size for Pop 1
    m12= proportion of pop1 made up of migrants from pop 2
    m21= proportion of pop2 made up of migrants from pop 1
    Tae= time of ancestral pop expansion 
    T1= time of split
    T2= time of secondary contact
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 
    """
    nu_ae,m12,m21,nu1,nu2,nu1b,nu2b,Tae,T1,T2,hrf,Q,O = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def AM_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2 

def AM_ae_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
   """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def AM_b_2M(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,T2,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2


def AM_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of Pop 1)
	s=proportion of the Pop 1 pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

	
def AM_NeC_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)

# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2


def AM_ae_NeC_2M(params, ns):
    """
	nu1= pop size for Pop1 
	nu2= pop size for the Pop2
	nu1b= pop size for Pop1 at time of SC
	nu2b=pop size for Pop2 at time of SC
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from Pop1 to Pop2
	m21= migration rate from Pop2 to Pop1
	i1= reduction in migration rate for Islands from the Pop2 to Pop1 
	i2= reduction in migration rate for Islands from the Pop1 to Pop2
	P= proportion of the genome made up of "islands"
    O= Proportion of miss-polirized SNPs  
   """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def AM_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3	
	
def AM_ae_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    Tae= pop size of ancestral pop after size change
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3    

def AM_b_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu1,nu2,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def AM_ae_b_2M2N(params, ns):
    """
    nu_ae= pop size of ancestral pop after size change
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1,nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf,nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def AM_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    nu1b= pop size for Pop1 at time T2
    nu2b=pop size for Pop2 at time T2
    T2= time of NeC
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def AM_ae_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    O= Proportion of oriented SNPs
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3
	
######################
####   2EP + 2EP+ ####
######################
def two_ep(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	T1= time of population split
	T2= time of second migration epoch
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
    """
    nu1,nu2,T1,T2,m12_0,m21_0,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def two_ep_ae(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second migration epoch
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs
	
def two_ep_b(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of secondary contact
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
    """
    nu1,nu2,s,T1,T2,m12_0,m21_0,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def two_ep_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for the Pop 1 population after the split)
	s=proportion of the Pop 1 pop which invaded Pop 2
	nu2= final size of Pop 2
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of second epoch of migration and start of population growth in Pop 2 
	m12_0= ancient migration rate from Pop 2 to Pop1
	m21_0= ancient migration rate from Pop 1 to Pop 2 
	m12= migration rate from Pop 2 to Pop 1
	m21= migration rate from Pop 1 to Pop 2
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12_0,m21_0,m12,m21,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def two_ep_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	nu1b= pop size for Pop1 at T2
	nu2b=pop size for Pop2 at T2
	T1= time of population split
	T2= time of second migration epoch
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12_0,m21_0,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def two_ep_ae_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	nu1b= pop size for Pop1 at T2
	nu2b=pop size for Pop2 at T2	Tae= timing of ancestral population expansion
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second migration epoch
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12_0,m21_0,m12,m21,O = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def two_ep_2N(params, ns):
    """
	nu1= pop size for Pop 1 
	nu2= pop size for Pop 2
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 	
    """
    nu1,nu2,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def two_ep_ae_2N(params, ns):
    """
	nu1= pop size after ancestral expansion
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO) 
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)
    return fs2

def two_ep_b_2N(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 
    """
    nu1,nu2,s,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)	
    return fs2

def two_ep_ae_b_2N(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of Pop 1)
	s=proportion of the Pop 1 pop which invaded Pop 2
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	Tae= time of ancestral pop expansion
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate teh spectrum for loer ne regions	
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu1*hrf], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate(nuhrf_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)	 	
    return fs2

def two_ep_NeC_2N(params, ns):
    """
	nu1= pop size after pop expansion (and pop size of Pop 1)
	nu2= final size of Pop 2
	nu1b= pop size of Pop1 at NeC
	nu2b= pop szie of Pop2 at NeC
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([nu_ae], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)	
    return fs2    

def two_ep_ae_NeC_2N(params, ns):
    """
	nu1= pop size after pop expansion (and pop size of Pop 1)
	nu2= final size of Pop 2
	nu1b= pop size of Pop1 at NeC
	nu2b= pop szie of Pop2 at NeC
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)	
        Q= proportion of the genome under linked selection 
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12_0,m21_0,m12,m21,hrf,Q,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2O = moments.Spectrum(stsn2)
    fsn2O.integrate([n1*hrf], Tae, dt_fac=0.01)
    fsn2O = moments.Manips.split_1D_to_2D(fsn2O, ns[0], ns[1])
    fsn2O.integrate([nu1*hrf, nu2*hrf], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsn2O.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsnM = moments.Numerics.reverse_array(fsn2O)
    
    fs2= Q*(O*fsn2O+(1-O)*fsnM)+(1-Q)*(O*fsO+(1-O)*fsM)	
    return fs2
	
def two_ep_2M(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def two_ep_ae_2M(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def two_ep_b_2M(params, ns):
    """
    nu1= pop size for Pop 1
	s=proportion of the Pop 1 pop which invaded Pop 2
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)	
    return fs2

def two_ep_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of Pop 1)
	s=proportion of the Pop 1 pop which invaded Pop 2 (i.e. original bottleneck)
	nu2= final size of Pop 2
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO) 
    
# calculate the spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu1*s], T1, dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2

def two_ep_NeC_2M(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	nu1b= pop size for Pop1 at T2
	nu2b=pop size for Pop2 at T2
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2
    
def two_ep_ae_NeC_2M(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2
	nu1b= pop size for Pop1 at T2
	nu2b=pop size for Pop2 at T2
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of second epoch of migration
	m12= proportion of Pop 1 made up of migrants from Pop 2
	m21= proportion of Pop 2 made up of migrants from Pop 1
	m12_0= ancient migration rate from Pop 2 to Pop 1
	m21_0= ancient migration rate from Pop 1 to Pop 2
	i1= reduction in migration rate m12 for islands
	i2= reduction in migration rate m21 for islands
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12_0,m21_0,m12,m21,i1,i2,P,O = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsO = moments.Spectrum(sts)
    fsO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsO = moments.Manips.split_1D_to_2D(fsO, ns[0], ns[1])
    fsO.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fsO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fsM = moments.Numerics.reverse_array(fsO)
    
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsiO = moments.Spectrum(stsi)
    fsiO.integrate([nu_ae], Tae, dt_fac=0.01)
    fsiO = moments.Manips.split_1D_to_2D(fsiO, ns[0], ns[1])
    fsiO.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fsiO.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsiM = moments.Numerics.reverse_array(fsiO)
    
    fs2= P*(O*fsiO+(1-O)*fsiM)+(1-P)*(O*fsO+(1-O)*fsM)
    return fs2
	
def two_ep_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,T1,T2,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3	

def two_ep_ae_2M2N(params, ns):
    """
    nu_ae= pop size at ancestral expansion
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def two_ep_b_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu1,nu2,T1,T2,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2*s], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def two_ep_ae_b_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu1*s], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu1*hrf*s], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate(nuhrf_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def two_ep_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3

def two_ep_ae_NeC_2M2N(params, ns):
    """
    nu1= pop size for Pop1
    nu2=pop size for Pop2
    T1= time of split
    T2= time of secondary contact
    m12= migration rate from Pop2 to Pop1
    m21= migration rate from Pop1 to Pop2
    m12_0= ancient migration rate from Pop 2 to Pop 1
    m21_0= ancient migration rate from Pop 1 to Pop 2
    i1= reduction in migration rate for Islands from pop2 to pop1 
    i2= reduction in migration rate for Islands from pop1 to pop2
    P= proportion of the genome made up of "islands"
    hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
    Q= proportion of the genome under linked selection 	
    R= proportion of genome under both linked selection and heterogenous migration
    """
    nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,m12_0,m21_0,i1,i2,hrf,P,Q,R,O = params
	
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae, dt_fac=0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum from genomic islands
    sts2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2M = moments.Spectrum(sts2M)
    fs2M.integrate([nu_ae], Tae, dt_fac=0.01)
    fs2M = moments.Manips.split_1D_to_2D(fs2M, ns[0], ns[1])
    fs2M.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2M.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
# calculate the spectrum for lower Ne regions of the genome
    sts2N = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N = moments.Spectrum(sts2N)
    fs2N.integrate([nu1*hrf], Tae)
    fs2N = moments.Manips.split_1D_to_2D(fs2N, ns[0], ns[1])
    fs2N.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0], [m21_0, 0]]))
    fs2N.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
	
# calculate the spectrum for lower Ne and M regions of the genome
    sts2N2M = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs2N2M = moments.Spectrum(sts2N2M)
    fs2N2M.integrate([nu1*hrf], Tae)
    fs2N2M = moments.Manips.split_1D_to_2D(fs2N2M, ns[0], ns[1])
    fs2N2M.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12_0*i2], [m21_0*i1, 0]]))
    fs2N2M.integrate([nu1b*hrf, nu2b*hrf], T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=(1-P-Q-R)*fs+P*fs2M+Q*fs2N+R*fs2N2M
    fsM= moments.Numerics.reverse_array(fs2)
    fs3= O*fs2+(1-O)*fsM
    return fs3