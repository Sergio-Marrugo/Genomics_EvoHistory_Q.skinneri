# Library importation
import matplotlib
matplotlib.use('PDF')
import sys
import os
import getopt
import numpy
from numpy import array
import moments
import pylab
import matplotlib.pyplot as plt
from datetime import datetime
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D,Plotting
import Model_2pop_unfolded
import Model_2pop_folded
import Model_3pop_folded
import Optimize_Functions


#Help function
def usage():
    """ Function for help """
    print("# This script allow you to test different demographic models on your genomic data using moments on python 3 \n"+
    "# It uses the optimisation routine propose by Portick et al, 2017 for dadi and adapted to moments by Momigliano et al., 2020 \n\n"+
    "# This is an exemple of the most complete command line :\n"+
    "# python Script_moments_model_optimisaton.py -f SNP_file -m IM_b -s True -r 1 -z -n 6  False \n\n"+
    "# that will performed the optimisaton for the IM model with a folded SFS (-u False) with singletons masked (-z) and masked centre (-s True) obtained from the file named SNP_file.data \n\n"+
    "# -h --help : Display the help you are looking at.\n"+
    "# -f SNPs input data file with .data as extention \n"+
    "# -m --model : model to test, which are found in Model_2pop and comprise five basic (SI, IM, AM, SC, 2EP) and their extention:\n"+
    "# SC_b for SC with bottleneck SC_ae for ancestral expension SC_2N for Hill-Robertson effect SC_2M for heterogeneous gene flow and combinaton of thereof  \n"+
    "# For more information on models see paper by  Momigliano et al. 2020\n"+
    "# -r --replicat : replicat number \n" +
    "# -z : mask the singletons (no argument to provide) \n"+
    "# -n --rounds:  number of rounds of optimization it must be 2 to 6 \n"+
    "# -s --maskmid: mask center of the spectrum to avoid paralogs \n"+
    "# -o --orientation: unfolded input SFS \n"+
    "########################## Enjoy ###########################")
    return()
    

def takearg(argv):
   fs_fil_name = ''
   model=''
   replicat='1'
   rounds='3'
   masked=False
   maskmid=False
   orientation=''
   try:
      opts, args = getopt.getopt(argv,"hf:m:r:n:z:s:o",["fs_fil_name=","model=","replicat=","rounds=","masked","maskdmid", "orientation="])
   except getopt.GetoptError:
    print("ERROR : INCORRECT INPUT PROVIDED.\nDescription of the script usage bellow: \n\n")
    print(usage())
    sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-f", "--fs_fil_name"):
         fs_fil_name = arg
      elif opt in ("-m", "--model"):
         model = arg
      elif opt in ("-r", "--replicat"):
         replicat = arg
      elif opt in ("-n", "--rounds"):
         rounds = arg
      elif opt in ("-z", "--masked"):
         masked = True
      elif opt in ("-s", "--maskmid"):
         maskmid = True
      elif opt in ("-o", "--orientation"):
         orientation = arg
   print ('Input file is ', fs_fil_name)
   print ('Mask singletons? ', masked)
   print ('Input SFS is  ', orientation)
   print ('Mask central entry? ', maskmid)
   print ('Replicate of ', model,' number ',replicat, 'nrounds', rounds)

   return(fs_fil_name, model,replicat, rounds,masked,maskmid,orientation)
if __name__ == "__takearg__":
   takearg(sys.argv[1:])

###########################
###      actual script         #####
###########################
# Load parameters
fs_fil_name, model,replicat, rounds, masked, maskmid, orientation = takearg(sys.argv[1:]) 

print ("nrounds",rounds,'\n')
# Load the data
#Create python dictionary from snps file

input_file=fs_fil_name+'.sfs'
data = Spectrum.from_file(input_file) 
ns = data.sample_sizes
numpy.set_printoptions(precision=3)   

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
#**************
#projection sizes, in ALLELES not individuals

#fig = pylab.figure(1)
#moments.Plotting.plot_single_2d_sfs(data, vmin=0.1)
#fig.savefig(fs_fil_name+'.svg')



meta=array([2,2,2])
maskc=numpy.divide(ns, meta)
maskc=maskc.astype(int)

if masked and not maskmid:
    data.mask[1,0] = True
    data.mask[0,1] = True
    data.mask[ns[0]-1,ns[1]] = True
    data.mask[ns[0],ns[1]-1] = True

if maskmid:
    data.mask[1,0] = True
    data.mask[0,1] = True
    data.mask[ns[0]-1,ns[1]] = True
    data.mask[ns[0],ns[1]-1] = True
    data.mask[maskc[0],maskc[1]] = True


print ("nrounds",rounds,'\n')

#Optimisaton set up

prefix = fs_fil_name+"_"+replicat

if rounds == '2':
    reps = [30,30] 
    maxiters = [20,30]
    folds = [3,2]
    
if rounds == '3':
    reps = [10,10,10] 
    maxiters = [10,20,30]
    folds = [3,2,1]
    
if rounds == '4':
    reps = [30,20,20,20]
    maxiters = [10,20,20,30]
    folds = [3,2,2,1]

if rounds == '5':
    reps = [30,20,20,20,30]
    maxiters = [10,20,20,20,30]
    folds = [3,2,2,2,1]

if rounds == '6':
    reps = [30,20,20,20,20,30]
    maxiters = [10,20,20,20,30,30]
    folds = [3,2,2,2,2,1]

if rounds == '7':
    reps = [30,20,20,20,20,20,30]
    maxiters = [20,20,20,20,20,30,30]
    folds = [3,2,2,2,2,2,1]

if rounds == '8':
    reps = [30,30,20,20,20,20,20,30]
    maxiters = [20,20,20,20,20,20,30,30]
    folds = [3,3,2,2,2,2,2,1]

if rounds == '9':
    reps = [30,20,20,20,20,20,20,30,30]
    maxiters = [20,20,20,20,20,20,30,30,30,]
    folds = [3,3,2,2,2,2,2,1,1]

if rounds == '10':
    reps = [30,30,30,30,30,30,30,30,30,30]
    maxiters = [30,30,30,30,30,30,30,30,30,30]
    folds = [3,3,2,2,2,2,2,1,1,1]
 


###########################
###  3 Population Model ###
###########################

   
if model == 'out_of_Africa':
    #if orientation == 'unfolded':
       # p_labels = "nuA,TA,nuB,TB,nuEu0,nuEuF,nuAs0,nuAsF,TF,mAfB,mAfEu,mAfAs,mEuAs, O"
       # upper = [100,10,100,10,100,100,100,100,10,200,200,200,200,0.999]
       # lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        #Optimize_Functions.Optimize_Routine(data, prefix, "SI", Model_2pop_unfolded.SI,rounds, 4, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        #if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    #if orientation == 'folded':
    p_labels = "nuA,TA,nuB,TB,nuEu0,nuEuF,nuAs0,nuAsF,TF,mAfB,mAfEu,mAfAs,mEuAs"
    upper = [100,10,100,10,100,100,100,100,10,200,200,200,200]
    lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
    Optimize_Functions.Optimize_Routine(data, prefix, "out_of_Africa", Model_3pop_folded.out_of_Africa ,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
    if done: print ("\n" + model + "folded = "+folded+ " is done\n")
 
if model == 'split_symmig_adjacent':
    p_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2"
    upper = [100,100,100,100,200,200,200,10,10]
    lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-3,1e-3]
    Optimize_Functions.Optimize_Routine(data, prefix, "split_symmig_adjacent", Model_3pop_folded.split_symmig_adjacent ,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
    if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    
if model == 'sim_split_sym_mig_all':
    p_labels = "nu1, nu2, nu3, m1, m2, T1"
    upper = [100,100,100,200,200,10]
    lower = [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
    Optimize_Functions.Optimize_Routine(data, prefix, "sim_split_sym_mig_all", Model_3pop_folded.sim_split_sym_mig_all ,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
    if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    
  
#################
###  Model SI ###
#################

if model == 'SI':
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, T1, O"
        upper = [100,100,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI", Model_2pop_unfolded.SI,rounds, 4, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, T1"
        upper = [100,100,10]
        lower = [1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI", Model_2pop_folded.SI,rounds, 3, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    
if model == 'SI_b':
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,O"
        upper = [0.999,100,100,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_b", Model_2pop_unfolded.SI_b,rounds, 5, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1, nu2, T1"
        upper = [0.999,100,100,10]
        lower = [1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_b", Model_2pop_folded.SI_b,rounds, 4, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")


if model == 'SI_ae':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,O"
        upper = [100,100,100,10,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae", Model_2pop_unfolded.SI_ae,rounds, 6, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1"
        upper = [100,100,100,10,10]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae", Model_2pop_folded.SI_ae,rounds, 5, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SI_ae_b':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,O"
        upper = [0.999,100,100,100,10,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_b", Model_2pop_unfolded.SI_ae_b,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1"
        upper = [0.999,100,100,100,10,10]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_b", Model_2pop_folded.SI_ae_b,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SI_NeC':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,O"
        upper = [100,100,100,100,10,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC", Model_2pop_unfolded.SI_NeC,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2"
        upper = [100,100,100,100,10,10]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC", Model_2pop_folded.SI_NeC,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")


if model == 'SI_ae_NeC':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,O"
        upper = [100,100,100,100,100,10,10,10,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_NeC", Model_2pop_unfolded.SI_ae_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2"
        upper = [100,100,100,100,100,10,10,10]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_NeC", Model_2pop_folded.SI_ae_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SI_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,hrf,Q,O"
        upper = [100,100,10,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_2N", Model_2pop_unfolded.SI_2N,rounds, 6, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,hrf,Q"
        upper = [100,100,10,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_2N", Model_2pop_folded.SI_2N,rounds, 5, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SI_ae_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,hrf,Q,O"
        upper = [100,100,100,10,10,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_2N", Model_2pop_unfolded.SI_ae_2N,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,hrf,Q"
        upper = [100,100,100,10,10,0.999,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_2N", Model_2pop_folded.SI_ae_2N,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SI_b_2N":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,hrf,Q,O"
        upper = [0.999,100,100,10,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_b_2N", Model_2pop_unfolded.SI_b_2N,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,hrf,Q"
        upper = [0.999,100,100,10,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_b_2N", Model_2pop_folded.SI_b_2N,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "SI_ae_b_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,hrf,Q,O"
        upper = [100,0.999,100,100,10,10,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_b_2N", Model_2pop_unfolded.SI_ae_b_2N,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,hrf,Q"
        upper = [100,0.999,100,100,10,10,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_b_2N", Model_2pop_folded.SI_ae_b_2N,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
	    
if model == 'SI_NeC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,hrf,Q,O"
        upper = [100,100,100,100,10,10,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC_2N", Model_2pop_unfolded.SI_NeC_2N,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,hrf,Q"
        upper = [100,100,100,100,10,10,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC_2N", Model_2pop_folded.SI_NeC_2N,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SI_ae_NeC_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,hrf,Q,O"
        upper = [100,100,100,100,100,10,10,10,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_NeC_2N", Model_2pop_unfolded.SI_ae_NeC_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,hrf,Q"
        upper = [100,100,100,100,100,10,10,10,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_ae_NeC_2N", Model_2pop_folded.SI_ae_NeC_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

#################
###   Model IM   ###
#################

if model == 'IM':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,m12,m21,O"
        upper = [100,100,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM", Model_2pop_unfolded.IM,rounds, 6, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':  
        p_labels = "nu1,nu2,T1,m12,m21"
        upper = [100,100,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM", Model_2pop_folded.IM,rounds, 5, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
       	 
if model == 'IM_b':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,O"
        upper = [100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b", Model_2pop_unfolded.IM_b,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21"
        upper = [100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b", Model_2pop_folded.IM_b,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_ae':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,O"
        upper = [100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae", Model_2pop_unfolded.IM_ae,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21"
        upper = [100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae", Model_2pop_folded.IM_ae,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_ae_b':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,m12,m21,O"
        upper = [0.999,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b", Model_2pop_unfolded.IM_ae_b,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,m12,m21"
        upper = [0.999,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b", Model_2pop_folded.IM_ae_b,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_NeC':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC", Model_2pop_unfolded.IM_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21"
        upper = [100,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC", Model_2pop_folded.IM_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_ae_NeC':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC", Model_2pop_unfolded.IM_ae_NeC,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21"
        upper = [100,100,100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC", Model_2pop_folded.IM_ae_NeC,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
	
if model == 'IM_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,m12,m21,hrf,Q,O"
        upper = [100,100,10,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2N", Model_2pop_unfolded.IM_2N,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,m12,m21,hrf,Q"
        upper = [100,100,10,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2N", Model_2pop_folded.IM_2N,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "IM_ae_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,hrf,Q,O"
        upper = [100,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2N", Model_2pop_unfolded.IM_ae_2N,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,hrf,Q"
        upper = [100,100,100,10,10,200,200,0.999,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2N", Model_2pop_folded.IM_ae_2N,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "IM_b_2N":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,m12,m21,hrf,Q,O"
        upper = [0.999,100,100,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2N", Model_2pop_unfolded.IM_b_2N,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,m12,m21,hrf,Q"
        upper = [0.999,100,100,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2N", Model_2pop_folded.IM_b_2N,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "IM_ae_b_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,m12,m21,hrf,Q,O"
        upper = [100,0.999,100,100,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2N", Model_2pop_unfolded.IM_ae_b_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,m12,m21,hrf,Q"
        upper = [100,0.999,100,100,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2N", Model_2pop_folded.IM_ae_b_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_NeC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2N", Model_2pop_unfolded.IM_NeC_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2N", Model_2pop_folded.IM_NeC_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "IM_ae_NeC_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2N", Model_2pop_unfolded.IM_ae_NeC_2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2N", Model_2pop_folded.IM_ae_NeC_2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_2M':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,m12,m21,i1,i2,P,O"
        upper = [100,100,10,200,200,0.999,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2M", Model_2pop_unfolded.IM_2M,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,m12,m21,i1,i2,P"
        upper = [100,100,10,200,200,0.999,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2M", Model_2pop_folded.IM_2M,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "IM_ae_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,i1,i2,P,O"
        upper = [100,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2M", Model_2pop_unfolded.IM_ae_2M,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,m12,m21,i1,i2,P"
        upper = [100,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2M", Model_2pop_folded.IM_ae_2M,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "IM_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,O"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2M", Model_2pop_unfolded.IM_b_2M,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2M", Model_2pop_folded.IM_b_2M,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "IM_ae_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,O"
        upper = [0.999,100,100,100,10,10,10, 200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2M", Model_2pop_unfolded.IM_ae_b_2M,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P"
        upper = [0.999,100,100,100,10,10,10, 200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2M", Model_2pop_folded.IM_ae_b_2M,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'IM_NeC_2M':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2M", Model_2pop_unfolded.IM_NeC_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2M", Model_2pop_folded.IM_NeC_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "IM_ae_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2M", Model_2pop_unfolded.IM_ae_NeC_2M,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':   
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2M", Model_2pop_folded.IM_ae_NeC_2M,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'IM_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2M2N", Model_2pop_unfolded.IM_2M2N,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_2M2N", Model_2pop_folded.IM_2M2N,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
		
if model == 'IM_ae_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,T1,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2M2N", Model_2pop_unfolded.IM_ae_2M2N,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,T1,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_2M2N", Model_2pop_folded.IM_ae_2M2N,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'IM_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2M2N", Model_2pop_unfolded.IM_b_2M2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_b_2M2N", Model_2pop_folded.IM_b_2M2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'IM_ae_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,T1,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2M2N", Model_2pop_unfolded.IM_ae_b_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,T1,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_b_2M2N", Model_2pop_folded.IM_ae_b_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'IM_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100, 100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2M2N", Model_2pop_unfolded.IM_NeC_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC_2M2N", Model_2pop_folded.IM_NeC_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'IM_ae_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2M2N", Model_2pop_unfolded.IM_ae_NeC_2M2N,rounds, 17, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_ae_NeC_2M2N", Model_2pop_folded.IM_ae_NeC_2M2N,rounds, 16, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
       
#################
### Model SC  ###
#################

if model =="SC":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,O"
        upper = [100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC", Model_2pop_unfolded.SC,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21"
        upper = [100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC", Model_2pop_folded.SC,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="SC_ae":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O"
        upper = [100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae", Model_2pop_unfolded.SC_ae,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21"
        upper = [100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae", Model_2pop_folded.SC_ae,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="SC_b":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,O"
        upper = [0.999,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b", Model_2pop_unfolded.SC_b,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21"
        upper = [0.999,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b", Model_2pop_folded.SC_b,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="SC_ae_b":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b", Model_2pop_unfolded.SC_ae_b,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21"
        upper = [0.999,100,100,100,10,10,10,200,2009]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b", Model_2pop_folded.SC_ae_b,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="SC_NeC":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC", Model_2pop_unfolded.SC_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21"
        upper = [100,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC", Model_2pop_folded.SC_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model =="SC_ae_NeC":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC", Model_2pop_unfolded.SC_ae_NeC,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1 nu2, nu1b, nu2b, Tae,T1,T2, m12, m21"
        upper = [100,100,100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC", Model_2pop_folded.SC_ae_NeC,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,10,10,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2N", Model_2pop_unfolded.SC_2N,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,10,10200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2N", Model_2pop_folded.SC_2N,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == "SC_ae_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2N", Model_2pop_unfolded.SC_ae_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2N", Model_2pop_folded.SC_ae_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_b_2N":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,hrf,Q,O"
        upper = [0.999,100,100,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2N", Model_2pop_unfolded.SC_b_2N,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,hrf,Q"
        upper = [0.999,100,100,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2N", Model_2pop_folded.SC_b_2N,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_ae_b_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,0.999,100,100,10,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2N", Model_2pop_unfolded.SC_ae_b_2N,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,0.999,100,100,10,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2N", Model_2pop_folded.SC_ae_b_2N,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SC_NeC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2N", Model_2pop_unfolded.SC_NeC_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2N", Model_2pop_folded.SC_NeC_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_ae_NeC_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2N", Model_2pop_unfolded.SC_ae_NeC_2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2N", Model_2pop_folded.SC_ae_NeC_2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_2M":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2M", Model_2pop_unfolded.SC_2M,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2M", Model_2pop_folded.SC_2M,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_ae_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2M", Model_2pop_unfolded.SC_ae_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2M", Model_2pop_folded.SC_ae_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

	    
if model == "SC_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2M", Model_2pop_unfolded.SC_b_2M,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P"
        upper= [0.999,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2M", Model_2pop_folded.SC_b_2M,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_ae_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2M", Model_2pop_unfolded.SC_ae_b_2M,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':  
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2M", Model_2pop_folded.SC_ae_b_2M,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
        
if model == "SC_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2M", Model_2pop_unfolded.SC_NeC_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2M", Model_2pop_folded.SC_NeC_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "SC_ae_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2M", Model_2pop_unfolded.SC_ae_NeC_2M,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2M", Model_2pop_folded.SC_ae_NeC_2M,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2M2N", Model_2pop_unfolded.SC_2M2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_2M2N", Model_2pop_folded.SC_2M2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
		
if model == 'SC_ae_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2M2N", Model_2pop_unfolded.SC_ae_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_2M2N", Model_2pop_folded.SC_ae_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SC_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2M2N", Model_2pop_unfolded.SC_b_2M2N,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1, nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_b_2M2N", Model_2pop_folded.SC_b_2M2N,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'SC_ae_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2M2N", Model_2pop_unfolded.SC_ae_b_2M2N,rounds, 16, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_b_2M2N", Model_2pop_folded.SC_ae_b_2M2N,rounds, 15, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'SC_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100, 100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2M2N", Model_2pop_unfolded.SC_NeC_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC_2M2N", Model_2pop_folded.SC_NeC_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'SC_ae_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2M2N", Model_2pop_unfolded.SC_ae_NeC_2M2N,rounds, 17, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_ae_NeC_2M2N", Model_2pop_folded.SC_ae_NeC_2M2N,rounds, 16, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
     		
#################
### Model AM  ###
#################

if model =="AM":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,O"
        upper = [100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM", Model_2pop_unfolded.AM,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21"
        upper = [100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM", Model_2pop_folded.AM,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

	    
if model =="AM_ae":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O"
        upper = [100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae", Model_2pop_unfolded.AM_ae,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21"
        upper = [100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae", Model_2pop_folded.AM_ae,rounds,8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="AM_b":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,O"
        upper = [0.999,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b", Model_2pop_unfolded.AM_b,rounds, 8, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21"
        upper = [0.999,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b", Model_2pop_folded.AM_b,rounds, 7, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="AM_ae_b":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,O"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b", Model_2pop_unfolded.AM_ae_b,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21"
        upper = [0.999,100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b", Model_2pop_folded.AM_ae_b,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="AM_NeC":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC", Model_2pop_unfolded.AM_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21"
        upper = [100,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC", Model_2pop_folded.AM_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="AM_ae_NeC":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC", Model_2pop_unfolded.AM_ae_NeC,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21"
        upper = [100,100,100,100,100,10,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC", Model_2pop_folded.AM_ae_NeC,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'AM_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,10,10200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2N", Model_2pop_unfolded.AM_2N,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,10,10200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2N", Model_2pop_folded.AM_2N,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2N", Model_2pop_unfolded.AM_ae_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2N", Model_2pop_folded.AM_ae_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_b_2N":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,hrf,Q,O"
        upper = [0.999,100,100,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2N", Model_2pop_unfolded.AM_b_2N,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,hrf,Q"
        upper = [0.999,100,100,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2N", Model_2pop_folded.AM_b_2N,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_b_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,s,nu1,n2,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,0.999,100,100,10,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2N", Model_2pop_unfolded.AM_ae_b_2N,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,s,nu1,nu2,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,0.999,100,100,10,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2N", Model_2pop_folded.AM_ae_b_2N,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'AM_NeC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2N", Model_2pop_unfolded.AM_NeC_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,10,10,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2N", Model_2pop_folded.AM_NeC_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_NeC_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2N", Model_2pop_unfolded.AM_ae_NeC_2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,hrf,Q"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2N", Model_2pop_folded.AM_ae_NeC_2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_2M":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2M", Model_2pop_unfolded.AM_2M,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2M", Model_2pop_folded.AM_2M,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2M", Model_2pop_unfolded.AM_ae_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2M", Model_2pop_folded.AM_ae_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2M", Model_2pop_unfolded.AM_b_2M,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P"
        upper= [0.999,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2M", Model_2pop_folded.AM_b_2M,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2M", Model_2pop_unfolded.AM_ae_b_2M,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':  
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2M", Model_2pop_folded.AM_ae_b_2M,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    
if model == "AM_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,O"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2M", Model_2pop_unfolded.AM_NeC_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2M", Model_2pop_folded.AM_NeC_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "AM_ae_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P,O"
        upper= [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2M", Model_2pop_unfolded.AM_ae_NeC_2M,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12,m21,i1,i2,P"
        upper= [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2M", Model_2pop_folded.AM_ae_NeC_2M,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'AM_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2M2N", Model_2pop_unfolded.AM_2M2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_2M2N", Model_2pop_folded.AM_2M2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
		
if model == 'AM_ae_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2M2N", Model_2pop_unfolded.AM_ae_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_2M2N", Model_2pop_folded.AM_ae_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'AM_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2M2N", Model_2pop_unfolded.AM_b_2M2N,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_b_2M2N", Model_2pop_folded.AM_b_2M2N,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'AM_ae_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2M2N", Model_2pop_unfolded.AM_ae_b_2M2N,rounds, 16, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_b_2M2N", Model_2pop_folded.AM_ae_b_2M2N,rounds, 15, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'AM_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100, 100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2M2N", Model_2pop_unfolded.AM_NeC_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1,nu2,nu1b,nu2b,T1,T2,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC_2M2N", Model_2pop_folded.AM_NeC_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'AM_ae_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2M2N", Model_2pop_unfolded.AM_ae_NeC_2M2N,rounds, 17, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,T1,T2,Tae,m12,m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,100,10,10,10,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_ae_NeC_2M2N", Model_2pop_folded.AM_ae_NeC_2M2N,rounds, 16, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
 

#################
###  Model two_ep  ###
#################

if model =="two_ep":
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,O"
        upper = [100,100,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep", Model_2pop_unfolded.two_ep,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21"
        upper = [100,100,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep", Model_2pop_folded.two_ep,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_ae":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0, m12, m21,O"
        upper = [100,100,100,10,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae", Model_2pop_unfolded.two_ep_ae,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0, m12, m21"
        upper = [100,100,100,10,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae", Model_2pop_folded.two_ep_ae,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_ae_b":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0, m12, m21,O"
        upper = [0.999,100,100,100,10,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b", Model_2pop_unfolded.two_ep_ae_b,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1,nu2,Tae,T1,T2,m12_0,m21_0, m12, m21"
        upper = [0.999,100,100,100,10,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b", Model_2pop_folded.two_ep_ae_b,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_b":
    if orientation == 'unfolded':
        p_labels = "s,nu1,nu2,T1,T2, m12_0,m21_0,m12, m21,O"
        upper = [0.999,100,100,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae", Model_2pop_unfolded.two_ep_ae,rounds, 10, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1,nu2,T1,T2, m12_0,m21_0,m12, m21"
        upper = [0.999,100,100,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae", Model_2pop_folded.two_ep_ae,rounds, 9, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_NeC":
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21,O"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC", Model_2pop_unfolded.two_ep_NeC,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21"
        upper = [100,100,100,100,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC", Model_2pop_folded.two_ep_NeC,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_ae_NeC":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12_0,m21_0,m12,m12, m21,O"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC", Model_2pop_unfolded.two_ep_ae_NeC,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1,nu2,nu1b,nu2b,Tae,T1,T2,m12_0,m21_0,m12,m12, m21"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC", Model_2pop_folded.two_ep_ae_NeC,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    	
if model == 'two_ep_2N':
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12, m21, hrf, Q, O"
        upper = [100,100,10,10,200,200,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2N", Model_2pop_unfolded.two_ep_2N,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12, m21, hrf, Q"
        upper = [100,100,10,10,200,200,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2N", Model_2pop_folded.two_ep_2N,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae, nu1, nu2, Tae, T1, T2, m12_0,m21_0,m12, m21, hrf, Q, O"
        upper = [100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2N", Model_2pop_unfolded.two_ep_ae_2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae, nu1, nu2, Tae, T1, T2, m12_0,m21_0,m12, m21, hrf, Q"
        upper = [100,100,100,10,10,10,200,200,200,200,0.999,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2N", Model_2pop_folded.two_ep_ae_2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_b_2N":
    if orientation == 'unfolded':
        p_labels = "s, nu1, nu2, T1, T2, m12_0,m21_0,m12, m21, hrf, Q, O"
        upper = [0.999,100,100,10,10,200,200,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2N", Model_2pop_unfolded.two_ep_b_2N,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s, nu1, nu2, T1, T2, m12_0,m21_0,m12, m21, hrf, Q"
        upper = [0.999,100,100,10,10,200,200,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2N", Model_2pop_folded.two_ep_b_2N,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_b_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae,s,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12, m21, hrf, Q,O"
        upper = [100,0.999,100,100,10,10,10,200,200,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2N", Model_2pop_unfolded.two_ep_ae_b_2N,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,s,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12, m21, hrf, Q"
        upper = [100,0.999,100,100,10,10,10,200,200,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2N", Model_2pop_folded.two_ep_ae_b_2N,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'two_ep_NeC_2N':
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2, m12_0,m21_0,m12, m21, hrf, Q, O"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999,0.499,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2N", Model_2pop_unfolded.two_ep_NeC_2N,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2, m12_0,m21_0,m12, m21, hrf, Q"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999,0.499]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2N", Model_2pop_folded.two_ep_NeC_2N,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_NeC_2N":
    if orientation == 'unfolded':
        p_labels = "nu_ae, nu1, nu2, nu1b, nu2b, Tae, T1, T2,m12_0,m21_0, m12, m21, hrf, Q, O"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2N", Model_2pop_unfolded.two_ep_ae_NeC_2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae, nu1, nu2, nu1b, nu2b, Tae, T1, T2,m12_0,m21_0, m12, m21, hrf, Q"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2N", Model_2pop_folded.two_ep_ae_NeC_2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_2M":
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper = [100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2M", Model_2pop_unfolded.two_ep_2M,rounds, 12, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper = [100,100,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2M", Model_2pop_folded.two_ep_2M,rounds, 11, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper= [100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2M", Model_2pop_unfolded.two_ep_ae_2M,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper= [100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2M", Model_2pop_folded.two_ep_ae_2M,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2M", Model_2pop_unfolded.two_ep_b_2M,rounds, 13, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper= [0.999,100,100,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2M", Model_2pop_folded.two_ep_b_2M,rounds, 12, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_b_2M":
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper= [0.999,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2M", Model_2pop_unfolded.two_ep_ae_b_2M,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1, nu2, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper= [0.999,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2M", Model_2pop_folded.two_ep_ae_b_2M,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
				
if model == "two_ep_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2M", Model_2pop_unfolded.two_ep_NeC_2M,rounds, 14, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2M", Model_2pop_folded.two_ep_NeC_2M,rounds, 13, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == "two_ep_ae_NeC_2M":
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1, nu2, nu1b, nu2b, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P,O"
        upper= [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2M", Model_2pop_unfolded.two_ep_ae_NeC_2M,rounds, 16, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1, nu2, nu1b, nu2b, Tae,T1, T2,m12_0,m21_0,m12,m21,i1,i2,P"
        upper= [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2M", Model_2pop_folded.two_ep_ae_NeC_2M,rounds, 15, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
 
if model == 'two_ep_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2M2N", Model_2pop_unfolded.two_ep_2M2N,rounds, 15, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_2M2N", Model_2pop_folded.two_ep_2M2N,rounds, 14, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
if model == 'two_ep_ae_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1, nu2, T1, T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2M2N", Model_2pop_unfolded.two_ep_ae_2M2N,rounds, 17, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1, nu2, T1, T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_2M2N", Model_2pop_folded.two_ep_ae_2M2N,rounds, 16, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'two_ep_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu1, nu2, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2M2N", Model_2pop_unfolded.two_ep_b_2M2N,rounds, 16, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu1, nu2, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_b_2M2N", Model_2pop_folded.two_ep_b_2M2N,rounds, 15, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'two_ep_ae_b_2M2N':
    if orientation == 'unfolded':
        p_labels = "s,nu_ae,nu1, nu2, T1, T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [0.999,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2M2N", Model_2pop_unfolded.two_ep_ae_b_2M2N,rounds, 18, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "s,nu_ae,nu1, nu2, T1, T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [0.999,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_b_2M2N", Model_2pop_folded.two_ep_ae_b_2M2N,rounds, 17, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")	

if model == 'two_ep_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100, 100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2M2N", Model_2pop_unfolded.two_ep_NeC_2M2N,rounds, 17, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC_2M2N", Model_2pop_folded.two_ep_NeC_2M2N,rounds, 16, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'two_ep_ae_NeC_2M2N':
    if orientation == 'unfolded':
        p_labels = "nu_ae,nu1, nu2, nu1b,nu2b, T1,T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R,O"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249,0.999]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2M2N", Model_2pop_unfolded.two_ep_ae_NeC_2M2N,rounds, 19, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if orientation == 'folded':
        p_labels = "nu_ae,nu1, nu2,nu1b,nu2b, T1, T2,Tae,m12_0,m21_0, m12, m21,i1,i2,P,hrf,Q,R"
        upper = [100,100,100,100,100,10,10,10,200,200,200,200,0.999,0.999,0.499,0.999,0.499,0.249]
        lower= [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_ae_NeC_2M2N", Model_2pop_folded.two_ep_ae_NeC_2M2N,rounds, 18, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")		
