#FIM uncertainty 2 pop models
import moments
import Model_2pop_folded
import sys

data = open(sys.argv[1], "r") # 1st sys.argv with location of sfs
pop1 = sys.argv[2]
pop2= sys.argv[3]
data = moments.Spectrum.from_file(data, mask_corners=False, return_comments=False)
params = [1.7664,1.6217,0.3179,0.6892,0.2369,0.1145] # paste comma separated list
ns = [18,14] # comma separated list with population sizes
model = Model_2pop_folded.SC #change model accordingly no parameter space needed

# FIM uncertainty function in moments - no bootstrap needed if data is unlinked
FIM = moments.Godambe.FIM_uncert(model, params, data, log=False, multinom=True, eps=0.01)
FIM = repr(FIM)
uncert = open(pop1 + pop2 + ".uncert.txt", "w") # the last entry in the array written in this file is uncertainty for theta
uncert.write(FIM)
uncert.close
# storing output in a text file
# run this script on bash terminal: python3 FIM_2pop_basic.py > POP1-POP2.uncert.txt
