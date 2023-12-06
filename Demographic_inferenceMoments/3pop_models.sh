#!/bin/sh

while read model run
do 
  nohup python3 Script_moments_model_optimisation_3pop.py -f NIQ-SBH-TPC -m ${model} -r ${run} -n 4 --orientation folded &

done < model_run_3pop
