#!/bin/bash
while IFS= read -r line; do
    python /home/script/python/run_sim.py "$line"
    Rscript /home/script/Rscript/run_sim/run_sim_hic.r "$line"
#Rscript /home/script/Rscript/run_sim/按均值范围计算灵敏度.r "$line"
done < /home/script/Rscript/run_sim/input_simdata.txt