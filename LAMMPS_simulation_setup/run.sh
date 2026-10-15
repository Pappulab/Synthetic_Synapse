#!/bin/bash
hmDir=$PWD

for N in 1 500 1000 1500 2000 2500 3000 3500 4000; do
    for rep in 1 2 3 4 5; do
        SEED=$((RANDOM * RANDOM))
        DIR=N_${N}/rep${rep}
        mkdir -p ${DIR}
        cp lammps_sim.in ${DIR}
        cp submit_N.sbatch ${DIR}
        cp prod.restart1 ${DIR}
        cd ${DIR}
        sed -i "s/NUMBER/${N}/g" "lammps_sim.in"
        sed -i "s/SEED/${SEED}/g" "lammps_sim.in"
        sed -i "s/NUMBER/${N}/g" "submit_N.sbatch"
        sbatch submit_N.sbatch
        cd ${hmDir}
        echo $PWD
    done
done

