#!/bin/bash

julia requirements.jl
foldername=$(julia heatmap_setup.jl)



max_processes=$(nproc --all)
current_processes=0
for mupv in $(seq 0 0.04 1); do 
    Rvv=0.1
        for sigmaV in $(seq 0 0.01 0.35); do 
            julia makeheatmap.jl $Rvv $mupv $sigmaV >> sims/$foldername/MUpv_phase.csv & 
            ((current_processes++))
            if [ $current_processes -ge $max_processes ]; then
                # Wait for background processes to finish before starting more
                wait
                current_processes=0
            fi
    done 
done

wait
for Rvv in $(seq 0 0.005 0.2); do
    mupv=0.1
    for sigmaV in $(seq 0 0.003 0.15); do 
        julia makeheatmap.jl $Rvv $mupv $sigmaV >> sims/$foldername/Rvv_phase.csv & 
        ((current_processes++))
        if [ $current_processes -ge $max_processes ]; then
            # Wait for background processes to finish before starting more
            wait
            current_processes=0
        fi
    done
done 

wait
