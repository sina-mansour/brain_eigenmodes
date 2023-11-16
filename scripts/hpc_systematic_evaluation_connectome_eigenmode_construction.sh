#!/bin/bash

# some colors for fancy logging :D
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# pang connectome eigenmodes at different densities
for density in 0.0005 0.001 0.005 0.01 0.05 0.1
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 1"
    echo "Parameters:" "${density}" "binary" "without_gyral_regression" "pang" "0" "1"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "${density}" "binary" "without_gyral_regression" "pang" "0" "1"
done

# pang connectome eigenmodes with gyral bias regression at different densities
for density in 0.0005 0.001 0.005 0.01 0.05 0.1
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 2"
    echo "Parameters:" "${density}" "binary" "with_gyral_regression" "pang" "0" "1"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "${density}" "binary" "with_gyral_regression" "pang" "0" "1"
done

# optimized connectome eigenmodes at different densities (binary)
for density in 0.0005 0.001 0.005 0.01 0.05 0.1
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 3"
    echo "Parameters:" "${density}" "binary" "without_gyral_regression" "optimized" "0" "1"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "${density}" "binary" "without_gyral_regression" "optimized" "0" "1"
done

# optimized connectome eigenmodes at different smoothing levels (binary)
for fwhm in 0 2 4 6 8 10
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 4"
    echo "Parameters:" "0.01" "binary" "without_gyral_regression" "optimized" "${fwhm}" "1"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "0.01" "binary" "without_gyral_regression" "optimized" "${fwhm}" "1"
done

# optimized connectome eigenmodes at different smoothing levels (weighted)
for fwhm in 0 2 4 6 8 10
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 5"
    echo "Parameters:" "0.1" "weighted" "without_gyral_regression" "optimized" "${fwhm}" "1e-6"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "0.1" "weighted" "without_gyral_regression" "optimized" "${fwhm}" "1e-6"
done

# optimized connectome eigenmodes at different density thresholds (weighted, smoothed)
for density in 0.0005 0.001 0.005 0.01 0.05 0.1
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 6"
    echo "Parameters:" "${density}" "weighted" "without_gyral_regression" "optimized" "8" "1e-6"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "${density}" "weighted" "without_gyral_regression" "optimized" "8" "1e-6"
done

# optimized connectome eigenmodes at different global-local integrations (weighted, smoothed)
for epsilon in "1e-6" "1e-5" "1e-4" "1e-3" "1e-2" "1e-1" "1"
do
    echo -e "${GREEN}[INFO]${NC} `date`: Part 7"
    echo "Parameters:" "0.1" "weighted" "without_gyral_regression" "optimized" "8" "${epsilon}"
    python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
                    "0.1" "weighted" "without_gyral_regression" "optimized" "8" "${epsilon}"
done

