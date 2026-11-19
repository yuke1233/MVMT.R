#!/bin/bash

# Run from the directory containing MVMT_power.R.
# Usage: bash submit_power.sh
# Submit nine power simulations with the requested parameter grid.

n_traits_values=(2 5 7)
shared_values=(0.2)
rg_values=(0.2)
re_values=(0.3 0.5 0.7)

for n_traits in "${n_traits_values[@]}"; do
  for shared in "${shared_values[@]}"; do
    for rg in "${rg_values[@]}"; do
      for re in "${re_values[@]}"; do
        echo "Submitting power: nt${n_traits}_sh${shared}_rg${rg}_re${re}"
        qsubshcom "Rscript ./MVMT_power.R $n_traits $shared $rg $re" \
                  32 100G power 80:00:00 "" || exit 1
      done
    done
  done
done
