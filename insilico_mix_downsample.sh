#!/bin/bash

for perc in 0.001 0.0025 0.005 0.01 0.025 0.05
do
    for cov in 1500 1250 1000 750 500 250 100 50
    do
            wgbstools mix_pat \
            --bed_cov $cov \
            -c $cov \
            -l \
            --shuffle \
            --lbeta \
            --reps 3 \
            -p spikein_PHep_${perc}_${cov} \
            -T $SCRATCH \
            --rates $perc \
            -v \
            -L /data/users/nkueng/Project_folders/twist_panel/1.5x_tiling_MTE-94463566_hg19/panel_info/submission3_methyl_final_CpG_info.bed \
            PHep_gDNA_S2.cons.sort.pat.gz \
            buffycoat_gDNA_S3.cons.sort.pat.gz
    done
done
