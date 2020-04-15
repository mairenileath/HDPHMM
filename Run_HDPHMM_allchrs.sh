#don't learn haplotype switching
R CMD BATCH '--no-restore-data --args PD4120 100 5 NA 200 100 10 F' RunHDPHMM_farm_allchrs.R 'HDPHMM_fixedHapBlocks_gamma100_alpha100_'$LSB_JOBINDEX'.Rout'
#do learn haplotype switching
#R CMD BATCH '--no-restore-data --args PD4120 100 5 NA 200 100 100 T' RunHDPHMM_farm_allchrs.R 'HDPHMM_learntHapBlocks_gamma100_alpha100_'$LSB_JOBINDEX'.Rout'
exit $?