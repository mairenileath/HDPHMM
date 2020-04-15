#$LSB_JOBINDEX
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 10000 50 5' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 10000 25 5' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 10000 50 10' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 1000 10000 50 10 100' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#params based on chr size and expected haplotype block size
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 NA 200 10' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 995 NA 200 10' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#relax concentration parameter
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 NA 200 100' RunHDPHMM_farm.R 'HDPHMM'$LSB_JOBINDEX'.Rout'
#set both concentration params
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 100 NA 200 10 25' RunHDPHMM_farm.R 'HDPHMM_alpha25_'$LSB_JOBINDEX'.Rout'
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 1000 NA 200 10 100 100' RunHDPHMM_farm.R 'HDPHMM_alpha100_'$LSB_JOBINDEX'.Rout'
#don't learn haplotype switching
R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 500 NA 200 100 10 F 95' RunHDPHMM_farm.R 'HDPHMM_fixedHapBlocks_gamma100_alpha100_'$LSB_JOBINDEX'.Rout'
#do learn haplotype switching
#R CMD BATCH '--no-restore-data --args PD4120 '$LSB_JOBINDEX' 30 95 NA 200 100 100 T' RunHDPHMM_farm.R 'HDPHMM_learntHapBlocks_gamma100_alpha100_'$LSB_JOBINDEX'.Rout'
exit $?