source("liftOver.function.R")

#script that loops over a list of three files
#pass these files into liftOver.subset.func
#files should be sorted

#list for conversion files
#list for unmapped file
#list for posttest on directory above

bio_posttest <- list.files("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif",
                           pattern="posttest_nodif$",
                           full.names=T)
unMapped <- list.files("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif/liftover_redo_male/",
                       pattern="^unMappedmale",
                       full.names = T)
convert <- list.files("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif/liftover_redo_male/",
                         pattern="^conversionsmale",
                         full.names = T)

for (i in 1:24){
  liftOver.subset.func(bio_posttest[i], unMapped[i], convert[i])
}


