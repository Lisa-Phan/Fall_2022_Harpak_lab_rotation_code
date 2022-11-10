library(stringr)

source('modify_functions.R')

#create a list for input files
merged_files <- list.files('/scratch/09069/dhp563/Summstats_redo_2/merge_AF_AC_AN',
			   pattern = "merged$", full.names=T)

#separate list
fem_merged <- merged_files[grep('/female', merged_files)]
male_merged <- merged_files[grep('/male', merged_files)]

#get names for outputfiles
nofullpath <- list.files('/scratch/09069/dhp563/Summstats_redo_2/merge_AF_AC_AN',
			 pattern = 'merged$', full.names=F)

#get path for outputfiles
outname <- nofullpath[grep('^female', nofullpath)]
path <- '/scratch/09069/dhp563/Summstats_redo_2/flattened/'

#string splicing to get a full name for outputfiles 
outfilenamelist <- list()
for (file in 1:length(outname)){
	outfilenamelist[file] <- str_replace(outname[file], "female_", "both_")
	outfilenamelist[file] <- paste0(path, outfilenamelist[file])
}


#loop that takes input arguments from list called fem_merged and male_merged
for (file in 1:24){
	flat_data <- consolidate_mf(fem_merged[file], male_merged[file])
	write.table(flat_data, as.character(outfilenamelist[file]))
}
