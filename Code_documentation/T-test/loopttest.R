source("ttest.R")

femalelist <- list.files(path = '/scratch/09069/dhp563/Summstats_redo_2/all_file', pattern="^female", all.files=F, full.names=T)
malelist <- list.files(path = '/scratch/09069/dhp563/Summstats_redo_2/all_file', pattern="^male", all.files=F, full.names=T)

for (i in 1:length(femalelist)){
  TtestFunc(femalelist[i], malelist[i])
}

