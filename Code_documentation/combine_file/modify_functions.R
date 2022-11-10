modify_AF <- function(tag, AF){
   if (!is.na(AF)){
        return(ifelse(tag == "match_ref", 1 - AF, AF))
    }
   else {
        return(NA)
    }
}

#Function
modify_AC <- function(tag, AC, AN){
    if (!is.na(AC) & !is.na(AN)){
        return(ifelse(tag == "match_ref", AN - AC, AC))
    }
    else {
        return(NA)
    }
}

#function that takes in two file paths, and then combine them 
#into one after making appropriate modifications 
#on AF, AC
#return a file

consolidate_mf <- function(fpath, mpath){
  ftab <- read.table(fpath, header=T)
  mtab <- read.table(mpath, header=T)
  
  #Make sure the files are formatted properly before binding
  #All colnames should be similar
  
  combine <- rbind(ftab, mtab)
  
  #reordering the columns
  combine <- combine[c("CHROM", "POS", "REF", "ALT", "A1",
                       "FEMALE_BETA", "FEMALE_SE", "FEMALE_COUNT",
                       "MALE_BETA", "MALE_SE", "MALE_COUNT", "T_TEST", "tag",
                       "AF_afr", "AF_asj", "AF_eas", "AF_sas", "AF_mid",
                       "AF_fin", "AF_nfe", "AF_oth", "AF_ami", "AF_amr",
                       "AN_afr", "AN_asj", "AN_eas", "AN_sas", "AN_mid",
                       "AN_fin", "AN_nfe", "AN_oth", "AN_ami", "AN_amr",
                       "AC_afr", "AC_asj", "AC_eas", "AC_sas", "AC_mid",
                       "AC_fin", "AC_nfe", "AC_oth", "AC_ami", "AC_amr")]

   
 for (a in 14:43){
    combine[[a]] <- as.numeric(combine[[a]])
  }
  #combine <- combine %>%  mutate_at(c(14:43), as.numeric)

  #remove all unessentials
  #filtering out rows with missing values(no AC, AN, AF data)
  subsetcombine <- combine[rowSums(is.na(combine)) < 30,]
 
  #fix AF
  combine <- combine %>% mutate_at(c(14:23),function(x) mapply(FUN=modify_AF, combine$tag, x))
  
  #these are AC columns 
  for (j in 34:43){
    subsetcombine[[j]] <- mapply(FUN=modify_AC, subsetcombine$tag,
                           subsetcombine[[j]], subsetcombine[[j - 10]])
  }

  return (subsetcombine)
}
