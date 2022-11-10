#' TtestFunc
#' @param 
#' f, m: directory to the summstat tsv file for female and male, the raw 
#' files that comes from PLINK generated GWAS result
#' name is the output file's name
#'
#' @return two_data_frames
#' 
#' @export
#'
#' @examples testFunc("test_female_all.diastolicBP", "test_male_all.diastolicBP")
TtestFunc <- function(f, m, name){
  library(dplyr)
  
  mtab <- read.delim(m)
  ftab <- read.delim(f)
  
  #making sure the file is as expected
  conf <- c("X.CHROM","POS","ID", "REF", "ALT",
            "A1", "AX", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
  
  if (all((colnames(mtab) != conf))){
      print("Please check file format")
      return (2)
  }
  #Welch ttest function
  welch_t_test <- function(beta1, beta2, se1, se2){
    return ((beta1-beta2)/ sqrt(se1^2 + se2^2))
  }
  
  ttest <- welch_t_test(as.numeric(mtab$BETA), as.numeric(ftab$BETA), as.numeric(mtab$SE), as.numeric(ftab$SE))
  

  #calculating df
  #R vectorize operation
  df <- mtab$OBS_CT + ftab$OBS_CT - 2

  #p_value. Null hypothesis: beta1 and 2 are similar
  p_value <- 2*pt(ttest, df=df, lower.tail = F)
  
  #creating a new column for p_value 
  mtab$P_VAL_TWO_SAMPLE <- p_value
  ftab$P_VAL_TWO_SAMPLE <- p_value
 
  #renaming the columns
  names(mtab)[9:13] <- c("MALE_COUNT", "MALE_BETA", "MALE_SE", "MALE_T_STAT", "MALE_P")
  names(ftab)[9:13] <- c("FEMALE_COUNT", "FEMALE_BETA", "FEMALE_SE", "FEMALE_T_STAT", "FEMALE_P")
 
  #names(mtab)[6] <- c("MALE_A1")
  #names(ftab)[6] <- c("FEMALE_A1")

  #merge the results into one dataframe
  #only keep SNPs that are present in both datasets
  merged <- merge(mtab, ftab)
  #write a table, using paste0 to keep naming consistent
  write.table(merged, file= paste0(name, ".posttest"), row.names=F)
  return (1)

}

