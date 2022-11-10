#' TtestFunc
#' @param 
#' f, m: directory to the summstat tsv file for female and male, enclosed in
#' double quote
#' if smallout = T, then the resulting table will only contain
#' c("CHROM", "POS", "REF", "ALT", "A1")
#' print out two tables 
#' @return two_data_frames
#' 
#' @export
#'
#' @examples testFunc("test_female_all.diastolicBP", "test_male_all.diastolicBP")
TtestFunc <- function(f, m, alpha=0.05){
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
  ttest <- welch_t_test(mtab$BETA, ftab$BETA, mtab$SE, ftab$SE)
  
  #conservative est of df
  df <- pmin(mtab$OBS_CT, ftab$OBS_CT) - 1
  
  #p_value. Null hypothesis: beta1 and 2 are similar
  p_value <- 2*pt(ttest, df=df, lower.tail = F)
  
  #creating a new column for p_value 
  mtab$p_value <- p_value
  ftab$p_value <- p_value
 
  #renaming the columns
  names(mtab)[10] <- c("MALE_BETA")
  names(mtab)[11] <- c("MALE_SE")

  names(ftab)[10] <- c("FEMALE_BETA")
  names(ftab)[11] <- c("FEMALE_SE")

  #This is part of the original attempt to make an arbitrary break using alpha
  #Not a very good idea
  sub_mtable <- mtab[mtab$p_value <= alpha, ]
  sub_mtable2 <- mtab[mtab$p_value > alpha, ]

  sub_ftable <- ftab[ftab$p_value <= alpha, ]
  sub_ftable2 <- ftab[ftab$p_value > alpha, ]

  sub_mtable <- data.frame(sub_mtable, sub_ftable$FEMALE_BETA, sub_ftable$FEMALE_SE, sub_ftable$OBS_CT)
  sub_mtable2 <- data.frame(sub_mtable2, sub_ftable2$FEMALE_BETA, sub_ftable2$FEMALE_SE, sub_ftable2$OBS_CT)

  #f and m are arbitrary labels and can cause confusion
  #not a good idea by design
  write.table(sub_mtable, file= paste0(as.character(f), ".posttest_diff"))
  write.table(sub_mtable2, file= paste0(as.character(m), ".posttest_nodif"))
  return (1)

}

