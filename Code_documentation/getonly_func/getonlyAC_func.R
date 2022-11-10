


#' getACinfo
#'
#' @param fileread: a phenotype file with coordinates lifted to hg38
#' achromfile: a chromosome file from the list of all chromosomes
#' i is the index, for creating the subfile with only values from 1 chromosome at a time
#' function aims to extract AC values that are found in achromfile, 
#' which is a file generated using VCF tools
#' @return a table with AC columns added
#' @export
#'
#' @examples
getACinfo <- function(fileread, achromfile, i){
  library(stringr)
  #only read the AC columns of the gnomAD AC AF AN table
  colstr <- c(rep("character", 4), rep("NULL", 20), rep("character",10))
  chromfile <- read.table(achromfile, header=T, colClasses = colstr)
  
  biofile <- fileread
  #removing the chr from both files 
 
  #biallelic filter the chrom file
  chromfile <- chromfile[!(duplicated(chromfile$POS) | duplicated(chromfile$POS, fromLast=T)), ]
  
  #convert chromfile AC numeric
  for (col in 5:14){
   chromfile[[col]] <- as.numeric(chromfile[[col]])
   }
  
  #limit the search to one chrom at a time
  subfile <- biofile[biofile$CHROM == as.character(i), ]
 
  #make keys
  #subfile is the phenotype file, so the key would be chrom_pos_A1
  subfile$key <- paste0(subfile$CHROM, "_", subfile$POS, "_", subfile$A1)
  chromfile$refkey <- paste0(chromfile$CHROM, "_", chromfile$POS, "_", chromfile$REF)
  chromfile$altkey <- paste0(chromfile$CHROM, "_", chromfile$POS, "_", chromfile$ALT)

  #match ref and match alt
  matchalt <- chromfile[chromfile$altkey %in% subfile$key, ]
  matchref <- chromfile[chromfile$refkey %in% subfile$key, ]
  
  #creating a tag by making a new column
  matchref$tag <- c(rep("match_ref", nrow(matchref)))
  matchalt$tag <- c(rep("match_alt", nrow(matchalt)))
  
  #creating a final dataframe by binding and merging values together
  #remove the key tag since this is no longer necessary
  fin <- rbind(matchref, matchalt)
  fin <- subset(fin, select = -c(refkey, altkey))
  mergedata <- merge(subset(subfile, select = -c(key)), fin, ID=c('POS'), all=T)
  
  return(mergedata)
}

#Example run would be something like 
#getAFinfo("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif/liftover_newcoord/female_all.albumin.glm.linear.posttest_diff_newcoord.comb",
         # "/scratch/09069/dhp563/gnomAD_data_test/allele_frequency/chr1VCFoutput.INFO", "albumin", 1)
