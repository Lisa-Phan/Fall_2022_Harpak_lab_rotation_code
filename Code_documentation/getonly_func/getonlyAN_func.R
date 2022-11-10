


#' getANinfo
#'
#' @param fileread: a phenotype file with coordinates lifted to hg38
#' achromfile: a chromosome file from the list of all chromosomes
#' i is the index
#' @return a table with AN columns added
#' @export
#'
#' @examples
getANinfo <- function(fileread, achromfile, i){
  library(stringr)

  #only read the AN columns of the gnomAD AF AC AN table
  colstr <- c(rep("character", 4), rep("NULL", 10), rep("character",10), rep("NULL", 10))
  chromfile <- read.table(achromfile, header=T, colClasses = colstr)
  
  biofile <- fileread
  #removing the chr from both files 
  #biofile$CHROM <- str_replace(biofile$CHROM, pattern="chr", replacement="")
  #chromfile$CHROM <- str_replace(chromfile$CHROM, pattern="chr", replacement="")
 
  #biallelic filter the chrom file
  chromfile <- chromfile[!(duplicated(chromfile$POS) | duplicated(chromfile$POS, fromLast=T)), ]
  
  head(chromfile)
  head(biofile)
  #convert chromfile AN numeric
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
  head(chromfile)
  head(subfile)

  #match ref and match alt
  matchalt <- chromfile[chromfile$altkey %in% subfile$key, ]
  matchref <- chromfile[chromfile$refkey %in% subfile$key, ]
  
  #replace AF
  #for (af in 5:14){
   #matchref[[af]] <- 1 - matchref[[af]]
   #}
  #ac -10 is corresponding an
  #for (ac in 25:34){
  #  matchref[[ac]] <- matchref[[ac-10]] - matchref[[ac]] 
  #}
  #add a tag just in case for latter
  
  matchref$tag <- c(rep("match_ref", nrow(matchref)))
  matchalt$tag <- c(rep("match_alt", nrow(matchalt)))
  
  fin <- rbind(matchref, matchalt)
  fin <- subset(fin, select = -c(refkey, altkey))
  mergedata <- merge(subset(subfile, select = -c(key)), fin, ID=c('POS'), all=T)
  
  return(mergedata)
}
  
#getAFinfo("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif/liftover_newcoord/female_all.albumin.glm.linear.posttest_diff_newcoord.comb",
         # "/scratch/09069/dhp563/gnomAD_data_test/allele_frequency/chr1VCFoutput.INFO", "albumin", 1)
