
getAFinfo <- function(fileread, achromfile, i){
  library(stringr)
  colstr <- c(rep("character", 14), rep("NULL", 20))
  chromfile <- read.table(achromfile, header=T, colClasses = colstr)
  
  biofile <- fileread
  
  #biallelic filter the chrom file
  chromfile <- chromfile[!(duplicated(chromfile$POS) | duplicated(chromfile$POS, fromLast=T)), ]
  
  #convert chromfile AF numeric
  for (col in 5:14){
   chromfile[[col]] <- as.numeric(chromfile[[col]])
   }
  
  #limit the search to one chrom at a time
  subfile <- biofile[biofile$CHROM == as.character(i), ]
 
  #make keys
  #subfile is the biometric file, so the key would be chrom_pos_A1
  subfile$key <- paste0(subfile$CHROM, "_", subfile$POS, "_", subfile$A1)
  chromfile$refkey <- paste0(chromfile$CHROM, "_", chromfile$POS, "_", chromfile$REF)
  chromfile$altkey <- paste0(chromfile$CHROM, "_", chromfile$POS, "_", chromfile$ALT)
  
  #match ref and match alt
  #if match alt, keep the same
  #if match ref, change ref AC, AF, keep AN
  matchalt <- chromfile[chromfile$altkey %in% subfile$key, ]
  matchref <- chromfile[chromfile$refkey %in% subfile$key, ]
  
  
  #ac -10 is corresponding an
  #for (ac in 25:34){
  #  matchref[[ac]] <- matchref[[ac-10]] - matchref[[ac]] 
  #}
  #add a tag just in case for latter
  matchref$tag <- c(rep("match_ref", nrow(matchref)))
  matchalt$tag <- c(rep("match_alt", nrow(matchalt)))
 	
  print("checkpoint 1")
  fin <- rbind(matchref, matchalt)
  fin <- subset(fin, select = -c(refkey, altkey))
  return(fin)
}
  
#getAFinfo("/scratch/09069/dhp563/Summstats_cp_full/diff_nodif/liftover_newcoord/female_all.albumin.glm.linear.posttest_diff_newcoord.comb",
         # "/scratch/09069/dhp563/gnomAD_data_test/allele_frequency/chr1VCFoutput.INFO", "albumin", 1)
