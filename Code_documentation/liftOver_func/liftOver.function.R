#TODO: function takes two inputs as two tables, big and small
# write a function that give the output as a subset of the big table, 
# containing all columns of big table that are not in small table


#big <- read.table("Test_data/testfemale_all.albumin.glm.linear.posttest_diff") 
#small <- read.table("Test_data/testunMappedmale_all.albumin.glm.linear.posttest_nodif.BED")
#conversion <- read.table("Test_data/testconversionsmale_all.albumin.glm.linear.posttest_nodif.BED")


#' liftOver.subset.func
#'
#' @param big: the path to the big file, the files that ends in .posttest  
#' @param small: the path to the small file, unMap outputs from liftOver
#'
#' @return
#' 
#' another file that contains the new coordinate, created from 
#' combining info from conv file with the other
#' @export
#'
#' @examples
liftOver.subset.func <- function(b, s, con){
  library(operators)
  big <- read.table(b)
  small <- read.table(s)
  #create keys
  small$key <- paste0(small[[1]], "_", small[[2]])
  big$key <- paste0("chr", big[[1]], "_", big[[2]])
  
  #take subset that are not part of the unMapped file
  #these are SNPs that has converted coordinates 
  subsetbig <- big[(big$key %!in% small$key), ]
  
  #these are SNPs that are unchanged, keep the same
  subsetbig2 <- big[(big$key %in% small$key), ]
  subsetbig2 <- subsetbig2[, c(1:5)]

  conversion <- read.table(con, header=T)
  
  #merge the data together, now with the new coordinates replacing the old ones
  #This should always execute, put condition just in case things go wrong
  if (nrow(conversion) == nrow(subsetbig)){
    final <- data.frame(conversion$NewChr, conversion$NewPos, subsetbig$REF, subsetbig$ALT, subsetbig$A1)
    appfinal <- as.data.frame(mapply(c, subsetbig2,final))
    write.table(appfinal, paste0(as.character(b), "_newcoord.comb"))
  } 
  else{
    print("Unequal number of rows!")
    sprintf("Subset row and conv row :%i %i", nrow(subsetbig), nrow(conversion))
  }
}






