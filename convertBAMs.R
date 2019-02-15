##Functions to easily read & show BAMs in R using ggplot2.
##Renske van Raaphorst, Lausanne, dec 18 - feb 19

#to install Rsamtools:
#install.packages("biocInstaller")
#biocLite("Rsamtools")

#to install ggplot2:
#install.packages("ggplot2")

##To get the BAM files into R and convert them to ggplot2
getBAMS <- function(IPname, inputname){ #IPname & inputname are the paths to the BAM files
  
  ##this requires Rsamtools, code below is to check & install this package if needed.
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){
      if (!requireNamespace("biocInstaller", quietly = TRUE)){install.packages("biocInstaller")}
      biocInstaller::biocLite("Rsamtools")}else{stop("Canceled")}
  }
  
  #read in the Bam Files, input & IP
  bamIP <- Rsamtools::scanBam(IPname)
  baminput <- Rsamtools::scanBam(inputname)
  
  elts <- setNames(names(bamIP[[1]]), names(bamIP[[1]])) #get the names of each component of the BAM
  lstbamIP <- lapply(elts, function(elt) .unlist(lapply(bamIP, "[[", elt))) #turn the bam into a list
  rm(bamIP)
  
  lstbaminput <- lapply(elts, function(elt) .unlist(lapply(baminput, "[[", elt))) #same for the input BAM 
  rm(baminput)
  
  baminputframe <- do.call("DataFrame", lstbaminput)
  
  bamIPframe <- do.call("DataFrame", lstbamIP) #turn the list into a data frame
  bamIPframe <- as.data.frame(bamIPframe)
  bamIPframe <- binBAMframes(bamIPframe, 5000)
  
  
  baminputframe <- as.data.frame(baminputframe)
  baminputframe <- binBAMframes(baminputframe, 5000)
  

  bamIPframe <- relBAMS(bamIPframe, baminputframe)
  return(list(bamIPframe, baminputframe))
}


.unlist <- function(x)   #<-- this function comes from https://gist.github.com/davetang/6460320 - I use it in the function above.
{
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

## bin the dataframes, makes them considirably easier to work with. standard bins = 5000 bps, but if you want more details in your plot, put it lower.
binBAMframes <- function(bamFRAME, bin=5000){
  maxpos <- max(bamFRAME$mpos,na.rm=TRUE) # find end of the sequence.
  
  binBAM <- hist(bamFRAME$mpos[!is.na(bamFRAME$mpos)], breaks= maxpos%/%bin, #use histogram function to get the read count of each sequence point, use cut()
                 plot=FALSE)                                                     #to bin these reads per 5000 basepairs.
  
  binBAM <- data.frame(breaks = binBAM$breaks[-1] - bin/2, count = binBAM$counts, density=binBAM$density) #turn the histogram into a data frame.
  
  totalcounts <- sum(binBAM$count)
  binBAM$perc_reads <- binBAM$count/totalcounts*100 #turn the count into a "percentage of total reads" so its possible to compare IP to input.
  return(binBAM)
}

relBAMS <- function(IPframe, inputframe){ #finally, remove background and correct for what's happening in the input.
  #relative
  IPframe$corrected <- IPframe$perc_reads - inputframe$perc_reads
  IPframe$corrected <- IPframe$corrected - min(IPframe$corrected)
  return(IPframe)
}

