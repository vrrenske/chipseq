##Functions to easily read & show BAMs in R using ggplot2.
##Renske van Raaphorst, Lausanne, dec 18 - feb 19

#to install Rsamtools:
#install.packages("biocInstaller")
#biocLite("Rsamtools")

#to install ggplot2:
#install.packages("ggplot2")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}

if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

if (!requireNamespace("biofiles", quietly=TRUE)){
  BiocManager::install("biofiles")
}

if (!requireNamespace("htmlwidgets", quietly=TRUE)){
  install.packages("htmlwidgets")
}

library(plotly)
library(ggplot2)

##To get the BAM files into R and convert them to ggplot2
getBAMS <- function(IPname, inputname){ #IPname & inputname are the paths to the BAM files

  #read in the Bam Files, input & IP
  bamIP <- read.delim(paste("Output/depth_files/", IPname,sep=""), header=FALSE)
  baminput <- read.delim(paste("Output/depth_files/", inputname, sep=""), header=FALSE)
  colnames(bamIP) <- c("genome", "bp", "reads")
  colnames(baminput) <- c("genome", "bp" , "reads")
  bamIPframe <- relBAMS(bamIP, baminput)
  bamIPframe$sample <- strsplit(inputname, "-")[[1]][2]
  bamIPframe$replicate <- strsplit(strsplit(inputname, "-")[[1]][3], ".txt")[[1]][1]
  return(bamIPframe)
}


relBAMS <- function(IPframe, inputframe){ #finally, remove background and correct for what's happening in the input.
  #relative
  IPframe$corrected <- IPframe$reads / inputframe$reads
  return(IPframe)
}

plotBAMS_1 <- function(bamIPframe, title, bins=200){
  bamIPframe$bins <- cut(bamIPframe$bp, max(bamIPframe$bp)/bins, 1:(max(bamIPframe$bp)/bins)*bins)
  bamIPframe$bins <- as.numeric(as.character(bamIPframe$bins))
  bamIPframe <- aggregate(list("reads" = bamIPframe$reads, "corrected" = bamIPframe$corrected), by=list("sample" = bamIPframe$sample, "replicate" = bamIPframe$replicate, "bp" = bamIPframe$bins), FUN=mean)
  return(
    ggplot(bamIPframe, aes(x=bp, y=corrected, color=sample, group=replicate), alpha=0.6) +
      geom_line() +
      theme_minimal() +
      ggtitle(title) +
      xlab("basepairs") +
      ylab("depth IP/input")
  )
}

###################read files & save
files_to_read <- list.files("Output/depth_files")
files_to_read <- files_to_read[order(files_to_read)]
files_to_read <- data.frame("input" = files_to_read[1:(length(files_to_read)/2)], "IP" = files_to_read[(length(files_to_read)/2+1):length(files_to_read)])
files_to_read$input <- as.character(files_to_read$input)
files_to_read$IP <- as.character(files_to_read$IP)

combined_dataset <- do.call("rbind", lapply(1:nrow(files_to_read), function(x) getBAMS(files_to_read$IP[x], files_to_read$input[x])))

save(combined_dataset, file="Output/plots/combined_dataframe.Rda")


##############plot files
ggsave(plotBAMS_1(combined_dataset, "overlay of samples", bins=200), filename="Output/plots/Overlay.PDF", width=10, height=4)
ggsave(plotBAMS_1(combined_dataset, "", bins=200)+facet_grid(sample~.), filename="Output/plots/facets.PDF", width=10, height=10)

################plot more elaborate using plotly

#get corresponding genome
D39 <- "CP027540 (D39V).gb"
prokka <- biofiles::gbRecord(D39)
genelist <- biofiles::select(prokka, .cols = c("gene", "key", "locus_tag", "start", "end"))
genelist <- genelist[genelist$key=="gene",]


#function for plotly plot

plotlyplot <- function(combined_dataset, genelist, bins=200){

  combined_dataset$bins <- cut(combined_dataset$bp, max(combined_dataset$bp)/bins, 1:(max(combined_dataset$bp)/bins)*bins)
  combined_dataset$bins <- as.numeric(as.character(combined_dataset$bins))
  combined_dataset <- aggregate(list("reads" = combined_dataset$reads, "corrected" = combined_dataset$corrected), by=list("sample" = combined_dataset$sample, "replicate" = combined_dataset$replicate, "bp" = combined_dataset$bins), FUN=mean)

  plot1 <- plotly::ggplotly(
    ggplot(data=genelist) +
      geom_line(data=combined_dataset, aes(x=bp, y=corrected, color=sample, group=replicate), alpha=0.5) +
      geom_rect(aes(xmin=start, xmax=end, ymin=-1, ymax=0,fill=gene), size=.1) +
      geom_rect(data=genelist[is.na(genelist$gene),], aes(xmin=start, xmax=end, ymin=-1, ymax=0, fill=locus_tag), size=.1) +
      theme_minimal() +
      scale_fill_manual(values=rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7"), 500)) +
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")) +
      xlab("basepairs (distance from origin of replication)") +
      ylab("relative enrichment (reads IP/reads input)") +
      theme(legend.position="none", panel.grid=element_blank()) +
      scale_x_continuous(breaks = seq(0, 2040000, 10000)) +
      ylim(-1, NA),
    hoverinfo="color"
  )
  plot1<- style(plot1%>%ggplotly(tooltip="color"), hoverinfo="none", traces=1:4)
  return(plot1)
}



htmlwidgets::saveWidget(plotlyplot(combined_dataset, genelist),file = "overlay_plotly.html", selfcontained=FALSE)
