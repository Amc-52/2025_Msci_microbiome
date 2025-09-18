#' # #  @title dadaPipe() 
#' @description performs the standard dada2 pipeline with optional batch correction with output as either a phyloseq object, TreeSummarisedExperiment (TSE) object or both.
#' @description currently only works with paired end reads, and is tailored to Illumina Miseq reads.
#' 
#' @author A. J McCarry <https://github.com/Amc-52>
#'
#' @references -DADA2- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). “DADA2:High-resolution sample inference from Illumina amplicon data.” Nature Methods, 13,581-583. doi:10.1038/nmeth.3869 <https://doi.org/10.1038/nmeth.3869>.
#' @references  -Tutorial referenced in creating the dada2 pipline- https://benjjneb.github.io/dada2/tutorial.html
#` # # 
#'
#' @importFrom dada2 plotQualityProfile filterAndTrim learnError plotErrors dada makeSequenceTable getSequences mergePairs removeBimeraDenovo getUniques assignTaxonomy
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq phyloseq
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom mia convertFromPhyloseq
#'
#' @param fastqdir FILEPATH | path to directory containing paired end fastq reads
#' @param metadata FILEPATH | path to csv file containing metadata 
#' @param R1patrn  STRING | file extension pattern for the forward reads e.g., "R1.fastq.gz" or "R1.fq.gz"
#' @param R1patrn  STRING | file extension pattern for the reverse reads e.g., "R2.fastq.gz" or "R2.fq.gz"
#' @param taxfasta FILEPATH | path to fq.gz reference for doing taxonomic assignment, e.g. silva_nr99_v138.2_toGenus_trainset.fa.gz from <https://zenodo.org/records/14169026>
#' @param out  STRING | ACCEPTED OPTIONS: ["phyloseq", "tse", "both"] | user selects whether the pipeline outputs a phyloseq object, TreeSummarisedExperiment (TSE) object or both.
#' 
#' @return a list containing:
#' @return - a phyloseq object if specified
#' @return - a TreeSummarisedExperiment object if specified
#' @return - list of 3 plots generated during quality control. The plots show the quality profiles and error rates of the reads
#' # #
# 
# install.packages(c("ggplot2", "devtools", "doParallel", "vegan", "dplyr"))
# BiocManager::install(c("dada2", "phyloseq", "Biostrings", "decontam", "MBEC", "microbiome"))
# 
# # check if required packaged are installed and loaded
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#   library(BiocManager)
# }
# # required packages
# req<- c("dada2", "phyloseq","ggplot2", "mia")
# 
# for (r in req){
#   if (!require(r, quietly = TRUE)) {
#     eval(bquote(BiocManager::install(.(r))))
#   }
# }

dadaPipe<- function(fastqdir, metadata, R1patrn, R2patrn, taxfasta, out= c("phyloseq", "tse", "both")) {

# Check if a valid output option has been selected (error if not) #
  print("Checking if input is valid...")
  out_option<- match.arg(out)
  
# Processing input files #
  
  # get matched read-pair names 
  F_reads <- sort(list.files(fastqdir, pattern = R1patrn, full.names = TRUE))
  R_reads <- sort(list.files(fastqdir, pattern = R2patrn, full.names = TRUE))
    
  # extract sample names
  sample.names <- sapply(strsplit(basename(F_reads), "R"), `[`, 1)
    
# Quality control #
  print("Starting quality control")
    
  # plot quality profiles
  F_qualityProfile<- dada2::plotQualityProfile(F_reads[[1]])
  R_qualityProfile<- dada2::plotQualityProfile(R_reads[[1]]) 
   
  # make sub-directory to store the filtered/trimmed reads
  filtF <- file.path(fastqdir, "dada_filtered_reads", paste0("filt_F", sample.names))
  filtR <- file.path(fastqdir, "dada_filtered_reads", paste0("filt_R", sample.names))
  
  # truncate!
  print('Performing quality trimming and filtering')
  # `rm.phix` discards reads that match against the phiX genome.The phiX genome
  #is a bacteriophage with a tiny genome used as a quality and calibration control for sequencing runs.
  QCd_reads <- dada2::filterAndTrim(F_reads, filtF, R_reads, filtR,
    maxN = 0, maxEE = c(2, 2), minQ=1, maxLen = 1000, rm.phix = TRUE,
    compress = TRUE)

  # make the sample.names match the rownames of the metadata table
  sample.names <- substr(sample.names, start = 1, stop = 5)
  sample.names <- gsub("_", "", sample.names)
  sample.names <- gsub("P*", "", sample.names)
  
  
# Model error rates #
  print("finished quality control. Starting error modelling")
  errF<- dada2::learnErrors(filtF)
  errR <- dada2::learnErrors(filtR)
  
  # plot errors
  errPlot<- dada2::plotErrors(errF, nominalQ = TRUE)

# Sample Inference #
  print("finished error modelling. Starting sample inference")
  dadaF <- dada2::dada(filtF, err = errF)
  dadaR <- dada2::dada(filtR, err = errR)
  
# Merging the reads #
  mergers <- dada2::mergePairs(dadaF, filtF, dadaR, filtR, verbose = TRUE)
  
  # inspect merged reads
  head(mergers[[1]])
  
# Construct sequence table #
  seqtbl <- dada2::makeSequenceTable(mergers)
  
  # inspect distribution of sequence lengths
  table(nchar(dada2::getSequences(seqtbl)))
  
# Remove chimeras #
  print("Finished sample inference. Removing chimeras...")
  seqtbl.nochim <- dada2::removeBimeraDenovo(seqtbl, method = "consensus", verbose = TRUE)
  dim(seqtbl.nochim)
  sum(seqtbl.nochim) / sum(seqtbl)
  
  # check the number of reads filtered during each stage
  getN <- function(x) sum(getUniques(x))
  print('H')
  track <- cbind(QCd_reads, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtbl.nochim))
  print('H')
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  print('H')
  rownames(track) <- sample.names
  print('H')
  track<- as.data.frame(track)
  
# Taxonomic Assignment #
  print("Assigning taxonomies")
  taxa <- dada2::assignTaxonomy(seqtbl.nochim, "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/Resources/silva_nr99_v138.2_toGenus_trainset.fa.gz")
  

  # make the sampleIDs the rownames in the metadata table
  rownames(metadata) <- metadata$sampleID
  
# Output Results #
  print("Creating output objects...")
  
  # adding plots to a list
  plots<- list( F_qualityProfile, R_qualityProfile)
  
  # change the seq table rownames to match metadata
  rownames(seqtbl.nochim) <- gsub("filt_F", "", rownames(seqtbl.nochim))
  rownames(seqtbl.nochim) <- gsub("_community_run", "", rownames(seqtbl.nochim))
  rownames(seqtbl.nochim) <- gsub("_feces_L001_", "", rownames(seqtbl.nochim))
  rownames(seqtbl.nochim) <- substr(rownames(seqtbl.nochim), start = 1, stop = 5)
  rownames(seqtbl.nochim) <- gsub("_.*", "", rownames(seqtbl.nochim))
  
  # make phyloseq object
  if (out == "both" || "phyloseq") {
  phy <- phyloseq::phyloseq(
    otu_table(seqtbl.nochim, taxa_are_rows = FALSE),
    sample_data(metadata),
    tax_table(taxa)
  )
  }
  
  # make TreeSummarisedExperiment object
  if (out == "both" || "tse") {
  tse<- mia::convertFromPhyloseq(phy)
  }
  
  # add objects to output list
  if (out == "both") {
  outlist<- list(phy, tse, plots)
  }
  
  if (out == "tse") {
    outlist<- list(tse, plots)
    
  }
  
  if (out == "phyloseq") {
    outlist<- list(phy, plots)

  }
  
  print("dada2 pipeline complete :]")
  return(outlist)

}

dadaResults<- dadaPipe(
  fastqdir= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/RawData/testdata",
  metadata= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/Tables/testMetadata.csv",
  R1patrn= "R1.fastq.gz",
  R2patrn= "R2.fastq.gz",
  taxfasta= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/Resources/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  out = "both"
)

#` example implementation of function`
#' @examples
#' \dontrun{
#' 
#'dadaResults<- dadaPipe(
#'  fastqdir= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/RawData/testdata",
#'  metadata= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/Tables/testMetadata.csv",
#'  R1patrn= "R1.fastq.gz",
#'  R2patrn= "R2.fastq.gz",
#'  taxfasta= "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/Resources/silva_nr99_v138.2_toGenus_trainset.fa.gz",
#'  out = "both"
#'  )
#'}
