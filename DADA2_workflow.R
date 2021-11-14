# Sherif Hamed Atif
# DADA2 Workflow tutorial
#Source ==> https://benjjneb.github.io/dada2/tutorial.html

getwd()
#setwd("Add your path")
setwd("E:/sherif/Bioinformatics/H3ABioNet/Int_16S/Module 5/Session 1/MiSeq_SOP")


# 1- DADA2  installation -----------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")


library(dada2); packageVersion("dada2")

# home <-Sys.getenv("HOME")

# CHANGE ME to the directory containing the fastq files after unzipping.
data <- paste("E:/sherif/Bioinformatics/H3ABioNet/Int_16S/Module 5/Session 1/MiSeq_SOP", sep = '')
list.files(data)


# 2- Sort the forward and Reverse reads to be corresponding to each o --------

# From the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(data, pattern="_R1_001.fastq", full.names = T)) #full.names = T; give you the absolute adress of the file while F gives the file name only.
fnRs <- sort(list.files(data, pattern="_R2_001.fastq", full.names = T))
#VIP:: By placing full.names = F, it gives the error (no input files found) with Quality control step
# 2.1- Extract sample names, assuming filenames have format: SAMPLENA --------
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # 1 ==> the first splitted part.
           # sapply() function takes list, vector or data frame as input and gives output in vector or matrix. 
# It is useful for operations on list objects and returns a list object of same length of original set. 
# Sapply function in R does the same job as lapply() function but returns a vector.

              # strsplit ==>  used to split the strings into substrings with split arguments.
# strsplit(df,split = '%') 
# Very imp Note: : the output of strsplit is list THAT'S WHY we used Sapply
# 3- Inspect read quality profiles(Quality Control) --------------------------

# We start by visualizing the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# 4- Filter and trim ---------------------------------------------------------
# # Place filtered files in filtered/ subdirectory
filtFs <- file.path(data, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # file.path to create a subdirectory in a certain path followed by its name and save it.
filtRs <- file.path(data, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
###Arguments###
#fnFs ==> Forward reads
#filtFs ==> Output filtered forward reads
# fnRs ==> Reverse reads 
# filtRs ==> Output filtered reverse reads
# truncLen=c(240,160) ==> trancates your reads at specific reads for (forward,reverse)
# maxN=0 ==> All reads with one or more ambiguous bases are discarded.
# maxEE=c(2,2) ==> sets the maximum number of "expected errors" allowed in a read i.e,forward & reverse reads that have at least 2 errors will be removed.
#Considerations for your own data: If you want to speed up downstream computation, consider tightening (lowering) maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.
# truncQ=2 ==> 2 not indicate specific errror rate but indicate final portion of the read must not be used in final analysis
                # Score 2 maens the probability of the base being incorrect is 63%
# rm.phix=TRUE ==> remove bacteriophage genome
# compress=TRUE ==> results are compressed .gz

# Upon running the command...Creating output directory: E:/sherif/Bioinformatics/H3ABioNet/Int_16S/Module 5/Session 1/MiSeq_SOP/filtered.

head(out)


# 5- Learn the Error Rates ---------------------------------------------------

# The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates.
# The learnErrors method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition
# until they converge on a jointly consistent solution.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# to visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
# The error rates for each possible transition (A???C, A???G, .) are shown. 
# Points are the observed error rates for each consensus quality score. 
# The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
# The red line shows the error rates expected under the nominal definition of the Q-score.
# Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected.



# 6- Sample Inference --------------------------------------------------------
# 
#dada(read_After_Filter, err= error from learnError, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #output is list
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)



# 7- Merge paired reads ------------------------------------------------------

# We now merge the forward and reverse reads together to obtain the full denoised sequences.
# Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged "contig" sequences.
# By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).
# To change min. no. overlapping nucleotides ==> minOverlap()
# To allow no. of changes in the overlapped region ==> maxMismatch()

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #Output is list of dataframes

# The output  data.frame(s) has a row for each unique pairing of forward/reverse denoised sequences, and the following columns:
#   
#   $abundance: Number of reads corresponding to this forward/reverse combination.
# 
# $sequence: The merged sequence.
# 
# $forward: The index of the forward denoised sequence.
# 
# $reverse: The index of the reverse denoised sequence.
# 
# $nmatch: Number of matches nts in the overlap region.
# 
# $nmismatch: Number of mismatches in the overlap region.
# 
# $nindel: Number of indels in the overlap region.
# 
# $prefer: The sequence used for the overlap region. 1=forward; 2=reverse.
# 
# $accept: TRUE if overlap between forward and reverse denoised sequences was at least minOverlap and had at most maxMismatch differences. FALSE otherwise.
# 
# $...: Additional columns specified in propagateCol.
# 
# A list of data.frames are returned if a list of input objects was provided.

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###   

The mergers object is a list of data.frames from each sample. 
### Each data.frame contains the merged $sequence, its $abundance, and the indices of the $forward and $reverse sequence variants that were merged. 
###VIP: Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.











# 8- Construct sequence table amplicon sequence variant table (ASV) table ---------------------------------------------

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab),type = "chars")) # get the frequencies i.e, 196 sequences their length is 253.
# getSequences ==> get the sequences
# nchar==>(Count the Number of Characters (or Bytes or Width); nchar takes a character vector as an argument and returns a vector whose elements contain the sizes of the corresponding elements of x.

# Considerations for your own data: Sequences that are much longer or shorter than expected may be the result of non-specific priming.
# You can remove non-target-length sequences from your sequence table 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] #seqtab2 table after non-target removal.
# This is analogous to "cutting a band" in-silico to get amplicons of the targeted length.
dim(seqtab2)
# Note; nchar(colnames(seqtab) == nchar(getSequences(seqtab)


# 9- Remove chimeras ------------------------------------------------------

# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences.


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#removeBimeraDenovo==> Sequence variants identified as bimeric are removed, and a bimera-free collection of unique sequences is returned.
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)






# 9.1- Track reads through the pipeline( check step) ----------------------

# we'll look at the number of reads that made it through each step in the pipeline.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# 10- Assign taxonomy -----------------------------------------------------

# The assignTaxonomy function takes as input 1- a set of sequences to be classified and 
# 2- a training set of reference sequences with known taxonomy.
#it outputs taxonomic assignments with at least minBoot bootstrap confidence.
# Download the training fasta file from: https://benjjneb.github.io/dada2/training.html

taxa <- assignTaxonomy(seqtab.nochim, "E:/sherif/Bioinformatics/H3ABioNet/Int_16S/Module 5/Session 1/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
head(taxa)

# to make species level assignments based on exact matching between ASVs and sequenced reference strains
# Download the silva_species_assignment_v132.fa.gz file, and place it in the directory with the fastq files.
taxa <- addSpecies(taxa, "E:/sherif/Bioinformatics/H3ABioNet/Int_16S/Module 5/Session 1/MiSeq_SOP/silva_species_assignment_v132.fa.gz")
# addSpcies==> add the taxonomic assignment of the speies to the rest of taxonomy. 

#  inspect the taxonomic assignments
taxa.print <- taxa # make new variable for the same data.
rownames(taxa.print) <- NULL # Removing sequence rownames for display only
head(taxa.print)
#Bacteroidetes are well represented among the most abundant taxa in these fecal samples.
unique(taxa.print[,"Phylum"])
# save the ASVs table in your working directory
write.csv(taxa, file="ASVs_taxonomy.csv")
saveRDS(taxa, "ASVs_taxonomy.rds")

# save the ASVs count table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
count.asv.tab <- t(seqtab.nochim)
row.names(count.asv.tab) <- sub(">", "", asv_headers)
write.csv(count.asv.tab, file="ASVs_counts.csv")
saveRDS(count.asv.tab, file="ASVs_counts.rds")

# 11- Alignment ---------------------------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

library(DECIPHER); packageVersion("DECIPHER")

seqs <- getSequences(seqtab.nochim)
class(seqs)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
#AlignSeqs ==> Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.

# 12- Construct Phylogenetic Tree ----------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

BiocManager::install("phangorn")

library(phangorn); packageVersion("phangorn")

# we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

# Change	sequence	alignment	output	into	phyDat	structure.
phang.align <- phyDat(as(alignment, "matrix"), type="DNA") 

# Create	distance	matrix	using	dist.ml.
dm <- dist.ml(phang.align)

#Perform	neighbor	joining.
treeNJ <- NJ(dm) # Note, tip order != sequence order

# Perform	internal	maximum	likelihood.
fit = pml(treeNJ, data=phang.align)
## warning: negative edges length changed to 0!

# change the negative edges length to 0. Then, we save the fitGTR file.
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, "phangorn.tree.RDS")
