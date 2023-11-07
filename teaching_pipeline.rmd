---
title: "Teaching Pipeline 2023"
author: "Lukas Dreyling & Henrique Valim"
date: "13. - 16.11.2023"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Program List

 * [BLAST+](https://github.com/ncbi/blast_plus_docs): sequence alignment tool

 * [cutadapt](https://cutadapt.readthedocs.io/en/stable/): demultiplexing and primer removal

 * [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): generating quality reports

 * [multiQC](https://multiqc.info/): aggregating multiple reports into a single report


 * [R(<https://www.r-project.org/>) (and [RStudio](https://posit.co/products/open-source/rstudio/)): statistical software environment
 ### Specific functions in R come with so-called packages so we need to install and load these as well. 
    + [dada2](https://benjjneb.github.io/dada2/): Denoising, filtering, clustering ASVs and assigning taxonomy
    + [ShortRead](https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html): Reading and examining raw sequence reads
    + [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html): Reading and examining sequencing data
    + [phyloseq](https://joey711.github.io/phyloseq/): Explore microbiome profiles using R (good for getting our data into a format that works with many other packages)
    + [tidyverse](https://www.tidyverse.org/): R packages for data science and data manipulation
    + [MicEco](https://github.com/Russel88/MicEco): R packge for analysis of microbial community data
    + [vegan](https://cran.r-project.org/web/packages/vegan/index.html): Methods for diversity analysis for community ecology
    + [fantaxtic](https://github.com/gmteunisse/fantaxtic): Fantastic plotting of microbiome data
    + [here](https://here.r-lib.org/): Package for easy file referencing for reproducible workflows
    + [microbiome](https://microbiome.github.io/tutorials/): Tools for microbiome analysis


Another important component of our program list is [conda](https://docs.conda.io/en/latest/), which is a package, dependency, and environment management tool for Unix systems (although it can be used on Windows, too). The two main advantages of conda is 1) to allow you to easily install new programs and tools, and 2) to keep your tools separated into "environments" that can avoid dependency-related problems between different versions of certain packages. This second advantage is something we will see in step 5 of our analysis.

If you want to know more about what this means and how to create and use different environments on conda, [you can go here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).

# 1. Adapter Trimming 

Normally you need to trim the Illumina sequencing adapters from the raw sequencing data in the beginning. We are lucky because the sequencing provider already did this for us. So we can skip this part in our pipeline. The most common tool used for Illumina Sequencing data is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) which was also used for our data. 

  1. What are sequencing adapters and why do we need to remove them? 
  
# 2. Quality Reports

```{bash, eval = F}
# Set the filepath to the location of the raw reads. 
cd /YOUR/FILEPATH/HERE/

# Create a subdirectory to store the FASTQC reports.
mkdir ./FASTQC

# Create a fastqc report of the read files, to examine the read quality.
fastqc 211103_SN234_A_L001_AUXV-4_AdapterTrimmed_R1.fastq.gz -o ./FASTQC &&
fastqc 211103_SN234_A_L001_AUXV-4_AdapterTrimmed_R2.fastq.gz -o ./FASTQC

# Create a multiqc report to have both quality reports in one html file.
multiqc FASTQC -o FASTQC --interactive

```

  2. What are the important variables in the Quality report? What differences do you notice between the two files?

# 3. Demultiplexing
  * we need the barcodes to distinguish the samples 
  * Maybe we need to include code to switch between directories here? 
```{bash, eval = F}
# The .txt file consisting of the octamer tags and primers (reffered to as barcodes from here on), 
# used for the study was previously uploaded to the computer. 

# Take each line of the .txt of the barcodes and adds a line containing info on which sample this barcode belongs to. 
awk '{print">fwd"NR"\n"$0}' ./BARCODES/algae_fwd_barcodes_big.txt > ./BARCODES/algae_barcodes_big.fwd.fasta 

# Take each line of the .txt of the barcodes and adds a line containing info on which sample this barcode belongs to. 
awk '{print">rev"NR"\n"$0}' ./BARCODES/algae_rev_barcodes_big.txt > ./BARCODES/algae_barcodes_big.rev.fasta 

```

## Demultiplexing with cutadapt

```{bash, eval = F}
# Create a subdirectory for the file that have the barcodes and primers removed. And enter it.
mkdir DEMULTIPLEXED

# Activate the cutadapt environment through conda. 
conda activate cutadaptenv 

# Increase the softlimit of the OS because cutadapt will open a lot of files. 
# One file for each forward and reverse combination. These will be filled with the reads containing the combination.
ulimit -S -n 3000


#Cutadapt Main Commands
# Everything needs to be run twice because the reads are in mixed orientation,
# because of the PCR free library preparation. 
# Because of the dual indexing approach we need to supply two barcode files. 
cutadapt \
-j 0 \
-e 0.15 --no-indels --minimum-length 50 \
-g file:./BARCODES/algae_barcodes_big.fwd.fasta \
-G file:./BARCODES/algae_barcodes_big.rev.fasta \
-o ./DEMULTIPLEXED/{name1}-{name2}.round1.1.fastq \
-p ./DEMULTIPLEXED/{name1}-{name2}.round1.2.fastq \
./211103_SN234_A_L001_AUXV-4_AdapterTrimmed_R1.fastq.gz \
./211103_SN234_A_L001_AUXV-4_AdapterTrimmed_R2.fastq.gz > algae_round1_cutadapt.txt &&

cutadapt \
-j 0 \
-e 0.15 --no-indels --minimum-length 50 \
-g file:./BARCODES/algae_barcodes_big.fwd.fasta \
-G file:./BARCODES/algae_barcodes_big.rev.fasta \
-o ./DEMULTIPLEXED/{name1}-{name2}.round2.1.fastq \
-p ./DEMULTIPLEXED/{name1}-{name2}.round2.2.fastq \
./DEMULTIPLEXED/unknown-unknown.round1.2.fastq \
./DEMULTIPLEXED/unknown-unknown.round1.1.fastq > algae_round2_cutadapt.txt

# Close the conda environment.
conda deactivate 
```

  3. Which differences do you notice between the two cutadapt commands? And why do we need to use two different commands? 
  
  * we had two passes through cutadapt so we need to merge the files 

First we rename the different samples. 
```{bash, eval = F}
# Before merging the files we need to rename them to make sorting easier and to only include reads that
# are real samples, multiplex controls, blanks and PCR negative controls. The text files were uploaded to the 
# RENAMING/ directory before. 

# We need to enter the directory where we stored the demultiplexed files. 
cd ./DEMULTIPLEXED/

# Then we rename the files. 
paste /YOUR/FILEPATH/HERE/RENAMING/renaming_R1_1_old_big.txt \
/YOUR/FILEPATH/HERE/RENAMING/renaming_R1_1_new_big.txt \
| while read n k; do rename -v $n $k * ; done > ./rename_logR1.1.txt

paste /YOUR/FILEPATH/HERE/RENAMING/renaming_R1_2_old_big.txt \
/YOUR/FILEPATH/HERE/RENAMING/renaming_R1_2_new_big.txt \
| while read n k; do rename -v $n $k * ; done > ./rename_logR1.2.txt

paste /YOUR/FILEPATH/HERE/RENAMING/renaming_R2_1_old_big.txt \
 /YOUR/FILEPATH/HERE/RENAMING/renaming_R2_1_new_big.txt \
| while read n k; do rename -v $n $k * ; done > ./rename_logR2.1.txt

paste /YOUR/FILEPATH/HERE/RENAMING/renaming_R2_2_old_big.txt \
 /YOUR/FILEPATH/HERE/RENAMING/renaming_R2_2_new_big.txt \
 | while read n k; do rename -v $n $k * ; done > ./rename_logR2.2.txt

```

For a later check of how many reads are passing the pipeline we need to create a FASTA file and count the reads. 

```{r}
# Transform the fastq file to a fasta file. 
cat A31_B1.round1.1.sample.fastq | \
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.round1.1.sample.fa

cat A31_B1.round2.1.sample.fastq | \
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.round2.1.sample.fa

# Grep the read numbers by counting the lines beginning with >.
grep -c '^>' *.fa | less 

```

  4. What is the difference between FASTQ and FASTA files? 
  
Now merge the two files. 

```{bash, eval = F}
# Create a subdirectory for the merged files. 
mkdir MERGED

# First create two lists of the filenames for the pairs from the cutadapt results. 
ls -1 *round1.1.sample.fastq | sed 's/round1.1.sample.fastq//' > listround1.1.

ls -1 *round2.1.sample.fastq | sed 's/round2.1.sample.fastq//' > listround2.1.

# Now we can merge the pairs. First for the files from R1.   
paste listround1.1. listround2.1. | while read n k; \
do cat $n"round1.1.sample.fastq" $k"round2.1.sample.fastq" > ./MERGED/$n"sample_demux.1.fastq"; done

# And again for the R2 reads. 
ls -1 *round1.2.sample.fastq | sed 's/round1.2.sample.fastq//'  > listround1.2.

ls -1 *round2.2.sample.fastq | sed 's/round2.2.sample.fastq//' > listround2.2.

paste listround1.2. listround2.2. | while read n k; \
do cat $n"round1.2.sample.fastq" $k"round2.2.sample.fastq" > ./MERGED/$n"sample_demux.2.fastq"; done

# To check if the merging has worked, we create a FASTA file of the merged sample 
# and check if the read numbers match the paired files. 

# Transform the fastq file to a fasta file. 
cat ./MERGED/A31_B1.sample_demux.1.fastq | \
 awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.sample_demux.1.fa

# Grep the read numbers by counting the lines beginning with >.
grep -c '^>' *sample*.fa | less 

# Exit the DEMULTIPLEXED subdirectory. 
cd ..

```

# 4. Checking for leftover primer sequences.

  * To run DADA2 we need to activate R 

Open R Studio

```{r, eval = F}
# Then we load all required packages. 

library(dada2) ; packageVersion('dada2')

library(ShortRead) ; packageVersion('ShortRead')

library(Biostrings) ; packageVersion('Biostrings')

# Create objects that contain the primer sequences.  
# For the algae that is:

FWD <- 'GTGARTCATCGAATCTTTG'
REV <- 'TCCTCCGCTTATTGATATGC'

# Make a custom function that creates all the possible orientations of the primers e.g. complement, reverse complement.
allOrients <- function(primer) {
  require(Biostrings)
  # Create all orientations of the input sequence
  dna     <- DNAString(primer)  # turn character to DNAString object
  orients <- c(Forward=dna, Complement=complement(dna), Reverse=reverse(dna),
               RevComp=reverseComplement(dna))
  return(sapply(orients, toString))  # back to character vector
}

# Make and save the orientation files.
FWD.orients <- allOrients(FWD)
FWD.orients

REV.orients <- allOrients(REV)
REV.orients

# Load in the demultiplexed files. 

fnFs <- sort(list.files(path = './DEMULTIPLEXED/MERGED', pattern = "sample_demux.1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path = './DEMULTIPLEXED/MERGED', pattern = "sample_demux.2.fastq", full.names = TRUE))

# Filter out ambiguous Ns with the filterAndTrim function setting maxN to zero.
# Place the N-filterd files into a filtN/ subdirectory.
fnFs.filtN <- file.path(path = './DEMULTIPLEXED/MERGED', "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path = './DEMULTIPLEXED/MERGED', "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Check for any leftover primers after the removal with Cutadapt.

# Create a function that counts the number of reads in which the primer is found.
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

# Search through all the reads and combine in a dataframe.
# If the samples come from the same library prep then it is enough to only process one of the files 
# (see the [1] at the end of the command). 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

```
  
  5. How many reads still contain primer sequences? 

Move out of R and into the terminal again. 
Go through a second pass of cutadapt to remove the remaining primers. 

```{bash, eval = F}
# Open the cutadapt environment.
# Remove leftover primers with Cutadapt
conda activate cutadaptenv

# Create a subdirectory for the first primer removal.
mkdir PRIMER_REMOVED1 

# Enter the MERGED subdirectory.
cd DEMULTIPLEXED/MERGED/

# The command below contains all possible orientations: 
# 1) fwd-rcrev + rev-rcfwd; 2) rcfwd-rev + rcrev-fwd; \
# 3) fwd + rcfwd; 4) rcfwd + fwd; 5) rev + rcrev; 6) rcrev + rev

# Use a for loop to run Cutadapt over all samples. 
ls *demu*.fastq | cut -f1 -d'.' > samples

for sample in $(cat samples); do

echo "On sample: $sample"

cutadapt --cores=0 --minimum-length 50 \
 -a ^GTGARTCATCGAATCTTTG...GCATATCAATAAGCGGAGGA -A ^TCCTCCGCTTATTGATATGC...CAAAGATTCGATGAYTCAC \
 -a ^CAAAGATTCGATGAYTCAC...TCCTCCGCTTATTGATATGC -A ^GCATATCAATAAGCGGAGGA...GTGARTCATCGAATCTTTG \
 -a GTGARTCATCGAATCTTTG -A CAAAGATTCGATGAYTCAC \
 -a CAAAGATTCGATGAYTCAC -A GTGARTCATCGAATCTTTG \
 -a TCCTCCGCTTATTGATATGC -A GCATATCAATAAGCGGAGGA \
 -a GCATATCAATAAGCGGAGGA -A TCCTCCGCTTATTGATATGC \
 -o ../../PRIMER_REMOVED1/${sample}.sample_demux_prirm.1.fastq \
 -p ../../PRIMER_REMOVED1/${sample}.sample_demux_prirm.2.fastq \
 ${sample}.sample_demux.1.fastq ${sample}.sample_demux.2.fastq \
 > cutadapt_primer_trimming_stats_{sample}.txt

done 

# Close cutadapt environment.
conda deactivate 

# Exit the subdirectory and go to our base directory for the fungal reads.
cd ../..

```
  6. How does the cutadapt call for the demultiplexing step differ from the primer removal? 
  
Now we check for primers again. 

```{r, eval = F}
# Then we load all required packages. 

library(dada2) ; packageVersion('dada2')

library(ShortRead) ; packageVersion('ShortRead')

library(Biostrings) ; packageVersion('Biostrings')

# Create objects that contain the primer sequences.  
# For the algae that is:

FWD <- 'GTGARTCATCGAATCTTTG'
REV <- 'TCCTCCGCTTATTGATATGC'

# Make a custom function that creates all the possible orientations of the primers e.g. complement, reverse complement.
allOrients <- function(primer) {
  require(Biostrings)
  # Create all orientations of the input sequence
  dna     <- DNAString(primer)  # turn character to DNAString object
  orients <- c(Forward=dna, Complement=complement(dna), Reverse=reverse(dna),
               RevComp=reverseComplement(dna))
  return(sapply(orients, toString))  # back to character vector
}

# Make and save the orientation files.
FWD.orients <- allOrients(FWD)
FWD.orients

REV.orients <- allOrients(REV)
REV.orients

# Load in the demultiplexed files. 

fnFs <- sort(list.files(path = './DEMULTIPLEXED/MERGED', pattern = "sample_demux.1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path = './DEMULTIPLEXED/MERGED', pattern = "sample_demux.2.fastq", full.names = TRUE))

# Filter out ambiguous Ns with the filterAndTrim function setting maxN to zero.
# Place the N-filterd files into a filtN/ subdirectory.
fnFs.filtN <- file.path(path = './DEMULTIPLEXED/MERGED', "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path = './DEMULTIPLEXED/MERGED', "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Check for any leftover primers after the removal with Cutadapt.

# Create a function that counts the number of reads in which the primer is found.
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

# Search through all the reads and combine in a dataframe.
# If the samples come from the same library prep then it is enough to only process one of the files 
# (see the [1] at the end of the command). 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

```

  7. Was the primer removal successful? Or do we need to go through cutadapt again? 

# 5. Sample inference with DADA2

Load in the files. 
```{r, eval = F}
# For the fungi we need to use the _prirm.fastq files. 

cutFs_1 <- sort(list.files("./PRIMER_REMOVED1", pattern = "1.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_1 <- sort(list.files("./PRIMER_REMOVED1", pattern = "1.sample_demux_prirm.2.fastq", full.names = TRUE))

cutFs_2 <- sort(list.files("./PRIMER_REMOVED1", pattern = "2.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_2 <- sort(list.files("./PRIMER_REMOVED1", pattern = "2.sample_demux_prirm.2.fastq", full.names = TRUE))

cutFs_3 <- sort(list.files("./PRIMER_REMOVED1", pattern = "3.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_3 <- sort(list.files("./PRIMER_REMOVED1", pattern = "3.sample_demux_prirm.2.fastq", full.names = TRUE))

# Create a function to obtain the sample names. "1.s" serves as the point at which the string will be split. 
get.sample.name <- function(fname) strsplit(basename(fname), "1.s")[[1]][1]

# Get the sample names.
sample.names <- unname(sapply(cutFs_1, get.sample.name))
head(sample.names)

```

Filter and trim the reads for quality.

```{r, eval = F}
# Assign filenames for the output of the filtered reads. 
filtFs_1 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutFs_1))
filtRs_1 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutRs_1))

filtFs_2 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutFs_2))
filtRs_2 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutRs_2))

filtFs_3 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutFs_3))
filtRs_3 <- file.path("./PRIMER_REMOVED1", "filtered", basename(cutRs_3))

# Apply the filtering parameters , no truncLen because these are ITS reads
# and therefore very variable in length, changed maxEE to 6,6 since the libraries are in mixed orientation.
out_1 <- filterAndTrim(cutFs_1, filtFs_1, cutRs_1, filtRs_1, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

out_2 <- filterAndTrim(cutFs_2, filtFs_2, cutRs_2, filtRs_2, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
	
out_3 <- filterAndTrim(cutFs_3, filtFs_3, cutRs_3, filtRs_3, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)	

# Check if a good number of reads passed the quality filtering. We are filtering out ~30% which is okay. 
head(out_1)
head(out_2)
head(out_3)

```

Learn the error profiles from the reads. 

```{r, eval = F}
# Learn the error rates for the R1 reads. 
errF_1 <- learnErrors(filtFs_1, multithread = TRUE)

errF_2 <- learnErrors(filtFs_2, multithread = TRUE)

errF_3 <- learnErrors(filtFs_3, multithread = TRUE)

# Learn the error rates for the R2 reads. 
errR_1 <- learnErrors(filtRs_1, multithread = TRUE)

errR_2 <- learnErrors(filtRs_2, multithread = TRUE)

errR_3 <- learnErrors(filtRs_3, multithread = TRUE)

```

De-replicate identical reads. 

```{r, eval = F}
# De-replicate identical reads. 

derepFs_1 <- derepFastq(filtFs_1, verbose = TRUE)
derepRs_1 <- derepFastq(filtRs_1, verbose = TRUE)

derepFs_2 <- derepFastq(filtFs_2, verbose = TRUE)
derepRs_2 <- derepFastq(filtRs_2, verbose = TRUE)

derepFs_3 <- derepFastq(filtFs_3, verbose = TRUE)
derepRs_3 <- derepFastq(filtRs_3, verbose = TRUE)

# Name the derep-class objects by the sample names.

names(derepFs_1) <- sample.names
names(derepRs_1) <- sample.names

names(derepFs_2) <- sample.names
names(derepRs_2) <- sample.names

names(derepFs_3) <- sample.names
names(derepRs_3) <- sample.names

```

  8. Why do we need to do everything three times? 

Run the DADA2 core for sample inference (i.e. call the ASVs). 
  
```{r, eval = F}
dadaFs_1 <- dada(derepFs_1, err = errF_1, multithread = TRUE)
dadaRs_1 <- dada(derepRs_1, err = errR_1, multithread = TRUE)

dadaFs_2 <- dada(derepFs_2, err = errF_2, multithread = TRUE)
dadaRs_2 <- dada(derepRs_2, err = errR_2, multithread = TRUE)

dadaFs_3 <- dada(derepFs_3, err = errF_3, multithread = TRUE)
dadaRs_3 <- dada(derepRs_3, err = errR_3, multithread = TRUE)

# Merge the paired reads within the replications.
mergers_1 <- mergePairs(dadaFs_1, derepFs_1, dadaRs_1, derepRs_1, verbose=TRUE)

mergers_2 <- mergePairs(dadaFs_2, derepFs_2, dadaRs_2, derepRs_2, verbose=TRUE)

mergers_3 <- mergePairs(dadaFs_3, derepFs_3, dadaRs_3, derepRs_3, verbose=TRUE)

# Construct the ASV table per replicate. 
seqtab_1 <- makeSequenceTable(mergers_1)
dim(seqtab_1)

seqtab_2 <- makeSequenceTable(mergers_2)
dim(seqtab_2)

seqtab_3 <- makeSequenceTable(mergers_3)
dim(seqtab_3)

```

Remove chimeric sequences created during PCR. 

```{r, eval = F}
# Chimera removal for the replicates. 
seqtab_1.nochim <- removeBimeraDenovo(seqtab_1, method="consensus", multithread=TRUE, verbose=TRUE)

seqtab_2.nochim <- removeBimeraDenovo(seqtab_2, method="consensus", multithread=TRUE, verbose=TRUE)

seqtab_3.nochim <- removeBimeraDenovo(seqtab_3, method="consensus", multithread=TRUE, verbose=TRUE)

```

  9. What are chimeric sequences? 

Because our reads are in mixed orientation we need to check for artificial diversity. 

```{r, eval = F}
# Check for reverse complement synthetic diversity. 
# Because the libraries are in mixed orientation we need to check for identical sequences, 
# which are read in reverse complement orientation.
sq_1 <- getSequences(seqtab_1.nochim)
sq.rc_1 <- dada2:::rc(sq_1)
rcdupes_1 <- sapply(seq_along(sq_1), function(i) {
    sq.rc_1[[i]] %in% sq_1[1:(i-1)]
})

sq_2 <- getSequences(seqtab_2.nochim)
sq.rc_2 <- dada2:::rc(sq_2)
rcdupes_2 <- sapply(seq_along(sq_2), function(i) {
    sq.rc_2[[i]] %in% sq_2[1:(i-1)]
})

sq_3 <- getSequences(seqtab_3.nochim)
sq.rc_3 <- dada2:::rc(sq_3)
rcdupes_3 <- sapply(seq_along(sq_3), function(i) {
    sq.rc_3[[i]] %in% sq_3[1:(i-1)]
})

# Merge the forward and reverse-complement reads.
colnames(seqtab_1.nochim)[rcdupes_1] <- dada2:::rc(colnames(seqtab_1.nochim)[rcdupes_1])
stm_1 <- collapseNoMismatch(seqtab_1.nochim)

colnames(seqtab_2.nochim)[rcdupes_2] <- dada2:::rc(colnames(seqtab_2.nochim)[rcdupes_2])
stm_2 <- collapseNoMismatch(seqtab_2.nochim)

colnames(seqtab_3.nochim)[rcdupes_3] <- dada2:::rc(colnames(seqtab_3.nochim)[rcdupes_3])
stm_3 <- collapseNoMismatch(seqtab_3.nochim)

```

Merge the resulting tables so we have one ASV table for all fungi from the three technical replicates. 

```{r, eval = F}
# Merge the two ASV tables.
# Put the tables in a list before merging them with mergeSequenceTables.
input_tables <- list(stm_1, stm_2, stm_3)

seqtab_merge <- mergeSequenceTables(tables = input_tables, repeats = 'sum', tryRC = TRUE)

# An additional run to remove chimeric sequences. 
seqtab_merge.nochim <- removeBimeraDenovo(seqtab_merge, method="consensus", multithread=TRUE, verbose=TRUE)

# Inspect the sequence lengths to look for anything that seems suspicious. Not the case here.
table(nchar(getSequences(seqtab_merge.nochim)))

# Give the sequence variants more manageable names, e.g. ASVn  
asv_seqs <- colnames(seqtab_merge.nochim)
asv_headers <- vector(dim(seqtab_merge.nochim)[2], mode="character")

for (i in 1:dim(seqtab_merge.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

```

  10. What is the range of the sequence length of our amplicons? 
  
Track our reads through the pipeline so we can see how many reads are filtered out. 

```{r, eval = F}
# Track the reads through the pipeline.
# This is a last sanity check to see if we are not loosing samples at an unexpected step. 

getN <- function(x) sum(getUniques(x))
track_1 <- cbind(out_1, sapply(dadaFs_1, getN), sapply(dadaRs_1, getN), sapply(mergers_1, getN), rowSums(seqtab_1.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_1) <- sample.names

track_2 <- cbind(out_2, sapply(dadaFs_2, getN), sapply(dadaRs_2, getN), sapply(mergers_2, getN), rowSums(seqtab_2.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_2) <- sample.names

track_3 <- cbind(out_3, sapply(dadaFs_3, getN), sapply(dadaRs_3, getN), sapply(mergers_3, getN), rowSums(seqtab_3.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_3) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_3) <- sample.names

head(track_1)
head(track_2)
head(track_3)

```

  11. How many sequences did we loose at each step? 
  
Now we go on to save our data. 

```{r, eval = F}
# Make and save a fasta of our final ASV sequences.
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_fungi.fa")

# Write our chimera screened and merge sequenc table to a text file. This will be our ASV table to work with from here. 
write.table(t(seqtab_merge.nochim), "asv_table_fungi.txt", sep="\t")

# Also save it as an .rds file to not corrupt any structure when we need to reload the table.
saveRDS(seqtab_merge.nochim, 'asv_table_fungi.rds')

```

# 6. Taxonomy assignment 

* DADA2 built in functions for assigning taxonomy 
* We work with the UNITE database of fungi

```{r, eval = F}
# Read in the UNITE database fasta.  
unite.ref <- './sh_general_release_dynamic_all_25.07.2023.fasta'  

# Run the taxonomy assignment on the ASV table. 
taxa <- assignTaxonomy(seqtab_merge.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# Save the taxonomy table as an R object.
saveRDS(taxa, 'tax_table_fungi.rds')

```

# 7. Match List 
* We will later run a post clustering algorithm and need to find similar sequences from our fasta file 

```{bash, eval = F}
# Open the conda environment containing blast. 
conda activate blastenv

# Create a database containing the algal ASVs in fasta format.
makeblastdb -in ASVs_fungi.fa -parse_seqids -dbtype nucl

# Compare all ASVs to each other and create a match list of sequences that are similar to each other.
blastn -db ASVs_fungi.fa \
-outfmt '6 qseqid sseqid pident' \
-out match_list_fungi.txt \
-qcov_hsp_perc 80 \
-perc_identity 84 \
-num_threads 20 \
-query ASVs_fungi.fa

conda deactivate

```

# 8. ASV curation and removal of potential contaminants

First we need to load additional packages. 

```{r, eval = F}
library(here)

library(decontam)

library(phyloseq)

library(lulu)

library(Biostrings)

library(tidyverse)

```

## Decontamination

Load in the data.

```{r, eval = F}
# Load ASV table for fungi (available as supplementary data).
fungi_asv <- utils::read.csv(here("Data", 'asv_table_fungi.txt'), header = T, sep = '\t')

# Load FASTA file to get ASV_IDs.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the ASV table and the representative sequences.
fungi_asv_IDs <- dplyr::left_join(fungi_rep_seqs, fungi_asv, by = 'sequence_fungi')

# Naming the Samples. 
fungi_asv_IDs <- fungi_asv_IDs %>% 
  dplyr::rename_with(~base::paste0("Sample_", .), -c(1:2))

# set ASV ID as the rownames.
base::rownames(fungi_asv_IDs) <- fungi_asv_IDs$seq_name_fungi
fungi_asv_IDs$seq_name_fungi <- NULL

# Load the DNA concentration necessary for decontam (available as supplementary data). 
fungi_conc <- utils::read.csv(here("Data", 'fungi_decontam_conc_big.csv'), header = T, sep =',')
base::rownames(fungi_conc) <- base::paste0("Sample_", fungi_conc$sample_ID)
base::rownames(fungi_conc) <- base::sub("[-]", "_", x = base::rownames(fungi_conc))

```

Combine them into a so-called phyloseq object. 

```{r, eval = F}
# Make the parts of the phyloseq object.
ASV_mat_decontam_fun <- base::data.matrix(fungi_asv_IDs)
ASV_fungi_decontam <- phyloseq::otu_table(ASV_mat_decontam_fun, taxa_are_rows = TRUE)
sampledata_fungi_decontam <- phyloseq::sample_data(fungi_conc)

# Combine with phyloseq
ps_fungi_decontam <- phyloseq::phyloseq(ASV_fungi_decontam, sampledata_fungi_decontam)
ps_fungi_decontam
```

  12. How many ASVs does our dataset contain before filtering out contaminants?
  
Run the decontamination algorithm. 

```{r, eval = F}
# Check the library sizes of the Samples and Controls.
df_fungi <- base::as.data.frame(phyloseq::sample_data(ps_fungi_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_fungi$LibrarySize <- phyloseq::sample_sums(ps_fungi_decontam)
df_fungi <- df_fungi[base::order(df_fungi$LibrarySize),]
df_fungi$Index <- base::seq(base::nrow(df_fungi))
ggplot2::ggplot(data=df_fungi, ggplot2::aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  ggplot2::geom_point()

# Decontam with combined frequency and prevalence approach. 
phyloseq::sample_data(ps_fungi_decontam)$is.neg <- phyloseq::sample_data(ps_fungi_decontam)$Sample_or_Control == "Control Sample"
contam_fungi_combi <- decontam::isContaminant(ps_fungi_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
base::table(contam_fungi_combi$contaminant)
base::which(contam_fungi_combi$contaminant)

# Trim the identified contaminants from the phyloseq object 
ps_fungi_noncontam <- phyloseq::prune_taxa(!contam_fungi_combi$contaminant,
                                 ps_fungi_decontam)
ps_fungi_noncontam

# Remove taxa without reads.
ps_fungi_noncontam_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_fungi_noncontam) > 0,
                                        ps_fungi_noncontam)
```

  13. Which taxa were removed from the original datasets as contaminants? 
  14. How many taxa remain after this step?
  15. What is the difference between the frequency, prevalence and combined approach in the decontam function?
  
## LULU curation

* LULU is an algorithm to cluster similar reads to obtain more accurate diversity estimates

```{r, eval = F}
ASV_table_fungi <- base::as.data.frame(phyloseq::otu_table(ps_fungi_noncontam_pruned))
fungi_matchlist <- utils::read.table(here::here("Data", 'match_list_fungi.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_fungi_cur <- lulu::lulu(ASV_table_fungi, fungi_matchlist)

ASV_table_fungi_cur$curated_count
ASV_table_fungi_cur$discarded_count

```

  16. How many reads were merged? 
  17. Which parameters does the LULU algorithm consider in it's merging? 
  18. What are the default parameters? 
  
# Diversity Analysis

# 9. Data initialization 

Load in the taxonomy table. 

```{r, eval = F}
# Load in taxonomy table. 
tax_fungi <- base::readRDS(here::here("tax_table_fungi.rds")) %>% 
  base::as.data.frame() %>%
  tibble::rownames_to_column('sequence') %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here('ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

```

Clean the taxonomy table so it looks nicer and is easier to handle. 

```{r, eval = F}
# Split the taxonomy into different columns of taxonomic levels.
fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__') %>% 
  tidyr::separate(Phylum, c(NA, 'Phylum') , sep = '__') %>% 
  tidyr::separate(Class, c(NA, 'Class') , sep = '__') %>% 
  tidyr::separate(Order, c(NA, 'Order') , sep = '__') %>% 
  tidyr::separate(Family, c(NA, 'Family') , sep = '__') %>% 
  tidyr::separate(Genus, c(NA, 'Genus') , sep = '__') %>% 
  tidyr::separate(Species, c(NA, 'Species') , sep = '__')

# Keep only Fungi 
fungi_tax_fin <- fungi_tax_fin %>% 
  dplyr::filter(Kingdom == "Fungi")

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)

# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL

```

Load in the ASV table we curated. 

```{r, eval = F}
# Load in ASV table previously curated using the LULU algorithm (https://doi.org/10.1038/s41467-017-01312-x).
# For code on how it was employed see https://github.com/LukDrey/beech_micro_communities. 
ASV_table_fungi_cur <- base::readRDS(here::here("ASV_table_fungi_cur.rds"))

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>% 
  dplyr::select(all_of(base::rownames(metadata)))
```

Load in some metadata. 

```{r, eval = F}
# Load in Metadata. 
metadata <- utils::read.csv(here::here('sample_data_3_AHS.csv'), sep = ';')

# Set the sample name as the rowname for the phyloseq creation.
base::row.names(metadata) <- metadata$sample
metadata$sample <- NULL

```

Combine the data into a phyloseq object.  

```{r, eval = F}
# First create a phyloseq object with all samples included.
# Create the matrices needed for the phyloseq functions.  
asv_mat <- base::as.matrix(asv_fungi)
taxmat <- base::as.matrix(fungi_tax_fin) 
  
# Create the ps object. 
ASV <- phyloseq::otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxmat)
sampledata <- phyloseq::sample_data(metadata)

full_physeq <- phyloseq::phyloseq(ASV, TAX, sampledata)

```

  19. How many taxa does our phyloseq object contain now? 
  
For our analysis we want to focus on the differences between soil and bark samples as well as differences between coniferous and deciduous trees. So we need to do a lot of splitting and filtering of the phyloseq objects. The coniferous tree species is not the same in the two regions of interest. It is *Picea abies* (Norway Spruce) in the Swabian Alb and *Pinus sylvestris* (Scots pine) in Schorfheide-Chorin, so we need to treat them differently. 

```{r, eval = F}
# Filter out samples that were sampled on trees that are not from our three target species. 

filtered_physeq <- phyloseq::subset_samples(full_physeq, 
                                            dominant_tree %in% c("Fagus_sylvatica",
                                                                 "Pinus_sylvestris", 
                                                                 "Picea_abies")) %>% 
  phyloseq::subset_samples(exploratory %in% c("Alb", "Schorfheide")) %>%
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories.

physeq_alb <- phyloseq::subset_samples(filtered_physeq, exploratory == "Alb") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch <- phyloseq::subset_samples(filtered_physeq, exploratory == "Schorfheide") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into bark and soil samples. 
physeq_bark <- phyloseq::subset_samples(filtered_physeq, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_soil <- phyloseq::subset_samples(filtered_physeq, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories & the substrate.
# Swabian Alb 
physeq_alb_bark <- phyloseq::subset_samples(physeq_alb, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil <- phyloseq::subset_samples(physeq_alb, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide-Chorin
physeq_sch_bark <- phyloseq::subset_samples(physeq_sch, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil <- phyloseq::subset_samples(physeq_sch, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseqobject by exploratories & the tree species. 
# Swabian Alb
physeq_alb_fagus <- phyloseq::subset_samples(physeq_alb, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_picea <- phyloseq::subset_samples(physeq_alb, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide Chorin
physeq_sch_fagus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_pinus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories & the substrate & the tree species.
# Swabian Alb
physeq_alb_bark_fagus <- phyloseq::subset_samples(physeq_alb_bark, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil_fagus <- phyloseq::subset_samples(physeq_alb_soil, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_bark_picea <- phyloseq::subset_samples(physeq_alb_bark, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil_picea <- phyloseq::subset_samples(physeq_alb_soil, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide-Chorin
physeq_sch_bark_fagus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_fagus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_bark_pinus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_pinus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

```

  20. How many taxa do the different splits contain? Make a little table. 
  
# 10. Analysis of library sizes

```{r; eval = F}
################Swabian Alb###########################

# Create a sample_data column containing a column that corresponds to tree species and substrate.
physeq_alb_curve <- physeq_alb
phyloseq::sample_data(physeq_alb_curve) <- phyloseq::sample_data(physeq_alb) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))

# Create rarefaction curve and color the lines by substrate (bark/soil) and host tree species. 
rare_alb <- ranacapa::ggrare(physeq_alb_curve, step = 50, color = "tree_substrate", se = FALSE) +
  theme(legend.text = element_text(size = 4),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))

################Schorfheide###########################

# Create a sample_data column containing a column that corresponds to tree species and substrate.
physeq_sch_curve <- physeq_sch
phyloseq::sample_data(physeq_sch_curve) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))

# Create rarefaction curve and color the lines by substrate (bark/soil) and host tree species. 
rare_sch <- ranacapa::ggrare(physeq_sch_curve, step = 50,
                             color = "tree_substrate", se = FALSE) +
  theme(legend.text = element_text(size = 4),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))

combined_rare_curves <- ggpubr::ggarrange(rare_alb, rare_sch,
                                          ncol = 2, nrow = 1)
combined_rare_curves

```

  21. What differences can we see from these rarefaction curves?
  
# 12. Community Composition

* Which taxa are in our dataset

```{r, eval = F}
############## Swabian Alb ###############
physeq_alb_barplot <- physeq_alb
phyloseq::sample_data(physeq_alb_barplot) <- phyloseq::sample_data(physeq_alb) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_alb_ord_top25 <- fantaxtic::top_taxa(physeq_alb_barplot,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 
phy_alb_ord_top25_named <- fantaxtic::name_na_taxa(phy_alb_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_alb_ord_top25_named_plot <-  phy_alb_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
phyloseq::taxa_names(phy_alb_ord_top25_named_plot) <- phyloseq::tax_table(phy_alb_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_alb_ord <- base::sort(phyloseq::taxa_names(phy_alb_ord_top25_named_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_alb <- base::data.frame(phyloseq::sample_data(phy_alb_ord_top25_named_plot))
sampledata_alb <- sampledata_alb %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-bark", "F. sylvatica bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-soil", "F. sylvatica soil")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Picea_abies-bark", "P. abies bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Picea_abies-soil", "P. abies soil")) 
sampledata_alb$tree_substrate <- factor(sampledata_alb$tree_substrate, 
                                       levels = c("F. sylvatica bark", "P. abies bark",
                                                  "F. sylvatica soil", "P. abies soil"))  

phyloseq::sample_data(phy_alb_ord_top25_named_plot) <- phyloseq::sample_data(sampledata_alb)


############## Schorfheide-Chorin ###############
physeq_sch_barplot <- physeq_sch
phyloseq::sample_data(physeq_sch_barplot) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_sch_ord_top25 <- fantaxtic::top_taxa(physeq_sch_barplot,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 
phy_sch_ord_top25_named <- fantaxtic::name_na_taxa(phy_sch_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_sch_ord_top25_named_plot <-  phy_sch_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
phyloseq::taxa_names(phy_sch_ord_top25_named_plot) <- phyloseq::tax_table(phy_sch_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_sch_ord <- sort(phyloseq::taxa_names(phy_sch_ord_top25_named_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_sch <- data.frame(phyloseq::sample_data(phy_sch_ord_top25_named_plot))
sampledata_sch <- sampledata_sch %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-bark", "F. sylvatica bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-soil", "F. sylvatica soil")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-bark", "P. sylvestris bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-soil", "P. sylvestris soil"))
          
sampledata_sch$tree_substrate <- factor(sampledata_sch$tree_substrate, 
                                        levels = c("F. sylvatica bark", "P. sylvestris bark",
                                                   "F. sylvatica soil", "P. sylvestris soil"))  

phyloseq::sample_data(phy_sch_ord_top25_named_plot) <- phyloseq::sample_data(sampledata_sch)


#################################################################
##                          Section 8.2                        ##
##              Create great looking plots                     ##
#################################################################

my_cols <- Polychrome::createPalette(33, seedcolors = c("#ff0000", "#00ff00", "#0000ff"))

unique_orders <- sort(unique(c(taxa_names_alb_ord, taxa_names_sch_ord))) 

custom_sort <- c("Agaricales", "Archaeorhizomycetales",
                 "Atheliales", "Boletales",
                 "Caliciales", "Cantharellales",
                 "Capnodiales", "Chaetothyriales",
                 "Eurotiales", "Filobasidiales",
                 "Helotiales", "Hypocreales",
                 "Lecanorales", "Mortierellales",
                 "Mycosphaerellales", "Mytilinidales",
                 "Orbiliales", "Pezizales", "Phaeothecales",
                 "Pleosporales", "Russulales",
                 "Sebacinales", "Thelebolales",
                 "Thelephorales", "Trapeliales", 
                 "Tremellales", "Umbelopsidales", "Verrucariales",
                 "Unknown Ascomycota (Phylum)", "Unknown Dothideomycetes (Class)",
                 "Unknown Fungi (Kingdom)", "Unknown Rozellomycota (Phylum)",
                              "Other")

full_cols <- data.frame(order = custom_sort, color = my_cols)

alb_cols <- full_cols %>% 
  dplyr::filter(order %in% taxa_names_alb_ord)

sch_cols <- full_cols %>% 
  dplyr::filter(order %in% taxa_names_sch_ord)


scale_fill_alb <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(alb_cols$color, alb_cols$order), 
    ...
  )
}

scale_fill_sch <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(sch_cols$color, sch_cols$order), 
    ...
  )
}

############## Swabian Alb ###############

# Custom plotting to make a nice stacked barplot. 
alb_ord_soil_plots <- phyloseq::subset_samples(phy_alb_ord_top25_named_plot, substrate == "soil") %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = alb_cols$order) +
  scale_fill_alb() +
  guides(fill = guide_legend(title.position = 'top', nrow = 3)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(B)")  
alb_ord_soil_plots

# Custom plotting to make a nice stacked barplot. 
alb_ord_bark_plots <- phyloseq::subset_samples(phy_alb_ord_top25_named_plot, substrate == "bark") %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = alb_cols$order) +
  scale_fill_alb() +
  guides(fill = guide_legend(title = 'Order',title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = 'bottom',
        text = element_text(colour = 'black', size = 20), 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 5),
        legend.title =  element_text("Order",size = 8),
        legend.key.size = unit(1, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(A) Swabian Alb")
alb_ord_bark_plots

############## Schorfheide-Chorin ###############

# Custom plotting to make a nice stacked barplot. 
sch_ord_soil_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "soil") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = sch_cols$order) +
  scale_fill_sch() +
  guides(fill = guide_legend(title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(B)")   
sch_ord_soil_plots

sch_ord_bark_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "bark") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = sch_cols$order) +
  scale_fill_sch() +
  guides(fill = guide_legend(title = 'Order',title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = 'bottom',
        text = element_text(colour = 'black', size = 20), 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 5),
        legend.title =  element_text("Order",size = 8),
        legend.key.size = unit(1, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(A) Schorfheide-Chorin") 
sch_ord_bark_plots

######## Create final arranged plot #######

final_alb_community_barplot <- ggpubr::ggarrange(alb_ord_bark_plots, alb_ord_soil_plots,
                                                 ncol = 1, nrow = 2, hjust = 5,
                                                 legend = "bottom", common.legend = TRUE)
final_alb_community_barplot

ggsave('final_alb_community_barplot.jpeg', device = 'jpeg',
       final_alb_community_barplot, width = 180, height = 140,
       units = 'mm', dpi = 300)

final_sch_community_barplot <- ggpubr::ggarrange(sch_ord_bark_plots, sch_ord_soil_plots,
                                                 ncol = 1, nrow = 2,
                                                 legend = "bottom", common.legend = TRUE)
final_sch_community_barplot

ggsave('final_sch_community_barplot.jpeg', device = 'jpeg',
       final_sch_community_barplot, width = 180, height = 140,
       units = 'mm', dpi = 300)


combined_community_barplot <- ggpubr::ggarrange(final_alb_community_barplot,
                                                final_sch_community_barplot, 
                                                ncol = 1, nrow = 2)
combined_community_barplot
```

# 13. Differences in community composition/beta diversity between bark and soil

* We can see that there are differences between trees but are they really there? 

```{r, eval = F}
#########################NMDS-Ordination############################
# Ordinate the Swabian Alb phyloseq using an NMDS with Bray-Curtis distance.
nmds_alb <- phyloseq::ordinate(physeq_alb, method = "NMDS", distance = "bray")

# Plot the ordination. 
ordination_alb <- phyloseq::plot_ordination(physeq_alb, nmds_alb, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Picea abies")) +
  ggplot2::labs(subtitle = "(A) Swabian Alb") +
  ggplot2::theme(legend.position = "none",
        axis.text = ggplot2::element_text(size = 15),
        axis.title = ggplot2::element_text(size = 15))

ordination_alb

# Ordinate the Schorfheide-Chorin phyloseq using an NMDS with Bray-Curtis distance.
nmds_sch <- phyloseq::ordinate(physeq_sch, method = "NMDS", distance = "bray")

# Plot the ordination.
ordination_sch <- phyloseq::plot_ordination(physeq_sch, nmds_sch, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkolivegreen4"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  ggplot2::labs(subtitle = "(B) Schorfheide-Chorin") +
  ggplot2::theme(legend.position = "none",
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15))

ordination_sch

# Run the ordination on the full dataset to grab a nice legend. 
full_ordination_ps <- filtered_physeq
sample_data(full_ordination_ps) <- data.frame(sample_data(full_ordination_ps)) %>% 
  dplyr::rename(habitat = substrate)

nmds_full <- phyloseq::ordinate(full_ordination_ps, method = "NMDS", distance = "bray")

ordination_full <- phyloseq::plot_ordination(full_ordination_ps, nmds_full,
                                             type="samples", color="dominant_tree", shape="habitat") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green", "darkgreen", "darkolivegreen4"), name = "tree species",
                               labels = c("Fagus sylvatica","Picea abies", "Pinus sylvestris")) +
  ggplot2::theme(legend.text = ggplot2::element_text(face = "italic", size = 9),
        legend.title = ggplot2::element_text(size = 8, face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.spacing = ggplot2::unit(2,"cm"),
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15)) 

ordination_full

# Extract the legend and store it as a ggplot object. 
ordination_legend <- ggpubr::get_legend(ordination_full)

# Combine the figures into one plot. 
ordination_final <- ggpubr::ggarrange(ordination_alb, ordination_sch,
                                      ncol = 1, nrow = 2,
                                      legend = "right", legend.grob = ordination_legend)
ordination_final

ggsave('ordination_final.jpeg', device = 'jpeg',
       ordination_final, width = 180, height = 150,
        units = 'mm', dpi = 300)
#########################PERMANOVA and Betadisper############################
# Swabian Alb

bray_mat_alb <- phyloseq::distance(physeq_alb, method = "bray")
# Run the PERMANOVA analysis thorugh vegans adonis2() including both effects of habitat and tree species.
vegan::adonis2(bray_mat_alb ~ factor(phyloseq::sample_data(physeq_alb)$substrate) +
                 factor(phyloseq::sample_data(physeq_alb)$dominant_tree), by = 'margin')

# Test the within group dispersion for the substrate.
dispr_substrate_alb <- vegan::betadisper(bray_mat_alb, 
                                          factor(phyloseq::sample_data(physeq_alb)$substrate))
dispr_substrate_alb

permutest(dispr_substrate_alb)

# Test the within group dispersion for the tree Species.
dispr_tree_alb <- vegan::betadisper(bray_mat_alb, 
                                         factor(phyloseq::sample_data(physeq_alb)$dominant_tree))
dispr_tree_alb

permutest(dispr_tree_alb)

# Schorfheide-Chorin 

bray_mat_sch <- phyloseq::distance(physeq_sch, method = "bray")
# Run the PERMANOVA analysis thorugh vegans adonis2() including both effects of habitat and tree species.
vegan::adonis2(bray_mat_sch ~ factor(phyloseq::sample_data(physeq_sch)$substrate) +
                 factor(phyloseq::sample_data(physeq_sch)$dominant_tree), by = 'margin')

# Test the within group dispersion for the substrate.
dispr_substrate_sch <- vegan::betadisper(bray_mat_sch, 
                                         factor(phyloseq::sample_data(physeq_sch)$substrate))
dispr_substrate_sch

vegan::permutest(dispr_substrate_sch)

# Test the within group dispersion for the tree Species.
dispr_tree_sch <- vegan::betadisper(bray_mat_sch, 
                                    factor(phyloseq::sample_data(physeq_sch)$dominant_tree))
dispr_tree_sch

vegan::permutest(dispr_tree_sch)

```

  22. What other options can you find for distance measures between samples? 
  23. Are there other ordination methods apart from NMDS? 
  24. What about the compositionality of the sequencing data? 
  25. Are there statistically meaningful differences between substrates? Between tree species? 

  
  
# Questions 

 * what can you not do with this data 