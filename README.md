# Bisulfite-Seq Data Analysis
Authors: 
<br>
https://github.com/varshapk14 (varsh.pk@gmail.com)
</br>
<br>
https://github.com/ankitasks1 (ankitverma9079@gmail.com)
</br>
Commands for Bisulfite Seq data analysis

# Step 0: Install the required softwares 

# Install FastQC

Go to 

FastQC v0.11.9 (Win/Linux zip file)

FastQC v0.11.9 (Mac DMG image)

Make sure the suitable java runtime environment (JRE) is installed 

<code>java -version</code>

Install fastqc 

Make the fastqc executable executable

<code>chmod 755 fastqc</code>

<code>./fastqc SRRaccession.fastq.gz</code>

# Check that cutadapt is installed: 
cutadapt --version

if not
https://cutadapt.readthedocs.io/en/stable/installation.html 
<br>
UNIX user use conda or pip
</br>
<br>
<code>conda create -n cutadaptenv cutadapt</code>

<code>conda activate cutadaptenv</code>

<code>cutadapt --version</code>
</br>
<br>
(mac user can use pip)
</br>
<code>pip install --user --upgrade cutadapt </code>

or

<code>sudo python3 -m venv /usr/local/cutadapt</code>

<code>sudo /usr/local/cutadapt/bin/pip install cutadapt</code>

<code>cd /usr/local/bin/</code>

<code>sudo ln -s ../cutadapt/bin/cutadapt</code>

# Check that FastQC is installed and in path

<code>fastqc -v</code>

#Install Trimgalore
#source installation

<code>wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz</code>

<code>tar xvzf trim_galore.tar.gz</code>

<code>~/TrimGalore-0.6.6/trim_galore</code>
	
# conda installation trimgalore

<code>conda install -c bioconda trim-galore</code>

# Bowtie2 installation from source

Go to https://bowtie-bio.sourceforge.net/bowtie2/index.shtml and download latest version (from SourceForge)

unzip and executable is ready to be run

<code>unzip bowtie2-2.5.0-linux-x86_64.zip</code>


# Bismark installation from source
Go to Babraham Institute page of Bismark and download latest version

unzip and executable is ready to be run

<code>tar xzf bismark_v0.X.Y.tar.gz</code>

# Samtools installation from source

Go to samtools webpage and download both samtools and htslib, 

install both

<code>tar xvjf samtools-1.1.tar.bz2</code>

<code>cd samtools-1.1</code>    # and similarly for  htslib

<code>./configure --prefix=/where/to/install</code>

<code>make</code>

<code>make install</code>


#Note: YOU CAN ALSO SET $PATH VARIABLE TO CALL THE SOFTWARES DIRECTLY WITHOUT GIVING PATH EVERYTIME (recommended)

#Example: Set path for Bowtie2 and Bismark

<code> export PATH=$PATH:/path_to_software_packages/bowtie2-2.5.0-linux-x86_64/ </code>

<code> export PATH=$PATH:/path_to_software_packages/Bismark-0.22.3/ </code>

# INSTALLATION OF R AND RSTUDIO
#Make sure that you installed R 4.2.2 and RStudio 

For Mac OSX : http://cran.r-project.org/bin/macosx/old/R-4.2.2.pkg 

For Windows :http://cran.r-project.org/bin/windows/base/old/4.2.2/ 

For Linux : http://cran.r-project.org/src/base/R-4/R-4.2.2.tar.gz 

For Rstudio any OS go to http://www.rstudio.com/ide/download/ and install as per developers instructions

## PROVIDED  R BASICS CHEATSHEET ##
#Ways to install packages

CRAN: http://cran.r-project.org/ 

Bioconductor: http://bioconductor.org/ 

Github: http://github.com/

#Install from CRAN

install.packages("devtools") 

#Install from the source: 

install.packages("devtools_1.1.tar.gz",repos=NULL, type="source") 

#Install from Bioconductor

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("methylKit")
	
	
#------------------------------------ Run the pipeline --------------------------------------#

# Step 1: Download Public data

sratoolkit (webversion) or CLI
<br>
fastq-dump --split-files SRR11207817</br>
<br>
fastq-dump --split-files SRR11207820</br>

:--skip this step if data is already available--:

# Step 2: Quality check data (You can also use GUI)

#Quality check data

<code> /Users/ankitverma/Documents/tutorial/dollar_education/fastqc *.fastq </code>

# Step 3: Trimming
#cutadapt and TrimGalore must be installed and in path, 

<code> /Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207817_1.fastq SRR11207817_2.fastq </code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207820_1.fastq SRR11207820_2.fastq</code>

#Filtered Reads

#SRR11207817_1_val_1.fastq SRR11207817_2_val_2.fastq

#SRR11207820_1_val_1.fastq SRR11207820_2_val_2.fastq


# Get reference genome

4 possible sources (UCSC, Gencode, NCBI, Ensembl)
<img src="Screenshot 2022-12-19 at 13.39.05.png" alt="Get reference">

<br>
<b>For Mouse and Human, I generally choose Gencode. For other species you can select any other assembly</b>
</br>


# Step 4: Genome Index Preparation
	
<code> path/to/bismark/Bismark-0.22.3/bismark_genome_preparation <path_to_genome_folder> </code>
	
Example: Generate indexes of reference genome  (only one time)

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --path_to_aligner /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/ . </code>

# Step 5: Alignment / Mapping

Current working directory: your current directory (./)

Software used: Bismark 

Reference assembly: Ensembl hg38 fasta (ref.fa and ref.fa.fai)

Alignment tool: Bismark (default) 

<code>path_to_bismark/Bismark-0.22.3/bismark --genome /path_to_genome_folder/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq -o ./ -score_min L,0,-0.2 -X 500 -I 0 </code>


Example

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/ </code>


<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207820_1_val_1.fq -2 SRR11207820_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/</code>


# Step 6: Sort aligned BAM

Software used: SAMtools

<code>samtools sort -n SRR11207817_1_val_1_bismark_bt2_pe.bam > SRR11207817_sortn.bam</code>

Example #BAM Sorting

<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam</code>

# Step 7:  Indexing

Example Index bam

<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207817_1_val_1_bismark_bt2_pe.sort.bam</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207820_1_val_1_bismark_bt2_pe.sort.bam</code>


# Step 8: Deduplicate aligned BAM (must for WGBS)

<code>path_to_bismark/Bismark-0.22.3/deduplicate_bismark -o SRR11207817 --bam SRR11207817_sortn.bam </code>

# Step 9: Extraction of methylation call

<code>path_to_bismark/Bismark-0.22.3/bismark_methylation_extractor --genome_folder /path_to_genome_folder/ -p --no_overlap --bedGraph --counts --buffer_size 10G --multicore 2 --cytosine_report --CX_context -o /path_to_output_folder/CX_report --gzip ./SRR11207817.deduplicated.bam
Example #Methylation call</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207817_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207820_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index</code>


Optional steps: 

Install Perl module for outputting M-bias plot

sudo apt install libgd-graph-perl

<code>bismark2report  --dir . -o SRR11207817_bismark_report.txt --alignment_report /path_to_file/SRR11207817_1_val_1_bismark_bt2_PE_report.txt --dedup_report /path_to_file/SRR11207817_sortn.deduplication_report.txt --splitting_report /path_to_file/SRR11207817.deduplicated_splitting_report.txt --mbias_report /path_to_file/SRR11207817.deduplicated.M-bias.txt</code>

# Step 10: Differential Methylation Analysis:
	
# methylkit
	
<b><u>Note: Specific format is required as input to methylkit</u></b>

<u>Extract lines corresponding to CG context from CX_report.txt</u>
<code>zcat Sample_CX_report.txt.gz | awk '{if($4=="CG") print $0}' -| gzip > Sample_CG_report.txt.gz</code>

<u>Extract chr21 and chr22 from PRACTICE DATASET provided</u>

<code>zcat GSM2971953_Test1.CpG_report.txt.gz | awk '{if($1=="chr21"||$2=="chr22") print $0}' -| gzip  > Test1_chr21_22_CpG_report.txt.gz</code



#Format CG_report.txt to methylkit compatible format in R

Run format_CG_report.R in RStudio or R environment

OR,

Use the CG report directly as input for Methylkit.

## Follow the methylkit_run.R and run the script in R environment.

In brief, follow this script
setwd("/mnt/home3/outfolder/methylkit")
Reading the methylation call files
library(methylKit)
Import deduplicated, sorted, BAM files
Store bam files in list
bam_files_list <- as.list(list.files(path = "/mnt/home3/outfolder/bismark_deduplicated/",
                                     pattern = "\\.deduplicated.sorted.bam$",
                                     full.names = TRUE))


List of sample IDs
sample_ids_list <- lapply(str_split(bam_files_list,"/"), function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam","",x[11]))


Set minimum CpG coverage desired
Used in processBismarkAln function
min_coverage <- 5

Set minimum quality of reads
Used in processBismarkAln function
min_qual <- 20


Set minimum methylation percentage difference between groups
Used in getMethylDiff function; 25 is the default value.
dml_diffs <- 25


#Import from BAM
Get methylation stats for CpGs with at least min_coverage coverage
<code>
meth_stats <- processBismarkAln(location = bam_files_list,
                                sample.id = sample_ids_list,
                                assembly = "GRCm39.fa ",
                                save.folder = NULL, 
                                save.context = c("CpG"), 
                                read.context = "CpG",
                                mincov = min_coverage,
                                minqual = min_qual, 
                                phred64 = FALSE, Your data is phred33,
                                treatment = rep(0, length(bam_files_list)),
                                save.db = FALSE)
</code>

<u>File count</u>
nFiles <- length(bam_files_list)
	
###OR you can import from CpG report produced by Bismark
	
#Import from CpG report
First convert .cov.gz file to .CpG.report.txt.gz using convert_cov_to_CpGreport.sh (custom)
test
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

test_report_list <- as.list(list.files(path = "/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage",
                                pattern = "\\.cov.CpG_test_report.txt.gz",
                                full.names = TRUE))

<code>
for (content1 in test_report_list){
  print(content1)
  temp1 <- gsub("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/", "", content1)
  temp1 <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_test_report.txt.gz", "", temp1)
  print(temp1)
  test_report <- fread(content1)
  colnames(test_report) <- c("chr", "base", "strand", "methylated", "unmethylated", "Ccontext", "trinu_context")
  test_report <- data.frame(test_report)
  test_report["chrBase"] <- paste0(test_report$chr,
                                   ".",
                                   test_report$base)

  test_report["coverage"] <- test_report$methylated + test_report$unmethylated
  test_report["freqC"] <- (test_report$methylated  *  100) / (test_report$methylated + test_report$unmethylated)
  test_report["freqT"] <- (test_report$unmethylated  *  100) / (test_report$methylated + test_report$unmethylated)
  test_report <- test_report[,c(8,1,2,3,9,10,11)]
  write.table(test_report, paste0("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/",temp1,".CpG_test_report.txt"), sep="\t", quote = F, append=F, row.names = F, col.names = T)
}
	</code>
Create a list of myCpG_report.txt
test_list <- lapply(test_report_list, function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov","",(gsub(".gz","",x))))
test.ids <- lapply(test_list, function(x) gsub("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/","",(gsub(".CpG_test_report.txt","",x))))

List of sample IDs
sample_ids_reports_list <- lapply(str_split(CpG.reports.list,"/"), function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_report.txt.gz","",x[11]))


read the files to a methylRawList object: myobj
myobj=methRead(test_list,
               sample.id=test.ids,
               assembly="GRCm39",
               treatment=rep(0, length(test_list)),
               context="CpG",
               mincov = 5
)

---Now run the script on real data----
samples
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

CpG.reports.list <- as.list(list.files(path = "/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage",
                                       pattern = "\\.cov.CpG_report.txt.gz",
                                       full.names = TRUE))

File count
nFiles <- length(CpG.reports.list)

for (report_content in CpG.reports.list){
  print(report_content)
  temp_content <- gsub("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/", "", report_content)
  temp_content <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_report.txt.gz", "", temp_content)
  print(temp_content)
  report <- fread(report_content)
  colnames(report) <- c("chr", "base", "strand", "methylated", "unmethylated", "Ccontext", "trinu_context")
  report <- data.frame(report)
  report["chrBase"] <- paste0(report$chr,
                              ".",
                              report$base)
  
  report["coverage"] <- report$methylated + report$unmethylated
  report["freqC"] <- (report$methylated  *  100) / (report$methylated + report$unmethylated)
  report["freqT"] <- (report$unmethylated  *  100) / (report$methylated + report$unmethylated)
  report <- report[,c(8,1,2,3,9,10,11)]
  write.table(report, paste0("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/",temp_content,".CpG_report.txt"), sep="\t", quote = F, append=F, row.names = F, col.names = T)
  rm(report)
  rm(temp_content)
}
Create a list of myCpG_report.txt
mysamples_list <- lapply(CpG.reports.list, function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov","",(gsub(".gz","",x))))
samples.id <- lapply(mysamples_list, function(x) gsub("/mnt/home3/outfolder/bismark_methylation_calls/methylation_coverage/","",(gsub(".CpG_report.txt","",x))))



read the files to a methylRawList object: myobj
myobjDB=methRead(mysamples_list,
               sample.id=samples.id,
               assembly="GRCm39",
               treatment=rep(0, length(mysamples_list)),
               context="CpG",
               mincov = 5,
               dbtype = "tabix",
               dbdir = "methylDB")

print(myobjDB[[1]]@dbpath)

Generate and save histograms showing Percent CpG Methylation
png("getMethylationStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getMethylationStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) Get %CpG methylation information
}
dev.off()

Generate and save histograms showing CpG Methylation Coverage
png("getCoverageStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getCoverageStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) Get %CpG methylation information
}
dev.off()


Filtering samples based on read coverage
filtered.myobjDB <- filterByCoverage(myobjDB, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

Normalize coverage
norm.filt.objDB <- normalizeCoverage(filtered.myobjDB,method="median")


Comparative analysis
Merging samples
Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.norm.filt.objDB <- unite(norm.filt.objDB, destrand=FALSE)
head(merge.norm.filt.objDB)

By default, unite function produces bases/regions covered in all samples. 
That requirement can be relaxed using “min.per.group” option in unite function.
This creates a methylBase object, where only CpGs covered with at least 1 sample per group will be returned
there were two groups defined by the treatment vector, 
given during the creation of myobj: treatment=c(1,1,0,0)
Not doing it meth.min=unite(myobj,min.per.group=1L)

Sample Correlation
png("getCorrelation.png", height = 2000, width = 2000)
getCorrelation(merge.norm.filt.objDB,plot=TRUE)
dev.off()

Clustering dendrogram
png("dendrogram_path.png", height = 600, width = 1000) 
clusterSamples(merge.norm.filt.objDB, dist="correlation", method="ward", plot=TRUE)
dev.off()

Run a PCA analysis on percent methylation for all samples
png("pca_path.png", height = 1000, width = 1000) 
PCASamples(merge.norm.filt.objDB)
dev.off()

Run the PCA analysis and plot variances against PC number in a screeplot
png("scree_path.png", height = 1000, width = 1000)
PCASamples(merge.norm.filt.objDB, screeplot = TRUE)
dev.off()


merge.norm.filt.objDB_data <- getData(merge.norm.filt.objDB)
colnames(merge.norm.filt.objDB_data) <- c(colnames(merge.norm.filt.objDB_data)[1:4], 
                                          paste0(unlist(lapply(getSampleID(merge.norm.filt.objDB), function(x) rep(x,3))), "_", 
                                                 colnames(getData(merge.norm.filt.objDB)[5:52])))

For some situations, it might be desirable to summarize methylation information over tiling windows rather than doing base-pair resolution analysis.
Tiling windows analysis
myobjDB.low <- methRead(mysamples_list,
                 sample.id=samples.id,
                 assembly="GRCm39",
                 treatment=rep(0, length(mysamples_list)),
                 context="CpG",
                 mincov = 2)

tiles = tileMethylCounts(myobjDB.low,win.size=1000,step.size=1000,cov.bases = 10)
head(tiles[[1]],3)


Finding differentially methylated bases or regions
myDiff=calculateDiffMeth(merge.norm.filt.objDB)


get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")


get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")


get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)



Finding differentially methylated bases using multiple-cores
myDiff=calculateDiffMeth(meth,mc.cores=2)

Annotating differentially methylated bases or regions
library(genomation)

read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                            package = "methylKit"))
annotate differentially methylated CpGs with 
promoter/exon/intron using annotation data

annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)


Read the CpG island annotation and annotate our differentially methylated bases/regions with them.
read the shores and flanking regions and name the flanks as shores 
and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))

convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

read the CpG island annotation and annotate our differentially methylated bases/regions with them.
read the shores and flanking regions and name the flanks as shores 
and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))

convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

Regional analysis
summarize methylation information over a set of defined regions such as promoters or CpG islands.
promoters=regionCounts(myobj,gene.obj$promoters)

head(promoters[[1]])

get the distance to TSS and nearest gene name
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

It is also desirable to get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")
percentage of differentially methylated bases are on CpG islands, CpG island shores and other regions.
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="differential methylation annotation")
percentage of intron/exon/promoters that overlap with differentially methylated bases.
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

#################   END OF ANALYSIS  ###################
