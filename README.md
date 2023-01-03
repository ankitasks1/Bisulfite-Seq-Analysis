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

fastq-dump --split-files SRR11207817

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

<br>
</br>
<br>
</br>



# Step 4: Genome Index Preparation
	
<code>path/to/bismark/Bismark-0.22.3/bismark_genome_preparation <path_to_genome_folder></code>
	
Example Generate indexes /enter index folder (only one time)

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --path_to_aligner /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/ ./ </code>

# Step 5: Alignment / Mapping

Current working directory: your current directory (./)

Software used: Bismark 

Reference assembly: Ensembl Hg38 fasta (ref.fa and ref.fa.fai)

Alignment tool: Bismark (default) and HISAT 

<code>path_to_bismark/Bismark-0.22.3/bismark --genome /path_to_genome_folder/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq -o ./ -score_min L,0,-0.2 -X 500 -I 0<c/ode>


Example  #Alignment

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/</code>

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
	
# specific format is required as input to methylkit

Extract lines corresponding to CG context from CX_report.txt
<code>zcat Sample_CX_report.txt.gz | awk '{if($4=="CG") print $0}' -| gzip > Sample_CG_report.txt.gz</code>

Extract chr21 and chr22 from PRACTICE DATASET provided

<code>zcat GSM2971953_Test1.CpG_report.txt.gz | awk '{if($1=="chr21"||$2=="chr22") print $0}' -| gzip  > Test1_chr21_22_CpG_report.txt.gz</code



#Format CG_report.txt to methylkit compatible format in R

Run format_CG_report.R in RStudio or R environment

OR,

Use the CG report directly as input for Methylkit.

## Follow the methylkit_run.R and run the script in R environment.
	
#####################################################################################
