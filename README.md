# Bisulfite-Seq-Analysis
Commands for Bisulfite Seq data analysis
#Installations
Go to 
FastQC v0.11.9 (Win/Linux zip file)
FastQC v0.11.9 (Mac DMG image)
Make sure the suitable java runtime environment (JRE) is installed 
java -version
Install fastqcn executable
Make the fastqc executable executable
chmod 755 fastqc
./fastqc SRRaccession.fastq.gz
# Check that cutadapt is installed: cutadapt --version, if not
https://cutadapt.readthedocs.io/en/stable/installation.html 
UNIX user use conda or pip
conda create -n cutadaptenv cutadapt
conda activate cutadaptenv
cutadapt --version

(mac use pip)
pip install --user --upgrade cutadapt 
or
sudo python3 -m venv /usr/local/cutadapt
sudo /usr/local/cutadapt/bin/pip install cutadapt
cd /usr/local/bin/
sudo ln -s ../cutadapt/bin/cutadapt

# Check that FastQC is installed fastqc -v
wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz
tar xvzf trim_galore.tar.gz
~/TrimGalore-0.6.6/trim_galore
conda install -c bioconda trim-galore

# Install necessary softwares:
tar xzf bismark_v0.X.Y.tar.gz
tar xvjf samtools-1.1.tar.bz2
cd /path_to_samtools/samtools-1.1
make

#SET $PATH VARIABLE TO CALL THE SOFTWARES DIRECTLY
#Example: Set path for Bowtie2 and Bismark
export PATH=$PATH:/path_to_software_packages/Bismark-0.22.3


#Run the pipeline
#Download Public data
sratoolkit (webversion) or CLI
#QC (You can also use GUI)
Fastq-dump –split-files SRR11207817
#Quality check data
/Users/ankitverma/Documents/tutorial/dollar_education/fastqc *.fastq
#Trimming
#cutadapt and TrimGalore must be installed and in path, 
/Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207817_1.fastq SRR11207817_2.fastq
/Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207820_1.fastq SRR11207820_2.fastq


#Filtered Reads
#SRR11207817_1_val_1.fastq SRR11207817_2_val_2.fastq
#SRR11207820_1_val_1.fastq SRR11207820_2_val_2.fastq


#Alignment
#Current working directory: alignment
#Software used: Bismark 
#Reference assembly: Ensembl Hg38 fasta (ref.fa and ref.fa.fai)
#Alignment tool: Bismark (default) and HISAT 
Step1: Genome Index Preparation
path/to/bismark/Bismark-0.22.3/bismark_genome_preparation <path_to_genome_folder>
Example #Generate indexes /enter index folder (only one time)
/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --path_to_aligner /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/ ./


Step2: Aligning reads to modified reference assembly
path_to_bismark/Bismark-0.22.3/bismark --genome /path_to_genome_folder/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq -o ./ -score_min L,0,-0.2 -X 500 -I 0


Example  #Alignment
/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/
/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207820_1_val_1.fq -2 SRR11207820_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/




Step3a: Sort and deduplicate aligned BAM
#Software used: SAMtools
samtools sort -n SRR11207817_1_val_1_bismark_bt2_pe.bam > SRR11207817_sortn.bam

Example #BAM Sorting
/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam
/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam

#Step3b:  Indexing
Example #Index bam
/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207817_1_val_1_bismark_bt2_pe.sort.bam
/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207820_1_val_1_bismark_bt2_pe.sort.bam


Step4: Sort and deduplicate aligned BAM
#Continue with bismark…
path_to_bismark/Bismark-0.22.3/deduplicate_bismark -o SRR11207817 --bam SRR11207817_sortn.bam 

Step5: Extraction of methylation call
path_to_bismark/Bismark-0.22.3/bismark_methylation_extractor --genome_folder /path_to_genome_folder/ -p --no_overlap --bedGraph --counts --buffer_size 10G --multicore 2 --cytosine_report --CX_context -o /path_to_output_folder/CX_report --gzip ./SRR11207817.deduplicated.bam
Example #Methylation call
/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207817_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index
/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207820_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index



Optional steps: 
#Install Perl module for outputting M-bias plot
#sudo apt install libgd-graph-perl
bismark2report  --dir . -o SRR11207817_bismark_report.txt --alignment_report /path_to_file/SRR11207817_1_val_1_bismark_bt2_PE_report.txt --dedup_report /path_to_file/SRR11207817_sortn.deduplication_report.txt --splitting_report /path_to_file/SRR11207817.deduplicated_splitting_report.txt --mbias_report /path_to_file/SRR11207817.deduplicated.M-bias.txt

## DIFFERENTIAL METHYLATION ANALYSIS:
FILE FORMATTING ON BASH
### Extract lines corresponding to CG context from CX_report.txt
zcat Sample_CX_report.txt.gz | awk '{if($4=="CG") print $0}' -| gzip > Sample_CG_report.txt.gz

### Extract chr21 and chr22 from PRACTICE DATASET provided
zcat GSM2971953_Test1.CpG_report.txt.gz | awk '{if($1=="chr21"||$2=="chr22") print $0}' -| gzip  > Test1_chr21_22_CpG_report.txt.gz

INSTALLATION OF R AND RSTUDIO
#Make sure that you installed R 4.2.2 and RStudio for Mac OSX : http://cran.r-project.org/bin/macosx/old/R-4.2.2.pkg for windows :http://cran.r-project.org/bin/windows/base/old/4.2.2/ For linux : http://cran.r-project.org/src/base/R-4/R-4.2.2.tar.gz 
http://www.rstudio.com/ide/download/
## PROVIDE R BASICS CHEATSHEET ###

INSTALLATION OF PACKAGES
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

#####################################################################################



### Format CG_report.txt to methylkit compatible format in R
Run format_CG_report.R in RStudio or R environment
OR,
Use the CG report directly as input for Methylkit.
## Follow the methylkit_run.R and run the script in R environment.
