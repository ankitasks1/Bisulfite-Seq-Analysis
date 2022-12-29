# Bisulfite-Seq Data Analysis
Authors: 
<br>
https://github.com/varshapk14 (varsh.pk@gmail.com)
</br>
<br>
https://github.com/ankitasks1 (ankitverma9079@gmail.com)
</br>
Commands for Bisulfite Seq data analysis

#get reference
<img src="Screenshot 2022-12-19 at 13.39.05.png" alt="Get regerence">

<br>
<b>For Mouse and Human, I generally choose Gencode. For other species you can select any other assembly</b>
</br>

<br>
</br>
<br>
</br>
#------------------------------------ Install the required softwares --------------------------------------#

# Step 0: Installations
#sratoolkit : easy to install if your using command line version other wise webversion is also available 

# FastQC

Go to 

FastQC v0.11.9 (Win/Linux zip file)

FastQC v0.11.9 (Mac DMG image)

Make sure the suitable java runtime environment (JRE) is installed 

<code>java -version</code>

#Install fastqc 

Make the fastqc executable executable

<code>chmod 755 fastqc</code>

<code>./fastqc SRRaccession.fastq.gz</code>
# Cutadapt, required by TrimHalore

Check that cutadapt is installed: 
<code>cutadapt --version</code> if not
https://cutadapt.readthedocs.io/en/stable/installation.html 
UNIX user use conda or pip

<code>conda create -n cutadaptenv cutadapt</code>

<code>conda activate cutadaptenv</code>

<code>cutadapt --version</code>

(mac use pip)
<code>pip install --user --upgrade cutadapt </code>
or
<code>sudo python3 -m venv /usr/local/cutadapt</code>
<code>sudo /usr/local/cutadapt/bin/pip install cutadapt</code>
<code>cd /usr/local/bin/</code>
<code>sudo ln -s ../cutadapt/bin/cutadapt</code>

# Check that FastQC is installed 
<code>fastqc -v</code>

#Install TrimGalore

<code>wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz</code>

<code>tar xvzf trim_galore.tar.gz</code>

<code>~/TrimGalore-0.6.6/trim_galore

<code>conda install -c bioconda trim-galore</code>



# Bismark
<code>tar xzf bismark_v0.X.Y.tar.gz</code>

# Samtools

<code>tar xvjf samtools-1.1.tar.bz2</code>

<code>cd /path_to_samtools/samtools-1.1</code>

<code>make</code>

#SET $PATH VARIABLE TO CALL THE SOFTWARES DIRECTLY

#Example: Set path for Bowtie2 and Bismark

<code> export PATH=$PATH:/path_to_software_packages/Bismark-0.22.3 </code>

 #INSTALLATION OF R AND RSTUDIO
#Make sure that you installed R 4.2.2 and RStudio for Mac OSX : http://cran.r-project.org/bin/macosx/old/R-4.2.2.pkg for windows :http://cran.r-project.org/bin/windows/base/old/4.2.2/ For linux : http://cran.r-project.org/src/base/R-4/R-4.2.2.tar.gz 
http://www.rstudio.com/ide/download/
## PROVIDE R BASICS CHEATSHEET ###

#INSTALLATION OF PACKAGES
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

# Step 1: Download Public data (optional)

sratoolkit (webversion) or CLI



fastq-dump –split-files SRR11207817
Note: If you already have data move to Step 2 and Skip Step 1

#Quality check data : FastQC

<code> /Users/ankitverma/Documents/tutorial/dollar_education/fastqc *.fastq </code>

# Step 2: Trimming
#cutadapt and TrimGalore must be installed and in path, 
<code> /Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207817_1.fastq SRR11207817_2.fastq </code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207820_1.fastq SRR11207820_2.fastq</code>

#Filtered Reads
#SRR11207817_1_val_1.fastq SRR11207817_2_val_2.fastq
#SRR11207820_1_val_1.fastq SRR11207820_2_val_2.fastq

	
# Step3: Genome Index Preparation
	
<code>path/to/bismark/Bismark-0.22.3/bismark_genome_preparation <path_to_genome_folder></code>
	
Example #Generate indexes /enter index folder (only one time)

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --path_to_aligner /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/ ./</code>


# Step4: Alignment
#Current working directory: your current directory (./)
#Software used: Bismark 
#Reference assembly: Ensembl Hg38 fasta (ref.fa and ref.fa.fai)
#Alignment tool: Bismark (default) and HISAT 

<code>path_to_bismark/Bismark-0.22.3/bismark --genome /path_to_genome_folder/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq -o ./ -score_min L,0,-0.2 -X 500 -I 0<c/ode>


Example  #Alignment
<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207817_1_val_1.fq -2 SRR11207817_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark --genome index/ -1 SRR11207820_1_val_1.fq -2 SRR11207820_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/bowtie2-2.5.0-macos-arm64/</code>


# Step5: Sort aligned BAM
#Software used: SAMtools
<code>samtools sort -n SRR11207817_1_val_1_bismark_bt2_pe.bam > SRR11207817_sortn.bam</code>

Example #BAM Sorting
<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools sort -o SRR11207817_1_val_1_bismark_bt2_pe.sort.bam SRR11207817_1_val_1_bismark_bt2_pe.bam</code>

# Step6:  Indexing
Example #Index bam
<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207817_1_val_1_bismark_bt2_pe.sort.bam</code>
<code>/Users/ankitverma/Documents/tutorial/dollar_education/samtools index SRR11207820_1_val_1_bismark_bt2_pe.sort.bam</code>


# Step7: Sort and deduplicate aligned BAM
#Continue with bismark…
<code>path_to_bismark/Bismark-0.22.3/deduplicate_bismark -o SRR11207817 --bam SRR11207817_sortn.bam </code>

# Step8: Extraction of methylation call
<code>path_to_bismark/Bismark-0.22.3/bismark_methylation_extractor --genome_folder /path_to_genome_folder/ -p --no_overlap --bedGraph --counts --buffer_size 10G --multicore 2 --cytosine_report --CX_context -o /path_to_output_folder/CX_report --gzip ./SRR11207817.deduplicated.bam
Example #Methylation call</code>
<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207817_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index</code>

<code>/Users/ankitverma/Documents/tutorial/dollar_education/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207820_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index</code>


Optional steps: 
#Install Perl module for outputting M-bias plot
#sudo apt install libgd-graph-perl
<code>bismark2report  --dir . -o SRR11207817_bismark_report.txt --alignment_report /path_to_file/SRR11207817_1_val_1_bismark_bt2_PE_report.txt --dedup_report /path_to_file/SRR11207817_sortn.deduplication_report.txt --splitting_report /path_to_file/SRR11207817.deduplicated_splitting_report.txt --mbias_report /path_to_file/SRR11207817.deduplicated.M-bias.txt</code>

# Step 9: DIFFERENTIAL METHYLATION ANALYSIS:
FILE FORMATTING ON BASH
### Extract lines corresponding to CG context from CX_report.txt
<code>zcat Sample_CX_report.txt.gz | awk '{if($4=="CG") print $0}' -| gzip > Sample_CG_report.txt.gz</code>

### Extract chr21 and chr22 from PRACTICE DATASET provided
<code>zcat GSM2971953_Test1.CpG_report.txt.gz | awk '{if($1=="chr21"||$2=="chr22") print $0}' -| gzip  > Test1_chr21_22_CpG_report.txt.gz</code



### Format CG_report.txt to methylkit compatible format in R
Run format_CG_report.R in RStudio or R environment
OR,
Use the CG report directly as input for Methylkit.
## Follow the methylkit_run.R and run the script in R environment.
	
	

#####################################################################################
