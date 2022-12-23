echo "\n  ----------------------------------------- \n"
echo "    Description: Bisulfite-Sequencing Analysis Wrapper Script "
echo "    Script compiled by Ankit Verma"
echo "    Contact: ankitverma9079@gmail.com  "
echo "\n  ----------------------------------------- \n"

sampleid="SRR11207820" 

echo "Downloading data .... $sampleid"
fastq-dump SRR11207820 --split-files

echo "Performing Quality Check ..... $sampleid"

/Users/ankitverma/Documents/tutorial/dollar_education/fastqc *.fastq

echo "Trimming data ..... $sampleid"

/Users/ankitverma/Documents/tutorial/dollar_education/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /usr/local/bin/cutadapt --length 36 --paired SRR11207820_1.fastq SRR11207820_2.fastq

#Preparing genome (only once)
#/Users/ankitverma/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --path_to_aligner /Users/ankitverma/Documents/tutorial/dollar_education/softwares/bowtie2-2.5.0-macos-arm64/ ./

echo "Aligning data  to custom reference genome ..... $sampleid"
/Users/ankitverma/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark --genome index/ -1 SRR11207820_1_val_1.fq -2 SRR11207820_2_val_2.fq --bowtie2 --path_to_bowtie2 /Users/ankitverma/Documents/tutorial/dollar_education/softwares/bowtie2-2.5.0-macos-arm64/

echo "Sorting BAM..... $sampleid"
samtools sort -o SRR11207820_1_val_1_bismark_bt2_pe.sort.bam SRR11207820_1_val_1_bismark_bt2_pe.bam

echo "Indexing BAM ..... $sampleid"
samtools index SRR11207820_1_val_1_bismark_bt2_pe.sort.bam

echo "Calling methylation ..... $sampleid"

/Users/ankitverma/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark_methylation_extractor --cutoff 10  SRR11207820_1_val_1_bismark_bt2_pe.sort.bam --cytosine_report --bedGraph  --genome index

