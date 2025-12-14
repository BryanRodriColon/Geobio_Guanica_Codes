#!/bin/bash

# QIIME2 cutadapt script assuming the data is demultiplexed and converted to an artifact. Im also including the FastQC script here

# Load Conda environment and load QIIME2

module load conda
conda activate qiime2-amplicon-2024.10

# Parse command-line arguments

# Define variables from command-line arguments
qiimedatapath="/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019"
demuxfile="SciRep_Data_demux.qza"
cores=4
Fadapter="GGATTAGATACCCBDGTAGTC"
forwardprimer="CCTACGGGNGGCWGCAG" #515F-Y
Radapter="CTGCWGCCNCCCGTAGG"
reverseprimer="GACTACHVGGGTATCTAATCC" #806R-Y
trimmedfolder="02_trimmed_files"


echo "QIIME data path: $qiimedatalspath"
echo "Demux file: $demuxfile"
echo "Trimmed output folder: $trimmedfolder"
echo "Number of cores: $cores"
echo "Forward adapter: $Fadapter"
echo "Forward primer: $forwardprimer"
echo "Reverse adapter: $Radapter"
echo "Reverse primer: $reverseprimer"

# Change to the QIIME data path
cd "$qiimedatapath" || { echo "Error: Failed to change directory to $qiimedatapath"; exit 1; }
echo "Changed directory to: $qiimedatapath"


# Run Cutadapt to trim primers and adapters (16S)
echo "Running Cutadapt to trim primers and adapters..."
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/01_demux_files/SciRep_Data_demux.qza" \
  --p-adapter-f "$Fadapter" \
  --p-front-f "$forwardprimer" \
  --p-adapter-r "$Radapter" \
  --p-front-r "$reverseprimer" \
  --o-trimmed-sequences "$trimmedfolder/${demuxfile%.qza}_trimmed.qza" \
  --verbose
echo "Cutadapt trimming completed."

# Summarize the trimmed data
echo "Summarizing trimmed data..."
qiime demux summarize \
   --i-data "$trimmedfolder/${demuxfile%.qza}_trimmed.qza" \
  --o-visualization "$trimmedfolder/${demuxfile%.qza}_trimmed.qzv"
echo "Summary of trimmed data generated."

#Extract the trimmed fastq files our of the artifact so they can also be examined with FastQC
qiime tools export \
  --input-path 02_trimmed_files/SciRep-Data-demux_trimmed.qza \
  --output-path 02_trimmed_files/exported_trimmed_data




echo "Cutadapt DONE!"

#Deactiave qiime2 

conda deactivate

## Activate fastqc to generate the reports

module load fastqc

# Define directories
TRIMMED_DIR="/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/02_trimmed_files/exported_trimmed_data"
QC_DIR="/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/02_trimmed_files/fastQC_reports"

# Run FastQC on each trimmed file and save reports to QC_DIR
for FILE in $TRIMMED_DIR/*.fastq.gz
do
    echo "Running FastQC on $FILE"
    fastqc -o $QC_DIR $FILE
done

echo "FastQC quality check completed for all samples."


#Go back to QIIME2

fastqc deactivate

module load conda/latest

conda activate qiime2-amplicon-2024.10
